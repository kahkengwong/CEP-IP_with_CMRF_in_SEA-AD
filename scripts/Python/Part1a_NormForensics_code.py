# =============================================================================
# Part1a_NormForensics_code.py
# SEA-AD CEP-IP/CMRF, Part 1a, normalisation of the expression matrix
#   Identifies the transform behind adata.X from the UMI counts: library-size
#   scaling or per-cell size factors, the log base, the target sum, exact
#   reconstruction, zero concordance and the source of the denominator
#   (tests T1 to T9).
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT)
#   python3 Part1a_NormForensics_code.py
#
# Outputs (Part1_PreChecks/): S1a_norm_forensics.txt (report).
# =============================================================================

import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import scipy.sparse as sp
from scipy.stats import spearmanr, pearsonr
import scanpy as sc

np.random.seed(42)
LN2 = np.log(2.0)

os.chdir(os.environ.get("SEAAD_ROOT", os.getcwd()))
FULL_H5AD = ("SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
             "RNAseq_final-nuclei.2026-06-22.h5ad")
OUT_DIR = "Part1_PreChecks"   # created in the SEA-AD folder
os.makedirs(OUT_DIR, exist_ok=True)

N_CELLS = 4000          # sample for row-wise tests
N_DENSE = 200           # smaller sample for per-gene ratio tests

_rep = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); _rep.append(s)
def hdr(t):
    say(""); say("=" * 74); say("  " + t); say("=" * 74)


hdr("Loading")
ad = sc.read_h5ad(FULL_H5AD)
n, g = ad.shape
say(f"  shape          : {ad.shape}")
say(f"  uns['log1p']   : {ad.uns.get('log1p', '(absent)')}")
say(f"  layers         : {[k for k in (ad.layers.keys() if ad.layers else []) if k]}")

CL = None
for c in ["UMIs", "counts", "raw_counts"]:
    if ad.layers and c in ad.layers:
        CL = c; break
if CL is None:
    say("  [FATAL] no count layer found - cannot run forensics"); raise SystemExit

pick = np.sort(np.random.choice(n, min(N_CELLS, n), replace=False))
X = ad.X[pick]
C = ad.layers[CL][pick]
if not sp.issparse(X): X = sp.csr_matrix(X)
if not sp.issparse(C): C = sp.csr_matrix(C)
X = X.tocsr().astype(np.float64)
C = C.tocsr().astype(np.float64)

depth_all = np.asarray(ad.layers[CL].sum(axis=1)).ravel().astype(float)
del ad
depth = depth_all[pick]
say(f"  sampled cells  : {len(pick)}")
say(f"  median depth (this file, all cells) : {np.median(depth_all):,.2f}")


# == T1: row sums of expm1(X) =================================================
hdr("T1 - row-sum of expm1(X): constant  <=>  library-size scaling")
# expm1 on a sparse matrix preserves sparsity (expm1(0)=0), so this is cheap.
rs_e = np.asarray(X.expm1().sum(axis=1)).ravel()
cv_e = rs_e.std() / rs_e.mean()
say(f"  base-e  row-sum : min={rs_e.min():,.3f}  median={np.median(rs_e):,.3f}  "
    f"max={rs_e.max():,.3f}")
say(f"  base-e  CV      : {cv_e:.3e}")

X2 = X.copy(); X2.data = X2.data * LN2          # 2^X - 1 = expm1(X*ln2)
rs_2 = np.asarray(X2.expm1().sum(axis=1)).ravel()
cv_2 = rs_2.std() / rs_2.mean()
say(f"  base-2  row-sum : min={rs_2.min():,.3f}  median={np.median(rs_2):,.3f}  "
    f"max={rs_2.max():,.3f}")
say(f"  base-2  CV      : {cv_2:.3e}")

LIBSIZE = min(cv_e, cv_2) < 1e-6
BASE = "e" if cv_e <= cv_2 else "2"
say("")
if LIBSIZE:
    T_hat = np.median(rs_e if BASE == "e" else rs_2)
    say(f"  >>> CONSTANT under base-{BASE}. Library-size scaling, target_sum "
        f"= {T_hat:,.2f}")
else:
    T_hat = np.nan
    say("  >>> not constant under either base")


# == T2: log base =============================================================
hdr("T2 - log base: e (scanpy log1p) vs 2 (scater logNormCounts)")
# For a single cell, X = log_b(1 + c/sf)  =>  (b^X - 1)/c = 1/sf for EVERY gene.
# So the correct base gives near-zero CV of that ratio within each cell.
sub = np.sort(np.random.choice(len(pick), min(N_DENSE, len(pick)), replace=False))
Xd = np.asarray(X[sub].todense())
Cd = np.asarray(C[sub].todense())
cv_row_e, cv_row_2 = [], []
for i in range(Xd.shape[0]):
    m = (Cd[i] > 0) & (Xd[i] > 0)
    if m.sum() < 20:
        continue
    r_e = np.expm1(Xd[i][m]) / Cd[i][m]
    r_2 = (np.exp2(Xd[i][m]) - 1.0) / Cd[i][m]
    cv_row_e.append(r_e.std() / r_e.mean())
    cv_row_2.append(r_2.std() / r_2.mean())
cv_row_e = np.array(cv_row_e); cv_row_2 = np.array(cv_row_2)
say(f"  cells tested                     : {len(cv_row_e)}")
say(f"  median within-cell CV, base e    : {np.median(cv_row_e):.3e}")
say(f"  median within-cell CV, base 2    : {np.median(cv_row_2):.3e}")
say("")
if np.median(cv_row_e) < 1e-8:
    BASE = "e"
    say("  >>> base e")
elif np.median(cv_row_2) < 1e-8:
    BASE = "2"
    say("  >>> base 2")
else:
    say("  >>> neither base")


# == T3: size factors and reconstruction ======================================
hdr("T3 - implied per-cell size factor, and exact reconstruction")
def implied_sf(Xd, Cd, base):
    out = np.full(Xd.shape[0], np.nan)
    for i in range(Xd.shape[0]):
        m = (Cd[i] > 0) & (Xd[i] > 0)
        if m.sum() < 5:
            continue
        num = np.expm1(Xd[i][m]) if base == "e" else (np.exp2(Xd[i][m]) - 1.0)
        out[i] = np.median(Cd[i][m] / num)
    return out

sf = implied_sf(Xd, Cd, BASE)
ok = np.isfinite(sf)
say(f"  implied sf: n={ok.sum()}  min={sf[ok].min():.6f}  "
    f"median={np.median(sf[ok]):.6f}  max={sf[ok].max():.6f}")

d_sub = depth[sub]
if LIBSIZE and np.isfinite(T_hat):
    sf_pred = d_sub / T_hat
    rel = np.abs(sf[ok] - sf_pred[ok]) / sf_pred[ok]
    say(f"  |sf - depth/T| / (depth/T):  max={rel.max():.3e}  "
        f"mean={rel.mean():.3e}")
    if rel.max() < 1e-6:
        say("  >>> sf = depth / target_sum")

    Xhat = np.log1p(Cd / d_sub[:, None] * T_hat) if BASE == "e" \
        else np.log2(1.0 + Cd / d_sub[:, None] * T_hat)
    err = np.abs(Xd - Xhat)
    say(f"  max |X - reconstruction|  = {err.max():.3e}")
    say(f"  mean|X - reconstruction|  = {err.mean():.3e}")
    if err.max() < 1e-6:
        say("  >>> exact")
    else:
        say("  >>> not exact")


# == T4: target sum against median depth ======================================
hdr("T4 - is target_sum the median depth of this file? (normalize_total default)")
if LIBSIZE and np.isfinite(T_hat):
    med_here = np.median(depth_all)
    say(f"  recovered target_sum     : {T_hat:,.2f}")
    say(f"  median depth, this file  : {med_here:,.2f}")
    say(f"  ratio                    : {T_hat/med_here:.6f}")
    for cand, lab in [(1e4, "1e4 (CP10K)"), (1e5, "1e5 (CP100K)"),
                      (1e6, "1e6 (CPM)"), (med_here, "median depth")]:
        tag = "  <== MATCH" if abs(T_hat - cand) / cand < 1e-6 else ""
        say(f"    vs {lab:<16}: {cand:,.2f}{tag}")
    say("")
    if abs(T_hat - med_here) / med_here < 1e-6:
        say("  >>> target_sum = median depth")
    else:
        say("  >>> target_sum differs from the median depth")
else:
    say("  (skipped: T1 says not library-size)")


# == T5: size factor against depth ============================================
hdr("T5 - size factor vs depth: the scran discriminator")
if ok.sum() > 50:
    rs_sp = spearmanr(sf[ok], d_sub[ok])[0]
    rs_pe = pearsonr(sf[ok], d_sub[ok])[0]
    say(f"  spearman(sf, depth) = {rs_sp:+.8f}")
    say(f"  pearson (sf, depth) = {rs_pe:+.8f}")
    lx, ly = np.log(d_sub[ok]), np.log(sf[ok])
    b, a = np.polyfit(lx, ly, 1)
    resid = ly - (a + b * lx)
    say(f"  log-log slope       = {b:+.6f}   (1.0 => strictly proportional)")
    say(f"  residual SD         = {resid.std():.3e}")
    say("")
    if rs_sp > 0.9999999 and resid.std() < 1e-8:
        say("  >>> sf is a deterministic function of depth")
    else:
        say("  >>> sf departs from depth")


# == T6: zero concordance =====================================================
hdr("T6 - zero concordance (catches imputation / ambient rewriting)")
xz = (X != 0); cz = (C != 0)
only_x = int((xz > cz).sum())
only_c = int((cz > xz).sum())
say(f"  nonzero in X but zero in counts : {only_x}")
say(f"  nonzero in counts but zero in X : {only_c}")
if only_x == 0 and only_c == 0:
    say("  >>> zeros concordant")
else:
    say("  >>> zeros discordant")


# == T7: exclude_highly_expressed variant =====================================
hdr("T7 - exclude_highly_expressed variant")
EXCLUDE_HE = False
if not LIBSIZE:
    frac = C.multiply(sp.csr_matrix(1.0 / np.maximum(depth, 1)).T)
    hi = np.asarray(frac.max(axis=0).todense()).ravel() > 0.05
    say(f"  genes with max fraction > 0.05 : {int(hi.sum())}")
    keep = np.where(~hi)[0]
    rs_k = np.asarray(X[:, keep].expm1().sum(axis=1)).ravel()
    cvk = rs_k.std() / rs_k.mean()
    say(f"  row-sum over kept genes: median={np.median(rs_k):,.3f}  CV={cvk:.3e}")
    if cvk < 1e-6:
        EXCLUDE_HE = True
        T_HE = float(np.median(rs_k))
        say("  >>> constant after exclusion")
        say(f"      recovered target_sum = {T_HE:,.2f}")
    else:
        say("  >>> not constant after exclusion")
else:
    say("  (skipped)")


# == T8: source of the denominator ============================================
hdr("T8 - denominator provenance: was the total taken over THESE genes?")
if BASE in ("e", "2") and ok.sum() > 50:
    # implied denominator = sum over genes of expm1(X) * sf  == what was divided
    num = np.asarray(X[sub].expm1().sum(axis=1)).ravel() if BASE == "e" \
        else np.asarray(X2[sub].expm1().sum(axis=1)).ravel()
    implied_total = num * sf
    ratio = implied_total[ok] / d_sub[ok]
    say(f"  implied_total / observed_depth :  min={ratio.min():.6f}  "
        f"median={np.median(ratio):.6f}  max={ratio.max():.6f}")
    say(f"  CV of that ratio               :  {ratio.std()/ratio.mean():.3e}")
    if abs(np.median(ratio) - 1.0) < 1e-6 and ratio.std()/ratio.mean() < 1e-6:
        say("  >>> denominator = row total over these genes")
    else:
        say("  >>> denominator differs from the row total over these genes")


# == T9: worked example =======================================================
hdr("T9 - worked example (manual audit)")
say(f"  {'cell':>6} {'gene':>8} {'count':>8} {'depth':>10} {'X':>10} "
    f"{'pred':>10} {'diff':>10}")
shown = 0
for i in range(Xd.shape[0]):
    nzj = np.where((Cd[i] > 0) & (Xd[i] > 0))[0]
    if nzj.size == 0:
        continue
    for j in nzj[:2]:
        if LIBSIZE and np.isfinite(T_hat):
            pred = (np.log1p(Cd[i, j] / d_sub[i] * T_hat) if BASE == "e"
                    else np.log2(1 + Cd[i, j] / d_sub[i] * T_hat))
        else:
            pred = np.nan
        say(f"  {pick[sub[i]]:>6} {j:>8} {Cd[i,j]:>8.0f} {d_sub[i]:>10.0f} "
            f"{Xd[i,j]:>10.6f} {pred:>10.6f} {Xd[i,j]-pred:>10.2e}")
        shown += 1
    if shown >= 10:
        break


# == Verdict ==================================================================
hdr("VERDICT")
if LIBSIZE and BASE == "e":
    say(f"  X = log1p( counts / total_counts * {T_hat:,.0f} )")
    say("  Method : normalize_total + log1p (natural log)")
    if abs(T_hat - 1e4) / 1e4 < 1e-6:
        say("  >>> target_sum = 1e4 (CP10K)")
    elif abs(T_hat - 1e5) / 1e5 < 1e-6:
        say("  >>> target_sum = 1e5 (CP100K)")
    else:
        say(f"  >>> target_sum = {T_hat:,.0f}, neither CP10K nor CP100K")
elif LIBSIZE and BASE == "2":
    say(f"  X = log2( 1 + counts / total_counts * {T_hat:,.0f} )")
    say("  Library-size scaling with a base-2 log")
elif EXCLUDE_HE:
    say(f"  X = log_{BASE}( 1 + counts / S_nonHE * {T_HE:,.0f} )")
    say("  where S_nonHE excludes genes exceeding 5% of any library.")
    say("  Method : normalize_total(exclude_highly_expressed=True) + log1p")
elif not LIBSIZE:
    say("  X is not library-size scaled (per-cell size factors)")
else:
    say("  Inconclusive.")

with open(os.path.join(OUT_DIR, "S1a_norm_forensics.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(_rep))
say("")
say(f"  saved: {OUT_DIR}/S1a_norm_forensics.txt")
