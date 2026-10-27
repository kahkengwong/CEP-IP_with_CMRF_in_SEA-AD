# check_lock_pineda.R
# Run at the end of the package step of the Pineda Containerfile. Stops the build when:
#   - any package of renv.lock (the SEA-AD container) is missing or has another version;
#   - cepip is not 0.1.0 installed from the expected commit;
#   - brglm2 is not the expected version;
#   - any package of pineda_r_packages.tsv (the R02 packages beyond renv.lock) is missing or has
#     another version than in the R02 run on Windows.
# Usage: Rscript check_lock_pineda.R <renv.lock> <cepip commit> <brglm2 version> <pineda_r_packages.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L)
    stop("Usage: Rscript check_lock_pineda.R <renv.lock> <cepip commit> <brglm2 version> <pineda_r_packages.tsv>")
lock <- jsonlite::fromJSON(args[1], simplifyVector = FALSE)$Packages
ip <- installed.packages()[, "Version"]

check_set <- function(want, label) {
    have <- ip[names(want)]
    bad <- names(want)[is.na(have) | have != want]
    cat(sprintf("%s: %d | same version installed: %d\n", label, length(want), length(want) - length(bad)))
    if (length(bad)) {
        cat(sprintf("  %s: wanted %s, installed %s\n", bad, want[bad], have[bad]), sep = "")
        stop("A package of ", label, " has another version; the image would not match.")
    }
}
check_set(vapply(lock, function(p) p$Version, character(1)), "renv.lock packages")
extra <- utils::read.delim(args[4], stringsAsFactors = FALSE, colClasses = "character")
check_set(stats::setNames(extra$version, extra$package), "R02 packages (pineda_r_packages.tsv)")

d <- packageDescription("cepip")
sha_ok <- identical(d$Version, "0.1.0") && !is.null(d$RemoteSha) && identical(d$RemoteSha, args[2])
cat(sprintf("cepip %s | RemoteSha %s\n", d$Version, if (is.null(d$RemoteSha)) "not recorded" else d$RemoteSha))
if (!sha_ok) stop("cepip is not 0.1.0 from commit ", args[2])

b <- as.character(packageVersion("brglm2"))
cat(sprintf("brglm2 %s\n", b))
if (!identical(b, args[3])) stop("brglm2 is ", b, ", expected ", args[3])

for (pk in c("cepip", "brglm2", "DropletUtils", "scDblFinder"))
    suppressPackageStartupMessages(library(pk, character.only = TRUE))
cat("cepip, brglm2, DropletUtils, and scDblFinder load\n")
added <- setdiff(names(ip), c(names(lock), extra$package, rownames(installed.packages(priority = "base"))))
cat(sprintf("Added beyond renv.lock and the R02 list: %s\n",
            if (length(added)) paste(sprintf("%s %s", added, ip[added]), collapse = ", ") else "none"))
cat("check_lock_pineda.R: PASS\n")
