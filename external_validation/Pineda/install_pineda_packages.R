# install_pineda_packages.R
# Run by the Pineda Containerfile. Installs the R packages that R02 (scDblFinder) needs beyond
# renv.lock of the SEA-AD container, at the versions of the R02 run on Windows (its session info,
# 2026-10-27). The list is pineda_r_packages.tsv (package, version, repository).
#   1. install.packages() installs the listed packages and any missing dependency from Bioconductor
#      3.23 and the Posit Package Manager CRAN snapshot. Packages already installed (those of
#      renv.lock) are left as they are.
#   2. A listed package installed at another version is installed again from the exact source file:
#      the current file of the repository, or the file in its Archive folder.
# check_lock_pineda.R then stops the build if a version differs.
# Usage: Rscript install_pineda_packages.R <pineda_r_packages.tsv> <Bioconductor version> <CRAN repository>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
    stop("Usage: Rscript install_pineda_packages.R <pineda_r_packages.tsv> <Bioconductor version> <CRAN repository>")
want <- utils::read.delim(args[1], stringsAsFactors = FALSE, colClasses = "character")
bioc <- sprintf("https://bioconductor.org/packages/%s/%s", args[2], c("bioc", "data/annotation", "data/experiment"))
cran <- args[3]
repos <- c(BioCsoft = bioc[1], BioCann = bioc[2], BioCexp = bioc[3], CRAN = cran)
lib <- .libPaths()[1]
options(Ncpus = max(1L, parallel::detectCores() - 1L), timeout = 600)
installed_version <- function(p) {
    ip <- installed.packages(lib.loc = .libPaths())
    if (p %in% rownames(ip)) ip[p, "Version"] else NA_character_
}

cat("Library:", lib, "\nRepositories:", paste(repos, collapse = " | "), "\n")
todo <- want$package[is.na(vapply(want$package, installed_version, character(1)))]
cat(sprintf("Step 1: %d of %d listed packages to install: %s\n", length(todo), nrow(want), paste(todo, collapse = ", ")))
if (length(todo))
    install.packages(todo, lib = lib, repos = repos, dependencies = c("Depends", "Imports", "LinkingTo"))

# Step 2: exact versions, in the order of the list (dependencies first).
src_urls <- function(p, v, repo) {
    f <- sprintf("%s_%s.tar.gz", p, v)
    if (repo == "Bioconductor")
        c(sprintf("%s/src/contrib/%s", bioc[1], f), sprintf("%s/src/contrib/Archive/%s/%s", bioc[1], p, f))
    else
        c(sprintf("https://cloud.r-project.org/src/contrib/%s", f),
          sprintf("https://cloud.r-project.org/src/contrib/Archive/%s/%s", p, f))
}
for (i in seq_len(nrow(want))) {
    p <- want$package[i]; v <- want$version[i]
    have <- installed_version(p)
    if (identical(have, v)) next
    cat(sprintf("Step 2: %s is %s, wanted %s; installing the exact source file\n", p, have, v))
    tf <- file.path(tempdir(), sprintf("%s_%s.tar.gz", p, v))
    ok <- FALSE
    for (u in src_urls(p, v, want$repository[i])) {
        r <- tryCatch(utils::download.file(u, tf, quiet = TRUE, mode = "wb"), error = function(e) 1L)
        if (identical(r, 0L) && file.size(tf) > 0) { cat("  from", u, "\n"); ok <- TRUE; break }
    }
    if (!ok) stop("No source file found for ", p, " ", v)
    install.packages(tf, lib = lib, repos = NULL, type = "source")
    unlink(tf)
}
cat("install_pineda_packages.R: done\n")
