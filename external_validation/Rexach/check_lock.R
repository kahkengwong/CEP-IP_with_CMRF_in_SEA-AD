# check_lock.R
# Run at the end of the package step of the Rexach Containerfile. Stops the build when:
#   - any package of renv.lock (the SEA-AD container) is missing or has another version;
#   - cepip is not 0.1.0 installed from the expected commit;
#   - brglm2 is not the expected version.
# Usage: Rscript check_lock.R <renv.lock> <cepip commit> <brglm2 version>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript check_lock.R <renv.lock> <cepip commit> <brglm2 version>")
lock <- jsonlite::fromJSON(args[1], simplifyVector = FALSE)$Packages
ip <- installed.packages()[, "Version"]

want <- vapply(lock, function(p) p$Version, character(1))
have <- ip[names(want)]
bad <- names(want)[is.na(have) | have != want]
cat(sprintf("renv.lock packages: %d | same version installed: %d\n", length(want), length(want) - length(bad)))
if (length(bad)) {
    cat(sprintf("  %s: lock %s, installed %s\n", bad, want[bad], have[bad]), sep = "")
    stop("A package of renv.lock changed version; the image would not match the SEA-AD container.")
}

d <- packageDescription("cepip")
sha_ok <- identical(d$Version, "0.1.0") && !is.null(d$RemoteSha) && identical(d$RemoteSha, args[2])
cat(sprintf("cepip %s | RemoteSha %s\n", d$Version, if (is.null(d$RemoteSha)) "not recorded" else d$RemoteSha))
if (!sha_ok) stop("cepip is not 0.1.0 from commit ", args[2])

b <- as.character(packageVersion("brglm2"))
cat(sprintf("brglm2 %s\n", b))
if (!identical(b, args[3])) stop("brglm2 is ", b, ", expected ", args[3])

for (pk in c("cepip", "brglm2")) suppressPackageStartupMessages(library(pk, character.only = TRUE))
cat("cepip and brglm2 load\n")
added <- setdiff(names(ip), c(names(want), rownames(installed.packages(priority = "base"))))
cat(sprintf("Added beyond renv.lock: %s\n", paste(sprintf("%s %s", added, ip[added]), collapse = ", ")))
cat("check_lock.R: PASS\n")
