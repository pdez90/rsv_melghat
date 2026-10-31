# 99_freeze.R -- thin wrapper. The freeze tool is scripts/99_freeze.py so that it
# runs in any environment that can see the project, including ones without R.
# That portability is what caught the earlier fault: the manifest had been
# generated against one machine's files and shipped to another's.
args <- commandArgs(trailingOnly = TRUE)
py <- file.path("scripts", "99_freeze.py")
if (!file.exists(py)) stop("scripts/99_freeze.py not found")
quit(status = system2("python3", c(py, args)))
