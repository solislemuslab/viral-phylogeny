# Re-root every OB45 result tree and reference tree on the outgroup "MYROD"
# Originals are kept; rooted copies go to a parallel "<dir>_MYROD" folder with the same file names
# Rooting only changes drawing / rooted comparisons -- unrooted RF distances are unaffected

library(ape)

outgroup <- "MYROD"
in_dirs <- c(Sys.glob("../OB45_*_result_trees"), "../OB45_reference_trees")

summary <- list()
for (d in in_dirs) {
  out_dir <- paste0(d, "_MYROD")
  dir.create(out_dir, showWarnings = FALSE)
  for (f in list.files(d, pattern = "\\.(nwk|treefile|tre)$")) {
    t <- read.tree(file.path(d, f))
    if (!outgroup %in% t$tip.label) stop(outgroup, " not found in ", file.path(d, f))
    r <- root(t, outgroup = outgroup, resolve.root = TRUE)
    root_kids <- r$edge[r$edge[, 1] == Ntip(r) + 1, 2]
    ok <- is.rooted(r) && length(root_kids) == 2 && outgroup %in% r$tip.label[root_kids]
    write.tree(r, file.path(out_dir, f))
    summary[[length(summary) + 1]] <- data.frame(dir = basename(d), file = f, ntip = Ntip(r), rooted_ok = ok)
  }
}

summary <- do.call(rbind, summary)
print(summary, row.names = FALSE)
if (!all(summary$rooted_ok)) stop("some trees were not rooted on ", outgroup)
