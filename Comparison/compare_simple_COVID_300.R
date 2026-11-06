library(ape)
library(phangorn)
library(phytools)

# COVID 300 (+4 animal-CoV outgroups), whole-genome nucleotide alignment (304 taxa x 30032 sites).
# Adapted from compare_simple_OB_603.R: no reference trees, tips coloured by continent, rooted on RaTG13,
# near-zero branches collapsed, and quartets sampled (all 4-taxon subsets of 304 tips is ~3.5e8).
tree_dir <- "../covid_result_trees"
files <- c("IQ-TREE"  = "covid_300_IQTree.nwk",
           "CMAPLE"   = "covid_300_CMAPLE.treefile",
           "NeuralNJ" = "covid_300_NeuralNJ.nwk")
files <- files[file.exists(file.path(tree_dir, files))]
trees <- lapply(file.path(tree_dir, files), read.tree)
names(trees) <- names(files)

# NeuralNJ needs taxa named T1..T304; map them back to the GISAID names
nnj_map <- "../NeuralNJ/examples/covid_300_aligned/taxon_map.tsv"
if ("NeuralNJ" %in% names(trees)) {
  m <- read.delim(nnj_map, stringsAsFactors = FALSE)
  trees$NeuralNJ$tip.label <- setNames(m$taxon, m$neuralnj_id)[trees$NeuralNJ$tip.label]
  trees$NeuralNJ$edge.length <- NULL # Argmax (README command, no --branch_optimize) writes a placeholder 0.12345 on every branch
}
# GISAID names are long ("hCoV-19/England/.../2020|EPI_ISL_423259|2020-03-24"); keep the EPI_ISL accession
for (nm in names(trees)) trees[[nm]]$tip.label <- sub("^[^|]*\\|(EPI_ISL_[0-9]+)\\|.*$", "\\1", trees[[nm]]$tip.label)

common <- Reduce(intersect, lapply(trees, `[[`, "tip.label")) # Extract the common tips across all trees
has_bl <- sapply(trees, function(t) !is.null(t$edge.length) && !anyNA(t$edge.length)) # trees with real branch lengths

# Root on the bat coronavirus RaTG13, the most distant outgroup.
# One substitution is 1/30032 = 3.3e-5 subs/site, and IQ-TREE reports unsupported branches as 1e-6,
# so collapse branches shorter than 1e-5 into polytomies (identical sequences cannot be resolved)
root_tip <- "EPI_ISL_402131"
trees <- lapply(trees, function(t) {
  if (is.null(t$edge.length) || anyNA(t$edge.length)) t <- compute.brlen(t, 1)
  ladderize(root(di2multi(keep.tip(t, common), tol = 1e-5), root_tip, resolve.root = TRUE))
})
class(trees) <- "multiPhylo"
ref <- trees[["IQ-TREE"]] # Use IQ-TREE as reference for cophylo and consensus trees

cols <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")[seq_along(trees)], names(trees))
out <- "figures_R"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(out, "compare_simple_COVID_300.pdf"), width = 12, height = 10)

## 1. circular trees, same colour for each continent in every tree
# Continent of each sampled sequence (select_covid_300_stratified.py); the 4 animal CoVs are the outgroup
meta <- read.csv("../COVID_data_fasta/covid_300_stratified_metadata.csv", stringsAsFactors = FALSE)
outgroups <- c("EPI_ISL_402131" = "bat RaTG13", "EPI_ISL_410721" = "pangolin GD/1",
               "EPI_ISL_410544" = "pangolin GD/P2S", "EPI_ISL_410543" = "pangolin GX/P3B")
clade_of <- c(setNames(meta$continent, meta$accession), setNames(rep("Outgroup", length(outgroups)), names(outgroups)))
clade_of <- setNames(ifelse(common %in% names(clade_of), clade_of[common], "Other"), common)
clade_cols <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "red", "grey40"),
                       c("Europe", "North America", "Oceania", "Asia", "Africa", "South America", "Outgroup", "Other"))

# An edge takes a continent's colour if all tips below it are from that continent; otherwise grey
clade_edge_cols <- function(t) {
  node_grp <- sapply(Descendants(t, seq_len(Ntip(t) + t$Nnode), "tips"), function(d) {
    g <- unique(clade_of[t$tip.label[d]])
    if (length(g) == 1) g else NA
  })
  ec <- clade_cols[node_grp[t$edge[, 2]]]
  ec[is.na(ec)] <- "grey70"
  ec
}

# The outgroup branches are ~1000x longer than any branch inside SARS-CoV-2, so draw topology only (see 4b for lengths)
par(mfrow = n2mfrow(length(trees) + 1), mar = c(1, 1, 2, 1), fg = "grey70", col = "black") # ape draws mixed-clade polytomy arcs in fg (black by default)
for (nm in names(trees)) {
  t <- trees[[nm]]
  plot(t, type = "fan", use.edge.length = FALSE, show.tip.label = FALSE, edge.color = clade_edge_cols(t), edge.width = 1.2)
  tiplabels(pch = 19, col = clade_cols[clade_of[t$tip.label]], cex = 0.5)
  title(nm)
}
plot.new()
legend("center", legend = names(clade_cols), col = clade_cols, pch = 19, lwd = 2, bty = "n", cex = 0.9, title = "continent")

## 2. co-phylogeny: reference vs each other tree
par(mfrow = c(1, 1), fg = "black")
for (nm in setdiff(names(trees), "IQ-TREE")) {
  plot(cophylo(ref, trees[[nm]], rotate = TRUE, print = FALSE), fsize = 0.15, lwd = 0.5)
  title(paste("IQ-TREE vs", nm))
}

## 3. topological distances between methods (unrooted)
u <- lapply(trees, unroot)
class(u) <- "multiPhylo"
rf <- as.matrix(RF.dist(u, normalize = TRUE))
print(round(rf, 3))
rf_breaks <- seq(0, 1, by = 0.1) # fixed 0-1 scale so colours are comparable across datasets
rf_cols <- hcl.colors(length(rf_breaks) - 1, "YlOrRd", rev = TRUE) # light = similar, dark = different
heatmap(rf, symm = TRUE, margins = c(11, 11), main = "normalised RF", col = rf_cols, breaks = rf_breaks)
# Colour key: the RF range covered by each colour
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
legend("topleft", title = "normalised RF", legend = sprintf("%.1f - %.1f", head(rf_breaks, -1), rf_breaks[-1]),
       fill = rf_cols, bty = "n", cex = 0.9)

## PCoA of the methods on RF distance (needs at least 3 trees for 2 axes)
if (length(trees) >= 3) {
  pc <- cmdscale(as.dist(rf), k = 2)
  plot(pc, pch = 19, col = cols[rownames(pc)], cex = 2, xlim = range(pc[, 1]) * 1.3, ylim = range(pc[, 2]) * 1.2, xlab = "PCoA 1", ylab = "PCoA 2", main = "methods by RF distance")
  text(pc, labels = rownames(pc), pos = 3)
}

## 3b. quartet consistency
# For every 4 taxa {a,b,c,d}, a tree either resolves them (ab|cd, ac|bd or ad|bc) or leaves a star (polytomy).
# Four-point condition on unit branch lengths: the pairing with the strictly smallest sum of path lengths is
# the quartet's topology; if there is no unique minimum the quartet sits on a polytomy and counts as unresolved.
# 304 tips give choose(304, 4) ~ 3.5e8 quartets, so use a random sample of distinct 4-taxon subsets
taxa <- sort(common)
n_quartets <- 200000
set.seed(42)
qs <- unique(t(apply(replicate(n_quartets, sample(length(taxa), 4)), 2, sort)))
a <- qs[, 1]; b <- qs[, 2]; c4 <- qs[, 3]; d <- qs[, 4]
quartet_topo <- function(t) {
  D <- cophenetic(compute.brlen(t, 1))[taxa, taxa]
  s1 <- D[cbind(a, b)] + D[cbind(c4, d)] # ab|cd
  s2 <- D[cbind(a, c4)] + D[cbind(b, d)] # ac|bd
  s3 <- D[cbind(a, d)] + D[cbind(b, c4)] # ad|bc
  ifelse(s1 < s2 & s1 < s3, 1L, ifelse(s2 < s1 & s2 < s3, 2L, ifelse(s3 < s1 & s3 < s2, 3L, 0L))) # 0 = unresolved
}
qt <- sapply(u, quartet_topo) # rows = quartets, columns = trees
cat("\nQuartets sampled:", nrow(qt), "of", choose(length(taxa), 4), "\nShare of quartets resolved in each tree (polytomies lower this):\n")
print(round(colMeans(qt > 0), 3))

# Pairwise agreement: of the quartets resolved in both trees, the share with the same topology
nt <- ncol(qt)
qagree <- matrix(1, nt, nt, dimnames = list(colnames(qt), colnames(qt)))
for (i in seq_len(nt)) for (j in seq_len(nt)) {
  both <- qt[, i] > 0 & qt[, j] > 0
  qagree[i, j] <- mean(qt[both, i] == qt[both, j])
}
qdist <- 1 - qagree
cat("\nQuartet distance (1 - agreement on quartets resolved in both trees):\n")
print(round(qdist, 3))
# Quartet distances are much smaller than RF, so use a finer fixed scale (0-0.2 in steps of 0.02, everything above in one bin)
q_breaks <- c(seq(0, 0.2, by = 0.02), 1)
q_cols <- hcl.colors(length(q_breaks) - 1, "YlOrRd", rev = TRUE) # light = similar, dark = different
heatmap(qdist, symm = TRUE, margins = c(11, 11), main = "quartet distance", col = q_cols, breaks = q_breaks)
# Colour key: the quartet-distance range covered by each colour
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
legend("topleft", title = "quartet distance", legend = sprintf("%.2f - %.2f", head(q_breaks, -1), q_breaks[-1]),
       fill = q_cols, bty = "n", cex = 0.9)

# Quartets that every tree resolves identically
all_same <- qt[, 1] > 0 & apply(qt == qt[, 1], 1, all)
cat("\nSampled quartets resolved identically by all trees:", sum(all_same), sprintf("(%.1f%%)\n", 100 * mean(all_same)))
# For each quartet: the largest number of trees sharing one resolved topology
n_agree <- apply(qt, 1, function(x) { x <- x[x > 0]; if (length(x)) max(tabulate(x, 3)) else 0L })
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
barplot(table(factor(n_agree, levels = 0:nt)) / length(n_agree), col = "#7570b3",
        xlab = "number of trees agreeing on the quartet topology", ylab = "share of sampled quartets",
        main = sprintf("quartet consistency (%d trees; %.1f%% shared by all)", nt, 100 * mean(all_same)))
common_q <- which(all_same) # topology column: "12|34" means taxon1 taxon2 | taxon3 taxon4

write.csv(data.frame(taxon1 = taxa[a[common_q]], taxon2 = taxa[b[common_q]], taxon3 = taxa[c4[common_q]], taxon4 = taxa[d[common_q]],
                     topology = c("12|34", "13|24", "14|23")[qt[common_q, 1]]),
          file.path(out, "quartets_common_COVID_300.csv"), row.names = FALSE)

## 4. densitree, one colour per method
g <- lapply(trees, compute.brlen, method = "Grafen"); class(g) <- "multiPhylo"
densiTree(g, consensus = ref, col = cols, alpha = 0.6, scaleX = TRUE, cex = 0.15, type = "cladogram")
legend("topleft", legend = names(cols), col = cols, lwd = 2, bty = "n", cex = 0.8)
title("densitree: topology only (branch lengths ignored)")

## 4b. real branch lengths, trees that have them, on one shared axis
# Draw each tree as a phylogram with the same x-axis (substitutions/site from the root).
# The outgroups sit far from the SARS-CoV-2 sequences, so also draw the trees without them
# to show the branch lengths inside SARS-CoV-2.
bl_trees <- trees[has_bl]
for (drop_out in c(FALSE, TRUE)) {
  bt <- if (drop_out) lapply(bl_trees, function(t) ladderize(drop.tip(t, names(outgroups)))) else bl_trees
  xmax <- max(sapply(bt, function(t) max(node.depth.edgelength(t))))
  par(mfrow = c(1, length(bt)), mar = c(4, 1, 3, 1), fg = "grey70", col = "black")
  for (nm in names(bt)) {
    t <- bt[[nm]]
    plot(t, x.lim = c(0, xmax * 1.12), cex = 0.15, edge.color = clade_edge_cols(t), edge.width = 0.8)
    axis(1, cex.axis = 0.8)
    mtext("substitutions/site from root", side = 1, line = 2.5, cex = 0.7)
    title(if (drop_out) paste(nm, "(outgroups removed)") else nm)
  }
}

## 5. which clades of the reference tree (IQ-TREE) are present in each other tree
# All trees are rooted on RaTG13, so a rooted clade here is the same as an unrooted split
others <- setdiff(names(trees), "IQ-TREE")
present <- sapply(others, function(nm) prop.clades(ref, trees[[nm]], rooted = TRUE))
present <- matrix(present, ncol = length(others), dimnames = list(NULL, others))
present[is.na(present)] <- 0
n_present <- rowSums(present)
cat("\nShare of IQ-TREE clades recovered by each tree:\n")
print(round(colMeans(present), 3))
count_cols <- hcl.colors(length(others) + 1, "YlGnBu", rev = TRUE)
par(mfrow = c(1, 1), mar = c(1, 1, 3, 1), fg = "grey70", col = "black")
plot(ref, cex = 0.15, no.margin = FALSE, edge.color = clade_edge_cols(ref))
nodelabels(n_present, bg = count_cols[n_present + 1], cex = 0.25, frame = "circle")
title(sprintf("IQ-TREE clades: number of the other %d trees that contain each clade", length(others)))
legend("bottomleft", legend = 0:length(others), fill = count_cols, title = "trees with clade", bty = "n", cex = 0.7)

## 6. consensus trees
par(mfrow = c(1, 2), fg = "black")
for (p in c(1, 0.5)) {
  ct <- consensus(trees, p = p, rooted = FALSE)
  plot(ct, cex = 0.15)
  title(if (p == 1) "strict" else "majority rule")
}

dev.off()
