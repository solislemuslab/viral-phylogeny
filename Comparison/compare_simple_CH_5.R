library(ape)
library(phangorn)
library(phytools)

outgroup <- "DANRE" # root of CH-astral.tre; used for every plot (RF/quartets stay unrooted)
tree_dir <- "../CH45_5_result_trees"
files <- c("IQ-TREE"         = "CH45_supermatrix_5_IQTree.treefile",
           "Phyloformer"     = "CH45_supermatrix_5_Phyloformer.nwk",
           "CMAPLE"          = "CH45_supermatrix_5_CMAPLE.treefile",
           "MSA-T attention" = "CH45_supermatrix_5_0_4_MSATransformer_AttentionTree.nwk",
           "MSA-T embedding" = "CH45_5_2_embedding_tree.nwk")
ref_dir <- "../CH45_reference_trees"
ref_files <- c("CH-astral"      = "CH-astral.tre",
               "CH-ncbi-astral" = "CH-ncbi-astral.tre",
               "CH-ncbi"        = "CH-ncbi.tre")
trees <- c(lapply(file.path(tree_dir, files), read.tree),
           lapply(file.path(ref_dir, ref_files), read.tree))
names(trees) <- c(names(files), names(ref_files))
common <- Reduce(intersect, lapply(trees, `[[`, "tip.label")) # Extract the common tips across all trees
has_bl <- sapply(trees, function(t) !is.null(t$edge.length) && !anyNA(t$edge.length)) # trees with real branch lengths

# Root every tree on the outgroup after pruning so all trees share the same root.
# Reference trees are topology-only, so give them unit branch lengths for plotting
trees <- lapply(trees, function(t) {
  if (is.null(t$edge.length) || anyNA(t$edge.length)) t <- compute.brlen(t, 1)
  ladderize(root(di2multi(keep.tip(t, common), tol = 1e-8), outgroup, resolve.root = TRUE))
})
class(trees) <- "multiPhylo"
ref <- trees[["IQ-TREE"]] # Use IQ-TREE as reference for cophylo and consensus trees

cols <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")[seq_along(trees)], names(trees))
out <- "figures_R"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(out, "compare_simple_CH_5.pdf"), width = 10, height = 8)

## 1. circular trees, same colour for each clade in every tree
# Clades follow the NCBI taxonomy (as in CH-ncbi.tre)
clades <- list(
  "Teleostei"           = c("MOLML", "TAKRU", "AMPCI", "ORYME", "XIPMA", "CYPVA", "KRYMA", "AMPPE", "ANATE", "SERDU",
                            "CYNSE", "ESOLU", "SALSA", "PYGNA", "ASTMX", "ICTPU"),
  "Shark + gar"         = c("CHIPU", "LEPOC"),
  "Amphibia"            = "XENLA",
  "Sauropsida"          = c("ALLSI", "PATFA", "ANAPP", "CHICK", "TAEGU"),
  "Laurasiatheria"      = c("BALAS", "PHYMC", "LIPVE", "TURTR", "DELLE", "CAPHI", "PIG", "PTEAL", "MYOLU", "VULVU",
                            "ODORO", "ENHLU", "URSAR", "LEPWE", "FELCA", "ERIEU", "HORSE"),
  "Euarchontoglires"    = c("OTOGA", "PROCO", "TARSY", "AOTNA", "CALJA", "CERAT", "HUMAN", "NOMLE", "RABIT", "DIPOR",
                            "MOUSE", "MESAU", "ICTTR", "HETGA", "CAVPO"),
  "Afrotheria + opossum" = c("TRIMA", "LOXAF", "MONDO"),
  "Outgroup (DANRE)"    = "DANRE")
clade_of <- setNames(rep(names(clades), lengths(clades)), unlist(clades))
clade_of <- setNames(ifelse(common %in% names(clade_of), clade_of[common], "Other"), common)
clade_cols <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "red", "grey40"),
                       c(names(clades), "Other"))

# An edge takes a clade's colour if all tips below it are in that clade; otherwise grey
clade_edge_cols <- function(t) {
  node_grp <- sapply(Descendants(t, seq_len(Ntip(t) + t$Nnode), "tips"), function(d) {
    g <- unique(clade_of[t$tip.label[d]])
    if (length(g) == 1) g else NA
  })
  ec <- clade_cols[node_grp[t$edge[, 2]]]
  ec[is.na(ec)] <- "grey70"
  ec
}

par(mfrow = n2mfrow(length(trees) + 1), mar = c(1, 1, 2, 1), fg = "grey70", col = "black") # ape draws mixed-clade polytomy arcs in fg (black by default)
for (nm in names(trees)) {
  t <- trees[[nm]]
  plot(t, type = "fan", show.tip.label = FALSE, edge.color = clade_edge_cols(t), edge.width = 1.5)
  tiplabels(pch = 19, col = clade_cols[clade_of[t$tip.label]], cex = 0.8)
  title(nm)
}
plot.new()
legend("center", legend = names(clade_cols), col = clade_cols, pch = 19, lwd = 2, bty = "n", cex = 0.9, title = "clade")

## 2. co-phylogeny: reference vs each other tree
par(mfrow = c(1, 1), fg = "black")
for (nm in setdiff(names(trees), "IQ-TREE")) {
  plot(cophylo(ref, trees[[nm]], rotate = TRUE, print = FALSE), fsize = 0.3)
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

## PCoA of the methods on RF distance
pc <- cmdscale(as.dist(rf), k = 2)
plot(pc, pch = 19, col = cols[rownames(pc)], cex = 2, xlim = range(pc[, 1]) * 1.3, ylim = range(pc[, 2]) * 1.2, xlab = "PCoA 1", ylab = "PCoA 2", main = "methods by RF distance")
text(pc, labels = rownames(pc), pos = 3)

## 3b. quartet consistency
# For every 4 taxa {a,b,c,d}, a tree either resolves them (ab|cd, ac|bd or ad|bc) or leaves a star (polytomy).
# Four-point condition on unit branch lengths: the pairing with the strictly smallest sum of path lengths is
# the quartet's topology; if there is no unique minimum the quartet sits on a polytomy and counts as unresolved.
taxa <- sort(common)
qs <- combn(length(taxa), 4) # all 4-taxon subsets, one per column
a <- qs[1, ]; b <- qs[2, ]; c4 <- qs[3, ]; d <- qs[4, ]
quartet_topo <- function(t) {
  D <- cophenetic(compute.brlen(t, 1))[taxa, taxa]
  s1 <- D[cbind(a, b)] + D[cbind(c4, d)] # ab|cd
  s2 <- D[cbind(a, c4)] + D[cbind(b, d)] # ac|bd
  s3 <- D[cbind(a, d)] + D[cbind(b, c4)] # ad|bc
  ifelse(s1 < s2 & s1 < s3, 1L, ifelse(s2 < s1 & s2 < s3, 2L, ifelse(s3 < s1 & s3 < s2, 3L, 0L))) # 0 = unresolved
}
qt <- sapply(u, quartet_topo) # rows = quartets, columns = trees
cat("\nQuartets:", nrow(qt), "\nShare of quartets resolved in each tree (polytomies lower this):\n")
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
cat("\nQuartets resolved identically by all trees:", sum(all_same), sprintf("(%.1f%%)\n", 100 * mean(all_same)))
# For each quartet: the largest number of trees sharing one resolved topology
n_agree <- apply(qt, 1, function(x) { x <- x[x > 0]; if (length(x)) max(tabulate(x, 3)) else 0L })
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
barplot(table(factor(n_agree, levels = 0:nt)) / length(n_agree), col = "#7570b3",
        xlab = "number of trees agreeing on the quartet topology", ylab = "share of quartets",
        main = sprintf("quartet consistency (%d trees; %.1f%% shared by all)", nt, 100 * mean(all_same)))
common_q <- which(all_same) # topology column: "12|34" means taxon1 taxon2 | taxon3 taxon4

write.csv(data.frame(taxon1 = taxa[a[common_q]], taxon2 = taxa[b[common_q]], taxon3 = taxa[c4[common_q]], taxon4 = taxa[d[common_q]],
                     topology = c("12|34", "13|24", "14|23")[qt[common_q, 1]]),
          file.path(out, "quartets_common_CH_5.csv"), row.names = FALSE)

## 4. densitree, one colour per method
g <- lapply(trees, compute.brlen, method = "Grafen"); class(g) <- "multiPhylo"
densiTree(g, consensus = ref, col = cols, alpha = 0.6, scaleX = TRUE, cex = 0.3, type = "cladogram")
legend("topleft", legend = names(cols), col = cols, lwd = 2, bty = "n", cex = 0.8)
title("densitree: topology only (branch lengths ignored)")

## 4b. real branch lengths, trees that have them, on one shared axis
# densiTree only draws ultrametric trees correctly, and these trees are far from clock-like
# (root-to-tip depth varies several-fold), so making them ultrametric (chronos / nnls) distorts the lengths badly.
# Instead draw each tree as a phylogram with the same x-axis (substitutions/site from the root).
bl_trees <- trees[has_bl]
xmax <- max(sapply(bl_trees, function(t) max(node.depth.edgelength(t))))
par(mfrow = c(1, length(bl_trees)), mar = c(4, 1, 3, 1), fg = "grey70", col = "black")
for (nm in names(bl_trees)) {
  t <- bl_trees[[nm]]
  plot(t, x.lim = c(0, xmax * 1.12), cex = 0.35, edge.color = clade_edge_cols(t), edge.width = 1.2)
  axis(1, cex.axis = 0.8)
  mtext("substitutions/site from root", side = 1, line = 2.5, cex = 0.7)
  title(nm)
}

## 5. which clades of the reference tree (IQ-TREE) are present in each other tree
# All trees are rooted on the outgroup, so a rooted clade here is the same as an unrooted split
others <- setdiff(names(trees), "IQ-TREE")
present <- sapply(others, function(nm) prop.clades(ref, trees[[nm]], rooted = TRUE))
present[is.na(present)] <- 0
n_present <- rowSums(present)
cat("\nShare of IQ-TREE clades recovered by each tree:\n")
print(round(colMeans(present), 3))
count_cols <- hcl.colors(length(others) + 1, "YlGnBu", rev = TRUE)
par(mfrow = c(1, 1), mar = c(1, 1, 3, 1), fg = "grey70", col = "black")
plot(ref, cex = 0.35, no.margin = FALSE, edge.color = clade_edge_cols(ref))
nodelabels(n_present, bg = count_cols[n_present + 1], cex = 0.45, frame = "circle")
title(sprintf("IQ-TREE clades: number of the other %d trees that contain each clade", length(others)))
legend("bottomleft", legend = 0:length(others), fill = count_cols, title = "trees with clade", bty = "n", cex = 0.7)

## 6. consensus trees
par(mfrow = c(1, 2), fg = "black")
for (p in c(1, 0.5)) {
  ct <- consensus(trees, p = p, rooted = FALSE)
  plot(ct, cex = 0.3)
  title(if (p == 1) "strict" else "majority rule")
}

dev.off()
