# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
.this_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]))
if (!file.exists(file.path(.this_dir, "_config.R"))) .this_dir <- "scripts"
source(file.path(.this_dir, "_config.R"))

#!/usr/bin/env Rscript
# ============================================================
# Export plotting data for revised Figures 2-4 (one small CSV per panel group).
# All values recomputed from cache/matrix.csv with the settings used in the revised text (scripts 10, 30, 85, 92).
# limma: Ritchie ME et al., Nucleic Acids Res 2015, PMID 25605792
# ============================================================
suppressMessages({library(limma); library(jsonlite)})
CACHE <- cfg$interim_dir
RES   <- cfg$results_dir
OUT   <- cfg$figdata_dir
dir.create(OUT, showWarnings = FALSE)
S <- cfg$samples
m <- read.csv(file.path(CACHE, "matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
pr <- m[[1]]
IP  <- as.matrix(m[, paste0(S, "(IP,Normalized)")]);  colnames(IP)  <- S; rownames(IP)  <- pr
SUP <- as.matrix(m[, paste0(S, "(Sup,Normalized)")]); colnames(SUP) <- S; rownames(SUP) <- pr
grp <- factor(cfg$groups, levels = c(cfg$reference, cfg$case)); design <- model.matrix(~ grp)
mp <- rbind(cbind(read.csv(file.path(CACHE, "probe_map_mrna.csv"))[, c("probeID","Transcript_ID","GeneSymbol")], cls = "mRNA"),
            cbind(read.csv(file.path(CACHE, "probe_map_lncrna.csv"))[, c("probeID","Transcript_ID","GeneSymbol")], cls = "lncRNA"),
            cbind(read.csv(file.path(CACHE, "probe_map_others.csv"))[, c("probeID","Transcript_ID","GeneSymbol")], cls = "other_ncRNA"))
mp <- mp[match(pr, mp$probeID), ]; stopifnot(identical(mp$probeID, pr))

# ---- Fig 2b: standard PCA (IP) ----
p <- prcomp(t(IP), center = TRUE, scale. = FALSE); ve <- 100 * p$sdev^2 / sum(p$sdev^2)
write.csv(data.frame(sample = S, group = as.character(grp), PC1 = p$x[, 1], PC2 = p$x[, 2]), file.path(OUT, "fig2b_pca.csv"), row.names = FALSE)
# ---- Fig 2c: correlation matrix (IP) ----
write.csv(round(cor(IP), 4), file.path(OUT, "fig2c_cor_ip.csv"))
# ---- Fig 2d: within-condition pairs (IP), all probes ----
pairs <- list(c("MOCK_1","MOCK_2"), c("MOCK_1","MOCK_3"), c("MOCK_2","MOCK_3"), c("PR8_1","PR8_2"), c("PR8_1","PR8_3"), c("PR8_2","PR8_3"))
write.csv(IP, file.path(OUT, "fig2d_ip_matrix.csv"))
write.csv(data.frame(x = sapply(pairs, `[`, 1), y = sapply(pairs, `[`, 2),
                     r = sapply(pairs, function(q) round(cor(IP[, q[1]], IP[, q[2]]), 3))), file.path(OUT, "fig2d_pairs_r.csv"), row.names = FALSE)
# ---- Fig 2e: per-sample distributions (IP and Sup) quantiles ----
qs <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
dist <- do.call(rbind, lapply(c("IP","Sup"), function(L) { X <- if (L == "IP") IP else SUP
  do.call(rbind, lapply(S, function(s) data.frame(layer = L, sample = s, q = qs, value = quantile(X[, s], qs)))) }))
write.csv(SUP, file.path(OUT, "fig2e_sup_matrix.csv"))
write.csv(dist, file.path(OUT, "fig2e_quantiles.csv"), row.names = FALSE)

# ---- Fig 3: limma IP and Sup (deposited) ----
fitL <- function(X) topTable(eBayes(lmFit(X, design)), coef = 2, number = Inf, sort.by = "none")
tIP <- fitL(IP); tSUP <- fitL(SUP)
f3 <- data.frame(probeID = pr, transcript = mp$Transcript_ID, gene = mp$GeneSymbol, rna_class = mp$cls,
                 ip_lfc = tIP$logFC, ip_fdr = tIP$adj.P.Val, sup_lfc = tSUP$logFC, sup_fdr = tSUP$adj.P.Val)
f3$ip_call <- ifelse(f3$ip_fdr < 0.05 & f3$ip_lfc >= 1, "higher_PR8", ifelse(f3$ip_fdr < 0.05 & f3$ip_lfc <= -1, "higher_MOCK", "ns"))
write.csv(f3, file.path(OUT, "fig3_limma_ip_sup.csv"), row.names = FALSE)
print(table(f3$rna_class, f3$ip_call))
mach <- c(METTL3 = "writer", METTL14 = "writer", WTAP = "writer", RBM15 = "writer", RBM15B = "writer", KIAA1429 = "writer", ZC3H13 = "writer",
          FTO = "eraser", ALKBH5 = "eraser", YTHDF1 = "reader", YTHDF2 = "reader", YTHDF3 = "reader", YTHDC1 = "reader", YTHDC2 = "reader",
          HNRNPC = "reader", HNRNPA2B1 = "reader", IGF2BP1 = "reader", IGF2BP2 = "reader", IGF2BP3 = "reader")
f3d <- f3[f3$rna_class == "mRNA" & f3$gene %in% names(mach), ]
f3d$role <- mach[f3d$gene]
write.csv(f3d, file.path(OUT, "fig3d_machinery.csv"), row.names = FALSE)
cat("machinery probes:", nrow(f3d), "genes:", length(unique(f3d$gene)), "\n")

# ---- Fig 4: %modified, detected probes, deposited and joint QN ----
det <- read.csv(file.path(RES, "probe_detection_flags.csv")); keep <- pr %in% det$probeID[det$sup_detected_n == 6]
fr <- function(a, b) { f <- 2^a / (2^a + 2^b); pmin(pmax(f, 1e-6), 1 - 1e-6) }
qn <- normalizeBetweenArrays(cbind(IP, SUP), method = "quantile")
one <- function(Fm, tag) {
  F <- Fm[keep, ]; tt <- topTable(eBayes(lmFit(log(F / (1 - F)), design)), coef = 2, number = Inf, sort.by = "none")
  d <- 100 * (rowMeans(F[, 4:6]) - rowMeans(F[, 1:3]))
  data.frame(probeID = pr[keep], dpp = d, fdr = tt$adj.P.Val, call = ifelse(tt$adj.P.Val < 0.05 & d >= 5, "higher_PR8", ifelse(tt$adj.P.Val < 0.05 & d <= -5, "higher_MOCK", "ns")))
}
a <- one(fr(IP, SUP), "dep"); b <- one(fr(qn[, 1:6], qn[, 7:12]), "qn")
supQN <- fitL(qn[, 7:12])$logFC
f4 <- data.frame(probeID = pr[keep], transcript = mp$Transcript_ID[keep], gene = mp$GeneSymbol[keep], rna_class = mp$cls[keep],
                 expr_lfc_dep = tSUP$logFC[keep], expr_lfc_qn = supQN[keep],
                 dpp_dep = a$dpp, fdr_dep = a$fdr, call_dep = a$call, dpp_qn = b$dpp, fdr_qn = b$fdr, call_qn = b$call)
f4$robust <- f4$call_dep != "ns" & f4$call_dep == f4$call_qn
write.csv(f4, file.path(OUT, "fig4_pctmod_detected.csv"), row.names = FALSE)
km <- f4$rna_class == "mRNA"
summ <- list(
  pca_var = round(ve[1:2], 2),
  ip_counts = as.list(table(f3$rna_class, f3$ip_call)["mRNA", ]),
  lnc_counts = as.list(table(f3$rna_class, f3$ip_call)["lncRNA", ]),
  sup_mrna = c(up = sum(f3$rna_class == "mRNA" & f3$sup_fdr < 0.05 & f3$sup_lfc >= 1), down = sum(f3$rna_class == "mRNA" & f3$sup_fdr < 0.05 & f3$sup_lfc <= -1)),
  pct = c(n = nrow(f4), dep_up = sum(f4$call_dep == "higher_PR8"), dep_down = sum(f4$call_dep == "higher_MOCK"),
          robust = sum(f4$robust), robust_up = sum(f4$robust & f4$call_dep == "higher_PR8"), robust_down = sum(f4$robust & f4$call_dep == "higher_MOCK"),
          median_dep = round(median(f4$dpp_dep), 3), median_qn = round(median(f4$dpp_qn), 3),
          n_mrna = sum(km), r_dep = round(cor(f4$expr_lfc_dep[km], f4$dpp_dep[km] / 100), 4), r_qn = round(cor(f4$expr_lfc_qn[km], f4$dpp_qn[km] / 100), 4)))
write_json(summ, file.path(OUT, "figdata_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(summ$pct)
