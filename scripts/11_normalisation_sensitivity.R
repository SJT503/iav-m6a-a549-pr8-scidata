# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
.this_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]))
if (!file.exists(file.path(.this_dir, "_config.R"))) .this_dir <- "scripts"
source(file.path(.this_dir, "_config.R"))

#!/usr/bin/env Rscript
# ============================================================
# Step 8: R2-4 global Delta-%modified shift -- decomposition + normalisation sensitivity
#   Normalisations: SI  = as deposited (provider calibration-spike-in normalisation)
#                   CEN = SI with per-array median log2(IP/Sup) removed (array-centred ratio)
#                   QNj = limma quantile normalisation of all 12 channels jointly
#                   QNs = quantile normalisation of IP channels and Sup channels separately
#   Layers: IP (m6A quantity), Sup (expression), %modified (logit scale for limma)
# limma: Ritchie ME et al., Nucleic Acids Res 2015, PMID 25605792 (casebook-verified)
# Input: cache/matrix.csv (deposited Matrix.xlsx) + cache/probe_map_{mrna,lncrna,others}.csv
# ============================================================
suppressMessages({library(limma); library(jsonlite)})

CACHE <- cfg$interim_dir
RES   <- cfg$results_dir
S <- cfg$samples

m <- read.csv(file.path(CACHE, "matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
probe <- m[[1]]
mat <- as.matrix(m[, -1]); rownames(mat) <- probe
IP  <- mat[, paste0(S, "(IP,Normalized)")];  colnames(IP)  <- S
SUP <- mat[, paste0(S, "(Sup,Normalized)")]; colnames(SUP) <- S
grp <- factor(cfg$groups, levels = c(cfg$reference, cfg$case))
design <- model.matrix(~ grp)

# ---- probe -> RNA class / gene (1:1 maps from the provider flag_pass tables) ----
rd_map <- function(f, cls) {
  x <- read.csv(file.path(CACHE, f), stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(probeID = x$probeID, Transcript_ID = x$Transcript_ID, GeneSymbol = x$GeneSymbol,
             rna_class = cls, region = if ("m6A_mRNA_region" %in% names(x)) x$m6A_mRNA_region else NA,
             stringsAsFactors = FALSE)
}
mp <- rbind(rd_map("probe_map_mrna.csv", "mRNA"), rd_map("probe_map_lncrna.csv", "lncRNA"),
            rd_map("probe_map_others.csv", "other_ncRNA"))
cat(sprintf("map rows=%d | unique probeID=%d | unique Transcript_ID=%d | matrix probes covered=%d/%d\n",
            nrow(mp), length(unique(mp$probeID)), length(unique(mp$Transcript_ID)),
            sum(probe %in% mp$probeID), length(probe)))
cls <- mp$rna_class[match(probe, mp$probeID)]
gene <- mp$GeneSymbol[match(probe, mp$probeID)]

# ---- A. decomposition of the array-level log2(IP/Sup) offset ----
# provider per-array factors, recovered exactly as (Normalised - log2 FE ProcessedSignal), sd < 1e-15
fIP  <- c(-10.65155, -10.17241, -10.85019, -10.54372, -10.31455, -10.12483)
fSUP <- c(-9.70681, -9.50237, -9.77775, -9.74515, -9.63554, -9.68526)
pos_in <- c(82.08, 78.26, 67.77, 85.03, 78.12, 87.75)   # positive m6A spike-in %Input (qPCR, MeRIP_QC.pdf)
neg_in <- c(1.65, 1.39, 1.18, 1.81, 2.15, 0.80)         # negative m6A spike-in %Input
LR <- IP - SUP
off_norm   <- apply(LR, 2, median)
off_factor <- fIP - fSUP
off_raw    <- off_norm - off_factor
med_pct    <- apply(100 * 2^IP / (2^IP + 2^SUP), 2, median)
offtab <- data.frame(sample = S, group = as.character(grp),
                     raw_median_log2ratio = round(off_raw, 4),
                     calibration_factor_IP_minus_Sup = round(off_factor, 4),
                     deposited_median_log2ratio = round(off_norm, 4),
                     deposited_median_pct_modified = round(med_pct, 2),
                     pos_spike_pct_input = pos_in, neg_spike_pct_input = neg_in)
write.csv(offtab, file.path(RES, "norm_offsets_table.csv"), row.names = FALSE)
cat("\n=== A. array-level offsets ===\n"); print(offtab, row.names = FALSE)
t_factor <- t.test(off_factor[4:6], off_factor[1:3])
t_raw    <- t.test(off_raw[4:6], off_raw[1:3])
t_pos    <- t.test(pos_in[4:6], pos_in[1:3])
c_pos    <- cor.test(off_norm, pos_in)
cat(sprintf("Welch: calibration-factor offset PR8-MOCK = %+.4f log2, t=%.3f df=%.2f p=%.4f\n",
            unname(diff(rev(t_factor$estimate))), t_factor$statistic, t_factor$parameter, t_factor$p.value))
cat(sprintf("Welch: raw (FE) median log ratio PR8-MOCK = %+.4f log2, p=%.4f\n",
            unname(diff(rev(t_raw$estimate))), t_raw$p.value))
cat(sprintf("Welch: positive spike %%Input PR8-MOCK = %+.2f, p=%.4f (MOCK %.1f+-%.1f; PR8 %.1f+-%.1f)\n",
            unname(diff(rev(t_pos$estimate))), t_pos$p.value, mean(pos_in[1:3]), sd(pos_in[1:3]),
            mean(pos_in[4:6]), sd(pos_in[4:6])))
cat(sprintf("Pearson: deposited median log ratio vs positive spike %%Input r=%.3f, p=%.4f (n=6)\n",
            c_pos$estimate, c_pos$p.value))

# ---- B. %modified under four normalisations ----
qn  <- normalizeBetweenArrays(cbind(IP, SUP), method = "quantile")
NS <- list(
  SI  = list(IP = IP, SUP = SUP),
  CEN = list(IP = sweep(IP, 2, off_norm), SUP = SUP),
  QNj = list(IP = qn[, 1:6], SUP = qn[, 7:12]),
  QNs = list(IP = normalizeBetweenArrays(IP, method = "quantile"),
             SUP = normalizeBetweenArrays(SUP, method = "quantile"))
)
pct_fit <- function(X) {
  fr <- 2^X$IP / (2^X$IP + 2^X$SUP); fr <- pmin(pmax(fr, 1e-6), 1 - 1e-6)
  lg <- log(fr / (1 - fr))
  tt <- topTable(eBayes(lmFit(lg, design)), coef = 2, number = Inf, sort.by = "none")
  welch_p <- apply(fr, 1, function(v) t.test(v[4:6], v[1:3])$p.value)
  data.frame(probeID = rownames(fr), rna_class = cls, GeneSymbol = gene,
             dfrac = rowMeans(fr[, 4:6]) - rowMeans(fr[, 1:3]), logitFC = tt$logFC,
             P = tt$P.Value, FDR = tt$adj.P.Val, welchP = welch_p, stringsAsFactors = FALSE)
}
PR <- lapply(NS, pct_fit)
cnt <- function(d, sel) c(up = sum(sel & d$dfrac > 0), down = sum(sel & d$dfrac < 0))
sens <- do.call(rbind, lapply(names(PR), function(nm) {
  d <- PR[[nm]]; E <- abs(d$dfrac) >= 0.05
  a <- cnt(d, E); b <- cnt(d, d$FDR < 0.05); cc <- cnt(d, d$FDR < 0.05 & E); w <- cnt(d, d$welchP < 0.05 & E)
  data.frame(norm = nm, median_dfrac_pp = round(100 * median(d$dfrac), 3), pct_probes_up = round(100 * mean(d$dfrac > 0), 1),
             eff_up = a[["up"]], eff_down = a[["down"]], fdr_up = b[["up"]], fdr_down = b[["down"]],
             fdr_eff_up = cc[["up"]], fdr_eff_down = cc[["down"]], welch_eff_up = w[["up"]], welch_eff_down = w[["down"]])
}))
cat("\n=== B. %modified sensitivity (probe level, all 40,955 probes) ===\n"); print(sens, row.names = FALSE)
write.csv(sens, file.path(RES, "pctmod_norm_sensitivity.csv"), row.names = FALSE)
for (nm in names(PR)) write.csv(PR[[nm]], file.path(RES, sprintf("limma_pct_%s.csv", nm)), row.names = FALSE)

callp <- function(d) ifelse(d$FDR < 0.05 & d$dfrac >= 0.05, 1L, ifelse(d$FDR < 0.05 & d$dfrac <= -0.05, -1L, 0L))
cSI <- callp(PR$SI); cQ <- callp(PR$QNj); cC <- callp(PR$CEN)
rob <- cSI != 0 & cSI == cQ
cat(sprintf("\nrobust %%modified set (FDR<0.05 & |dfrac|>=0.05, same direction in SI and QNj): %d (up %d / down %d); of these also in CEN: %d\n",
            sum(rob), sum(rob & cSI == 1), sum(rob & cSI == -1), sum(rob & cC == cSI)))
cat("   by RNA class:\n"); print(table(class = cls[rob], direction = ifelse(cSI[rob] == 1, "higher_PR8", "higher_MOCK")))
robdf <- data.frame(PR$SI[rob, c("probeID", "rna_class", "GeneSymbol")],
                    dfrac_SI = PR$SI$dfrac[rob], FDR_SI = PR$SI$FDR[rob],
                    dfrac_QNj = PR$QNj$dfrac[rob], FDR_QNj = PR$QNj$FDR[rob],
                    direction = ifelse(cSI[rob] == 1, "higher_PR8", "higher_MOCK"))
write.csv(robdf, file.path(RES, "pctmod_robust_set.csv"), row.names = FALSE)
cat("   SI-significant calls retained in QNj (same direction): ", sum(rob), "/", sum(cSI != 0), "\n")
cat("   SI x QNj call table:\n"); print(table(SI = cSI, QNj = cQ))

# ---- C. IP and Sup layers: SI vs single-channel QN ----
lay <- function(X) topTable(eBayes(lmFit(X, design)), coef = 2, number = Inf, sort.by = "none")
call2 <- function(t) ifelse(t$adj.P.Val < 0.05 & t$logFC >= 1, 1L, ifelse(t$adj.P.Val < 0.05 & t$logFC <= -1, -1L, 0L))
layer_out <- list()
for (L in c("IP", "SUP")) {
  a <- lay(NS$SI[[L]]); b <- lay(NS$QNs[[L]])
  ca <- call2(a); cb <- call2(b)
  keep <- sum(ca != 0 & ca == cb)
  bycls <- table(class = cls[ca != 0], dir = ifelse(ca[ca != 0] == 1, "higher_PR8", "higher_MOCK"))
  cat(sprintf("\n=== C. %s layer (FDR<0.05 & |log2FC|>=1) ===\n", L))
  cat(sprintf("   SI : up %d / down %d | QN : up %d / down %d | SI calls kept in QN (same dir) %d/%d (%.1f%%) | r(logFC SI,QN)=%.4f | median logFC SI=%+.3f QN=%+.3f\n",
              sum(ca == 1), sum(ca == -1), sum(cb == 1), sum(cb == -1), keep, sum(ca != 0),
              100 * keep / sum(ca != 0), cor(a$logFC, b$logFC), median(a$logFC), median(b$logFC)))
  cat("   SI calls by RNA class:\n"); print(bycls)
  write.csv(data.frame(probeID = probe, rna_class = cls, GeneSymbol = gene, logFC = a$logFC, P = a$P.Value, FDR = a$adj.P.Val),
            file.path(RES, sprintf("limma_%s_SI.csv", L)), row.names = FALSE)
  write.csv(data.frame(probeID = probe, rna_class = cls, GeneSymbol = gene, logFC = b$logFC, P = b$P.Value, FDR = b$adj.P.Val),
            file.path(RES, sprintf("limma_%s_QN.csv", L)), row.names = FALSE)
  layer_out[[L]] <- list(SI_up = sum(ca == 1), SI_down = sum(ca == -1), QN_up = sum(cb == 1), QN_down = sum(cb == -1),
                         kept = keep, r_logFC = round(cor(a$logFC, b$logFC), 4),
                         SI_by_class = as.data.frame.matrix(bycls))
  assign(paste0("lay_", L, "_SI"), a); assign(paste0("lay_", L, "_QN"), b)
}

# ---- D. decoupling (Figure 4 claim) under SI and QN, mRNA probes ----
isM <- cls == "mRNA"
dec <- function(x, y, tag) {
  de <- abs(x) >= 1; se <- abs(x) < 0.5
  o <- c(pearson = cor(x, y), spearman = cor(x, y, method = "spearman"),
         n_expr_changed = sum(de), expr_changed_frac_stable = sum(de & abs(y) < 0.05),
         expr_changed_frac_changed = sum(de & abs(y) >= 0.05),
         n_expr_stable = sum(se), expr_stable_frac_changed_0.1 = sum(se & abs(y) >= 0.1))
  cat(sprintf("   %-4s r=%+.3f rho=%+.3f | |expr lfc|>=1: n=%d, of which |dfrac|<0.05: %d (%.1f%%) | |lfc|<0.5 & |dfrac|>=0.1: %d\n",
              tag, o[["pearson"]], o[["spearman"]], o[["n_expr_changed"]], o[["expr_changed_frac_stable"]],
              100 * o[["expr_changed_frac_stable"]] / o[["n_expr_changed"]], o[["expr_stable_frac_changed_0.1"]]))
  o
}
cat("\n=== D. decoupling of expression change vs Delta-fraction (mRNA probes) ===\n")
d_si <- dec(lay_SUP_SI$logFC[isM], PR$SI$dfrac[isM], "SI")
d_qn <- dec(lay_SUP_QN$logFC[isM], PR$QNj$dfrac[isM], "QN")

# ---- E. ISG induction under SI vs QN (Sup layer, probe-level mean per gene) ----
isg <- c("IFNL1", "IFNB1", "RSAD2", "ISG15", "CCL5", "IFIT3", "IFNL3", "OAS1", "PML", "OASL", "OAS2", "IFIT1", "IRF7", "STAT1", "MX1")
ctl <- c("ACTB", "GAPDH", "TUBB", "RPLP0")
gtab <- do.call(rbind, lapply(c(isg, ctl), function(g) {
  i <- which(gene == g & isM)
  if (!length(i)) return(NULL)
  data.frame(gene = g, set = ifelse(g %in% isg, "ISG", "control"), n_probes = length(i),
             log2FC_SI = round(mean(lay_SUP_SI$logFC[i]), 3), log2FC_QN = round(mean(lay_SUP_QN$logFC[i]), 3),
             minFDR_SI = signif(min(lay_SUP_SI$adj.P.Val[i]), 3), minFDR_QN = signif(min(lay_SUP_QN$adj.P.Val[i]), 3))
}))
cat("\n=== E. ISG panel, Sup layer, SI vs QN ===\n"); print(gtab, row.names = FALSE)
write.csv(gtab, file.path(RES, "isg_panel_SI_vs_QN.csv"), row.names = FALSE)

out <- list(offsets = offtab,
            welch_calibration_offset = list(diff = unname(diff(rev(t_factor$estimate))), p = t_factor$p.value, df = unname(t_factor$parameter)),
            welch_raw_offset = list(diff = unname(diff(rev(t_raw$estimate))), p = t_raw$p.value),
            welch_pos_spike_input = list(diff = unname(diff(rev(t_pos$estimate))), p = t_pos$p.value),
            cor_offset_vs_pos_input = list(r = unname(c_pos$estimate), p = c_pos$p.value),
            pct_sensitivity = sens,
            robust_pct = list(n = sum(rob), up = sum(rob & cSI == 1), down = sum(rob & cSI == -1),
                              SI_calls = sum(cSI != 0), also_CEN = sum(rob & cC == cSI)),
            layers = layer_out, decoupling = list(SI = as.list(d_si), QN = as.list(d_qn)))
write_json(out, file.path(RES, "norm_sensitivity_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 6)
cat("\nSaved -> norm_sensitivity_summary.json and per-normalisation CSVs in", RES, "\n")
