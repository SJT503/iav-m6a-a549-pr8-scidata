# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
.this_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]))
if (!file.exists(file.path(.this_dir, "_config.R"))) .this_dir <- "scripts"
source(file.path(.this_dir, "_config.R"))

#!/usr/bin/env Rscript
# ============================================================
# Step 9: %modified analysis restricted to probes with a measurable denominator
#   keep = Sup (Cy3) signal well above background (FE gIsWellAboveBG == 1) on all six arrays
#   Normalisation is done on all probes; testing (limma, logit scale, BH FDR) on kept probes only.
# limma: Ritchie ME et al., Nucleic Acids Res 2015, PMID 25605792 (casebook-verified)
# ============================================================
suppressMessages({library(limma); library(jsonlite)})
CACHE <- cfg$interim_dir
RES   <- cfg$results_dir
S <- cfg$samples
m <- read.csv(file.path(CACHE, "matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
probe <- m[[1]]; mat <- as.matrix(m[, -1]); rownames(mat) <- probe
IP  <- mat[, paste0(S, "(IP,Normalized)")];  colnames(IP)  <- S
SUP <- mat[, paste0(S, "(Sup,Normalized)")]; colnames(SUP) <- S
grp <- factor(cfg$groups, levels = c(cfg$reference, cfg$case)); design <- model.matrix(~ grp)

detf <- read.csv(file.path(RES, "probe_detection_flags.csv"), stringsAsFactors = FALSE)
keep <- probe %in% detf$probeID[detf$sup_detected_n == 6]
mp <- rbind(read.csv(file.path(CACHE, "probe_map_mrna.csv"))[, c("probeID", "GeneSymbol")],
            read.csv(file.path(CACHE, "probe_map_lncrna.csv"))[, c("probeID", "GeneSymbol")],
            read.csv(file.path(CACHE, "probe_map_others.csv"))[, c("probeID", "GeneSymbol")])
cls <- ifelse(probe %in% read.csv(file.path(CACHE, "probe_map_mrna.csv"))$probeID, "mRNA",
              ifelse(probe %in% read.csv(file.path(CACHE, "probe_map_lncrna.csv"))$probeID, "lncRNA", "other_ncRNA"))
gene <- mp$GeneSymbol[match(probe, mp$probeID)]
cat(sprintf("probes kept (Sup detected on all 6 arrays): %d of %d | mRNA %d, lncRNA %d, other %d\n",
            sum(keep), length(keep), sum(keep & cls == "mRNA"), sum(keep & cls == "lncRNA"), sum(keep & cls == "other_ncRNA")))

off <- apply(IP - SUP, 2, median)
qn <- normalizeBetweenArrays(cbind(IP, SUP), method = "quantile")
NS <- list(SI = list(IP = IP, SUP = SUP), CEN = list(IP = sweep(IP, 2, off), SUP = SUP),
           QNj = list(IP = qn[, 1:6], SUP = qn[, 7:12]))
fitp <- function(X) {
  fr <- 2^X$IP / (2^X$IP + 2^X$SUP); fr <- pmin(pmax(fr, 1e-6), 1 - 1e-6)
  fr <- fr[keep, , drop = FALSE]; lg <- log(fr / (1 - fr))
  tt <- topTable(eBayes(lmFit(lg, design)), coef = 2, number = Inf, sort.by = "none")
  data.frame(probeID = rownames(fr), rna_class = cls[keep], GeneSymbol = gene[keep],
             dfrac = rowMeans(fr[, 4:6]) - rowMeans(fr[, 1:3]), FDR = tt$adj.P.Val, P = tt$P.Value,
             stringsAsFactors = FALSE)
}
PR <- lapply(NS, fitp)
callp <- function(d) ifelse(d$FDR < 0.05 & d$dfrac >= 0.05, 1L, ifelse(d$FDR < 0.05 & d$dfrac <= -0.05, -1L, 0L))
tab <- do.call(rbind, lapply(names(PR), function(nm) {
  d <- PR[[nm]]; c0 <- callp(d); E <- abs(d$dfrac) >= 0.05
  data.frame(norm = nm, n_tested = nrow(d), median_dfrac_pp = round(100 * median(d$dfrac), 3),
             pct_up = round(100 * mean(d$dfrac > 0), 1), eff_up = sum(E & d$dfrac > 0), eff_down = sum(E & d$dfrac < 0),
             fdr_eff_up = sum(c0 == 1), fdr_eff_down = sum(c0 == -1))
}))
cat("\n=== %modified, detection-filtered ===\n"); print(tab, row.names = FALSE)
cS <- callp(PR$SI); cQ <- callp(PR$QNj); cC <- callp(PR$CEN)
rob <- cS != 0 & cS == cQ
cat(sprintf("\nrobust (SI & QNj, same direction): %d (higher PR8 %d / higher mock %d); also in CEN: %d; SI calls kept in QNj: %d/%d\n",
            sum(rob), sum(rob & cS == 1), sum(rob & cS == -1), sum(rob & cC == cS), sum(rob), sum(cS != 0)))
print(table(class = PR$SI$rna_class[rob], direction = ifelse(cS[rob] == 1, "higher_PR8", "higher_MOCK")))
robdf <- data.frame(PR$SI[rob, c("probeID", "rna_class", "GeneSymbol")], dfrac_SI = PR$SI$dfrac[rob], FDR_SI = PR$SI$FDR[rob],
                    dfrac_QNj = PR$QNj$dfrac[rob], FDR_QNj = PR$QNj$FDR[rob],
                    direction = ifelse(cS[rob] == 1, "higher_PR8", "higher_MOCK"))
robdf <- robdf[order(-abs(robdf$dfrac_SI)), ]
write.csv(robdf, file.path(RES, "pctmod_robust_set_detected.csv"), row.names = FALSE)
cat("\ntop robust probes:\n"); print(head(robdf[, c("GeneSymbol", "rna_class", "dfrac_SI", "dfrac_QNj", "direction")], 16), row.names = FALSE)

# decoupling on detected mRNA probes: expression log2FC from limma on Sup (all probes), restricted to kept mRNA probes
supSI <- topTable(eBayes(lmFit(SUP, design)), coef = 2, number = Inf, sort.by = "none")$logFC
# expression for the quantile-normalised decoupling must come from the same joint QN as the %modified values (QNj),
# otherwise the pair mixes a Sup-only normalisation with a joint normalisation and disagrees with Figure 4b.
supQN <- topTable(eBayes(lmFit(qn[, 7:12], design)), coef = 2, number = Inf, sort.by = "none")$logFC
km <- cls[keep] == "mRNA"
dec <- function(x, y, tag) {
  de <- abs(x) >= 1; se <- abs(x) < 0.5
  cat(sprintf("   %-4s n=%d r=%+.3f | |expr lfc|>=1: %d, of which |dfrac|<0.05: %d (%.1f%%) | |lfc|<0.5 & |dfrac|>=0.1: %d\n",
              tag, length(x), cor(x, y), sum(de), sum(de & abs(y) < 0.05), 100 * mean(abs(y[de]) < 0.05), sum(se & abs(y) >= 0.1)))
  c(r = cor(x, y), n_de = sum(de), de_stable = sum(de & abs(y) < 0.05), pct = 100 * mean(abs(y[de]) < 0.05), se_changed = sum(se & abs(y) >= 0.1))
}
cat("\n=== decoupling, detected mRNA probes ===\n")
dSI <- dec(supSI[keep][km], PR$SI$dfrac[km], "SI")
dQN <- dec(supQN[keep][km], PR$QNj$dfrac[km], "QN")
write_json(list(n_kept = sum(keep), table = tab, robust = list(n = sum(rob), up = sum(rob & cS == 1), down = sum(rob & cS == -1),
                also_CEN = sum(rob & cC == cS), SI_calls = sum(cS != 0)),
                decoupling = list(SI = as.list(dSI), QN = as.list(dQN))),
           file.path(RES, "pctmod_detected_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 6)
cat("\nSaved -> pctmod_detected_summary.json, pctmod_robust_set_detected.csv\n")
