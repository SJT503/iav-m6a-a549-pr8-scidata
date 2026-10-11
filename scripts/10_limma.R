# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
.this_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]))
if (!file.exists(file.path(.this_dir, "_config.R"))) .this_dir <- "scripts"
source(file.path(.this_dir, "_config.R"))

#!/usr/bin/env Rscript
# ============================================================
# Step 1: limma + BH FDR 重分析（IP / Sup / %Modified 三层）
# 输入: data/interim/matrix.csv (40,955 probes x 12 channels)
# 输出: results/limma_*.csv + summary json
# 参考: R1-M1, R2-3
# ============================================================
suppressMessages({library(limma); library(jsonlite)})

CACHE <- cfg$interim_dir
RES   <- cfg$results_dir
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

m <- read.csv(file.path(CACHE, "matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cat("Matrix:", nrow(m), "probes x", ncol(m), "cols\n")
probe <- m[[1]]
mat <- as.matrix(m[, -1])
colnames(mat) <- colnames(m)[-1]
rownames(mat) <- probe   # keep probe IDs: the limma output and the t-test cross-check below both rely on them

IPcols  <- grep("(IP,Normalized)",  colnames(mat), value = TRUE, fixed = TRUE)
Supcols <- grep("(Sup,Normalized)", colnames(mat), value = TRUE, fixed = TRUE)
# 统一顺序 MOCK_1..3, PR8_1..3
ord <- function(x) {
  s <- sub("^([A-Z0-9]+_[0-9]).*$", "\1", x)
  x[order(match(s, c("MOCK_1","MOCK_2","MOCK_3","PR8_1","PR8_2","PR8_3")))]
}
IPcols <- ord(IPcols); Supcols <- ord(Supcols)
cat("IP cols :", IPcols, "\n")
cat("Sup cols:", Supcols, "\n")

IP  <- mat[, IPcols,  drop = FALSE]
Sup <- mat[, Supcols, drop = FALSE]
grp <- factor(cfg$groups, levels = c(cfg$reference, cfg$case))
design <- model.matrix(~ grp)

# ---- %Modified: IP/(IP+Sup) on linear (2^) scale ----
lin <- function(x) 2^x
Pct <- (lin(IP) / (lin(IP) + lin(Sup))) * 100
colnames(Pct) <- sub("(IP,Normalized)", "%Modified", IPcols, fixed = TRUE)

run_limma <- function(X, tag) {
  ok <- rowSums(is.na(X)) == 0
  fit <- lmFit(X[ok, , drop = FALSE], design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$probeID <- rownames(tt)
  tt$n_tested <- sum(ok)
  tt$n_dropped_NA <- sum(!ok)
  write.csv(tt, file.path(RES, sprintf("limma_%s.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s] tested=%d droppedNA=%d | p<0.05: %d | FDR<0.05: %d | p<0.05 & |log2FC|>=1: %d | FDR<0.05 & |log2FC|>=1: %d\n",
      tag, sum(ok), sum(!ok), sum(tt$P.Value < 0.05, na.rm = TRUE),
      sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
      sum(tt$P.Value < 0.05 & abs(tt$logFC) >= 1, na.rm = TRUE),
      sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1, na.rm = TRUE)))
  tt
}

res_IP  <- run_limma(IP,  "IP_m6Aquantity")
res_Sup <- run_limma(Sup, "Sup_expression")
res_Pct <- run_limma(Pct, "PctModified")

# ---- 复现旧口径（unpaired t-test）对照 ----
ttest_count <- function(X, tag) {
  ok <- rowSums(is.na(X)) == 0
  Xo <- X[ok, , drop = FALSE]
  p <- apply(Xo, 1, function(v) tryCatch(t.test(v[4:6], v[1:3])$p.value, error = function(e) NA))
  fc <- rowMeans(Xo[, 4:6], na.rm = TRUE) - rowMeans(Xo[, 1:3], na.rm = TRUE)  # log2 差值
  cat(sprintf("[%s] t-test p<0.05: %d | p<0.05 & |dLog2|>=1: %d\n",
      tag, sum(p < 0.05, na.rm = TRUE),
      sum(p < 0.05 & abs(fc) >= 1, na.rm = TRUE)))
  data.frame(probeID = rownames(Xo), t_p = p, t_dlog2 = fc)
}
tt_IP  <- ttest_count(IP,  "IP")
tt_Sup <- ttest_count(Sup, "Sup")

summary <- list(
  matrix_dim        = dim(mat),
  IP_cols           = IPcols,
  Sup_cols          = Supcols,
  probes_total      = nrow(m),
  limma = list(
    IP  = list(tested = sum(rowSums(is.na(IP)) == 0),  FDR05 = sum(res_IP$adj.P.Val < 0.05, na.rm=TRUE),
               p05 = sum(res_IP$P.Value < 0.05, na.rm=TRUE),
               FDR05_FC2 = sum(res_IP$adj.P.Val < 0.05 & abs(res_IP$logFC) >= 1, na.rm=TRUE),
               p05_FC2   = sum(res_IP$P.Value < 0.05 & abs(res_IP$logFC) >= 1, na.rm=TRUE)),
    Sup = list(tested = sum(rowSums(is.na(Sup)) == 0), FDR05 = sum(res_Sup$adj.P.Val < 0.05, na.rm=TRUE),
               p05 = sum(res_Sup$P.Value < 0.05, na.rm=TRUE),
               FDR05_FC2 = sum(res_Sup$adj.P.Val < 0.05 & abs(res_Sup$logFC) >= 1, na.rm=TRUE),
               p05_FC2   = sum(res_Sup$P.Value < 0.05 & abs(res_Sup$logFC) >= 1, na.rm=TRUE)),
    Pct = list(tested = sum(rowSums(is.na(Pct)) == 0), FDR05 = sum(res_Pct$adj.P.Val < 0.05, na.rm=TRUE),
               p05 = sum(res_Pct$P.Value < 0.05, na.rm=TRUE))
  )
)
write_json(summary, file.path(RES, "limma_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat("\nSaved to", RES, "\n")
