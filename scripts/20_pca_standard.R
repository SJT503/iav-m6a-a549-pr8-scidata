# ---------------------------------------------------------------------------
# Paths and settings come from config/config.ini; nothing is machine-specific.
# ---------------------------------------------------------------------------
.this_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]))
if (!file.exists(file.path(.this_dir, "_config.R"))) .this_dir <- "scripts"
source(file.path(.this_dir, "_config.R"))

#!/usr/bin/env Rscript
# ============================================================
# Step 2: 标准中心化 PCA 重算 (R2-5, R1-minor 10)
# 输入: cache/matrix.csv (40955 probes x 6 IP channels + 6 Sup channels)
# 输出: results/pca_standard.csv, pca_variance.csv, pca_samples.csv
# 要点: prcomp 默认 center=TRUE, scale.=FALSE；同时给出 scale 版本作稳健性
# ============================================================
suppressMessages(library(jsonlite))

CACHE <- cfg$interim_dir
RES   <- cfg$results_dir

m <- read.csv(file.path(CACHE, "matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
mat <- as.matrix(m[, -1])
ipc  <- grep("(IP,Normalized)",  colnames(mat), fixed = TRUE, value = TRUE)
supc <- grep("(Sup,Normalized)", colnames(mat), fixed = TRUE, value = TRUE)
ordc <- function(x) {
  s <- sub("^([A-Za-z0-9]+_[0-9]).*$", "\\1", x)
  x[order(match(s, c("MOCK_1","MOCK_2","MOCK_3","PR8_1","PR8_2","PR8_3")))]
}
ipc <- ordc(ipc); supc <- ordc(supc)
grp <- cfg$groups

do_pca <- function(X, tag) {
  X <- t(X)                                   # samples x probes
  # 只用无缺失探针（此处应为全集）
  keep <- colSums(is.na(X)) == 0
  Xk <- X[, keep, drop = FALSE]
  pca <- prcomp(Xk, center = TRUE, scale. = FALSE)
  ve <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  sc <- as.data.frame(pca$x[, 1:min(6, ncol(pca$x))])
  sc$sample <- rownames(sc)
  sc$group  <- grp
  cat(sprintf("\n[%s] probes used = %d | n samples = %d\n", tag, sum(keep), nrow(Xk)))
  cat("   variance explained (%):", paste(sprintf("PC%d=%.2f", seq_along(ve[1:6]), ve[1:6]), collapse = "  "), "\n")
  cat("   PC scores (PC1, PC2):\n")
  for (i in seq_len(nrow(sc))) {
    cat(sprintf("     %-9s %-5s PC1=%+8.4f  PC2=%+8.4f\n",
                sc$sample[i], sc$group[i], sc$PC1[i], sc$PC2[i]))
  }
  write.csv(sc, file.path(RES, sprintf("pca_%s_scores.csv", tag)), row.names = FALSE)
  write.csv(data.frame(PC = seq_along(ve), var_explained_pct = ve),
            file.path(RES, sprintf("pca_%s_variance.csv", tag)), row.names = FALSE)
  invisible(list(scores = sc, var = ve, npc_probes = sum(keep)))
}

cat("=========================================================\n")
cat("PCA from deposited Matrix.xlsx (mean-centred, prcomp)\n")
cat("=========================================================\n")
r_ip  <- do_pca(mat[, ipc,  drop = FALSE], "IP_standard")
r_sup <- do_pca(mat[, supc, drop = FALSE], "Sup_standard")

# 稳健性：scale.=TRUE 版本
do_pca_scaled <- function(X, tag) {
  X <- t(X); keep <- colSums(is.na(X)) == 0; Xk <- X[, keep, drop = FALSE]
  pca <- prcomp(Xk, center = TRUE, scale. = TRUE)
  ve <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  cat(sprintf("\n[%s scaled] variance: %s\n", tag,
              paste(sprintf("PC%d=%.2f", 1:4, ve[1:4]), collapse = "  ")))
  write.csv(data.frame(PC = seq_along(ve), var_explained_pct = ve),
            file.path(RES, sprintf("pca_%s_variance.csv", tag)), row.names = FALSE)
  invisible(ve)
}
do_pca_scaled(mat[, ipc, drop = FALSE], "IP_scaled")

# 记录旧版(provider) PCA 供对照
old <- list(PC1 = 50.3, PC2 = 48.0, note = "manuscript Figure 2a as submitted (provider-computed)")
out <- list(
  standard_IP  = list(var_pct = round(r_ip$var[1:6], 3),  probes = r_ip$npc_probes),
  standard_Sup = list(var_pct = round(r_sup$var[1:6], 3), probes = r_sup$npc_probes),
  provider_as_submitted = old
)
write_json(out, file.path(RES, "pca_summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat("\nSaved PCA results to", RES, "\n")
