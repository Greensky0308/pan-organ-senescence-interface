suppressMessages({library(limma); library(edgeR)})

COH <- list(
  kidney = list(ctrl = c("LD", "TN"), cases = c("DN", "FSGS", "MCD", "RPGN"), k = 2),
  colon  = list(ctrl = "control",      cases = c("ulcerative colitis", "CD"), k = 1),
  skin   = list(ctrl = "control",      cases = "case",                      k = 1),
  lung   = list(ctrl = "control",      cases = "ipf",                       k = 1),
  liver  = list(ctrl = "control",      cases = c("cirrhosis"), k = 1),
  heart  = list(ctrl = "Non-Failing Donor", cases = "Dilated cardiomyopathy (DCM)", k = 1),
  brain  = list(ctrl = "control", cases = "case", k = 1),
  artery = list(ctrl = "control", cases = "case", k = 1),
  esophagus = list(ctrl = "control", cases = "case", k = 1),
  adipose = list(ctrl = "control", cases = "case", k = 1),
  endometrium = list(ctrl = "control", cases = "case", k = 1)
)

is_counts <- function(m) all(m >= 0) && all(abs(m - round(m)) < 1e-6)

run_de <- function(expr, g, con, dis) {
  idx <- which(g %in% c(con, dis))
  e <- expr[, idx]; gr <- factor(ifelse(g[idx] %in% dis, "dis", "con"), levels = c("con", "dis"))
  if (any(table(gr) < 2)) return(NULL)
  if (is_counts(e)) {
    e <- e[rowSums(e) > 10, , drop = FALSE]
    dge <- DGEList(counts = round(e)); dge <- calcNormFactors(dge)
    v <- voom(dge, model.matrix(~ gr), plot = FALSE)
    fit <- eBayes(lmFit(v, model.matrix(~ gr)))
  } else {
    e <- e[rowMeans(e, na.rm = TRUE) > 0, , drop = FALSE]
    fit <- eBayes(lmFit(e, model.matrix(~ gr)))
  }
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, FDR = tt$adj.P.Val)
}

up_list <- list()
for (org in names(COH)) {
  cfg <- COH[[org]]
  expr <- as.matrix(read.delim(sprintf("out/cohorts/%s_expr.tsv", org), row.names = 1, check.names = FALSE))
  grp <- read.delim(sprintf("out/cohorts/%s_groups.tsv", org))
  g <- setNames(grp$stratum, grp$sample)
  g <- g[colnames(expr)]
  ups <- list()
  for (cs in cfg$cases) {
    de <- run_de(expr, g, cfg$ctrl, cs)
    if (is.null(de)) { cat(org, cs, "skip\n"); next }
    write.csv(de, sprintf("out/cohorts/%s_DE_%s.csv", org, gsub("[^A-Za-z0-9]", "_", cs)), row.names = FALSE)
    ups[[cs]] <- de$gene[de$FDR < 0.05 & de$logFC >= 1]
    cat(sprintf("%s | %s vs ctrl: %d up (of %d tested)\n", org, cs, length(ups[[cs]]), nrow(de)))
  }
  tab <- table(unlist(ups))
  cons <- names(tab)[tab >= cfg$k & names(tab) %in% unique(unlist(ups))]
  up_list[[org]] <- cons
  cat(sprintf(">>> %s consensus up (>=%d/%d): %d genes\n\n", org, cfg$k, length(cfg$cases), length(cons)))
}
saveRDS(up_list, "out/organ_disease_up.rds")
for (o in names(up_list)) writeLines(up_list[[o]], sprintf("out/organ_up_%s.txt", o))
