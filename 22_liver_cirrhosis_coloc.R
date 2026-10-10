suppressMessages({library(coloc); library(data.table)})
IN <- 'out/causal'
genes <- c('MMP7', 'MMP12')
res <- list()
fg <- fread(file.path(IN, 'fg_liver_cirrhosis.csv'))
fg <- fg[!duplicated(fg$snp)]
for (g in genes) {
  pq <- fread(file.path(IN, paste0('pqtl_', g, '_cirrhosis.csv')))
  pq <- pq[!is.na(beta) & !is.na(se) & se > 0]
  pq <- pq[!duplicated(pq$snp)]
  mg <- merge(pq, fg, by = 'snp', suffixes = c('.p', '.f'))
  mg <- mg[!is.na(beta.f) & !is.na(se.f) & se.f > 0]
  if (nrow(mg) < 50) next
  d1 <- list(snp = mg$snp, beta = mg$beta.p, varbeta = mg$se.p^2, MAF = mg$MAF.p, N = mg$n[1], type = 'quant')
  d2 <- list(snp = mg$snp, beta = mg$beta.f, varbeta = mg$se.f^2, MAF = mg$MAF.f, N = mg$N[1], s = mg$s[1], type = 'cc')
  cz <- suppressWarnings(coloc.abf(d1, d2))
  sm <- cz$summary
  top <- cz$results$snp[which.max(cz$results$SNP.PP.H4)]
  res[[g]] <- data.frame(gene = g, organ = 'liver', nsnp = as.integer(sm[['nsnps']]),
                         H3 = sm[['PP.H3.abf']], H4 = sm[['PP.H4.abf']], top_snp = top, stringsAsFactors = FALSE)
}
R <- do.call(rbind, res)
write.csv(R, 'out/coloc_full_cirrhosis.csv', row.names = FALSE)
cat('\n=== Coloc_full (pQTL x cirrhosis) ===\n'); print(R, row.names = FALSE, digits = 3)
