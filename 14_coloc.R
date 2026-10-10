suppressMessages({library(coloc); library(data.table)})
IN <- 'out/coloc_input'
genes  <- c('MMP7','MMP12')
organs <- c('kidney','liver','heart','lung','brain','colon','skin','artery','esophagus')

res <- list()
for (g in genes) {
  p <- fread(file.path(IN, paste0('pqtl_', g, '.csv')))
  p <- p[!is.na(beta) & !is.na(se) & se > 0]
  p <- p[!duplicated(p$snp)]
  for (o in organs) {
    f <- fread(file.path(IN, paste0('fg_', o, '.csv')))
    f <- f[!duplicated(f$snp)]
    m <- merge(p, f, by = 'snp', suffixes = c('.p', '.f'))
    m <- m[!is.na(beta.f) & !is.na(se.f) & se.f > 0]
    if (nrow(m) < 50) next
    d1 <- list(snp = m$snp, beta = m$beta.p, varbeta = m$se.p^2, MAF = m$MAF.p, N = m$n[1], type = 'quant')
    d2 <- list(snp = m$snp, beta = m$beta.f, varbeta = m$se.f^2, MAF = m$MAF.f, N = m$N[1], s = m$s[1], type = 'cc')
    r <- suppressWarnings(coloc.abf(d1, d2))
    s <- r$summary
    top <- r$results$snp[which.max(r$results$SNP.PP.H4)]
    res[[paste(g, o)]] <- data.frame(gene = g, organ = o, nsnp = as.integer(s[['nsnps']]),
                                     H3 = s[['PP.H3.abf']], H4 = s[['PP.H4.abf']], top_snp = top,
                                     stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, res)
R <- R[order(R$gene, -R$H4), ]
write.csv(R, 'out/coloc_results.csv', row.names = FALSE)
print(R, row.names = FALSE, digits = 3)
