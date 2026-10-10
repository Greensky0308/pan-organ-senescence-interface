suppressMessages({library(coloc); library(data.table)})
dn <- fread('out/coloc_input_scan/_disease_N.csv')
files <- list.files('cache/narrow', pattern = '\\.csv$', full.names = TRUE)
res <- list()
for (f in files) {
  m <- fread(f)
  m <- m[!is.na(beta.f) & !is.na(se.f) & se.f > 0 & !is.na(beta.p) & se.p > 0 & !is.na(eaf.p) & !is.na(eaf.f)]
  m <- m[!duplicated(m$snp)]
  if (nrow(m) < 30) next
  g <- m$gene[1]; o <- m$organ[1]
  p0 <- m$pos[which.min(m$p)]
  m <- m[abs(m$pos - p0) <= 300000]
  d <- dn[organ == o]
  r <- suppressWarnings(coloc.abf(
    list(snp = m$snp, beta = m$beta.p, varbeta = m$se.p^2, MAF = pmin(m$eaf.p, 1-m$eaf.p), N = m$n[1], type = 'quant'),
    list(snp = m$snp, beta = m$by,   varbeta = m$se.f^2, MAF = pmin(m$eaf.f, 1-m$eaf.f), N = d$N, s = d$s, type = 'cc')))
  s <- r$summary
  res[[paste(g, o)]] <- data.frame(gene = g, organ = o, n = nrow(m),
                                   H3 = s[['PP.H3.abf']], H4 = s[['PP.H4.abf']])
}
R <- do.call(rbind, res); R <- R[order(-R$H4), ]
write.csv(R, 'out/coloc_narrow_results.csv', row.names = FALSE)
print(R, row.names = FALSE, digits = 3)
