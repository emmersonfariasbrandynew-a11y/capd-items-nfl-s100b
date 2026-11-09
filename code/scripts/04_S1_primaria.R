# 04_S1_primaria.R
# Total-score correlations with NfL and S100B on days 1 and 5 (Holm adjustment).
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
library(boot)
set.seed(CFG$seed)

D   <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
VAR <- attr(D, "VAR")
A   <- D[D$avaliavel_capd, ]
nA  <- nrow(A)

cor_bca <- function(x, y, R = CFG$n_boot, seed = CFG$seed) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  est <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE,
                                   alternative = CFG$alternativa))

  set.seed(seed)
  bo <- boot::boot(data.frame(x = x, y = y), R = R,
                   statistic = function(d, i)
                     suppressWarnings(cor(d$x[i], d$y[i], method = "spearman")))

  ic <- c(NA, NA); tipo <- CFG$boot_type
  ci <- try(boot::boot.ci(bo, type = "bca", conf = 0.95), silent = TRUE)
  if (!inherits(ci, "try-error") && !is.null(ci$bca)) {
    ic <- ci$bca[4:5]
  } else {
    cp <- try(boot::boot.ci(bo, type = "perc", conf = 0.95), silent = TRUE)
    if (!inherits(cp, "try-error")) { ic <- cp$percent[4:5]; tipo <- "percentil" }
    warning("BCa indisponivel; IC por percentil.")
  }
  list(n = n, rho = unname(est$estimate), p = est$p.value,
       ic_inf = ic[1], ic_sup = ic[2], tipo_ic = tipo,
       n_empates_x = n - length(unique(x)), n_empates_y = n - length(unique(y)),
       boot = bo)
}

COMPARACOES <- list(
  list(rot = "CAPD total D1 x NfL D1",   vx = VAR$capd_d1, vy = VAR$nfl_d1),
  list(rot = "CAPD total D1 x S100B D1", vx = VAR$capd_d1, vy = VAR$s100b_d1),
  list(rot = "CAPD total D5 x NfL D5",   vx = VAR$capd_d5, vy = VAR$nfl_d5),
  list(rot = "CAPD total D5 x S100B D5", vx = VAR$capd_d5, vy = VAR$s100b_d5)
)

cat("\n=== H1 — ESCORE TOTAL x BIOMARCADORES (n =", nA, ") ===\n")
cat("Bootstrap BCa com", CFG$n_boot, "reamostragens. Pode levar 1-2 minutos.\n\n")

res  <- vector("list", length(COMPARACOES))
boots <- vector("list", length(COMPARACOES))
for (k in seq_along(COMPARACOES)) {
  cp <- COMPARACOES[[k]]
  cat("  [", k, "/", length(COMPARACOES), "] ", cp$rot, "\n", sep = "")
  r <- cor_bca(A[[cp$vx]], A[[cp$vy]])
  boots[[k]] <- r$boot; r$boot <- NULL
  res[[k]] <- data.frame(Comparacao = cp$rot, n = r$n,
                         rho = r$rho, p_bruto = r$p,
                         IC_inf = r$ic_inf, IC_sup = r$ic_sup,
                         Tipo_IC = r$tipo_ic,
                         Empates_CAPD = r$n_empates_x,
                         stringsAsFactors = FALSE)
}
T6 <- do.call(rbind, res)

T6$p_holm <- p.adjust(T6$p_bruto, method = CFG$ajuste_S1, n = 4)
T6$Significativo_Holm <- ifelse(T6$p_holm < 0.05, "sim", "nao")

T6$Magnitude <- cut(abs(T6$rho), breaks = c(-Inf, .10, .30, .50, Inf),
                    labels = c("desprezivel","pequena","moderada","grande"))

write.csv(T6, file.path(DIR_OUT, "T6_H1_numerico.csv"), row.names=FALSE)
saveRDS(T6, file.path(DIR_OUT, "T6_H1_numerico.rds"))
saida <- data.frame(
  Comparacao = T6$Comparacao,
  n          = T6$n,
  rho        = sprintf("%.3f", T6$rho),
  IC95       = sprintf("%.3f a %.3f", T6$IC_inf, T6$IC_sup),
  Tipo_IC    = T6$Tipo_IC,
  p_bruto    = fmt_p(T6$p_bruto),
  p_Holm     = fmt_p(T6$p_holm),
  Holm_sig   = T6$Significativo_Holm,
  Magnitude  = as.character(T6$Magnitude),
  stringsAsFactors = FALSE)

cat("\n=== RESULTADOS S1 ===\n")
print(saida, row.names = FALSE)
write.csv(saida, file.path(DIR_OUT, "T6_S1_CAPDtotal_biomarcadores.csv"), row.names = FALSE)
saveRDS(list(tabela = T6, boots = boots),
        file.path(DIR_OUT, "S1_boot_objetos.rds"))

cat("\n=== COMO LER ESTA TABELA ===\n")
cat("1. O IC BCa e a medida de incerteza destas correlacoes. Um IC que cruza zero indica que a\n")
cat("   direcao da associacao nao esta estabelecida, mesmo com p pequeno.\n")
cat(sprintf("2. Com n = %d, o estudo so tinha 80%% de poder para |rho| >= %.2f\n",
            nA, tanh((qnorm(0.975)+qnorm(0.80))*1.03/sqrt(nA-3))))
cat("   sob alfa 0.05. Resultados nulos abaixo disso NAO excluem associacao.\n")
cat("3. Achados nulos indicam que a associacao nao foi detectada,\n")
cat("   e nao que ela inexista\n")
cat("   (precisao limitada pelo tamanho amostral).\n")
cat("4. A coluna Empates_CAPD indica quantas observacoes compartilham posto.\n")
cat("   Muitos empates reduzem a precisao do IC bootstrap; se algum valor for\n")
cat("   alto, o IC deve ser lido com cautela.\n")

cat("\n=== SCRIPT 04 CONCLUIDO ===\n")
cat("Proximo passo: 05_S2_dominios.R\n\n")
