# 09_H4_sensibilidades.R
# Correlations between total-score changes and biomarker changes (Holm); sensitivity analyses S1
# (WAT-1), S2 (FLACC), and S3 (stable classification).
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
library(boot)
if (!exists("teste_delta_rho"))
  source(file.path("scripts", "00b_motor_deltarho.R"))

D   <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
VAR <- attr(D, "VAR"); ITENS_D1 <- attr(D, "ITENS_D1")
ITENS_D5 <- attr(D, "ITENS_D5"); ROTULO_ITEM <- attr(D, "ROTULO_ITEM")
A  <- D[D$avaliavel_capd, ]
nA <- nrow(A)

boot_bca <- function(dat, estat, R = CFG$n_boot, seed = CFG$seed) {
  obs <- estat(dat, seq_len(nrow(dat)))
  set.seed(seed)
  bo <- boot::boot(dat, statistic = estat, R = R)
  ic <- c(NA, NA); tipo <- "BCa"
  ci <- try(boot::boot.ci(bo, type = "bca", conf = 0.95), silent = TRUE)
  if (!inherits(ci, "try-error") && !is.null(ci$bca)) ic <- ci$bca[4:5]
  else {
    cp <- try(boot::boot.ci(bo, type = "perc", conf = 0.95), silent = TRUE)
    if (!inherits(cp, "try-error")) { ic <- cp$percent[4:5]; tipo <- "percentil" }
  }
  se <- stats::sd(bo$t[, 1], na.rm = TRUE)
  if (all(is.na(ic))) tipo <- "indisponivel"
  list(est = obs, ic_inf = ic[1], ic_sup = ic[2], tipo_ic = tipo, se = se)
}

cat("\n=== H4 — VARIACAO INTRAINDIVIDUAL (n =", nA, ") ===\n")
cat("Estimando: rho(delta log biomarcador, delta CAPD total), D5 - D1\n")
cat("Duas comparacoes, ajuste de Holm. Todos os testes bilaterais.\n\n")

d_capd <- A[[VAR$capd_d5]] - A[[VAR$capd_d1]]
linhas <- list()
for (bio in c("NfL", "S100B")) {
  v1 <- if (bio == "NfL") VAR$nfl_d1 else VAR$s100b_d1
  v5 <- if (bio == "NfL") VAR$nfl_d5 else VAR$s100b_d5
  d_bio <- log(A[[v5]]) - log(A[[v1]])
  dat <- data.frame(a = d_bio, b = d_capd)
  r  <- boot_bca(dat, function(d, i) sp(d$a[i], d$b[i]))
  ct <- suppressWarnings(cor.test(d_bio, d_capd, method = "spearman",
                                  exact = FALSE, alternative = "two.sided"))
  linhas[[length(linhas) + 1]] <- data.frame(
    Biomarcador = bio, n = sum(is.finite(d_bio) & is.finite(d_capd)),
    rho = r$est, IC_inf = r$ic_inf, IC_sup = r$ic_sup, Tipo_IC = r$tipo_ic,
    p_bruto = unname(ct$p.value), stringsAsFactors = FALSE)
  cat(sprintf("  %-6s rho = %+.3f  [%.3f a %.3f]\n", bio, r$est, r$ic_inf, r$ic_sup))
}
H4 <- do.call(rbind, linhas)
H4$p_holm <- p.adjust(H4$p_bruto, method = CFG$ajuste_S1, n = 2)

write.csv(H4, file.path(DIR_OUT, "T16_H4_numerico.csv"), row.names=FALSE)
saveRDS(H4, file.path(DIR_OUT, "T16_H4_numerico.rds"))
saida_H4 <- data.frame(
  Biomarcador = H4$Biomarcador, n = H4$n,
  rho = sprintf("%.3f", H4$rho),
  IC95_BCa = sprintf("%.3f a %.3f", H4$IC_inf, H4$IC_sup),
  p_bruto = fmt_p(H4$p_bruto), p_holm = fmt_p(H4$p_holm),
  stringsAsFactors = FALSE)
cat("\n=== RESULTADOS H4 ===\n"); print(saida_H4, row.names = FALSE)
write.csv(saida_H4, file.path(DIR_OUT, "T16_H4_deltas_totais.csv"), row.names = FALSE)

sensibilidade <- function(rotulo, dados, covar_d1, covar_d5, arquivo) {
  cat("\n\n=== ", rotulo, " ===\n", sep = "")
  parcial <- !is.null(covar_d1)
  nS <- nrow(dados)
  cat("n = ", nS, if (parcial) "  | ajuste por covariavel continua\n" else
                  "  | restricao de subamostra\n", sep = "")

  cat("\n-- H1 --\n")
  lin <- list()
  for (mom in c("D1", "D5")) {
    vcapd <- if (mom == "D1") VAR$capd_d1 else VAR$capd_d5
    vcov  <- if (parcial) (if (mom == "D1") covar_d1 else covar_d5) else NA
    for (bio in c("NfL", "S100B")) {
      vbio <- if (bio == "NfL") (if (mom == "D1") VAR$nfl_d1 else VAR$nfl_d5)
              else              (if (mom == "D1") VAR$s100b_d1 else VAR$s100b_d5)
      dat <- if (parcial) data.frame(x = dados[[vcapd]], y = dados[[vbio]], z = dados[[vcov]])
             else         data.frame(x = dados[[vcapd]], y = dados[[vbio]])
      dat <- dat[stats::complete.cases(dat), ]
      rs <- if (parcial) boot_bca(dat, function(d, i) sp_parcial(d$x[i], d$y[i], d$z[i]))
            else         boot_bca(dat, function(d, i) sp(d$x[i], d$y[i]))
      rb <- if (parcial) boot_bca(dat, function(d, i) sp(d$x[i], d$y[i])) else NULL
      est_na <- if (parcial) rb$est else NA_real_
      lin[[length(lin) + 1]] <- data.frame(
        Camada = "H1", Momento = mom, Alvo = bio, n = nrow(dat),
        estimativa_nao_ajustada = est_na, estimativa_sensibilidade = rs$est,
        IC_inf = rs$ic_inf, IC_sup = rs$ic_sup, SE = rs$se,
        tipo_IC = rs$tipo_ic,
        estado = if (rs$tipo_ic == "indisponivel") "IC indisponivel"
                 else if (rs$tipo_ic == "percentil") "IC percentil (BCa falhou)"
                 else "ok",
        stringsAsFactors = FALSE)
      cat(sprintf("  %s %-6s %s%+.3f  [%.3f a %.3f]  (%s)\n", mom, bio,
                  if (parcial) sprintf("nao ajustada %+.3f -> ajustada ", rb$est) else "",
                  rs$est, rs$ic_inf, rs$ic_sup, rs$tipo_ic))
    }
  }
  A1 <- do.call(rbind, lin)

  cat("\n-- H2: Delta_rho por dominio --\n")
  lin <- list(); k <- 0
  for (mom in c("D1", "D5")) {
    itens <- if (mom == "D1") ITENS_D1 else ITENS_D5
    vnfl  <- if (mom == "D1") VAR$nfl_d1   else VAR$nfl_d5
    vs100 <- if (mom == "D1") VAR$s100b_d1 else VAR$s100b_d5
    vcov  <- if (parcial) (if (mom == "D1") covar_d1 else covar_d5) else NA
    for (i in 1:8) {
      k <- k + 1
      cat(sprintf("  [%2d/16] Q%d %s ... ", k, i, mom)); flush.console()
      r <- if (parcial)
             teste_delta_parcial(dados[[itens[i]]], dados[[vnfl]], dados[[vs100]], dados[[vcov]])
           else
             teste_delta_rho(dados[[itens[i]]], dados[[vnfl]], dados[[vs100]], sobreposto = TRUE)
      lin[[k]] <- data.frame(
        Camada = "H2", Momento = mom, Alvo = paste0("Q", i, " ", ROTULO_ITEM[i]),
        n = r$n, estimativa_nao_ajustada = NA_real_,
        estimativa_sensibilidade = r$delta,
        IC_inf = r$ic_inf, IC_sup = r$ic_sup, SE = r$se_obs,
        tipo_IC = if (is.na(r$ic_inf)) "indisponivel" else "bootstrap-t",
        estado = r$estado, stringsAsFactors = FALSE)
      cat(r$estado, "\n")
    }
  }
  A2 <- do.call(rbind, lin)

  BRUTO <- rbind(A1, A2)
  write.csv(BRUTO, file.path(DIR_OUT, sub("\\.csv$", "_bruto.csv", arquivo)),
            row.names = FALSE)

  saida <- data.frame(
    Camada = BRUTO$Camada, Momento = BRUTO$Momento, Alvo = BRUTO$Alvo, n = BRUTO$n,
    estimativa_nao_ajustada = ifelse(is.na(BRUTO$estimativa_nao_ajustada), "",
                                     sprintf("%.3f", BRUTO$estimativa_nao_ajustada)),
    estimativa_sensibilidade = sprintf("%.3f", BRUTO$estimativa_sensibilidade),
    SE = sprintf("%.3f", BRUTO$SE),
    IC95 = sprintf("%.3f a %.3f", BRUTO$IC_inf, BRUTO$IC_sup),
    tipo_IC = BRUTO$tipo_IC, estado = BRUTO$estado, stringsAsFactors = FALSE)
  print(saida, row.names = FALSE)
  write.csv(saida, file.path(DIR_OUT, arquivo), row.names = FALSE)
  cat(sprintf("\n%d estimativas nesta sensibilidade. Nenhum valor-p, nenhum ajuste.\n",
              nrow(BRUTO)))
  invisible(BRUTO)
}

mde_spearman <- function(n, alfa = 0.05, poder = 0.80, k = 1.03)
  tanh((qnorm(1 - alfa/2) + qnorm(poder)) * k / sqrt(n - 3))

S1 <- sensibilidade("S1 — AJUSTE POR WAT-1 CONTINUO (abstinencia iatrogenica)",
                    A, VAR$wat1_d1, VAR$wat1_d5, "T17_S1_ajuste_WAT1.csv")
S2 <- sensibilidade("S2 — AJUSTE POR FLACC CONTINUO (dor)",
                    A, VAR$dor_d1, VAR$dor_d5, "T18_S2_ajuste_FLACC.csv")

d1_pos  <- A[[VAR$capd_d1]] >= CFG$capd_corte
d5_pos  <- A[[VAR$capd_d5]] >= CFG$capd_corte
E  <- A[d1_pos == d5_pos, ]
nE <- nrow(E)
cat("\n\n=== S3 — PRECISAO DA SUBAMOSTRA ===\n")
cat(sprintf("  Fenotipo estavel : %s   | transicoes excluidas: %s\n",
            fmt_contagem(nE, nA), fmt_contagem(nA - nE, nA)))
cat(sprintf("  |rho| minimo detectavel: %.2f (n = %d) contra %.2f (n = %d)\n",
            mde_spearman(nE), nE, mde_spearman(nA), nA))
cat("  Reportar esta perda de precisao ao lado do resultado.\n")

S3 <- sensibilidade("S3 — RESTRICAO AO FENOTIPO ESTAVEL",
                    E, NULL, NULL, "T19_S3_fenotipo_estavel.csv")

cat(sprintf("\n=== TOTAL DE ESTIMATIVAS DE SENSIBILIDADE: %d ===\n",
            nrow(S1) + nrow(S2) + nrow(S3)))
cat("Nenhuma entra no inventario de 102 testes: nao ha valor-p nem ajuste.\n")

leitura <- c(
"COMO LER AS SENSIBILIDADES",
"--------------------------",
"1. Elas NAO REPORTAM valor-p nem ajuste de multiplicidade: a saida traz",
"   estimativa, erro-padrao e intervalo. O valor-p calculado internamente",
"   pelo motor de bootstrap-t NAO e exportado nem usado.",
"   Nenhuma destas estimativas entra no inventario",
"   de 102 testes.",
"",
"2. A consistencia com a analise principal e julgada por ESTIMATIVA, DIRECAO",
"   e INCERTEZA. NAO por sobreposicao de intervalos, que nao e teste, nem por",
"   'sobreviveu ao ajuste numa e nao na outra', que e a falacia de comparar",
"   significancias.",
"",
"3. S1 e S2 sao ajustes SEPARADOS. Nao ha analise ajustada",
"   simultaneamente por WAT-1 e FLACC.",
"",
"4. A correlacao parcial foi recalculada DENTRO de cada reamostra bootstrap.",
"   O intervalo corresponde ao estimando ajustado, e nao ao bruto.",
"",
"5. LIMITE DA VERIFICACAO: a verificacao sintetica do teste de Delta-rho",
"   avaliou somente a correlacao de Spearman SIMPLES. A variante PARCIAL de",
"   S1 e S2 compartilha o pivo e o erro-padrao do procedimento principal,",
"   mas nao esta coberta pela verificacao sintetica (Apendice 2 do Material Suplementar).",
"",
"",
"6. S3 reduz a amostra e a precisao. O |rho| minimo detectavel da subamostra",
"   e reportado acima e deve acompanhar o resultado.",
"",
"7. Nao ha quarta sensibilidade. Em particular, nao foi feita",
"   restricao por exposicao farmacologica no Dia 5: deixaria cerca de 18",
"   pacientes e elevaria o |rho| minimo detectavel de 0,48 para 0,63.",
"",
"8. Os arquivos *_bruto.csv trazem os valores SEM arredondamento, incluindo",
"   erro-padrao e tipo de intervalo. As tabelas formatadas derivam deles.",
"",
"9. A coluna estimativa_nao_ajustada so e preenchida em S1 e S2, onde ha de",
"   fato duas estimativas na MESMA amostra: sem e com a covariavel. Em S3 ela",
"   fica vazia — a estimativa nao ajustada na subamostra E a propria",
"   estimativa de S3. A comparacao com a analise principal se faz contra a",
"   tabela T6 do script 04, nunca recalculando na subamostra.")
writeLines(leitura, file.path(DIR_OUT, "09_sensibilidades_leitura.txt"))
cat("\n"); cat(leitura, sep = "\n"); cat("\n")

cat("\n=== SCRIPT 09 CONCLUIDO ===\n")
cat("Ultimo script do pipeline analitico.\n\n")
