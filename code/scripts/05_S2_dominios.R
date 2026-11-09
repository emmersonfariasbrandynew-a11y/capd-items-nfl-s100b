# 05_S2_dominios.R
# Correlations of each CAPD item with each biomarker on days 1 and 5 (descriptive estimates with BCa
# intervals).
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
library(boot)
set.seed(CFG$seed)

D <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
VAR <- attr(D, "VAR"); ITENS_D1 <- attr(D, "ITENS_D1")
ITENS_D5 <- attr(D, "ITENS_D5"); ROTULO_ITEM <- attr(D, "ROTULO_ITEM")
A  <- D[D$avaliavel_capd, ]
nA <- nrow(A)

triagem <- function(itens, momento) {
  do.call(rbind, lapply(seq_along(itens), function(i) {
    x <- A[[itens[i]]]
    data.frame(Item = paste0("Q", i), Momento = momento,
               niveis = length(unique(x)), sd = round(sd(x), 3),
               prop_modal = round(max(table(x)) / length(x), 3),
               stringsAsFactors = FALSE)
  }))
}
TRI <- rbind(triagem(ITENS_D1, "D1"), triagem(ITENS_D5, "D5"))
TRI$baixa_variabilidade <- ifelse(
  TRI$niveis >= CFG$item_min_niveis & TRI$sd >= CFG$item_min_sd, "nao", "SIM")

cat("\n=== 5.1 TRIAGEM DE VARIABILIDADE (cega aos biomarcadores) ===\n")
cat("Criterio de baixa variabilidade: >=", CFG$item_min_niveis,
    "niveis distintos E desvio-padrao >=", CFG$item_min_sd, "\n")
cat("DIAGNOSTICO, nao filtro: nenhum item e removido de nenhuma familia.\n\n")
print(TRI, row.names = FALSE)

n_baixa <- sum(TRI$baixa_variabilidade == "SIM")
cat(sprintf("\nItens-momento sinalizados por baixa variabilidade : %d de %d\n",
            n_baixa, nrow(TRI)))
cat("TODOS permanecem nas familias de analise. A sinalizacao serve a\n")
cat("interpretacao: correlacoes com esses itens tem menor resolucao.\n")

n_desc_baixa <- n_baixa * 2
n_desc_tot   <- nrow(TRI) * 2

decl <- c(
 "CORRELACOES POR ITEM (SEM TESTE DE HIPOTESE)",
 "-----------------------------------------------------------",
 sprintf("Amostra analitica                      : n = %d", nA),
 sprintf("Itens do CAPD x momentos               : 8 x 2 = %d", nrow(TRI)),
 sprintf("Biomarcadores por item-momento         : 2 (NfL e S100B, contemporaneos)"),
 "",
 sprintf("TOTAL DE CORRELACOES DESCRITIVAS POR ITEM  : %d", n_desc_tot),
 sprintf("  das quais sinalizadas por baixa variabilidade : %d", n_desc_baixa),
 "  (sinalizacao NAO reduz o conjunto: todas sao reportadas)",
 "",
 "TESTES INFERENCIAIS NESTE SCRIPT           : 0",
 "",
 "Estas correlacoes NAO sao testadas. Reportam-se estimativa pontual e",
 "intervalo de confianca BCa. Nao ha valor-p e nao ha correcao de",
 "multiplicidade, porque nao ha teste.",
 "",
 "A inferencia no nivel de dominio ocorre no script 06 (H2), cujo estimando",
 "e a diferenca entre correlacoes, e no script 07 (H5).",
 "",
 "Nenhuma correlacao adicional foi conduzida e descartada."
)
writeLines(decl, file.path(DIR_OUT, "05_contagem_comparacoes.txt"))
cat("\n"); cat(decl, sep = "\n"); cat("\n")

cor_bca <- function(x, y, R = CFG$n_boot, seed = CFG$seed) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  est <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  set.seed(seed)
  bo <- boot::boot(data.frame(x = x, y = y), R = R,
                   statistic = function(d, i)
                     suppressWarnings(cor(d$x[i], d$y[i], method = "spearman")))
  ic <- c(NA, NA); tipo <- "BCa"
  ci <- try(boot::boot.ci(bo, type = "bca", conf = 0.95), silent = TRUE)
  if (!inherits(ci, "try-error") && !is.null(ci$bca)) ic <- ci$bca[4:5]
  else {
    cp <- try(boot::boot.ci(bo, type = "perc", conf = 0.95), silent = TRUE)
    if (!inherits(cp, "try-error")) { ic <- cp$percent[4:5]; tipo <- "percentil" }
  }
  list(n = n, rho = unname(est$estimate), p = est$p.value,
       ic_inf = ic[1], ic_sup = ic[2], tipo_ic = tipo)
}

cat("\n=== 5.4 CALCULANDO", n_desc_tot, "CORRELACOES DESCRITIVAS ===\n")
cat("Bootstrap BCa com", CFG$n_boot, "reamostragens cada. Aguarde.\n\n")

grade <- expand.grid(i = 1:8, momento = c("D1","D5"),
                     bio = c("NfL","S100B"), stringsAsFactors = FALSE)
linhas <- vector("list", nrow(grade))

for (k in seq_len(nrow(grade))) {
  g <- grade[k, ]
  item_col <- if (g$momento == "D1") ITENS_D1[g$i] else ITENS_D5[g$i]
  bio_col  <- if (g$bio == "NfL")
                (if (g$momento == "D1") VAR$nfl_d1   else VAR$nfl_d5)
              else
                (if (g$momento == "D1") VAR$s100b_d1 else VAR$s100b_d5)

  r <- cor_bca(A[[item_col]], A[[bio_col]])
  cam <- TRI$baixa_variabilidade[TRI$Item == paste0("Q", g$i) & TRI$Momento == g$momento]

  linhas[[k]] <- data.frame(
    Item = paste0("Q", g$i), Dominio = ROTULO_ITEM[g$i],
    Momento = g$momento, Biomarcador = g$bio, baixa_variabilidade = cam,
    n = r$n, rho = r$rho, IC_inf = r$ic_inf, IC_sup = r$ic_sup,
    Tipo_IC = r$tipo_ic, p_bruto = r$p, stringsAsFactors = FALSE)
  cat(sprintf("  [%2d/%d] Q%d %s x %s\n", k, nrow(grade), g$i, g$momento, g$bio))
}
T7 <- do.call(rbind, linhas)

T7$IC_exclui_zero <- ifelse(T7$IC_inf > 0 | T7$IC_sup < 0, "sim", "nao")
T7 <- T7[order(T7$baixa_variabilidade, -abs(T7$rho)), ]

write.csv(T7[, setdiff(names(T7), "p_bruto"), drop=FALSE], file.path(DIR_OUT, "T7_dominios_numerico.csv"), row.names=FALSE)
saveRDS(T7[, setdiff(names(T7), "p_bruto"), drop=FALSE], file.path(DIR_OUT, "T7_dominios_numerico.rds"))
fmt7 <- data.frame(
  Item = T7$Item, Dominio = T7$Dominio, Momento = T7$Momento,
  Biomarcador = T7$Biomarcador, baixa_variabilidade = T7$baixa_variabilidade, n = T7$n,
  rho  = sprintf("%.3f", T7$rho),
  IC95 = sprintf("%.3f a %.3f", T7$IC_inf, T7$IC_sup),
  IC_exclui_zero = T7$IC_exclui_zero, stringsAsFactors = FALSE)

write.csv(fmt7, file.path(DIR_OUT, "T7_S2_dominios_completo.csv"), row.names = FALSE)
write.csv(fmt7[fmt7$baixa_variabilidade == "nao", ],
          file.path(DIR_OUT, "T8_S2_itens_boa_variabilidade.csv"), row.names = FALSE)

cat("\n=== ITENS SEM SINALIZACAO DE BAIXA VARIABILIDADE ===\n")
print(fmt7[fmt7$baixa_variabilidade == "nao", ], row.names = FALSE)
cat("\n=== ITENS SINALIZADOS POR BAIXA VARIABILIDADE ===\n")
print(fmt7[fmt7$baixa_variabilidade == "SIM", ], row.names = FALSE)
cat("\nOs dois blocos tem o MESMO estatuto. A separacao e de leitura, nao de\n")
cat("camada: correlacoes do segundo bloco tem menor resolucao e mais empates,\n")
cat("e todas permanecem integralmente nas familias de H2 e H5.\n")

cat(sprintf("\nCorrelacoes cujo IC 95%% exclui zero: %d de %d.\n",
            sum(T7$IC_exclui_zero == "sim"), nrow(T7)))
cat("Este numero e DESCRITIVO. Nao constitui teste de hipotese, nao foi\n")
cat("corrigido para multiplicidade e NAO deve ser reportado como achado\n")
cat("significativo. A inferencia no nivel de dominio esta no script 06.\n")

cat("\n=== LIMITES DE INTERPRETACAO ===\n")
cat("1. Os 8 itens do CAPD NAO possuem propriedades psicometricas individuais\n")
cat("   publicadas. O instrumento foi validado como escore TOTAL. 'Dominio'\n")
cat("   aqui e agrupamento conceitual post-hoc, nao subescala validada. Isto\n")
cat("   e limitacao principal da analise por item.\n")
cat("2. Esta tabela nao sustenta afirmacao inferencial; descreve o\n")
cat("   padrao de associacao, sem implicar predicao.\n")
cat("3. Os itens do CAPD sobrepoem-se, em parte, a abstinencia\n")
cat("   e a dor. Nenhum achado desta tabela deve ser interpretado\n")
cat("   sem as analises de sensibilidade do script 09.\n")
cat("4. A tabela COMPLETA e reportada; apresentar apenas as associacoes\n")
cat("   sinalizadas seria relato seletivo.\n")

cat("\n=== SCRIPT 05 CONCLUIDO ===\n")
cat("Proximo passo: 06_H2H3_diferenca_rho.R\n\n")
