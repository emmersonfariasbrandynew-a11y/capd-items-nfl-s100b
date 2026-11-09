# 07_H5_comparadores.R
# Contrasts between item correlations and correlations with pain (FLACC) and withdrawal (WAT-1)
# scores (single Benjamini-Hochberg family of 64 tests).
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

N_FAMILIA_H5 <- 64

cor_bca <- function(x, y, R = CFG$n_boot, seed = CFG$seed) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
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
  list(n = length(x), rho = sp(x, y), ic_inf = ic[1], ic_sup = ic[2], tipo_ic = tipo)
}

cat("\n=== H5 — MATRIZ COMPARADORA, CAMADA DESCRITIVA (n =", nA, ") ===\n\n")

linhas <- list(); k <- 0
for (mom in c("D1", "D5")) {
  itens <- if (mom == "D1") ITENS_D1 else ITENS_D5
  vdor  <- if (mom == "D1") VAR$dor_d1  else VAR$dor_d5
  vwat  <- if (mom == "D1") VAR$wat1_d1 else VAR$wat1_d5
  vcapd <- if (mom == "D1") VAR$capd_d1 else VAR$capd_d5
  for (bio in c("NfL", "S100B")) {
    vbio <- if (bio == "NfL") (if (mom == "D1") VAR$nfl_d1 else VAR$nfl_d5)
            else              (if (mom == "D1") VAR$s100b_d1 else VAR$s100b_d5)
    alvos <- c(setNames(as.list(itens), paste0("Q", 1:8)),
               list("CAPD total" = vcapd, "FLACC" = vdor, "WAT-1" = vwat))
    for (nm in names(alvos)) {
      k <- k + 1
      r <- cor_bca(A[[vbio]], A[[alvos[[nm]]]])
      linhas[[k]] <- data.frame(
        Momento = mom, Biomarcador = bio, Medida = nm,
        Tipo = if (nm %in% paste0("Q", 1:8)) "dominio CAPD"
               else if (nm == "CAPD total") "CAPD total" else "comparador",
        Rotulo = ifelse(nm %in% paste0("Q", 1:8),
                        ROTULO_ITEM[as.integer(sub("Q", "", nm))], nm),
        n = r$n, rho = r$rho, IC_inf = r$ic_inf, IC_sup = r$ic_sup,
        Tipo_IC = r$tipo_ic, stringsAsFactors = FALSE)
    }
    cat(sprintf("  %s %s — 11 medidas\n", mom, bio))
  }
}
M <- do.call(rbind, linhas)
M <- M[order(M$Momento, M$Biomarcador, -abs(M$rho)), ]
write.csv(M[, setdiff(names(M), c("p_bruto", "p")), drop=FALSE], file.path(DIR_OUT, "T11_matriz_numerico.csv"), row.names=FALSE)
saveRDS(M[, setdiff(names(M), c("p_bruto", "p")), drop=FALSE], file.path(DIR_OUT, "T11_matriz_numerico.rds"))
saida_M <- data.frame(
  Momento = M$Momento, Biomarcador = M$Biomarcador, Medida = M$Medida,
  Rotulo = M$Rotulo, Tipo = M$Tipo, n = M$n,
  rho = sprintf("%.3f", M$rho),
  IC95_BCa = sprintf("%.3f a %.3f", M$IC_inf, M$IC_sup),
  Tipo_IC = M$Tipo_IC, stringsAsFactors = FALSE)
print(saida_M, row.names = FALSE)
write.csv(saida_M, file.path(DIR_OUT, "T11_H5_matriz_comparadora.csv"), row.names = FALSE)

cat("\n=== RANKING DE |rho| (descritivo, comparavel a Matsuishi 2022) ===\n")
cat("Ordenacao descritiva. O ranking NAO e teste e nao sustenta afirmacao de\n")
cat("que um dominio se associe mais que outro.\n")
for (mm in c("D1", "D5")) for (bb in c("NfL", "S100B")) {
  sub <- M[M$Momento == mm & M$Biomarcador == bb & M$Medida != "CAPD total", ]
  sub <- sub[order(-abs(sub$rho)), ]
  cat(sprintf("\n  %s / %s — cinco medidas mais associadas:\n", mm, bb))
  for (j in 1:min(5, nrow(sub)))
    cat(sprintf("    %d. %-11s %-32s rho = %+.3f  [%s]\n", j, sub$Medida[j],
                substr(sub$Rotulo[j], 1, 32), sub$rho[j], sub$Tipo[j]))
}

baixa_var <- function(x) length(unique(x)) < CFG$item_min_niveis ||
                        sd(x) < CFG$item_min_sd

decl <- c(
 "DECLARACAO DE MULTIPLICIDADE — H5",
 "----------------------------------",
 sprintf("Amostra analitica                 : n = %d", nA),
 "Desenho                           : FATORIAL COMPLETO, sem reducao",
 "Momentos                          : D1 e D5",
 "Dominios                          : Q1 a Q8, todos os oito",
 "Comparadores                      : FLACC e WAT-1 (continuos)",
 "Biomarcadores                     : NfL e S100B",
 "",
 sprintf("TOTAL DE CONTRASTES               : %d  (8 x 2 x 2 x 2)", N_FAMILIA_H5),
 "",
 sprintf("Controle de multiplicidade        : Benjamini-Hochberg, q < %.2f",
         CFG$limiar_q_S2),
 "Aplicado em FAMILIA UNICA de 64: todos os contrastes sustentam a mesma",
 "conclusao sobre H5. Contrastes sem valor-p permanecem no denominador — o",
 "tamanho integral da familia e preservado no ajuste.",
 "",
 "A triagem de variabilidade dos itens NAO reduz esta contagem. Ela e",
 "reportada como coluna diagnostica. Ser cega ao desfecho evita selecao",
 "dependente do desfecho, o que nao equivale a garantia de controle do erro",
 "tipo I, e nenhuma garantia desse tipo e alegada.",
 "",
 "Este ajuste por familia NAO constitui controle global dos 102 testes",
 "do estudo. Os IC de 95% sao individuais, nao simultaneos."
)
writeLines(decl, file.path(DIR_OUT, "07_contagem_H5.txt"))
cat("\n\n"); cat(decl, sep = "\n"); cat("\n")

cat("\n=== EXECUTANDO", N_FAMILIA_H5, "CONTRASTES ===\n\n")
linhas <- list(); k <- 0
for (mom in c("D1", "D5")) {
  itens <- if (mom == "D1") ITENS_D1 else ITENS_D5
  comps <- list(FLACC   = if (mom == "D1") VAR$dor_d1  else VAR$dor_d5,
                "WAT-1" = if (mom == "D1") VAR$wat1_d1 else VAR$wat1_d5)
  for (bio in c("NfL", "S100B")) {
    vbio <- if (bio == "NfL") (if (mom == "D1") VAR$nfl_d1 else VAR$nfl_d5)
            else              (if (mom == "D1") VAR$s100b_d1 else VAR$s100b_d5)
    for (cn in names(comps)) for (i in 1:8) {
      k <- k + 1
      cat(sprintf("  [%2d/%d] %s %s Q%d vs %s ... ", k, N_FAMILIA_H5, mom, bio, i, cn))
      flush.console()
      r <- teste_delta_rho(A[[vbio]], A[[itens[i]]], A[[comps[[cn]]]],
                           sobreposto = TRUE)
      linhas[[k]] <- linha_delta(r, Momento = mom, Biomarcador = bio,
                                 Item = paste0("Q", i), Dominio = ROTULO_ITEM[i],
                                 Comparador = cn)
      linhas[[k]]$baixa_variabilidade <- ifelse(baixa_var(A[[itens[i]]]), "sim", "nao")
      cat(r$estado, "\n")
    }
  }
}
H5 <- do.call(rbind, linhas)
names(H5)[names(H5) == "rho1"] <- "rho_dominio"
names(H5)[names(H5) == "rho2"] <- "rho_comparador"

H5$q_BH <- p.adjust(H5$p_bootstrap_t, method = CFG$ajuste_S2, n = N_FAMILIA_H5)
H5$IC_exclui_zero <- ifelse(is.na(H5$IC_inf), NA_character_,
                     ifelse(H5$IC_inf > 0 | H5$IC_sup < 0, "sim", "nao"))
H5 <- H5[order(H5$q_BH, na.last = TRUE), ]

write.csv(H5, file.path(DIR_OUT, "T12_H5_numerico.csv"), row.names=FALSE)
saveRDS(H5, file.path(DIR_OUT, "T12_H5_numerico.rds"))
saida_H5 <- data.frame(
  Momento = H5$Momento, Biomarcador = H5$Biomarcador, Item = H5$Item,
  Dominio = H5$Dominio, Comparador = H5$Comparador, n = H5$n,
  rho_dominio    = sprintf("%.3f", H5$rho_dominio),
  rho_comparador = sprintf("%.3f", H5$rho_comparador),
  Delta_rho      = sprintf("%.3f", H5$Delta_rho),
  SE_jack        = sprintf("%.3f", H5$SE_jack),
  IC95_bootstrap_t = sprintf("%.3f a %.3f", H5$IC_inf, H5$IC_sup),
  IC_exclui_zero = H5$IC_exclui_zero,
  p_bootstrap_t = fmt_p(H5$p_bootstrap_t), q_BH = fmt_p(H5$q_BH),
  IC95_BCa_diagnostico = sprintf("%.3f a %.3f", H5$BCa_inf, H5$BCa_sup),
  baixa_variabilidade = H5$baixa_variabilidade,
  R_valido = H5$R_valido, estado = H5$estado, stringsAsFactors = FALSE)

cat("\n=== RESULTADOS H5 ===\n"); print(saida_H5, row.names = FALSE)
write.csv(saida_H5, file.path(DIR_OUT, "T12_H5_delta_vs_comparadores.csv"), row.names = FALSE)

n_falha <- sum(H5$estado != "ok")
if (n_falha > 0) {
  cat("\nATENCAO:", n_falha, "de", N_FAMILIA_H5, "contrastes sem valor-p:\n")
  print(table(H5$estado[H5$estado != "ok"]))
  cat("Permanecem na tabela com estimativa pontual e estado, e no denominador\n")
  cat("do ajuste. Devem ser relatados com a justificativa.\n")
}

cat("\n=== COMO LER =================================================\n")
cat("1. Delta_rho e a diferenca COM SINAL entre a associacao com o dominio e a\n")
cat("   associacao com o comparador. Positivo: maior com o dominio do CAPD.\n")
cat("   Negativo: maior com dor ou abstinencia — se ocorrer, e o achado mais\n")
cat("   importante da tabela e deve ser reportado com a mesma enfase.\n\n")
cat("2. IC contendo zero significa que a diferenca NAO foi estabelecida. NAO\n")
cat("   significa que dominio e comparador sejam equivalentes ou que os\n")
cat("   construtos sejam indistinguiveis. H5 compara associacoes, nao\n")
cat("   estabelece separabilidade.\n\n")
cat("3. O IC principal e o bootstrap-t, do mesmo pivo do valor-p. A coluna BCa\n")
cat("   e diagnostica: a verificacao sintetica avaliou o teste, nao a cobertura\n")
cat("   de intervalos. Discordancia entre os dois deve ser relatada, nunca\n")
cat("   usada para escolher o resultado mais favoravel.\n\n")
cat("4. baixa_variabilidade e DIAGNOSTICO. Nao removeu contrastes e nao altera\n")
cat("   o denominador do ajuste. Serve para o leitor calibrar a confianca.\n\n")
cat("5. Todas as hipoteses deste estudo sao EXPLORATORIAS. Nenhum q < 0.10\n")
cat("   autoriza afirmacao confirmatoria ou causal.\n")

cat("\n=== SCRIPT 07 CONCLUIDO ===\n")
cat("Proximo passo: 08_transicoes.R, depois 09_H4_sensibilidades.R\n\n")
