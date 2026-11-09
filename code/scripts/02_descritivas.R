# 02_descritivas.R
# Descriptive tables (sample characteristics, item distributions, concurrent screens).
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
D <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
VAR <- attr(D, "VAR"); ITENS_D1 <- attr(D, "ITENS_D1")
ITENS_D5 <- attr(D, "ITENS_D5"); ROTULO_ITEM <- attr(D, "ROTULO_ITEM")

A  <- D[D$avaliavel_capd, ]
nA <- nrow(A)

cat("\n=== AMOSTRA ANALITICA: n =", nA, "===\n\n")

linha_cont <- function(rotulo, x, dig = 1)
  data.frame(Variavel = rotulo, Valor = fmt_mediana_iqr(x, dig),
             stringsAsFactors = FALSE)

linha_cat <- function(rotulo, x, valor = 1)
  data.frame(Variavel = rotulo,
             Valor = fmt_contagem(sum(x == valor, na.rm = TRUE), nA),
             stringsAsFactors = FALSE)

T1 <- rbind(
  data.frame(Variavel = "-- DEMOGRAFIA --", Valor = ""),
  linha_cat ("Age group <1 year",               A$age_group, "<1"),
  linha_cat ("Age group 1-5 years",             A$age_group, "1-5"),
  linha_cat ("Age group >=6 years",             A$age_group, ">=6"),
  linha_cat ("Sexo masculino",                  A$sexo),
  data.frame(Variavel = "-- GRAVIDADE --", Valor = ""),
  linha_cont("PRISM IV, risco de mortalidade",  A[["PRISMIV%"]], 1),
  linha_cont("PELOD-2",                         A$PELOD2, 0),
  linha_cat ("Sepse",                           A$sepsis),
  linha_cat ("Choque",                          A$SHOCK),
  linha_cat ("SDRA",                            A$sdra),
  data.frame(Variavel = "-- SEDACAO E ANALGESIA --", Valor = ""),
  linha_cont("Duracao da sedacao, dias",        A[["sedation time"]], 1),
  linha_cont("RASS Dia 1",                      A[[VAR$rass_d1]], 0),
  linha_cont("RASS Dia 5",                      A[[VAR$rass_d5]], 0),
  linha_cat ("Fentanil",                        A[["fentanyl use"]]),
  linha_cat ("Morfina",                         A[["morphine use"]]),
  linha_cat ("Benzodiazepinico",                A$benzo),
  linha_cat ("Dexmedetomidina",                 A$precedx),
  linha_cat ("Cetamina",                        A$cetamina),
  linha_cat ("Propofol",                        A$Propofol),
  data.frame(Variavel = "-- ABSTINENCIA E DOR --", Valor = ""),
  linha_cat ("Sindrome de abstinencia iatrogenica",
             A[["Iatrogenic withdrawal syndrome (IWS)"]]),
  linha_cont("WAT-1 Dia 1",                     A[[VAR$wat1_d1]], 0),
  linha_cont("WAT-1 Dia 5",                     A[[VAR$wat1_d5]], 0),
  linha_cont("FLACC Dia 1",                     A[[VAR$dor_d1]], 0),
  linha_cont("FLACC Dia 5",                     A[[VAR$dor_d5]], 0),
  data.frame(Variavel = "-- ESTADO NEUROLOGICO --", Valor = ""),
  linha_cont("CAPD total Dia 1",                A[[VAR$capd_d1]], 0),
  linha_cont("CAPD total Dia 5",                A[[VAR$capd_d5]], 0),
  data.frame(Variavel = "Delirium (CAPD >= 9) no Dia 1",
             Valor = fmt_contagem(sum(A[[VAR$capd_d1]] >= CFG$capd_corte), nA)),
  data.frame(Variavel = "Delirium (CAPD >= 9) no Dia 5",
             Valor = fmt_contagem(sum(A[[VAR$capd_d5]] >= CFG$capd_corte), nA)),
  data.frame(Variavel = "-- DESFECHO --", Valor = ""),
  linha_cat ("Obito hospitalar",                A$obitooverall)
)

com_delirium <- A$tipodelierium %in% c(1, 2, 3)
nD <- sum(com_delirium)
sub_lab <- c(`1` = "hiperativo", `2` = "hipoativo", `3` = "misto")

T1b <- rbind(
  data.frame(Variavel = "-- CARACTERIZACAO DO DELIRIUM --", Valor = ""),
  data.frame(Variavel = "Pacientes com delirium (denominador dos itens abaixo)",
             Valor = fmt_contagem(nD, nA)),
  do.call(rbind, lapply(names(sub_lab), function(k) data.frame(
    Variavel = paste0("Subtipo motor: ", sub_lab[[k]]),
    Valor = fmt_contagem(sum(A$tipodelierium == as.integer(k), na.rm = TRUE), nD)))),
  data.frame(Variavel = "Tratamento farmacologico do delirium",
             Valor = fmt_contagem(sum(A$ttodelirium == 1, na.rm = TRUE), nD)),
  linha_cont("Duracao do tratamento farmacologico, dias",
             A[["deliriumtto days"]][com_delirium], 0),
  linha_cont("Dia de inicio do delirium",
             A[["delirium onset"]][com_delirium], 0)
)

T1 <- rbind(T1, T1b)
print(T1, row.names = FALSE)
write.csv(T1, file.path(DIR_OUT, "T1_caracteristicas.csv"), row.names = FALSE)

fim_tto <- A[["delirium onset"]] + A[["deliriumtto days"]] - 1
ativo_d5 <- com_delirium & is.finite(fim_tto) & fim_tto >= 5
cat("\n=== TRATAMENTO FARMACOLOGICO ATIVO NA AVALIACAO DO DIA 5 ===\n")
cat("Criterio: (dia de inicio + duracao farmacologica - 1) >= 5\n")
cat(sprintf("  Sob tratamento farmacologico ativo no Dia 5: %s dos pacientes com delirium\n",
            fmt_contagem(sum(ativo_d5), nD)))
cat(sprintf("  Sobre a amostra analitica completa            : %s\n",
            fmt_contagem(sum(ativo_d5), nA)))
cat("DERIVACAO, nao medida direta: assume inicio do tratamento no dia de\n")
cat("inicio do delirium.\n")
cat("Nao foi feita analise de\n")
cat("sensibilidade adicional: a restricao aos nao expostos deixaria cerca de 18\n")
cat("pacientes e elevaria o |rho| minimo detectavel de 0,48 para 0,63.\n")
if (any(is.na(fim_tto[com_delirium])))
  cat("ATENCAO: ", sum(is.na(fim_tto[com_delirium])),
      " paciente(s) com delirium sem inicio ou duracao registrados.\n", sep = "")

write.csv(data.frame(
  criterio = "inicio + duracao farmacologica - 1 >= 5",
  com_delirium = nD, ativo_d5 = sum(ativo_d5),
  amostra_analitica = nA,
  pressuposto = "inicio do tratamento = dia de inicio do delirium"),
  file.path(DIR_OUT, "T1c_exposicao_farmacologica_D5.csv"), row.names = FALSE)

cat("\n=== ANCORA DESENVOLVIMENTAL DO CAPD ===\n")
anc <- table(A$capd_anchor)
for (k in names(anc)) cat(sprintf("  %-10s %s\n", k, fmt_contagem(anc[[k]], nA)))

T2 <- do.call(rbind, lapply(1:8, function(i) {
  x1 <- A[[ITENS_D1[i]]]; x5 <- A[[ITENS_D5[i]]]
  data.frame(
    Item      = paste0("Q", i),
    Dominio   = ROTULO_ITEM[i],
    D1_med    = fmt_mediana_iqr(x1, 1),
    D1_sd     = round(sd(x1), 2),
    D1_niveis = length(unique(x1)),
    D5_med    = fmt_mediana_iqr(x5, 1),
    D5_sd     = round(sd(x5), 2),
    D5_niveis = length(unique(x5)),
    stringsAsFactors = FALSE)
}))
T2$D1_baixa_variab <- ifelse(T2$D1_niveis >= CFG$item_min_niveis &
                             T2$D1_sd    >= CFG$item_min_sd, "nao", "SIM")
T2$D5_baixa_variab <- ifelse(T2$D5_niveis >= CFG$item_min_niveis &
                             T2$D5_sd    >= CFG$item_min_sd, "nao", "SIM")

cat("\n=== DISTRIBUICAO DOS DOMINIOS DO CAPD (n =", nA, ") ===\n")
print(T2, row.names = FALSE)
write.csv(T2, file.path(DIR_OUT, "T2_distribuicao_itens_CAPD.csv"), row.names = FALSE)

cat("\nItens sinalizados por BAIXA VARIABILIDADE (diagnostico, nao selecao):\n")
reb <- c(paste0(T2$Item[T2$D1_baixa_variab == "SIM"], "-D1"),
         paste0(T2$Item[T2$D5_baixa_variab == "SIM"], "-D5"))
cat("  ", ifelse(length(reb) == 0, "nenhum", paste(reb, collapse = ", ")), "\n")
cat("Estes itens PERMANECEM em todas as familias de analise. A sinalizacao\n")
cat("serve a interpretacao: correlacoes envolvendo itens de baixa\n")
cat("variabilidade tem menor resolucao e mais empates.\n")

T3 <- data.frame(
  Biomarcador = c("NfL Dia 1","NfL Dia 5","S100B Dia 1","S100B Dia 5"),
  Unidade     = "ng/mL",
  Mediana_IQR = c(fmt_mediana_iqr(A[[VAR$nfl_d1]],   3),
                  fmt_mediana_iqr(A[[VAR$nfl_d5]],   3),
                  fmt_mediana_iqr(A[[VAR$s100b_d1]], 3),
                  fmt_mediana_iqr(A[[VAR$s100b_d5]], 3)),
  Min_Max     = sapply(list(A[[VAR$nfl_d1]], A[[VAR$nfl_d5]],
                            A[[VAR$s100b_d1]], A[[VAR$s100b_d5]]),
                       function(x) sprintf("%.3f-%.3f", min(x), max(x))),
  Shapiro_p_bruto = sapply(list(A[[VAR$nfl_d1]], A[[VAR$nfl_d5]],
                                A[[VAR$s100b_d1]], A[[VAR$s100b_d5]]),
                           function(x) fmt_p(shapiro.test(x)$p.value)),
  Shapiro_p_log   = sapply(list(A[[VAR$nfl_d1]], A[[VAR$nfl_d5]],
                                A[[VAR$s100b_d1]], A[[VAR$s100b_d5]]),
                           function(x) fmt_p(shapiro.test(log(x))$p.value)),
  stringsAsFactors = FALSE)

cat("\n=== BIOMARCADORES (n =", nA, ") ===\n")
print(T3, row.names = FALSE)
write.csv(T3, file.path(DIR_OUT, "T3_biomarcadores.csv"), row.names = FALSE)

cat("\nNOTA: se a normalidade nao melhorar com log (compare as\n")
cat("duas ultimas colunas), a transformacao logaritmica se justifica por\n")
cat("coerencia entre analises e pela interpretacao multiplicativa dos\n")
cat("coeficientes nos modelos parametricos — NAO por adequacao distribucional.\n")
cat("As analises primarias sao de postos e independem dessa escolha.\n")

cat("\n=== SOBREPOSICAO DE CONSTRUTOS, POR MOMENTO (n =", nA, ") ===\n")
cat("WAT-1 positivo definido como >=", CFG$wat1_corte,
    "(limiar validado do instrumento).\n\n")

sobrep <- list()
for (mom in c("D1", "D5")) {
  vcapd <- if (mom == "D1") VAR$capd_d1 else VAR$capd_d5
  vwat  <- if (mom == "D1") VAR$wat1_d1 else VAR$wat1_d5
  cp <- A[[vcapd]] >= CFG$capd_corte
  wp <- A[[vwat]]  >= CFG$wat1_corte
  tb <- table(CAPD = ifelse(cp, "CAPD+", "CAPD-"), WAT1 = ifelse(wp, "WAT1+", "WAT1-"))
  cat("--- ", mom, " ---\n", sep = ""); print(addmargins(tb))
  cat(sprintf("  Ambos positivos NO MESMO momento: %s\n\n",
              fmt_contagem(sum(cp & wp), nA)))
  sobrep[[mom]] <- as.data.frame.matrix(addmargins(tb))
}
ambos_d1 <- (A[[VAR$capd_d1]] >= CFG$capd_corte) & (A[[VAR$wat1_d1]] >= CFG$wat1_corte)
ambos_d5 <- (A[[VAR$capd_d5]] >= CFG$capd_corte) & (A[[VAR$wat1_d5]] >= CFG$wat1_corte)
cat(sprintf("Participantes com CAPD+ e WAT-1+ simultaneos em pelo menos um\n"))
cat(sprintf("momento avaliado: %s\n", fmt_contagem(sum(ambos_d1 | ambos_d5), nA)))
write.csv(do.call(rbind, lapply(names(sobrep), function(m)
            cbind(Momento = m, sobrep[[m]]))),
          file.path(DIR_OUT, "T4_sobreposicao_CAPD_WAT1_por_momento.csv"), row.names = TRUE)

cat("\n=== WAT-1: COMPARACAO DE LIMIARES ===\n")
for (mom in c("D1", "D5")) {
  v <- if (mom == "D1") VAR$wat1_d1 else VAR$wat1_d5
  cat(sprintf("  %s  >= 3 (validado): %s   |   > 3 (comparacao historica): %s\n",
              mom, fmt_contagem(sum(A[[v]] >= CFG$wat1_corte), nA),
              fmt_contagem(sum(A[[v]] >= CFG$wat1_hist), nA)))
}
cat("\nNOTA: IWS (sindrome de abstinencia registrada) e variavel DISTINTA da\n")
cat("positividade pelo escore WAT-1. Diagnostico registrado e escore positivo\n")
cat("sao informacoes diferentes e nao devem ser usados de forma intercambiavel.\n")
cat(sprintf("IWS registrada: %s\n",
            fmt_contagem(sum(A[["Iatrogenic withdrawal syndrome (IWS)"]] == 1, na.rm = TRUE), nA)))
cat("\nNas correlacoes e nos ajustes das sensibilidades, WAT-1 e usado como\n")
cat("variavel NUMERICA, sem dicotomizacao. Os cortes acima servem apenas a\n")
cat("descricao da coorte.\n")

cat("\n=== SCRIPT 02 CONCLUIDO ===\n")
cat("Proximo passo: 03_precisao_MDE.R\n\n")
