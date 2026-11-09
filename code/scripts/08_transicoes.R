# 08_transicoes.R
# Aggregate transitions between days 1 and 5; correlations between item changes and biomarker changes
# (Benjamini-Hochberg).
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
library(boot); library(ggplot2)

D   <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
VAR <- attr(D, "VAR"); ITENS_D1 <- attr(D, "ITENS_D1")
ITENS_D5 <- attr(D, "ITENS_D5"); ROTULO_ITEM <- attr(D, "ROTULO_ITEM")
A  <- D[D$avaliavel_capd, ]
nA <- nrow(A)

cat("\n=== MATRIZ DE TRANSIÇÃO D1 -> D5 ===\n")
tab <- table(D1 = D$fenotipo_D1, D5 = D$fenotipo_D5)
print(addmargins(tab))
write.csv(as.data.frame.matrix(addmargins(tab)),
          file.path(DIR_OUT, "T13_matriz_transicao.csv"))

A$transicao <- ifelse(A$fenotipo_D1 == A$fenotipo_D5, "estavel",
               ifelse(A$fenotipo_D1 == "G1_free" & A$fenotipo_D5 == "G2_delirium",
                      "incidente", "resolucao"))
cat("\nEntre os", nA, "avaliáveis:\n")
for (k in names(table(A$transicao)))
  cat(sprintf("  %-12s %d\n", k, sum(A$transicao == k)))

TR <- A[A$transicao != "estavel", ]
cat(sprintf("\nParticipantes com mudança de fenótipo: %d\n", nrow(TR)))

TR <- TR[order(TR$transicao, TR$ID), ]
TR$rotulo <- LETTERS[seq_len(nrow(TR))]

write.csv(data.frame(rotulo = TR$rotulo, ID = TR$ID),
          file.path(DIR_OUT, "INTERNO_chave_rotulos_trajetorias.csv"), row.names = FALSE)
cat("Rótulos atribuídos:", paste(TR$rotulo, collapse = ", "), "\n")
cat("Chave rótulo->ID gravada em INTERNO_chave_rotulos_trajetorias.csv —\n")
cat("arquivo de uso local, contem identificadores; nao publicar.\n")

traj <- do.call(rbind, lapply(seq_len(nrow(TR)), function(r) {
  p <- TR[r, ]
  base <- data.frame(
    Participante = p$rotulo, transicao = p$transicao,
    fenotipo_D1 = p$fenotipo_D1, fenotipo_D5 = p$fenotipo_D5,
    stringsAsFactors = FALSE)
  itens <- setNames(
    as.data.frame(t(c(
      sapply(1:8, function(i) p[[ITENS_D1[i]]]),
      sapply(1:8, function(i) p[[ITENS_D5[i]]]),
      sapply(1:8, function(i) p[[ITENS_D5[i]]] - p[[ITENS_D1[i]]])))),
    c(paste0("Q", 1:8, "_D1"), paste0("Q", 1:8, "_D5"), paste0("Q", 1:8, "_delta")))
  bio <- data.frame(
    CAPD_D1 = p[[VAR$capd_d1]], CAPD_D5 = p[[VAR$capd_d5]],
    CAPD_delta = p[[VAR$capd_d5]] - p[[VAR$capd_d1]],
    NfL_D1 = p[[VAR$nfl_d1]], NfL_D5 = p[[VAR$nfl_d5]],
    NfL_delta_log = p$delta_log_nfl,
    S100B_D1 = p[[VAR$s100b_d1]], S100B_D5 = p[[VAR$s100b_d5]],
    S100B_delta_log = p$delta_log_s100b,
    FLACC_D1 = p[[VAR$dor_d1]], FLACC_D5 = p[[VAR$dor_d5]],
    WAT1_D1 = p[[VAR$wat1_d1]], WAT1_D5 = p[[VAR$wat1_d5]])
  cbind(base, bio, itens)
}))

cat("\n=== TRAJETÓRIAS INDIVIDUAIS ===\n")
print(traj[, c("Participante","transicao","CAPD_D1","CAPD_D5","CAPD_delta",
               "NfL_D1","NfL_D5","S100B_D1","S100B_D5",
               "FLACC_D1","FLACC_D5","WAT1_D1","WAT1_D5")], row.names = FALSE)
write.csv(traj, file.path(DIR_OUT, "T14_trajetorias_individuais.csv"), row.names = FALSE)

cat("\nATENÇÃO: nenhuma estatística de concordância será calculada sobre estas\n")
cat("linhas. Descreva-as individualmente, nomeando o número de casos.\n")

long_it <- do.call(rbind, lapply(seq_len(nrow(TR)), function(r) {
  p <- TR[r, ]
  data.frame(ID = factor(paste0(p$rotulo, " (", p$transicao, ")")),
             Item = factor(paste0("Q", 1:8), levels = paste0("Q", 1:8)),
             delta = sapply(1:8, function(i) p[[ITENS_D5[i]]] - p[[ITENS_D1[i]]]),
             stringsAsFactors = FALSE)
}))
long_bio <- do.call(rbind, lapply(seq_len(nrow(TR)), function(r) {
  p <- TR[r, ]
  data.frame(ID = factor(paste0(p$rotulo, " (", p$transicao, ")")),
             Marcador = factor(c("NfL", "S100B"), levels = c("NfL", "S100B")),
             delta = c(p$delta_log_nfl, p$delta_log_s100b),
             stringsAsFactors = FALSE)
}))

g1 <- ggplot(long_it, aes(Item, delta)) +
  geom_col(width = 0.65, fill = "grey35") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ ID, ncol = 2) +
  labs(x = NULL, y = "Variação do escore do domínio (D5 - D1)",
       title = "Domínios do CAPD") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

g2 <- ggplot(long_bio, aes(Marcador, delta)) +
  geom_col(width = 0.5, fill = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ ID, ncol = 2) +
  labs(x = NULL, y = "log(D5) - log(D1)", title = "Biomarcadores") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- g1 + g2 + patchwork::plot_layout(widths = c(2.4, 1)) +
    patchwork::plot_annotation(
      caption = paste0("Trajetórias individuais dos ", nrow(TR),
                       " participantes com mudança de fenótipo entre o Dia 1 e o Dia 5. ",
                       "Apresentação descritiva; nenhuma inferência estatística."))
  ggsave(file.path(DIR_FIG, "FIG_trajetorias_transicao.png"), fig,
         width = 11, height = 8, dpi = 300)
  ggsave(file.path(DIR_FIG, "FIG_trajetorias_transicao.tiff"), fig,
         width = 11, height = 8, dpi = 300, compression = "lzw")
  cat("\nFigura salva em", DIR_FIG, "\n")
} else {
  ggsave(file.path(DIR_FIG, "FIG_trajetorias_dominios.png"), g1, width = 8, height = 7, dpi = 300)
  ggsave(file.path(DIR_FIG, "FIG_trajetorias_biomarcadores.png"), g2, width = 5, height = 7, dpi = 300)
  cat("\npatchwork ausente: figuras salvas separadamente.\n")
}

cat("\n\n=== H4b — VARIAÇÃO POR DOMÍNIO x VARIAÇÃO DO BIOMARCADOR ===\n")
cat("Camada exploratória, n =", nA, ", BH q <", CFG$limiar_q_S2, "\n")
cat("Bootstrap BCa,", CFG$n_boot, "reamostragens.\n\n")

cor_bca <- function(x, y, R = CFG$n_boot, seed = CFG$seed) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
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
  list(n = length(x), rho = unname(est$estimate), p = est$p.value,
       ic_inf = ic[1], ic_sup = ic[2], tipo_ic = tipo)
}

linhas <- list(); k <- 0
for (bio in c("NfL", "S100B")) {
  vdelta <- if (bio == "NfL") "delta_log_nfl" else "delta_log_s100b"
  for (i in 1:8) {
    k <- k + 1
    cat(sprintf("  [%2d/16] delta Q%d x delta %s\n", k, i, bio))
    r <- cor_bca(A[[paste0("delta_Q", i)]], A[[vdelta]])
    linhas[[k]] <- data.frame(
      Item = paste0("Q", i), Dominio = ROTULO_ITEM[i], Biomarcador = bio,
      n = r$n, rho = r$rho, IC_inf = r$ic_inf, IC_sup = r$ic_sup,
      Tipo_IC = r$tipo_ic, p_bruto = r$p, stringsAsFactors = FALSE)
  }
}
H4b <- do.call(rbind, linhas)
H4b$q_BH <- p.adjust(H4b$p_bruto, method = CFG$ajuste_S2, n = 16)
H4b$Sinalizado <- ifelse(H4b$q_BH < CFG$limiar_q_S2, "sim", "nao")
H4b <- H4b[order(H4b$q_BH), ]

write.csv(H4b, file.path(DIR_OUT, "T15_H4b_numerico.csv"), row.names=FALSE)
saveRDS(H4b, file.path(DIR_OUT, "T15_H4b_numerico.rds"))
saida <- data.frame(
  Item = H4b$Item, Dominio = H4b$Dominio, Biomarcador = H4b$Biomarcador, n = H4b$n,
  rho = sprintf("%.3f", H4b$rho),
  IC95 = sprintf("%.3f a %.3f", H4b$IC_inf, H4b$IC_sup),
  p_bruto = fmt_p(H4b$p_bruto), q_BH = fmt_p(H4b$q_BH),
  Sinalizado = H4b$Sinalizado, stringsAsFactors = FALSE)

cat("\n=== RESULTADOS H4b ===\n"); print(saida, row.names = FALSE)
write.csv(saida, file.path(DIR_OUT, "T15_H4b_deltas_por_dominio.csv"), row.names = FALSE)

cat("\n=== TRANSIÇÕES: LEITURA ===\n")
cat("Seis casos; apresentação apenas descritiva:\n\n")
cat("  Cinco participantes passaram de sem delirium a delirium\n")
cat("  entre o Dia 1 e o Dia 5, e em um o delirium\n")
cat("  se resolveu. Nenhuma inferência estatística é feita\n")
cat("  sobre estes seis casos.\n\n")
cat("Não sustentado por seis casos:\n")
cat("  - \"in patients who developed delirium, NfL increased\"\n")
cat("  - \"biomarker changes tracked domain changes\"\n")
cat("  - qualquer proporção, percentual ou valor-p sobre estes seis\n\n")
cat("Manter a linguagem descritiva e numerada.\n")
cat("\n")

cat("\n=== SCRIPT 08 CONCLUÍDO ===\n")
cat("Proximo passo: 09_H4_sensibilidades.R\n\n")
