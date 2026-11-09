# 06_H2H3_diferenca_rho.R
# Main comparison: item-level differences between NfL and S100B correlations (Benjamini-Hochberg);
# temporal differences between day-1 and day-5 correlations.
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
if (!exists("teste_delta_rho"))
  source(file.path("scripts", "00b_motor_deltarho.R"))

D   <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
VAR <- attr(D, "VAR"); ITENS_D1 <- attr(D, "ITENS_D1")
ITENS_D5 <- attr(D, "ITENS_D5"); ROTULO_ITEM <- attr(D, "ROTULO_ITEM")
A  <- D[D$avaliavel_capd, ]
nA <- nrow(A)

cat("\n=== H2 — DIFERENCA ENTRE NfL E S100B POR DOMINIO (n =", nA, ") ===\n")
cat("Estimando: Delta_rho = rho(Qi, NfL) - rho(Qi, S100B), com sinal\n")
cat("Bootstrap-t estudentizado,", CFG$n_boot, "replicas por contraste.\n")
cat("Familia de 16 contrastes; BH q <", CFG$limiar_q_S2, "\n\n")

linhas <- list(); k <- 0
for (mom in c("D1", "D5")) {
  itens <- if (mom == "D1") ITENS_D1 else ITENS_D5
  vnfl  <- if (mom == "D1") VAR$nfl_d1   else VAR$nfl_d5
  vs100 <- if (mom == "D1") VAR$s100b_d1 else VAR$s100b_d5
  for (i in 1:8) {
    k <- k + 1
    cat(sprintf("  [%2d/16] Q%d %s ... ", k, i, mom)); flush.console()
    r  <- teste_delta_rho(A[[itens[i]]], A[[vnfl]], A[[vs100]], sobreposto = TRUE)
    linhas[[k]] <- linha_delta(r, Item = paste0("Q", i), Dominio = ROTULO_ITEM[i],
                               Momento = mom)
    cat(r$estado, "\n")
  }
}
H2 <- do.call(rbind, linhas)
names(H2)[names(H2) == "rho1"] <- "rho_NfL"
names(H2)[names(H2) == "rho2"] <- "rho_S100B"

N_FAMILIA_H2 <- 16
H2$q_BH <- p.adjust(H2$p_bootstrap_t, method = CFG$ajuste_S2, n = N_FAMILIA_H2)
H2$IC_exclui_zero <- ifelse(is.na(H2$IC_inf), NA_character_,
                     ifelse(H2$IC_inf > 0 | H2$IC_sup < 0, "sim", "nao"))
H2 <- H2[order(H2$q_BH, na.last = TRUE), ]

write.csv(H2, file.path(DIR_OUT, "T9_H2_numerico.csv"), row.names=FALSE)
saveRDS(H2, file.path(DIR_OUT, "T9_H2_numerico.rds"))
saida_H2 <- data.frame(
  Item = H2$Item, Dominio = H2$Dominio, Momento = H2$Momento, n = H2$n,
  rho_NfL   = sprintf("%.3f", H2$rho_NfL),
  rho_S100B = sprintf("%.3f", H2$rho_S100B),
  Delta_rho = sprintf("%.3f", H2$Delta_rho),
  SE_jack   = sprintf("%.3f", H2$SE_jack),
  IC95_bootstrap_t = sprintf("%.3f a %.3f", H2$IC_inf, H2$IC_sup),
  IC_exclui_zero = H2$IC_exclui_zero,
  p_bootstrap_t = fmt_p(H2$p_bootstrap_t), q_BH = fmt_p(H2$q_BH),
  IC95_BCa_diagnostico = sprintf("%.3f a %.3f", H2$BCa_inf, H2$BCa_sup),
  R_valido = H2$R_valido, estado = H2$estado, stringsAsFactors = FALSE)

cat("\n=== RESULTADOS H2 ===\n"); print(saida_H2, row.names = FALSE)
write.csv(saida_H2, file.path(DIR_OUT, "T9_H2_diferenca_NfL_vs_S100B.csv"), row.names = FALSE)

n_falha <- sum(H2$estado != "ok")
if (n_falha > 0) {
  cat("\nATENCAO:", n_falha, "de 16 contrastes sem valor-p. Estados observados:\n")
  print(table(H2$estado[H2$estado != "ok"]))
  cat("Estes contrastes constam da tabela com estimativa pontual e estado, e\n")
  cat("permanecem no denominador do ajuste. Devem ser relatados.\n")
}

cat("\n\n=== H3 — CONSISTENCIA TEMPORAL (EXPLORATORIO, n =", nA, ") ===\n")
cat("Estimando: Delta_rho = rho(Qi_D1, bio_D1) - rho(Qi_D5, bio_D5)\n")
cat("16 contrastes: 8 dominios x 2 biomarcadores.\n")
cat("Contraste NAO SOBREPOSTO: fora da cobertura da verificacao sintetica (Material Suplementar, Apendice 2).\n")
cat("SEM valor-p, SEM ajuste e SEM declaracao de equivalencia.\n\n")

delta_temporal <- function(x1, y1, x5, y5, R = CFG$n_boot, seed = CFG$seed) {
  ok <- !is.na(x1) & !is.na(y1) & !is.na(x5) & !is.na(y5)
  x1 <- x1[ok]; y1 <- y1[ok]; x5 <- x5[ok]; y5 <- y5[ok]; n <- length(x1)
  dstat <- function(i) sp(x1[i], y1[i]) - sp(x5[i], y5[i])
  D <- dstat(seq_len(n))
  dj <- vapply(seq_len(n), function(i) dstat(seq_len(n)[-i]), numeric(1))
  SE <- if (anyNA(dj)) NA_real_ else sqrt((n - 1) / n * sum((dj - mean(dj))^2))
  set.seed(seed)
  Ts <- rep(NA_real_, R)
  for (b in seq_len(R)) {
    i  <- sample.int(n, n, replace = TRUE)
    db <- dstat(i)
    dbj <- vapply(seq_len(n), function(j) dstat(i[-j]), numeric(1))
    sb <- if (anyNA(dbj)) NA_real_ else sqrt((n - 1) / n * sum((dbj - mean(dbj))^2))
    if (is.finite(db) && is.finite(sb) && sb > 0) Ts[b] <- (db - D) / sb
  }
  val <- is.finite(Ts)
  ic <- c(NA_real_, NA_real_)
  if (sum(val) >= 0.95 * R && is.finite(SE) && SE > 0) {
    qt <- stats::quantile(Ts[val], c(0.975, 0.025), names = FALSE, type = 7)
    ic <- c(D - qt[1] * SE, D - qt[2] * SE)
  }
  list(n = n, delta = D, rD1 = sp(x1, y1), rD5 = sp(x5, y5),
       se = SE, ic_inf = ic[1], ic_sup = ic[2], R_valido = sum(val))
}

linhas <- list(); k <- 0
for (bio in c("NfL", "S100B")) {
  v1 <- if (bio == "NfL") VAR$nfl_d1 else VAR$s100b_d1
  v5 <- if (bio == "NfL") VAR$nfl_d5 else VAR$s100b_d5
  for (i in 1:8) {
    k <- k + 1
    cat(sprintf("  [%2d/16] Q%d %s\n", k, i, bio)); flush.console()
    r <- delta_temporal(A[[ITENS_D1[i]]], A[[v1]], A[[ITENS_D5[i]]], A[[v5]])
    linhas[[k]] <- data.frame(
      Item = paste0("Q", i), Dominio = ROTULO_ITEM[i], Biomarcador = bio, n = r$n,
      rho_D1 = r$rD1, rho_D5 = r$rD5, Delta_rho = r$delta, SE_jack = r$se,
      IC_inf = r$ic_inf, IC_sup = r$ic_sup, R_valido = r$R_valido,
      stringsAsFactors = FALSE)
  }
}
H3 <- do.call(rbind, linhas)

write.csv(H3, file.path(DIR_OUT, "T10_H3_numerico.csv"), row.names=FALSE)
saveRDS(H3, file.path(DIR_OUT, "T10_H3_numerico.rds"))
saida_H3 <- data.frame(
  Item = H3$Item, Dominio = H3$Dominio, Biomarcador = H3$Biomarcador, n = H3$n,
  rho_D1 = sprintf("%.3f", H3$rho_D1), rho_D5 = sprintf("%.3f", H3$rho_D5),
  Delta_rho = sprintf("%.3f", H3$Delta_rho),
  IC95_bootstrap_t = sprintf("%.3f a %.3f", H3$IC_inf, H3$IC_sup),
  R_valido = H3$R_valido, stringsAsFactors = FALSE)

cat("\n=== RESULTADOS H3 ===\n"); print(saida_H3, row.names = FALSE)
write.csv(saida_H3, file.path(DIR_OUT, "T10_H3_consistencia_D1_vs_D5.csv"), row.names = FALSE)

cat("\n=== COMO LER =================================================\n")
cat("1. O estimando e Delta_rho COM SINAL. Negativo = associacao mais forte\n")
cat("   com S100B; positivo = mais forte com NfL. Nao e diferenca de modulos.\n\n")
cat("2. Os biomarcadores NAO diferem num item apenas porque\n")
cat("   rho_NfL foi significativo e rho_S100B nao foi. Essa inferencia e\n")
cat("   invalida (Gelman & Stern 2006).\n\n")
cat("3. O IC principal e o bootstrap-t, do MESMO pivo que produziu o valor-p.\n")
cat("   A coluna BCa e DIAGNOSTICO: a verificacao sintetica avaliou somente o\n")
cat("   teste, nao a cobertura de intervalos. Discordancia entre os dois deve\n")
cat("   ser relatada, nunca usada para escolher o mais favoravel.\n\n")
cat("4. H3 e EXPLORATORIA e NAO SOBREPOSTA: sem valor-p, sem ajuste, sem\n")
cat("   equivalencia. A verificacao sintetica nao cobre este caso, e isso\n")
cat("   constitui limitacao.\n\n")

cat("\n=== SCRIPT 06 CONCLUIDO ===\n")
cat("Proximo passo: 07_H5_comparadores.R\n\n")
