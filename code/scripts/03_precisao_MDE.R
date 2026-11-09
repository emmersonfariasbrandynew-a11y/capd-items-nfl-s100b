# 03_precisao_MDE.R
# Precision: minimum detectable Spearman correlation for the available sample sizes.
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
D  <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
nA <- sum(D$avaliavel_capd)

mde_spearman <- function(n, alfa = 0.05, poder = 0.80, k = 1.03) {
  z <- (qnorm(1 - alfa/2) + qnorm(poder)) * k / sqrt(n - 3)
  tanh(z)
}

grade <- expand.grid(
  n     = c(nA, 45),
  alfa  = c(0.05, 0.0125, 0.10),
  poder = c(0.80, 0.90)
)
grade$MDE_rho <- round(mapply(mde_spearman, grade$n, grade$alfa, grade$poder), 3)
grade <- grade[order(grade$n, grade$poder, -grade$alfa), ]

cat("\n=== EFEITO MINIMO DETECTAVEL — correlacao de Spearman ===\n")
cat("Amostra analitica: n =", nA, "\n\n")
print(grade, row.names = FALSE)
write.csv(grade, file.path(DIR_OUT, "T5_precisao_MDE.csv"), row.names = FALSE)

m80_05   <- mde_spearman(nA, 0.05,   0.80)
m80_holm <- mde_spearman(nA, 0.0125, 0.80)

cat("\n=== EFEITO MINIMO DETECTAVEL: RESUMO ===\n")
cat(sprintf("  n = %d | poder 80%% | alfa 0.05 bilateral  -> |rho| >= %.2f\n", nA, m80_05))
cat(sprintf("  n = %d | poder 80%% | alfa 0.0125 (1o passo Holm)-> |rho| >= %.2f\n", nA, m80_holm))

cat("\n--- RESUMO EM INGLES ---\n")
cat(sprintf(
"No formal a priori sample size calculation was performed for this secondary
analysis, which was constrained by the size of the parent cohort. Rather than
reporting uninformative observed power, we quantified the precision boundaries
of the design. With %d CAPD-evaluable participants, the study had 80%% power to
detect a Spearman correlation of |rho| >= %.2f at alpha = 0.05 (two-sided), and
|rho| >= %.2f under the most conservative step of the Holm procedure (alpha/4 for the
smallest of four p-values; subsequent steps face alpha/3, alpha/2 and alpha). These thresholds apply to a single correlation and do not transfer to the
difference between two dependent correlations. Associations weaker than these
magnitudes could not be reliably excluded; accordingly, null findings are
reported as 'we did not detect an association' rather than as evidence of
absence.\n", nA, m80_05, m80_holm))

cat("\n--- ADVERTENCIA SOBRE Delta-rho ---\n")
cat("Os valores acima aplicam-se a UMA correlacao. NAO se aplicam aos\n")
cat("estimandos de diferenca de H2 e H5: Var(rho1-rho2) = Var+Var-2Cov, e a\n")
cat("covariancia entre correlacoes que compartilham variavel costuma ser\n")
cat("positiva e grande. A precisao de cada contraste e o seu EP de jackknife,\n")
cat("na coluna SE_jack das tabelas T9 e T12.\n")

cat("\nConstante de correcao para postos: k = 1.03 (Bonett & Wright 2000).\n")
cat("Sem ela o valor seria o de Pearson, inadequado a correlacoes de postos.\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ns <- 10:60
  cur <- data.frame(
    n    = rep(ns, 2),
    alfa = factor(rep(c("0.05", "0.0125 (1o passo de Holm)"), each = length(ns)),
                  levels = c("0.05", "0.0125 (1o passo de Holm)")),
    mde  = c(mde_spearman(ns, 0.05, 0.80), mde_spearman(ns, 0.0125, 0.80)))
  g <- ggplot2::ggplot(cur, ggplot2::aes(n, mde, linetype = alfa)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = nA, colour = "grey40", linetype = "dotted") +
    ggplot2::annotate("text", x = nA + 1, y = 0.95, hjust = 0, size = 3.2,
                      label = paste0("n = ", nA, " (amostra analitica)")) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Tamanho amostral",
                  y = expression(paste("|", rho, "| minimo detectavel (poder 80%)")),
                  linetype = "Limiar de significancia") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(DIR_FIG, "SDC_curva_precisao.png"), g,
                  width = 6, height = 4.2, dpi = 300)
  cat("\nFigura suplementar salva em:", file.path(DIR_FIG, "SDC_curva_precisao.png"), "\n")
}

cat("\n=== SCRIPT 03 CONCLUIDO ===\n")
cat("Proximo passo: 04_S1_primaria.R\n\n")
