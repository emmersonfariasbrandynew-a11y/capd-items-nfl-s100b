# 00_setup.R
# Configuration: required packages, fixed analysis parameters (seed, resamples, thresholds), dataset
# path, and structural check of the dataset.
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

pkgs <- c(
  "boot",
  "digest",
  "ggplot2",
  "patchwork"
)

faltando_pkg <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(faltando_pkg) > 0) {
  stop("Pacotes ausentes: ", paste(faltando_pkg, collapse = ", "),
       "\nInstale-os e reexecute. Nao ha instalacao automatica, para manter o\n",
       "ambiente declarado em docs/sessionInfo.txt.\n",
       "Comando sugerido:  install.packages(c(",
       paste(sprintf('"%s"', faltando_pkg), collapse = ", "), "))")
}
invisible(lapply(pkgs, library, character.only = TRUE))

VERSOES_PKG <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), "")
cat("Pacotes carregados:\n")
for (p in pkgs) cat(sprintf("  %-10s %s\n", p, VERSOES_PKG[[p]]))

SEED_GLOBAL <- 20260818
set.seed(SEED_GLOBAL)

CFG <- list(
  seed         = SEED_GLOBAL,
  n_boot       = 5000,
  boot_type    = "bca",
  metodo_cor   = "spearman",
  ajuste_S1    = "holm",
  ajuste_S2    = "BH",
  limiar_q_S2  = 0.10,
  capd_corte   = 9,
  rass_coma    = -4,
  alternativa  = "two.sided",
  wat1_corte   = 3,
  wat1_hist    = 3.0001,
  item_min_niveis = 4,
  item_min_sd     = 0.75
)

DIR_OUT <- "output"
DIR_FIG <- file.path(DIR_OUT, "figuras")
for (d in c(DIR_OUT, DIR_FIG)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

CAMINHO_DADOS <- file.path("..", "data_restricted", "capd_items_biomarkers_anonymized.csv")
ABA_DADOS     <- "dados"

ASSINATURA_BANCO <- list(
  rotulo = "anonymized repository dataset (2026-10-14)",
  linhas = 45L, colunas = 50L,
  tipodelierium_na = 0L,
  tipodelierium_tab = c(`0` = 25L, `1` = 8L, `2` = 7L, `3` = 5L),
  deliriumtto_na = 25L,
  antipsicoticos = c("haldol", "quetiapine", "respiridona", "clorpromazine"),
  antipsicotico_na_delirium = 0L
)

conferir_banco <- function(caminho = CAMINHO_DADOS, aba = ABA_DADOS,
                           esp = ASSINATURA_BANCO, parar = TRUE) {
  d <- read.csv(caminho, check.names = FALSE, na.strings = "", stringsAsFactors = FALSE)
  falhas <- character(0)
  if (nrow(d) != esp$linhas || ncol(d) != esp$colunas)
    falhas <- c(falhas, sprintf("dimensoes %d x %d, esperado %d x %d",
                                nrow(d), ncol(d), esp$linhas, esp$colunas))
  if (!"tipodelierium" %in% names(d)) {
    falhas <- c(falhas, "coluna tipodelierium ausente")
  } else {
    na_t <- sum(is.na(d$tipodelierium))
    if (na_t != esp$tipodelierium_na)
      falhas <- c(falhas, sprintf(
        "tipodelierium com %d celulas vazias, esperado %d",
        na_t, esp$tipodelierium_na))
    tab <- table(factor(d$tipodelierium, levels = names(esp$tipodelierium_tab)))
    if (!identical(as.integer(tab), as.integer(esp$tipodelierium_tab)))
      falhas <- c(falhas, sprintf("distribuicao de tipodelierium %s, esperado %s",
                                  paste(as.integer(tab), collapse = "/"),
                                  paste(as.integer(esp$tipodelierium_tab), collapse = "/")))
  }
  if ("deliriumtto days" %in% names(d)) {
    na_d <- sum(is.na(d[["deliriumtto days"]]))
    if (na_d != esp$deliriumtto_na)
      falhas <- c(falhas, sprintf("deliriumtto days com %d vazios, esperado %d",
                                  na_d, esp$deliriumtto_na))
  } else falhas <- c(falhas, "coluna 'deliriumtto days' ausente")

  if (all(c(esp$antipsicoticos, "delirium") %in% names(d))) {
    dl <- d[d$delirium == 1, esp$antipsicoticos, drop = FALSE]
    na_a <- sum(is.na(dl))
    if (na_a != esp$antipsicotico_na_delirium)
      falhas <- c(falhas, sprintf(
        "%d celula(s) vazia(s) nos antipsicoticos entre os pacientes com delirium, esperado %d",
        na_a, esp$antipsicotico_na_delirium))
  }

  if (length(falhas)) {
    msg <- paste0(paste0("BANCO NAO CONFERE COM A ", esp$rotulo, ":\n  - "),
                  paste(falhas, collapse = "\n  - "),
                  "\nNAO prossiga: verifique o arquivo de dados.")
    if (parar) stop(msg) else warning(msg)
  } else {
    cat("Banco conferido: dimensoes e valores-chave verificaveis correspondem a", esp$rotulo, "\n")
  }
  invisible(length(falhas) == 0)
}

if (!file.exists(CAMINHO_DADOS)) {
  message("Arquivo nao encontrado no diretorio de trabalho.")
  message("Coloque o arquivo de dados em data_restricted/ ou ajuste CAMINHO_DADOS.")
  message("O arquivo e de acesso restrito; veja o README.")
}

fmt_contagem <- function(x, n) {
  sprintf("%d/%d", as.integer(x), as.integer(n))
}

fmt_mediana_iqr <- function(x, dig = 2) {
  x <- x[!is.na(x)]
  sprintf("%.*f [%.*f-%.*f]", dig, median(x), dig, quantile(x, .25),
          dig, quantile(x, .75))
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

cat("\n=========================================================\n")
cat("AMBIENTE CONFIGURADO — CAPD ITEMS, NfL AND S100B\n")
cat("=========================================================\n")
cat("R           :", R.version.string, "\n")
cat("Semente     :", CFG$seed, "\n")
cat("Bootstrap   :", CFG$boot_type, "(correlacoes isoladas; Delta-rho por bootstrap-t) com", CFG$n_boot, "iteracoes\n")
cat("Correlacao  :", CFG$metodo_cor, "\n")
cat("Ajuste total:", CFG$ajuste_S1, "| Ajuste itens:", CFG$ajuste_S2,
    "(q <", CFG$limiar_q_S2, ")\n")
cat("Alternativa :", CFG$alternativa, "(todos os testes)\n")
cat("Estatuto    : TODAS as hipoteses exploratorias; H2 e a principal\n")
cat("H5          : 64 contrastes, BH em familia unica\n")
cat("WAT-1       : positivo >= 3 (limiar validado); > 3 apenas descritivo\n")
cat("AVISO       : ajuste por familia NAO e controle global dos 102 testes\n")
cat("Saida em    :", normalizePath(DIR_OUT), "\n")
cat("Pacotes     :", paste(sprintf("%s %s", names(VERSOES_PKG), VERSOES_PKG), collapse = " | "), "\n")
cat("Dados       :", ifelse(file.exists(CAMINHO_DADOS),
                            "encontrados", "AUSENTES"), "\n")
cat("=========================================================\n\n")
cat("Proximo passo: 01_data_prep.R\n\n")
