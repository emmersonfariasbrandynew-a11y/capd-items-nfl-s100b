# 01_data_prep.R
# Data input, SHA-256 integrity check, structural validation, neurological classification on days 1
# and 5, analytic sample, and day 5 - day 1 changes.
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))
log_con <- file(file.path(DIR_OUT, "01_log_validacao.txt"), open = "wt")
registrar <- function(...) { txt <- sprintf(...); cat(txt, "\n"); writeLines(txt, log_con) }

if (exists("conferir_banco")) conferir_banco()

D <- read.csv(CAMINHO_DADOS, check.names = FALSE, na.strings = "", stringsAsFactors = FALSE)
names(D)[names(D) == "study_id"] <- "ID"

D <- D[, !(names(D) == "" | is.na(names(D))), drop = FALSE]
dups <- names(D)[duplicated(names(D))]
if (length(dups)) registrar("Colunas duplicadas removidas: %s", paste(unique(dups), collapse = ", "))
D <- D[, !duplicated(names(D)), drop = FALSE]
D <- D[!is.na(D$ID), ]
D$ID <- as.integer(D$ID)

num_cols <- setdiff(names(D), c("capd_anchor", "flacc_flaccR", "PCCC_categoria"))
for (cn in num_cols) if (!is.numeric(D[[cn]])) {
  conv <- suppressWarnings(as.numeric(D[[cn]]))
  if (sum(is.na(conv)) == sum(is.na(D[[cn]]))) D[[cn]] <- conv
  else registrar("ATENCAO: coluna '%s' nao e numerica e nao pode ser convertida.", cn)
}

HASH_BASELINE <- "763e2673267b58128d0de888f3be893170b158847b4d9fac376b1fe7ceb6e830"

if (!requireNamespace("digest", quietly = TRUE))
  stop("Pacote 'digest' ausente: a integridade do banco NAO pode ser conferida.\n",
       "Instale-o e reexecute. Execucao interrompida.")
hash_arquivo <- digest::digest(file = CAMINHO_DADOS, algo = "sha256")

registrar("=== 1.0 INTEGRIDADE DO ARQUIVO ===")
registrar("Arquivo : %s", CAMINHO_DADOS)
registrar("SHA-256 : %s", hash_arquivo)
if (!identical(hash_arquivo, HASH_BASELINE))
  stop("O arquivo de dados difere do SHA-256 esperado.\n",
       "  esperado: ", HASH_BASELINE, "\n",
       "  obtido  : ", hash_arquivo, "\n",
       "Execucao interrompida. Verifique o arquivo de dados.")
registrar("Integridade do arquivo conferida (SHA-256 esperado).")

registrar("=== 1.1 LEITURA ===")
registrar("Linhas lidas: %d | Colunas: %d", nrow(D), ncol(D))
stopifnot(nrow(D) == 45)

VAR <- list(
  nfl_d1   = "NFLD1",     nfl_d5   = "NFLD5",
  s100b_d1 = "SB100D1",   s100b_d5 = "SB100D5",
  capd_d1  = "CAPDTOTALD1", capd_d5 = "CAPDTOTALD5",
  rass_d1  = "rassd1",    rass_d5  = "rassd5",
  wat1_d1  = "wat-1 d1",  wat1_d5  = "wat-1 d5",
  dor_d1   = "dor_d1",    dor_d5   = "dor_d5"
)
ITENS_D1 <- paste0("capd_Q", 1:8, " d1")
ITENS_D5 <- paste0("capd_Q", 1:8, " d5")

faltando <- setdiff(c(unlist(VAR), ITENS_D1, ITENS_D5), names(D))
if (length(faltando)) stop("Colunas ausentes no banco: ", paste(faltando, collapse = ", "))
registrar("Todas as colunas esperadas estao presentes.")

ROTULO_ITEM <- c(
  Q1 = "Contato visual com o cuidador",
  Q2 = "Movimentos com proposito",
  Q3 = "Percepcao do ambiente",
  Q4 = "Comunicacao de necessidades",
  Q5 = "Inquietacao",
  Q6 = "Inconsolabilidade",
  Q7 = "Hipoatividade / interacao reduzida",
  Q8 = "Latencia de resposta"
)

registrar("\n=== 1.3 VALIDACAO ESTRUTURAL ===")

chk_soma <- function(itens, total, rotulo) {
  m <- rowSums(!is.na(D[, itens])) == 8 & !is.na(D[[total]])
  dif <- rowSums(D[m, itens]) - D[[total]][m]
  n_dif <- sum(dif != 0)
  registrar("(a) %s: n avaliavel = %d | discrepancias item x total = %d",
            rotulo, sum(m), n_dif)
  if (n_dif > 0) {
    print(data.frame(ID = D$ID[m][dif != 0],
                     total = D[[total]][m][dif != 0],
                     soma  = rowSums(D[m, itens])[dif != 0]))
    stop("Banco inconsistente: a soma dos itens difere do escore total.")
  }
}
chk_soma(ITENS_D1, VAR$capd_d1, "D1")
chk_soma(ITENS_D5, VAR$capd_d5, "D5")

stopifnot(all(unlist(D[, c(ITENS_D1, ITENS_D5)]) %in% c(NA, 0:4)))
stopifnot(all(D[[VAR$rass_d1]] >= -5 & D[[VAR$rass_d1]] <= 4, na.rm = TRUE))
stopifnot(all(D[[VAR$rass_d5]] >= -5 & D[[VAR$rass_d5]] <= 4, na.rm = TRUE))
registrar("(b) Itens do CAPD em 0-4 e RASS em -5..+4: OK")

for (v in c(VAR$dor_d1, VAR$dor_d5)) {
  x <- D[[v]][!is.na(D[[v]])]
  stopifnot(all(x >= 0 & x <= 10))
  registrar("(c) %s: n = %d | amplitude observada %g-%g", v, length(x), min(x), max(x))
}

classificar <- function(rass, capd) {
  ifelse(is.na(rass), NA_character_,
    ifelse(rass <= CFG$rass_coma, "G3_coma",
      ifelse(is.na(capd), NA_character_,
        ifelse(capd >= CFG$capd_corte, "G2_delirium", "G1_free"))))
}
D$fenotipo_D1 <- classificar(D[[VAR$rass_d1]], D[[VAR$capd_d1]])
D$fenotipo_D5 <- classificar(D[[VAR$rass_d5]], D[[VAR$capd_d5]])

registrar("\n=== 1.4 FENOTIPO CONTEMPORANEO ===")
registrar("D1: %s", paste(names(table(D$fenotipo_D1)), table(D$fenotipo_D1),
                          sep = "=", collapse = " | "))
registrar("D5: %s", paste(names(table(D$fenotipo_D5)), table(D$fenotipo_D5),
                          sep = "=", collapse = " | "))
registrar("\nMatriz de transicao D1 x D5 (contagens absolutas):")
tab_trans <- table(D1 = D$fenotipo_D1, D5 = D$fenotipo_D5)
print(addmargins(tab_trans))
write.csv(as.data.frame.matrix(addmargins(tab_trans)),
          file.path(DIR_OUT, "01_matriz_transicao_D1xD5.csv"))

D$avaliavel_capd <- !is.na(D[[VAR$capd_d1]]) & !is.na(D[[VAR$capd_d5]])
n_av <- sum(D$avaliavel_capd)

registrar("\n=== 1.5 AMOSTRA ANALITICA ===")
registrar("Coorte total          : %d", nrow(D))
registrar("Avaliaveis pelo CAPD  : %d", n_av)
registrar("Nao avaliaveis (coma) : %d", nrow(D) - n_av)
stopifnot(n_av == 34)

stopifnot(all(D[[VAR$rass_d1]][D$avaliavel_capd] > CFG$rass_coma))
stopifnot(all(D[[VAR$rass_d5]][D$avaliavel_capd] > CFG$rass_coma))
registrar("Todos os %d avaliaveis tem RASS >= -3 em D1 e D5: OK", n_av)

if (!"capd_anchor" %in% names(D)) stop("Coluna capd_anchor ausente no banco.")
if ("idm" %in% names(D)) {
esperada <- cut(D$idm, breaks = c(-Inf, 6, 12, 24, Inf),
                labels = c("0-6m", ">6-12m", ">12-24m", ">24m"), right = TRUE)
divergentes <- which(as.character(esperada) != as.character(D$capd_anchor))
registrar("\n=== 1.6 ANCORA DESENVOLVIMENTAL ===")
registrar("Preenchidas: %d/%d | divergentes da faixa etaria: %d",
          sum(!is.na(D$capd_anchor)), nrow(D), length(divergentes))
if (length(divergentes))
  print(data.frame(ID = D$ID[divergentes], idade = D$idm[divergentes],
                   anchor = D$capd_anchor[divergentes],
                   esperada = esperada[divergentes]))
} else registrar("Anchor check skipped: age in months not available in the repository dataset.")
D$capd_anchor <- factor(D$capd_anchor, levels = c("0-6m", ">6-12m", ">12-24m", ">24m"))
print(table(D$capd_anchor[D$avaliavel_capd]))

D$delta_log_nfl   <- log(D[[VAR$nfl_d5]])   - log(D[[VAR$nfl_d1]])
D$delta_log_s100b <- log(D[[VAR$s100b_d5]]) - log(D[[VAR$s100b_d1]])
D$delta_capd      <- D[[VAR$capd_d5]] - D[[VAR$capd_d1]]
for (i in 1:8) D[[paste0("delta_Q", i)]] <- D[[ITENS_D5[i]]] - D[[ITENS_D1[i]]]

registrar("\n=== 1.9 COMPLETUDE (n = %d avaliaveis) ===", n_av)
chave <- c(unlist(VAR), ITENS_D1, ITENS_D5)
falt <- sapply(chave, function(v) sum(is.na(D[[v]][D$avaliavel_capd])))
if (all(falt == 0)) {
  registrar("Nenhum valor faltante nas variaveis-chave. Analise de casos")
  registrar("completos e a analise principal; imputacao nao se aplica.")
} else {
  registrar("ATENCAO: valores faltantes nas variaveis-chave:")
  print(falt[falt > 0])
}

attr(D, "VAR") <- VAR; attr(D, "ITENS_D1") <- ITENS_D1
attr(D, "ITENS_D5") <- ITENS_D5; attr(D, "ROTULO_ITEM") <- ROTULO_ITEM
saveRDS(D, file.path(DIR_OUT, "dados_artigo2.rds"))

.chk <- readRDS(file.path(DIR_OUT, "dados_artigo2.rds"))
.at  <- c("VAR", "ITENS_D1", "ITENS_D5", "ROTULO_ITEM")
.sem <- .at[vapply(.at, function(a) is.null(attr(.chk, a)), TRUE)]
if (nrow(.chk) != nrow(D) || length(.sem))
  stop("Falha ao gravar dados_artigo2.rds",
       if (length(.sem)) paste0(": atributos ausentes — ", paste(.sem, collapse = ", ")) else "")
registrar("dados_artigo2.rds gravado e relido: %d linhas, %d atributos.",
          nrow(.chk), length(.at))
rm(.chk, .at, .sem)

close(log_con)

cat("\n=== SCRIPT 01 CONCLUIDO ===\n")
cat("Objeto salvo: output/dados_artigo2.rds (formato largo)\n")
cat("Log de validacao: output/01_log_validacao.txt\n")
cat("Proximo passo: 02_descritivas.R\n\n")
