# 00b_motor_deltarho.R
# Resampling engine: Spearman correlations, signed differences between dependent correlations
# (studentized bootstrap-t with jackknife standard errors), and BCa intervals.
# Part of the analysis code for the CAPD item-level NfL/S100B analysis (MIT License).

stopifnot(exists("CFG"))

MOTOR_VERSAO <- "motor-deltarho-2026-09-08"
FRAC_MIN_REPLICAS <- 0.95

sp <- function(a, b) {
  ra <- rank(a); rb <- rank(b)
  ra <- ra - mean(ra); rb <- rb - mean(rb)
  s <- sqrt(sum(ra * ra) * sum(rb * rb))
  if (s <= 0) return(NA_real_)
  sum(ra * rb) / s
}

delta_stat <- function(x, y1, y2) sp(x, y1) - sp(x, y2)

jack_deltas <- function(x, y1, y2) {
  n <- length(x); d <- numeric(n)
  for (i in seq_len(n)) d[i] <- delta_stat(x[-i], y1[-i], y2[-i])
  d
}
jack_se <- function(x, y1, y2) {
  d <- jack_deltas(x, y1, y2)
  if (anyNA(d)) return(NA_real_)
  n <- length(d)
  sqrt((n - 1) / n * sum((d - mean(d))^2))
}

teste_delta_rho <- function(x, y1, y2, R = CFG$n_boot, seed = CFG$seed,
                            sobreposto = TRUE) {
  ok  <- !is.na(x) & !is.na(y1) & !is.na(y2)
  x <- x[ok]; y1 <- y1[ok]; y2 <- y2[ok]; n <- length(x)

  out <- list(n = n, delta = NA_real_, rho1 = NA_real_, rho2 = NA_real_,
              se_obs = NA_real_, t_obs = NA_real_, p = NA_real_,
              ic_inf = NA_real_, ic_sup = NA_real_, tipo_ic = NA_character_,
              bca_inf = NA_real_, bca_sup = NA_real_,
              R_valido = 0L, estado = "ok", sobreposto = sobreposto,
              motor = MOTOR_VERSAO)

  if (n < 4) { out$estado <- "n insuficiente"; return(out) }
  if (length(unique(x)) < 2 || length(unique(y1)) < 2 || length(unique(y2)) < 2) {
    out$estado <- "nao aplicavel: variavel constante"; return(out)
  }

  r1 <- sp(x, y1); r2 <- sp(x, y2); D <- r1 - r2
  out$rho1 <- r1; out$rho2 <- r2
  if (!is.finite(D)) { out$estado <- "estatistica observada nao finita"; return(out) }
  out$delta <- D

  dj <- jack_deltas(x, y1, y2)
  if (anyNA(dj)) { out$estado <- "jackknife nao finito na amostra original"; return(out) }
  SE <- sqrt((n - 1) / n * sum((dj - mean(dj))^2))
  out$se_obs <- SE
  if (!is.finite(SE) || SE <= 0) {
    out$estado <- "SE da amostra original nula ou nao finita"; return(out)
  }
  Tobs <- D / SE; out$t_obs <- Tobs

  set.seed(seed)
  Ts <- rep(NA_real_, R); Ds <- rep(NA_real_, R)
  for (b in seq_len(R)) {
    i  <- sample.int(n, n, replace = TRUE)
    db <- delta_stat(x[i], y1[i], y2[i])
    sb <- jack_se(x[i], y1[i], y2[i])
    Ds[b] <- db
    if (is.finite(db) && is.finite(sb) && sb > 0) Ts[b] <- (db - D) / sb
  }
  val <- is.finite(Ts); Rv <- sum(val); out$R_valido <- Rv
  if (Rv < FRAC_MIN_REPLICAS * R) {
    out$estado <- sprintf("R_valido insuficiente (%d/%d)", Rv, R); return(out)
  }
  if (stats::var(Ts[val]) <= .Machine$double.eps) {
    out$estado <- "distribuicao de T* degenerada"; return(out)
  }

  out$p <- (1 + sum(abs(Ts[val]) >= abs(Tobs))) / (Rv + 1)

  qt <- stats::quantile(Ts[val], c(0.975, 0.025), names = FALSE, type = 7)
  out$ic_inf <- D - qt[1] * SE
  out$ic_sup <- D - qt[2] * SE
  out$tipo_ic <- "bootstrap-t"

  okD <- is.finite(Ds)
  if (sum(okD) >= 0.95 * R) {
    z0v <- mean(Ds[okD] < D) + 0.5 * mean(Ds[okD] == D)
    if (z0v > 0 && z0v < 1) {
      z0 <- stats::qnorm(z0v)
      dbar <- mean(dj); u <- dbar - dj
      den <- 6 * (sum(u^2))^1.5
      a <- if (den > 0) sum(u^3) / den else 0
      zc <- stats::qnorm(c(0.025, 0.975))
      alf <- stats::pnorm(z0 + (z0 + zc) / (1 - a * (z0 + zc)))
      if (all(is.finite(alf))) {
        qb <- stats::quantile(Ds[okD], alf, names = FALSE, type = 7)
        out$bca_inf <- qb[1]; out$bca_sup <- qb[2]
      }
    }
  }
  out
}

sp_parcial <- function(x, y, z) {
  rxy <- sp(x, y); rxz <- sp(x, z); ryz <- sp(y, z)
  if (!is.finite(rxy) || !is.finite(rxz) || !is.finite(ryz)) return(NA_real_)
  den <- sqrt((1 - rxz^2) * (1 - ryz^2))
  if (!is.finite(den) || den <= 0) return(NA_real_)
  (rxy - rxz * ryz) / den
}
delta_parcial <- function(x, y1, y2, z) sp_parcial(x, y1, z) - sp_parcial(x, y2, z)

jack_se_parcial <- function(x, y1, y2, z) {
  n <- length(x); d <- numeric(n)
  for (i in seq_len(n)) d[i] <- delta_parcial(x[-i], y1[-i], y2[-i], z[-i])
  if (anyNA(d)) return(NA_real_)
  sqrt((n - 1) / n * sum((d - mean(d))^2))
}

teste_delta_parcial <- function(x, y1, y2, z, R = CFG$n_boot, seed = CFG$seed) {
  ok <- !is.na(x) & !is.na(y1) & !is.na(y2) & !is.na(z)
  x <- x[ok]; y1 <- y1[ok]; y2 <- y2[ok]; z <- z[ok]; n <- length(x)
  out <- list(n = n, delta = NA_real_, rho1 = NA_real_, rho2 = NA_real_,
              se_obs = NA_real_, t_obs = NA_real_, p = NA_real_,
              ic_inf = NA_real_, ic_sup = NA_real_, tipo_ic = "bootstrap-t",
              bca_inf = NA_real_, bca_sup = NA_real_,
              R_valido = 0L, estado = "ok", sobreposto = TRUE,
              motor = paste0(MOTOR_VERSAO, "-parcial"))
  if (n < 5) { out$estado <- "n insuficiente"; return(out) }
  if (length(unique(z)) < 2) { out$estado <- "covariavel constante"; return(out) }

  r1 <- sp_parcial(x, y1, z); r2 <- sp_parcial(x, y2, z); D <- r1 - r2
  out$rho1 <- r1; out$rho2 <- r2
  if (!is.finite(D)) { out$estado <- "parcial observada nao finita"; return(out) }
  out$delta <- D
  SE <- jack_se_parcial(x, y1, y2, z); out$se_obs <- SE
  if (!is.finite(SE) || SE <= 0) {
    out$estado <- "SE da amostra original nula ou nao finita"; return(out)
  }
  Tobs <- D / SE; out$t_obs <- Tobs

  set.seed(seed)
  Ts <- rep(NA_real_, R)
  for (b in seq_len(R)) {
    i  <- sample.int(n, n, replace = TRUE)
    db <- delta_parcial(x[i], y1[i], y2[i], z[i])
    sb <- jack_se_parcial(x[i], y1[i], y2[i], z[i])
    if (is.finite(db) && is.finite(sb) && sb > 0) Ts[b] <- (db - D) / sb
  }
  val <- is.finite(Ts); Rv <- sum(val); out$R_valido <- Rv
  if (Rv < FRAC_MIN_REPLICAS * R) {
    out$estado <- sprintf("R_valido insuficiente (%d/%d)", Rv, R); return(out)
  }
  out$p <- (1 + sum(abs(Ts[val]) >= abs(Tobs))) / (Rv + 1)
  qt <- stats::quantile(Ts[val], c(0.975, 0.025), names = FALSE, type = 7)
  out$ic_inf <- D - qt[1] * SE; out$ic_sup <- D - qt[2] * SE
  out
}

linha_delta <- function(r, ...) {
  data.frame(..., n = r$n,
             rho1 = r$rho1, rho2 = r$rho2, Delta_rho = r$delta,
             SE_jack = r$se_obs, t_obs = r$t_obs,
             IC_inf = r$ic_inf, IC_sup = r$ic_sup, Tipo_IC = r$tipo_ic,
             BCa_inf = r$bca_inf, BCa_sup = r$bca_sup,
             p_bootstrap_t = r$p, R_valido = r$R_valido,
             estado = r$estado, stringsAsFactors = FALSE)
}

cat("Motor de Delta-rho carregado:", MOTOR_VERSAO, "\n")
cat("Bootstrap-t estudentizado, EP por jackknife delete-one com postos\n")
cat("recalculados; pivo na estatistica observada; R =", CFG$n_boot, "\n")
cat("IC principal: bootstrap-t. BCa reportado como diagnostico secundario,\n")
cat("com cobertura NAO verificada.\n\n")
