#' Regressao polinomial com erro experimental do modelo fatorial completo
#' @export
regressao <- function(analise, grau = 1L, fator = NULL, dentro = list()) {
  e <- selecionar_efeito(analise, fator, dentro)
  x <- e$original
  if (!is.numeric(x) || any(!is.finite(x)))
    stop("Regressao requer a coluna original do fator numerica e finita.")
  if (length(grau) != 1L || !is.numeric(grau) || !is.finite(grau) ||
      grau != as.integer(grau) || grau < 1L || grau >= length(unique(x)))
    stop("grau deve ser inteiro positivo menor que o numero de doses.")
  doses <- as.numeric(names(e$medias))
  P <- outer(doses, 0:grau, `^`)
  if (qr(P)$rank < ncol(P)) stop("Polinomio singular: reescale as doses ou reduza grau.")
  # GLS das medias estimadas; a covariancia usa sempre o erro puro completo.
  W <- solve(e$vcov)
  covbeta <- solve(t(P) %*% W %*% P)
  beta <- as.vector(covbeta %*% t(P) %*% W %*% e$medias)
  se <- sqrt(diag(covbeta))
  t <- beta / se
  coeficientes <- cbind(Estimate = beta, `Std. Error` = se, `t value` = t,
                        `Pr(>|t|)` = 2 * stats::pt(-abs(t), e$gl))
  rownames(coeficientes) <- c("(Intercept)", paste0("dose^", seq_len(grau)))
  residuos <- e$medias - as.vector(P %*% beta)
  gl_falta <- length(doses) - grau - 1L
  falta <- NULL
  if (gl_falta > 0L) {
    F <- as.numeric(t(residuos) %*% W %*% residuos) / gl_falta
    falta <- data.frame(Res.Df = c(e$gl + gl_falta, e$gl),
      Df = c(NA, gl_falta), F = c(NA, F),
      `Pr(>F)` = c(NA, stats::pf(F, gl_falta, e$gl, lower.tail = FALSE)),
      check.names = FALSE)
    rownames(falta) <- c("Polinomio", "Medias dos tratamentos")
  }
  # Modelo global reduzido: manter as interacoes com todos os demais fatores.
  d <- analise$dados
  d$dose <- analise$originais[[e$fator]]
  polinomio <- paste(sprintf("I(dose^%d)", seq_len(grau)), collapse = " + ")
  outros <- setdiff(unname(analise$fatores), e$interno)
  estrutura <- if (length(outros)) paste0("(", polinomio, ") * ", paste(outros, collapse = " * ")) else polinomio
  modelo <- stats::lm(stats::reformulate(c(analise$ajuste, estrutura), "resposta"), d)
  list(modelo = modelo, coeficientes = coeficientes, vcov = covbeta,
       anova = stats::anova(modelo), falta_de_ajuste = falta, grau = grau,
       fator = e$fator, dentro = dentro, gl_erro = e$gl, qm_erro = e$qm,
       medias = data.frame(dose = doses, media = unname(e$medias)))
}
