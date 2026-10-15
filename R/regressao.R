#' Regressao polinomial de tratamentos quantitativos
#' @export
regressao <- function(analise, grau = 1L) {
  validar_analise(analise)
  x <- analise$valores_tratamento
  if (!is.numeric(x) || any(!is.finite(x)))
    stop("Regressao requer a coluna original de tratamentos numerica e finita.")
  if (length(grau) != 1L || !is.numeric(grau) || !is.finite(grau) ||
      grau != as.integer(grau) || grau < 1L || grau >= length(unique(x)))
    stop("grau deve ser inteiro positivo menor que o numero de doses.")
  d <- analise$dados
  d$dose <- x
  termos <- c(if (analise$delineamento == "DBC") "bloco",
              if (analise$delineamento == "DQL") c("linha", "coluna"),
              sprintf("I(dose^%d)", seq_len(grau)))
  f <- stats::reformulate(termos, response = "resposta")
  modelo <- stats::lm(f, data = d)
  if (modelo$rank < length(stats::coef(modelo)))
    stop("Modelo singular: verifique a escala das doses e o grau.")
  # Compara o polinomio com o modelo de medias para detectar falta de ajuste.
  falta <- if (grau < length(unique(x)) - 1L)
    stats::anova(modelo, analise$modelo) else NULL
  list(modelo = modelo, coeficientes = summary(modelo)$coefficients,
       anova = stats::anova(modelo), falta_de_ajuste = falta, grau = grau)
}
