#' Analise de variancia de experimentos balanceados
#' @export
dic <- function(dados, resposta, tratamento) {
  ajustar(dados, resposta, tratamento, delineamento = "DIC")
}

#' @export
dbc <- function(dados, resposta, tratamento, bloco) {
  ajustar(dados, resposta, tratamento, bloco = bloco, delineamento = "DBC")
}

#' @export
dql <- function(dados, resposta, tratamento, linha, coluna) {
  ajustar(dados, resposta, tratamento, linha = linha, coluna = coluna,
          delineamento = "DQL")
}

ajustar <- function(dados, resposta, tratamento, bloco = NULL,
                    linha = NULL, coluna = NULL, delineamento) {
  nomes <- c(resposta, tratamento, bloco, linha, coluna)
  if (!is.data.frame(dados) || any(lengths(list(resposta, tratamento)) != 1L) ||
      anyDuplicated(nomes) || !all(nomes %in% names(dados)))
    stop("Informe um data.frame e nomes distintos de colunas existentes.")
  d <- dados[, nomes, drop = FALSE]
  if (anyNA(d)) stop("Dados ausentes: resolva os valores NA antes da analise.")
  if (!is.numeric(d[[resposta]]) || any(!is.finite(d[[resposta]])))
    stop("A resposta deve ser numerica e finita.")
  z <- data.frame(resposta = d[[resposta]], tratamento = factor(d[[tratamento]]))
  if (nlevels(z$tratamento) < 2L) stop("Sao necessarios ao menos dois tratamentos.")
  if (delineamento == "DIC") {
    n <- table(z$tratamento)
    if (length(unique(n)) != 1L) stop("Esta versao requer DIC balanceado.")
    f <- resposta ~ tratamento
  } else if (delineamento == "DBC") {
    z$bloco <- factor(d[[bloco]])
    if (nlevels(z$bloco) < 2L || any(table(z$tratamento, z$bloco) != 1L))
      stop("DBC requer uma observacao por tratamento em cada bloco e dois blocos.")
    f <- resposta ~ bloco + tratamento
  } else {
    z$linha <- factor(d[[linha]])
    z$coluna <- factor(d[[coluna]])
    k <- nlevels(z$tratamento)
    if (k < 3L || nlevels(z$linha) != k || nlevels(z$coluna) != k ||
        any(table(z$linha, z$coluna) != 1L) ||
        any(table(z$tratamento, z$linha) != 1L) ||
        any(table(z$tratamento, z$coluna) != 1L))
      stop("DQL requer quadrado latino completo, de ordem ao menos 3.")
    f <- resposta ~ linha + coluna + tratamento
  }
  modelo <- stats::aov(f, data = z)
  if (stats::df.residual(modelo) < 1L) stop("Nao ha graus de liberdade residuais.")
  structure(list(delineamento = delineamento, modelo = modelo,
    anova = stats::anova(modelo), dados = z,
    valores_tratamento = d[[tratamento]],
    medias = stats::aggregate(resposta ~ tratamento, z, mean)), class = "exp_anova")
}

#' @export
print.exp_anova <- function(x, ...) {
  cat("ANOVA -", x$delineamento, "\n")
  print(x$anova, ...)
  invisible(x)
}

validar_analise <- function(x) {
  if (!inherits(x, "exp_anova")) stop("Use um resultado de dic(), dbc() ou dql().")
}
