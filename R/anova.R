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
  colunas <- list(resposta, tratamento, bloco, linha, coluna)
  nomes <- unlist(colunas, use.names = FALSE)
  if (!is.data.frame(dados) || !is.character(nomes) ||
      length(resposta) != 1L || length(tratamento) < 1L || length(tratamento) > 3L ||
      any(lengths(colunas[c(1, 3:5)]) > 1L) || anyDuplicated(nomes) ||
      !all(nomes %in% names(dados)))
    stop("Informe colunas distintas existentes e de um a tres fatores em tratamento.")
  d <- dados[, nomes, drop = FALSE]
  if (anyNA(d)) stop("Dados ausentes: resolva os valores NA antes da analise.")
  if (!is.numeric(d[[resposta]]) || any(!is.finite(d[[resposta]])))
    stop("A resposta deve ser numerica e finita.")
  fatores <- if (length(tratamento) == 1L) "tratamento" else paste0("f", seq_along(tratamento))
  z <- data.frame(resposta = d[[resposta]])
  for (i in seq_along(fatores)) {
    z[[fatores[i]]] <- factor(d[[tratamento[i]]])
    if (nlevels(z[[fatores[i]]]) < 2L) stop("Cada fator requer ao menos dois niveis.")
  }
  celulas <- do.call(interaction, c(z[fatores], list(drop = FALSE, lex.order = TRUE)))
  if (any(table(celulas) == 0L)) stop("O fatorial deve conter todas as combinacoes.")
  estrutura <- paste(fatores, collapse = " * ")
  ajuste <- character()
  if (delineamento == "DIC") {
    if (length(unique(table(celulas))) != 1L) stop("DIC requer repeticoes iguais por combinacao.")
  } else if (delineamento == "DBC") {
    z$bloco <- factor(d[[bloco]])
    if (nlevels(z$bloco) < 2L || any(table(celulas, z$bloco) != 1L))
      stop("DBC requer uma observacao de cada combinacao por bloco e dois blocos.")
    ajuste <- "bloco"
  } else {
    z$linha <- factor(d[[linha]])
    z$coluna <- factor(d[[coluna]])
    k <- nlevels(celulas)
    if (k < 3L || nlevels(z$linha) != k || nlevels(z$coluna) != k ||
        any(table(z$linha, z$coluna) != 1L) ||
        any(table(celulas, z$linha) != 1L) || any(table(celulas, z$coluna) != 1L))
      stop("DQL requer quadrado latino completo: cada combinacao uma vez por linha e coluna.")
    ajuste <- c("linha", "coluna")
  }
  f <- stats::reformulate(c(ajuste, estrutura), response = "resposta")
  modelo <- stats::aov(f, data = z)
  if (stats::df.residual(modelo) < 1L) stop("Nao ha graus de liberdade residuais.")
  tabela <- stats::anova(modelo)
  # Manter nomes internos seguros para formulas; exibir os nomes fornecidos.
  rownames(tabela) <- vapply(strsplit(rownames(tabela), ":", fixed = TRUE), function(x) {
    idx <- match(x, fatores)
    x[!is.na(idx)] <- tratamento[idx[!is.na(idx)]]
    paste(x, collapse = ":")
  }, character(1))
  medias <- stats::aggregate(z$resposta, z[fatores], mean)
  names(medias) <- c(tratamento, "media")
  structure(list(delineamento = delineamento, modelo = modelo,
    anova = tabela, dados = z, fatores = stats::setNames(fatores, tratamento),
    originais = d[tratamento], ajuste = ajuste,
    valores_tratamento = if (length(tratamento) == 1L) d[[tratamento]] else NULL,
    medias = medias), class = "exp_anova")
}

#' @export
print.exp_anova <- function(x, ...) {
  cat("ANOVA -", x$delineamento, "-", length(x$fatores), "fator(es)\n")
  print(x$anova, ...)
  invisible(x)
}

validar_analise <- function(x) {
  if (!inherits(x, "exp_anova")) stop("Use um resultado de dic(), dbc() ou dql().")
}
