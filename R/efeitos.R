# Medias marginais ou efeitos simples, sempre com erro do modelo completo.
selecionar_efeito <- function(analise, fator, dentro) {
  validar_analise(analise)
  if (is.null(fator)) {
    if (length(analise$fatores) > 1L)
      stop("Em fatoriais, informe fator e, para efeitos simples, dentro.")
    fator <- names(analise$fatores)[1]
  }
  if (length(fator) != 1L || !fator %in% names(analise$fatores))
    stop("fator deve ser o nome de um dos fatores da analise.")
  if (!is.list(dentro) || (length(dentro) &&
      (is.null(names(dentro)) || any(!nzchar(names(dentro))) ||
       anyDuplicated(names(dentro)) || fator %in% names(dentro) ||
       !all(names(dentro) %in% names(analise$fatores)) || any(lengths(dentro) != 1L))))
    stop("dentro deve ser lista nomeada de niveis dos outros fatores.")
  d <- analise$dados
  selecionar <- rep(TRUE, nrow(d))
  for (nm in names(dentro)) {
    if (is.na(dentro[[nm]]) || !as.character(dentro[[nm]]) %in%
        levels(d[[analise$fatores[[nm]]]])) stop("Nivel inexistente em dentro.")
    selecionar <- selecionar & d[[analise$fatores[[nm]]]] == as.character(dentro[[nm]])
  }
  # Nao ocultar interacoes relevantes ao interpretar medias marginais.
  restantes <- setdiff(names(analise$fatores), c(fator, names(dentro)))
  if (length(restantes)) {
    tab <- stats::anova(analise$modelo)
    interacoes <- strsplit(rownames(tab), ":", fixed = TRUE)
    relevantes <- vapply(interacoes, function(x)
      length(x) > 1L && analise$fatores[[fator]] %in% x &&
        any(unname(analise$fatores[restantes]) %in% x), logical(1))
    if (any(tab$`Pr(>F)`[relevantes] < 0.05, na.rm = TRUE))
      warning("Ha interacao significativa (5%) com fator nao fixado; considere desdobrar com dentro.")
  }
  interno <- analise$fatores[[fator]]
  m <- analise$modelo
  X <- stats::model.matrix(m)
  niveis <- levels(d[[interno]])
  L <- do.call(rbind, lapply(niveis, function(nv)
    colMeans(X[selecionar & d[[interno]] == nv, , drop = FALSE])))
  rownames(L) <- niveis
  media <- as.vector(L %*% stats::coef(m))
  V <- L %*% stats::vcov(m) %*% t(L)
  gl <- stats::df.residual(m)
  qm <- stats::deviance(m) / gl
  list(fator = fator, dentro = dentro, interno = interno, L = L,
       medias = stats::setNames(media, niveis), vcov = V,
       gl = gl, qm = qm, dados = d[selecionar, , drop = FALSE],
       original = analise$originais[[fator]][selecionar])
}

#' Teste F de um efeito marginal ou simples
#' @export
desdobrar <- function(analise, fator = NULL, dentro = list()) {
  e <- selecionar_efeito(analise, fator, dentro)
  k <- length(e$medias)
  C <- cbind(-1, diag(k - 1L))
  h <- as.vector(C %*% e$medias)
  V <- C %*% e$vcov %*% t(C)
  F <- as.numeric(crossprod(h, solve(V, h))) / (k - 1L)
  list(fator = e$fator, dentro = dentro, medias = e$medias,
       teste = data.frame(GL = k - 1L, GL_erro = e$gl, QM_erro = e$qm,
                          F = F, p_valor = stats::pf(F, k - 1L, e$gl, lower.tail = FALSE)))
}
