#' Comparacoes de medias marginais ou efeitos simples
#' @export
comparar <- function(analise, metodo = c("tukey", "scheffe", "scott_knott",
                                        "dunnett", "duncan"),
                     alfa = 0.05, controle = NULL, fator = NULL, dentro = list()) {
  metodo <- match.arg(metodo)
  if (length(alfa) != 1L || !is.numeric(alfa) || !is.finite(alfa) ||
      alfa <= 0 || alfa >= 1) stop("alfa deve estar entre zero e um.")
  e <- selecionar_efeito(analise, fator, dentro)
  m <- analise$modelo
  if (metodo == "tukey") {
    pares <- utils::combn(seq_along(e$medias), 2)
    C <- matrix(0, ncol(pares), length(e$medias))
    for (i in seq_len(nrow(C))) {C[i, pares[1, i]] <- -1; C[i, pares[2, i]] <- 1}
    dif <- as.vector(C %*% e$medias)
    se <- sqrt(diag(C %*% e$vcov %*% t(C)))
    limite <- stats::qtukey(1 - alfa, length(e$medias), e$gl) * se / sqrt(2)
    tab <- cbind(diff = dif, lwr = dif - limite, upr = dif + limite,
                 `p adj` = stats::ptukey(abs(dif) * sqrt(2) / se,
                                        length(e$medias), e$gl, lower.tail = FALSE))
    rownames(tab) <- paste(names(e$medias)[pares[2, ]], names(e$medias)[pares[1, ]], sep = "-")
    # Conservar o formato de TukeyHSD da primeira versao.
    return(structure(stats::setNames(list(tab), e$interno), class = c("TukeyHSD", "multicomp"),
                     conf.level = 1 - alfa, ordered = FALSE, orig.call = match.call(),
                     fator = e$fator, dentro = dentro))
  }
  pacote <- switch(metodo, scheffe = "agricolae", duncan = "agricolae",
                    scott_knott = "ExpDes.pt", dunnett = "multcomp")
  if (!requireNamespace(pacote, quietly = TRUE)) stop("Instale o pacote '", pacote, "'.")
  y <- e$dados$resposta
  trt <- e$dados[[e$interno]]
  if (metodo == "scheffe")
    return(agricolae::scheffe.test(y, trt, DFerror = e$gl, MSerror = e$qm,
                                 alpha = alfa, console = FALSE))
  if (metodo == "duncan")
    return(agricolae::duncan.test(y, trt, DFerror = e$gl, MSerror = e$qm,
                                alpha = alfa, console = FALSE))
  if (metodo == "scott_knott") {
    anterior <- getwd()
    temporario <- tempfile("exp-sk-")
    dir.create(temporario)
    on.exit({setwd(anterior); unlink(temporario, recursive = TRUE)}, add = TRUE)
    setwd(temporario)
    resultado <- NULL
    saida <- utils::capture.output(resultado <- ExpDes.pt::scottknott(
      y, trt, e$gl, e$qm * e$gl, alpha = alfa))
    return(list(resultado = resultado, saida = saida, alfa = alfa,
                fator = e$fator, dentro = dentro, gl_erro = e$gl, qm_erro = e$qm))
  }
  if (is.null(controle) || length(controle) != 1L || is.na(controle) ||
      !as.character(controle) %in% names(e$medias))
    stop("Informe em controle o nivel testemunha do fator selecionado.")
  base <- match(as.character(controle), names(e$medias))
  outros <- setdiff(seq_along(e$medias), base)
  K <- e$L[outros, , drop = FALSE] - e$L[rep(base, length(outros)), , drop = FALSE]
  rownames(K) <- paste(names(e$medias)[outros], controle, sep = " - ")
  g <- multcomp::glht(m, linfct = K)
  list(teste = summary(g), intervalos = stats::confint(g, level = 1 - alfa),
       alfa = alfa, controle = controle, fator = e$fator, dentro = dentro)
}
