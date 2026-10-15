#' Comparacoes de medias
#' @export
comparar <- function(analise, metodo = c("tukey", "scheffe", "scott_knott",
                                        "dunnett", "duncan"),
                     alfa = 0.05, controle = NULL) {
  validar_analise(analise)
  metodo <- match.arg(metodo)
  if (length(alfa) != 1L || !is.numeric(alfa) || !is.finite(alfa) ||
      alfa <= 0 || alfa >= 1) stop("alfa deve estar entre zero e um.")
  m <- analise$modelo
  if (metodo == "tukey")
    return(stats::TukeyHSD(m, "tratamento", conf.level = 1 - alfa))
  pacote <- switch(metodo, scheffe = "agricolae", duncan = "agricolae",
                    scott_knott = "ExpDes.pt", dunnett = "multcomp")
  if (!requireNamespace(pacote, quietly = TRUE))
    stop("Instale o pacote '", pacote, "' para este teste.")
  if (metodo == "scheffe")
    return(agricolae::scheffe.test(m, "tratamento", alpha = alfa, console = FALSE))
  if (metodo == "duncan")
    return(agricolae::duncan.test(m, "tratamento", alpha = alfa, console = FALSE))
  if (metodo == "scott_knott") {
    # A dependencia usa um arquivo temporario de nome fixo.
    anterior <- getwd()
    temporario <- tempfile("exp-sk-")
    dir.create(temporario)
    on.exit({setwd(anterior); unlink(temporario, recursive = TRUE)}, add = TRUE)
    setwd(temporario)
    # ExpDes.pt imprime a tabela; preservar tambem a saida para consulta.
    resultado <- NULL
    saida <- utils::capture.output(resultado <- ExpDes.pt::scottknott(
      analise$dados$resposta, analise$dados$tratamento,
      stats::df.residual(m), stats::deviance(m), alpha = alfa))
    return(list(resultado = resultado, saida = saida, alfa = alfa))
  }
  if (is.null(controle) || length(controle) != 1L ||
      !as.character(controle) %in% levels(analise$dados$tratamento))
    stop("Informe em controle o nome de um tratamento testemunha existente.")
  d <- analise$dados
  d$tratamento <- stats::relevel(d$tratamento, ref = as.character(controle))
  m <- stats::aov(stats::formula(m), data = d)
  g <- multcomp::glht(m, linfct = multcomp::mcp(tratamento = "Dunnett"))
  list(teste = summary(g), intervalos = stats::confint(g, level = 1 - alfa),
       alfa = alfa, controle = controle)
}
