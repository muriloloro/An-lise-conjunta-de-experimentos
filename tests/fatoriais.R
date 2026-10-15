library(exp)
set.seed(2026)
falha <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
perto <- function(x, y) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-7)))
casos <- 0L
for (nf in 1:3) {
  fatores <- c("dose", "cultivar", "manejo")[seq_len(nf)]
  grade <- expand.grid(c(list(dose = c(0, 1, 2, 3)),
                        if (nf >= 2) list(cultivar = c("A", "B")),
                        if (nf == 3) list(manejo = c("M1", "M2"))))
  k <- nrow(grade)
  for (del in c("DIC", "DBC", "DQL")) {
    if (del == "DQL") {
      pos <- expand.grid(linha = seq_len(k), coluna = seq_len(k))
      indice <- (pos$linha + pos$coluna - 2L) %% k + 1L
      d <- cbind(grade[indice, , drop = FALSE], pos)
    } else {
      d <- grade[rep(seq_len(k), 4), , drop = FALSE]
      d$bloco <- rep(seq_len(4), each = k)
    }
    d$y <- 10 + 3 * d$dose + 0.5 * d$dose^2 + rnorm(nrow(d), sd = 0.3)
    if (nf >= 2) d$y <- d$y + 4 * d$dose * (d$cultivar == "B")
    if (nf == 3) d$y <- d$y + 5 * d$dose * (d$cultivar == "B") * (d$manejo == "M2")
    a <- switch(del, DIC = dic(d, "y", fatores),
                DBC = dbc(d, "y", fatores, "bloco"),
                DQL = dql(d, "y", fatores, "linha", "coluna"))
    rd <- d
    rd[fatores] <- lapply(rd[fatores], factor)
    nuisance <- switch(del, DIC = character(), DBC = "factor(bloco)",
                       DQL = c("factor(linha)", "factor(coluna)"))
    referencia <- aov(reformulate(c(nuisance, paste(fatores, collapse = "*")), "y"), rd)
    perto(a$anova$`Sum Sq`, anova(referencia)$`Sum Sq`)
    stopifnot(nrow(a$anova) == 2^nf + length(nuisance),
              identical(names(a$medias), c(fatores, "media")),
              paste(fatores, collapse = ":") %in% rownames(a$anova))
    gl <- df.residual(referencia)
    qm <- deviance(referencia) / gl
    niveis <- list(list())
    if (nf >= 2) niveis <- c(niveis, list(list(cultivar = "A")))
    if (nf == 3) niveis <- c(niveis, list(list(cultivar = "B", manejo = "M2")))
    for (dentro in niveis) {
      sel <- rep(TRUE, nrow(d))
      for (nm in names(dentro)) sel <- sel & d[[nm]] == dentro[[nm]]
      medias <- tapply(d$y[sel], d$dose[sel], mean)
      n <- sum(sel) / 4
      efeito <- suppressWarnings(desdobrar(a, fator = "dose", dentro = dentro))
      perto(efeito$teste$F, n * sum((medias - mean(medias))^2) / (3 * qm))
      perto(efeito$teste$GL_erro, gl)
      tukey <- suppressWarnings(comparar(a, fator = "dose", dentro = dentro))
      stopifnot(length(capture.output(print(tukey))) > 0L)
      t <- tukey[[1]]
      perto(t["1-0", "diff"], medias[2] - medias[1])
      perto(t["1-0", "upr"] - t["1-0", "diff"], qtukey(.95, 4, gl) * sqrt(qm / n))
      for (metodo in c("scheffe", "duncan", "scott_knott", "dunnett")) {
        r <- suppressWarnings(comparar(a, metodo, fator = "dose", dentro = dentro, controle = "0"))
        if (metodo %in% c("scheffe", "duncan")) {
          perto(r$statistics$MSerror, qm)
          perto(r$statistics$Df, gl)
          stopifnot(nrow(r$groups) == 4L)
        }
        if (metodo == "scott_knott") {
          perto(r$gl_erro, gl)
          stopifnot(any(grepl("Tratamentos", r$saida)))
        }
        if (metodo == "dunnett") {
          perto(r$teste$test$coefficients[1], medias[2] - medias[1])
          perto(r$teste$test$sigma[1], sqrt(2 * qm / n))
          perto(r$teste$df, gl)
        }
      }
      reg <- suppressWarnings(regressao(a, 2, fator = "dose", dentro = dentro))
      P <- cbind(1, 0:3, (0:3)^2)
      beta <- solve(crossprod(P), crossprod(P, medias))
      perto(reg$coeficientes[, 1], beta)
      perto(reg$coeficientes[, 2], sqrt(diag(solve(crossprod(P))) * qm / n))
      esperado <- sum((as.numeric(medias) - as.vector(P %*% beta))^2) * n / qm
      perto(reg$falta_de_ajuste$F[2], esperado)
      perto(reg$gl_erro, gl)
    }
    if (nf > 1) {
      falha(comparar(a))
      falha(regressao(a))
      falha(suppressWarnings(regressao(a, fator = "cultivar")))
      falha(comparar(a, fator = "dose", dentro = list(cultivar = "inexistente")))
      avisou <- FALSE
      withCallingHandlers(comparar(a, fator = "dose"), warning = function(w) {
        avisou <<- TRUE; invokeRestart("muffleWarning")
      })
      stopifnot(avisou)
      # Comparacao do segundo fator dentro de uma dose (e do terceiro fator).
      dentro <- c(list(dose = 1), if (nf == 3) list(manejo = "M1"))
      r <- comparar(a, "dunnett", controle = "A", fator = "cultivar", dentro = dentro)
      stopifnot(length(r$teste$test$pvalues) == 1L)
      for (alvo in fatores[-1]) {
        cond <- lapply(d[1, setdiff(fatores, alvo), drop = FALSE], function(x) x[1])
        for (metodo in c("tukey", "scheffe", "duncan", "scott_knott", "dunnett")) {
          rr <- comparar(a, metodo, fator = alvo, dentro = cond, controle = as.character(d[[alvo]][1]))
          stopifnot(!is.null(rr))
        }
      }
    }
    casos <- casos + 1L
  }
}
# Fatorial incompleto e DQL de ordem incorreta nao sao aceitos.
d <- expand.grid(A = letters[1:2], B = 1:3, rep = 1:3)
d$y <- rnorm(nrow(d))
falha(dic(d[!(d$A == "a" & d$B == 1), ], "y", c("A", "B")))
falha(dic(d, "y", c("A", "B", "rep", "y")))
falha(comparar(dic(d, "y", c("A", "B")), fator = "A", dentro = list(A = "a")))
cat("Fatoriais: ", casos, " configuracoes DIC/DBC/DQL x uni/bi/tri validadas.\n", sep = "")
