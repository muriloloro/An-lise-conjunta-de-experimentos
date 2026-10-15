library(exp)
falha <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
set.seed(123)
d <- expand.grid(trat = c(0, 10, 20, 30), bloco = letters[1:5])
d$y <- 10 + 0.7 * d$trat + as.integer(factor(d$bloco)) + rnorm(nrow(d), sd = 2)
a <- dic(d, "y", "trat")
b <- dbc(d, "y", "trat", "bloco")
stopifnot(isTRUE(all.equal(unname(a$anova[, "Sum Sq"]),
                          unname(anova(aov(y ~ factor(trat), d))[, "Sum Sq"]))),
          isTRUE(all.equal(unname(b$anova[, "Sum Sq"]),
                          unname(anova(aov(y ~ factor(bloco) + factor(trat), d))[, "Sum Sq"]))),
          df.residual(b$modelo) == 12,
          isTRUE(all.equal(comparar(b)$tratamento,
                           TukeyHSD(b$modelo, "tratamento")$tratamento)))
q <- expand.grid(linha = 1:4, coluna = 1:4)
q$trat <- (q$linha + q$coluna - 2) %% 4
q$y <- q$trat * 2 + q$linha + q$coluna + rnorm(16)
c <- dql(q, "y", "trat", "linha", "coluna")
stopifnot(df.residual(c$modelo) == 6)
for (z in list(a, b, c)) {
  for (metodo in c("scheffe", "duncan", "scott_knott", "dunnett")) {
    r <- comparar(z, metodo, controle = "0")
    stopifnot(!is.null(r))
    if (metodo == "scott_knott") stopifnot(any(nzchar(r$saida)))
    if (metodo == "dunnett") stopifnot(length(r$teste$test$pvalues) == 3L)
    if (metodo %in% c("scheffe", "duncan")) stopifnot(nrow(r$groups) == 4L)
  }
  r <- regressao(z, 2)
  stopifnot(nrow(r$falta_de_ajuste) == 2L, r$falta_de_ajuste$Df[2] == 1L,
            is.null(regressao(z, 3)$falta_de_ajuste))
}
r <- regressao(b, 1)
ref <- lm(y ~ factor(bloco) + trat, d)
stopifnot(isTRUE(all.equal(unname(coef(r$modelo)), unname(coef(ref)))))
falha(dic(d[-1, ], "y", "trat"))
falha(dbc(d[-1, ], "y", "trat", "bloco"))
q$trat[1] <- q$trat[2]
falha(dql(q, "y", "trat", "linha", "coluna"))
d$y[1] <- NA
falha(dic(d, "y", "trat"))
falha(comparar(a, "dunnett"))
falha(comparar(a, alfa = 1))
falha(regressao(a, 4))
cat("Validacao: tres delineamentos, cinco testes de medias, regressao e entradas invalidas OK\n")
