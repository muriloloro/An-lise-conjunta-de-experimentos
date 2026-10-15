# exp

Pacote R para ANOVA de experimentos **unifatoriais, bifatoriais e trifatoriais**
em DIC, DBC e DQL, com efeitos fixos, dados completos e balanceados.

## Instalação

```r
install.packages("remotes")
remotes::install_github("muriloloro/An-lise-conjunta-de-experimentos",
                        dependencies = TRUE)
library(exp)
```

O pacote se chama `exp`; o repositório mantém o nome existente.
A função matemática `exp()` de R continua disponível normalmente.

## ANOVA com um, dois ou três fatores

Os nomes das colunas são passados entre aspas. Use um vetor em `tratamento`
para os fatoriais. Todas as interações entram na ANOVA.

```r
set.seed(42)
dados <- expand.grid(dose = 0:3, cultivar = c("A", "B"),
                     manejo = c("M1", "M2"), bloco = 1:4)
dados$y <- 10 + 2 * dados$dose +
  3 * dados$dose * (dados$cultivar == "B") +
  4 * dados$dose * (dados$cultivar == "B") * (dados$manejo == "M2") +
  rnorm(nrow(dados))

a <- dic(dados, "y", c("dose", "cultivar", "manejo"))
b <- dbc(dados, "y", c("dose", "cultivar", "manejo"), "bloco")
b
plot(b$modelo)

# Unifatorial e bifatorial seguem a mesma interface:
uni <- subset(dados, cultivar == "A" & manejo == "M1")
bi <- subset(dados, manejo == "M1")
dic(uni, "y", "dose")
dbc(bi, "y", c("dose", "cultivar"), "bloco")
```

## Desdobramentos e testes de médias

Em fatoriais, informe `fator`. `dentro` fixa níveis dos demais fatores.
Sem fixar um fator, a análise calcula a média sobre seus níveis.
Uma interação significativa a 5% com um fator não fixado gera aviso;
a escolha do desdobramento permanece explícita.

```r
# Dose dentro de cultivar A e manejo M1: efeito simples no trifatorial.
condicao <- list(cultivar = "A", manejo = "M1")
desdobrar(b, fator = "dose", dentro = condicao) # Teste F do efeito simples

comparar(b, "tukey", fator = "dose", dentro = condicao)
comparar(b, "scheffe", fator = "dose", dentro = condicao)
comparar(b, "scott_knott", fator = "dose", dentro = condicao)$saida
comparar(b, "dunnett", fator = "dose", dentro = condicao, controle = "0")
comparar(b, "duncan", fator = "dose", dentro = condicao)

# Cultivares dentro da dose 1 e manejo M1:
comparar(b, "tukey", fator = "cultivar", dentro = list(dose = 1, manejo = "M1"))

# Efeito marginal: use quando a interpretação das interações permitir.
comparar(b, "tukey", fator = "cultivar")
```

Os testes usam **o erro e os graus de liberdade do modelo completo**, inclusive
nos desdobramentos. Dunnett exige um nível testemunha do fator selecionado.
A correção de multiplicidade vale por chamada, sem ajuste conjunto entre
várias famílias de desdobramentos. Duncan não oferece o mesmo controle de erro
familiar de Tukey. Os testes não são condicionados automaticamente ao F.

Tukey usa a amplitude studentizada de `stats`; Scheffé e Duncan usam
`agricolae`; Scott-Knott usa `ExpDes.pt`; Dunnett usa contrastes de `multcomp`.

## Regressão para fatores quantitativos

```r
r <- regressao(b, grau = 2, fator = "dose", dentro = condicao)
r$coeficientes
r$falta_de_ajuste
r$medias
```

O fator original precisa ser numérico. Grau 1 corresponde à reta, 2 à quadrática
e 3 à cúbica; o grau deve ser menor que o número de doses distintas.
Coeficientes e falta de ajuste usam o erro experimental completo, também para
efeitos simples. Se o polinômio saturar as médias, a falta de ajuste é `NULL`.

`r$modelo` é o modelo **global reduzido**, com o polinômio interagindo com os
demais fatores e ajuste por blocos ou linhas/colunas. `summary(r$modelo)` usa o
resíduo desse modelo reduzido: para inferência do efeito selecionado com o erro
puro do experimento, consulte `r$coeficientes` e `r$falta_de_ajuste`.

## DQL fatorial

Cada combinação é um tratamento do quadrado. Um fatorial 2 × 3 exige ordem 6;
um 2 × 2 × 2 exige ordem 8. Cada combinação deve aparecer exatamente uma vez
por linha e coluna. Não é um quadrado greco-latino.

```r
combinacoes <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"))
k <- nrow(combinacoes)
pos <- expand.grid(linha = seq_len(k), coluna = seq_len(k))
indice <- (pos$linha + pos$coluna - 2) %% k + 1
q <- cbind(combinacoes[indice, ], pos)
q$y <- 10 + as.integer(q$A) + as.integer(q$B) + rnorm(nrow(q))
c <- dql(q, "y", c("A", "B"), "linha", "coluna")
c
comparar(c, "tukey", fator = "A", dentro = list(B = "b1"))
```

## Requisitos e validação

- DIC: repetições iguais por combinação; DBC: uma observação de cada combinação
  por bloco; DQL: quadrado completo de ordem ao menos 3.
- Sem dados ausentes, combinações faltantes, parcelas subdivididas ou efeitos
  aleatórios. Regressão não linear e análise conjunta não estão implementadas.
- Verifique independência, normalidade e homogeneidade dos erros.
- Os testes validam as nove configurações (três delineamentos × um, dois e três
  fatores), efeitos marginais e simples, os cinco métodos de médias e regressões
  contra cálculos de referência, além das entradas inválidas.

Com R >= 4.0 e os pacotes sugeridos instalados:

```sh
R CMD build .
R CMD check --no-manual exp_0.2.0.tar.gz
```

O GitHub Actions repete a checagem a cada push e pull request.
Licença: GPL-3.
