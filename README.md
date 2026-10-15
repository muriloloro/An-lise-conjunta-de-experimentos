# exp

Pacote R para análise de experimentos com efeitos fixos em delineamento
inteiramente casualizado (DIC), blocos casualizados (DBC) e quadrado latino (DQL).

## Instalação

```r
install.packages("remotes")
remotes::install_github("muriloloro/An-lise-conjunta-de-experimentos",
                        dependencies = TRUE)
library(exp)
```

O nome do pacote é `exp`; o nome do repositório permanece o existente.
A função matemática `exp()` de R continua disponível normalmente.

## Exemplo reproduzível

```r
set.seed(42)
dados <- expand.grid(dose = c(0, 10, 20, 30), bloco = letters[1:5])
dados$producao <- 10 + 0.7 * dados$dose + rnorm(20, sd = 2)

# Os argumentos de colunas são nomes entre aspas.
a <- dic(dados, "producao", "dose")
b <- dbc(dados, "producao", "dose", "bloco")
b                         # Tabela ANOVA
b$medias                  # Médias dos tratamentos
plot(b$modelo)            # Diagnósticos dos resíduos

comparar(b, "tukey")
comparar(b, "scheffe")
comparar(b, "scott_knott")$saida
comparar(b, "dunnett", controle = "0")
comparar(b, "duncan")

r <- regressao(b, grau = 2) # Linear: 1; quadrática: 2; cúbica: 3
summary(r$modelo)
r$falta_de_ajuste

quadrado <- expand.grid(linha = 1:4, coluna = 1:4)
quadrado$tratamento <- (quadrado$linha + quadrado$coluna - 2) %% 4
quadrado$y <- 5 + quadrado$tratamento + rnorm(16)
c <- dql(quadrado, "y", "tratamento", "linha", "coluna")
c
```

## Escopo e interpretação

- Versão inicial para delineamentos completos e balanceados, sem valores ausentes.
- DBC: exatamente uma observação por tratamento por bloco.
- DQL: quadrado completo de ordem pelo menos 3; cada tratamento ocorre uma vez
  por linha e coluna.
- Tukey usa `stats`; Scheffé e Duncan usam `agricolae`; Scott-Knott usa
  `ExpDes.pt`; Dunnett usa `multcomp` e exige testemunha explícita.
- Os testes são solicitados pelo usuário, sem condicionamento automático ao F.
  Duncan não oferece o mesmo controle de erro familiar de Tukey.
- A regressão exige doses numéricas e preserva o ajuste para blocos ou para
  linhas e colunas. O grau precisa ser menor que o número de doses distintas.
  A falta de ajuste usa o erro do modelo de médias de tratamentos. Não há
  teste de falta de ajuste quando o grau satura essas médias.
- Verifique independência, normalidade e homogeneidade dos erros antes da
  interpretação. O pacote não determina automaticamente se as hipóteses valem.
- Regressões não lineares e análise conjunta de experimentos não fazem parte
  desta primeira versão.

## Desenvolvimento

Com R >= 4.0 e os pacotes sugeridos instalados:

```sh
R CMD build .
R CMD check --no-manual exp_0.1.0.tar.gz
```

O check executa exemplos e os testes em `tests/validacao.R`. O GitHub Actions
repete a checagem a cada push e pull request.

Licença: GPL-3.
