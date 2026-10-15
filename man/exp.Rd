\name{exp}
\alias{dic}
\alias{dbc}
\alias{dql}
\alias{comparar}
\alias{regressao}
\title{Analise de experimentos balanceados}
\description{ANOVA, comparacoes de medias e regressao polinomial em DIC, DBC e DQL.}
\usage{
dic(dados, resposta, tratamento)
dbc(dados, resposta, tratamento, bloco)
dql(dados, resposta, tratamento, linha, coluna)
comparar(analise, metodo = c("tukey", "scheffe", "scott_knott", "dunnett", "duncan"),
         alfa = 0.05, controle = NULL)
regressao(analise, grau = 1L)
}
\arguments{
\item{dados}{Data frame completo e balanceado, sem valores ausentes.}
\item{resposta}{Nome da coluna numerica de respostas, como texto.}
\item{tratamento}{Nome da coluna de tratamentos. Para regressao, deve ser numerica.}
\item{bloco}{Nome da coluna de blocos.}
\item{linha, coluna}{Nomes das colunas de linhas e colunas do quadrado latino.}
\item{analise}{Objeto retornado por dic, dbc ou dql.}
\item{metodo}{Teste de comparacao de medias.}
\item{alfa}{Nivel de significancia entre zero e um.}
\item{controle}{Rotulo da testemunha, obrigatorio para Dunnett.}
\item{grau}{Grau inteiro positivo menor que o numero de doses distintas.}
}
\details{
Modelos aditivos com efeitos fixos, erros independentes, normais e de variancia
constante. Inspecione os residuos com \code{plot(analise$modelo)}.
DBC requer uma observacao por tratamento por bloco. DQL requer quadrado latino
completo de ordem ao menos 3. As comparacoes nao sao condicionadas ao teste F.
Duncan nao garante o mesmo controle do erro familiar de Tukey.
Tukey usa stats; Scheffe e Duncan requerem agricolae; Dunnett requer multcomp;
Scott-Knott usa ExpDes.pt e devolve tambem sua tabela textual.
A regressao ajusta blocos ou linhas e colunas. A falta de ajuste compara o
polinomio com o modelo de medias usando seu erro experimental; retorna NULL
quando o polinomio satura as medias dos tratamentos.
}
\value{
Os delineamentos retornam objeto exp_anova com modelo aov, ANOVA, dados e medias.
comparar retorna o objeto do teste; Dunnett retorna teste, intervalos, alfa e
controle; Scott-Knott retorna resultado, saida textual e alfa. regressao retorna
modelo lm, coeficientes, ANOVA sequencial, falta de ajuste e grau.
}
\examples{
set.seed(42)
d <- data.frame(dose = rep(c(0, 10, 20, 30), each = 5))
d$y <- 5 + 0.8 * d$dose + rnorm(nrow(d))
a <- dic(d, "y", "dose")
a
comparar(a, "tukey")
regressao(a, grau = 2)
}
