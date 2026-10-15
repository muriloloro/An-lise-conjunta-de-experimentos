\name{exp}
\alias{dic}
\alias{dbc}
\alias{dql}
\alias{comparar}
\alias{regressao}
\alias{desdobrar}
\title{Experimentos uni, bi e trifatoriais balanceados}
\description{ANOVA com todas as interacoes, testes F de efeitos simples,
comparacoes de medias e regressao polinomial em DIC, DBC e DQL.}
\usage{
dic(dados, resposta, tratamento)
dbc(dados, resposta, tratamento, bloco)
dql(dados, resposta, tratamento, linha, coluna)
comparar(analise, metodo = c("tukey", "scheffe", "scott_knott", "dunnett", "duncan"),
         alfa = 0.05, controle = NULL, fator = NULL, dentro = list())
regressao(analise, grau = 1L, fator = NULL, dentro = list())
desdobrar(analise, fator = NULL, dentro = list())
}
\arguments{
\item{dados}{Data frame completo e balanceado, sem valores ausentes.}
\item{resposta}{Nome da coluna numerica de respostas, como texto.}
\item{tratamento}{Vetor com um, dois ou tres nomes de colunas de fatores. Todas as combinacoes devem existir.}
\item{bloco}{Nome da coluna de blocos.}
\item{linha, coluna}{Nomes das colunas de linhas e colunas do quadrado latino.}
\item{analise}{Objeto retornado por dic, dbc ou dql.}
\item{metodo}{Teste de comparacao de medias.}
\item{alfa}{Nivel de significancia entre zero e um.}
\item{controle}{Nivel testemunha do fator selecionado, obrigatorio para Dunnett.}
\item{grau}{Grau inteiro positivo menor que o numero de doses distintas.}
\item{fator}{Nome original do fator de interesse. Obrigatorio em bi e trifatoriais.}
\item{dentro}{Lista nomeada que fixa niveis dos outros fatores, por exemplo
\code{list(cultivar = "A")}. Fatores nao fixados sao promediados com pesos iguais.}
}
\details{
Modelos com efeitos fixos e todas as interacoes fatoriais; blocos e linhas/colunas
entram como efeitos aditivos. Requerem erros independentes, normais e de variancia
constante. Inspecione os residuos com \code{plot(analise$modelo)}.
DIC requer repeticoes iguais por combinacao. DBC requer uma observacao de cada
combinacao por bloco. DQL requer quadrado latino completo de ordem ao menos 3:
a ordem e o produto dos numeros de niveis dos fatores, e cada combinacao ocorre
uma vez por linha e coluna. Nao se trata de delineamento greco-latino.

Sem \code{dentro}, o efeito e marginal. Fixar outros fatores permite efeitos
simples, inclusive um fator dentro de uma combinacao dos outros dois.
Uma interacao significativa a 5 por cento envolvendo fator nao fixado gera
aviso. O aviso nao seleciona automaticamente uma analise; interprete a ANOVA
antes de escolher os desdobramentos. \code{desdobrar} testa por F a igualdade das
medias selecionadas com o erro experimental do modelo completo.

Os cinco testes e os desdobramentos usam os graus de liberdade e o quadrado
medio residual do modelo fatorial completo, sem reestimar o erro em subconjuntos.
O ajuste por multiplicidade se refere a cada familia selecionada em uma chamada;
nao ha ajuste conjunto entre multiplos desdobramentos. As comparacoes nao sao
condicionadas ao teste F. Duncan nao garante o controle de erro familiar de Tukey.
Scheffe e Duncan requerem agricolae; Dunnett requer multcomp; Scott-Knott usa
ExpDes.pt e devolve tambem sua tabela textual. Tukey usa a amplitude studentizada
sobre medias estimadas pelo modelo. Nos delineamentos completos balanceados
aceitos, as medias marginais e simples coincidem com as medias aritmeticas.

A regressao exige o fator original numerico. Seus coeficientes selecionados sao
estimados das medias marginais ou simples, com covariancia e testes t usando o
erro experimental completo. A falta de ajuste testa esse polinomio contra as
medias selecionadas; retorna NULL quando nao ha graus de liberdade de falta de
ajuste. O componente modelo e um lm global reduzido, com o polinomio interagindo
com os demais fatores e mantendo blocos ou linhas/colunas. Sua ANOVA e o seu
summary usam o residuo reduzido e nao substituem a tabela coeficientes e o teste
falta_de_ajuste, que sao a inferencia do efeito selecionado com erro puro completo.
}
\value{
Os delineamentos retornam objeto exp_anova com modelo aov, ANOVA, dados, mapa de
fatores e medias por combinacao. comparar retorna o objeto do teste; Tukey retorna
lista com matriz de diferencas, limites e p ajustados; Dunnett retorna teste e
intervalos; Scott-Knott retorna resultado e saida textual. desdobrar retorna medias
e tabela de F com GL, GL_erro, QM_erro, F e p_valor. regressao retorna modelo lm
global, coeficientes do efeito selecionado, covariancia, ANOVA global reduzida,
falta de ajuste, grau, fator, dentro, erro experimental e medias utilizadas.
}
\examples{
set.seed(42)
d <- expand.grid(dose = 0:3, cultivar = c("A", "B"), bloco = 1:4)
d$y <- 5 + 2 * d$dose + 3 * d$dose * (d$cultivar == "B") + rnorm(nrow(d))
a <- dbc(d, "y", c("dose", "cultivar"), "bloco")
a
desdobrar(a, fator = "dose", dentro = list(cultivar = "A"))
comparar(a, "tukey", fator = "dose", dentro = list(cultivar = "A"))
regressao(a, grau = 2, fator = "dose", dentro = list(cultivar = "A"))$coeficientes
}
