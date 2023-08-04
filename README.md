# Bayesian-Inference

---
params: 
    set_title: "Introduction to Bayesian Inference"
    set_subtitle: "Group 13"
title: "`r params$set_title`"
subtitle: "`r params$set_subtitle`"
author: |
  | Group members:
  | Calo, Theresse Joy
  | Colman, Roos
  | Merezhko, Maria
  | Noe, Sebastian
  | Ormando Aramburu, Santiago
date: "Due date: 01.06.2022"
output:
  bookdown::pdf_document2:
    author_block: TRUE
    toc: FALSE
    number_sections: FALSE
    fig_caption: yes
    tab_caption: yes
header-includes:
  - \usepackage{authblk, mdframed}
  - \usepackage{fancyhdr}
  - \usepackage[official]{eurosym}
  - \usepackage{titling}
  - \usepackage{xcolor}
  - \pretitle{\begin{flushleft}\LARGE}
  - \posttitle{\end{flushleft}\LARGE}
  - \preauthor{\begin{flushleft}}
  - \postauthor{\end{flushleft}} 
  - \predate{\begin{flushleft}\tiny}
  - \postdate{\end{flushleft}\tiny}
  - \usepackage{helvet}
  - \renewcommand{\familydefault}{\sfdefault}
  - \pagestyle{fancy}
  - \fancyhead[LO,LE]{`r params$set_subtitle`}
  - \usepackage{setspace}\doublespacing

---

(ref:rho) $\rho$

\newpage

```{r setup, include=FALSE}
knitr::opts_chunk$set(include=TRUE, echo = FALSE)

library(coda)

library(DescTools)
library(dplyr)

library(ggplot2)
library(gridExtra)

library(kableExtra)

library(nimble)

library(rootSolve)

```

```{r functions, inlcude=F}

pairwise.correlation = function(Arg1, Arg2){
  result = cor.test(Arg1, Arg2, method="spearman")
  est = result[[4]]
  if(est < 0.01){est = "<0.01"} else {est = paste0("=", round(est,2))}
  p = result$p.value
  if(p < 0.001){p = "<0.001"} else {p = paste0("=", round(p,3))}
  paste0("(ref:rho)",est," (p", p, ")")
}

```

# Exercise No. 4

## 4.1 Model the data as using a nonlinear logistic growth function. Take flat uniform priors for the parameters $\alpha$, $\beta$, $\gamma$ and a noninformative prior for $\sigma$. Write out the Bayesian model (likelihood and priors).

### Likelihood

$$ lik(..) = \prod_i^n \frac{\alpha}{1+\beta^{\gamma \cdot x_i}} + \epsilon_i $$

### Priors and prior assumptions

$$ \alpha = \frac{1}{b-a} $$

$$ \beta = \frac{1}{b-a} $$

$$ \gamma = \frac{1}{b-a} $$

$$ \sigma \sim N \{0; 100\} $$

## 4.2 Run 5000 iterations, then look at history plots and autocorrelation plots of the sample traces and calculate the Gelman and Rubin convergence diagnostic for each of the parameters you have monitored. Do the simulations look like they have converged? If not, carry out some more updates and check again, until you are happy with convergence.

```{r, include=F}
set.seed(123)

library('nimble')


x = c( 1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0,
       7.0, 8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
       13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5)
N = 27

Y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
      2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43, 2.47,
      2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57)

model4.data <- list('x' = x, 'Y'=Y)
model4.constant <- list('N' = N)

# DEFINE INITIAL VALUES

model4.inits <- list(sigma=1, alpha=0, beta=0, gamma=0 )


# MODEL SPECIFICATION RC changed this to original to obtain curve fit

Model4.code <- nimbleCode(
  {
    # likelihood specification
    for (i in 1:N) {
      Y[i] ~  dnorm( alpha / (1 + beta*exp(gamma *x[i])),tau) 
    }
    # prior information
    alpha ~ dunif(-10,10)
    beta ~ dunif(-10,10)
    gamma ~ dunif(-10,10)
    sigma ~ dflat()
    tau <- 1/(sigma^2)
  })

# specify model, data, number of parallel chains
Model4 <- nimbleModel(Model4.code, 
                     constants=model4.constant,
                     data=model4.data,
                     inits=model4.inits)


Model4.compiled <- compileNimble(Model4)

Model4a.output = nimbleMCMC(Model4, 
                           niter=6000,
                           nburnin=1000,
                           thin=10, 
                           summary=T)
Model4a.output$summary

#code to check observed vs predicted
pred <- Model4a.output$summary[1,1] / (1 + Model4a.output$summary[2,1] *exp(Model4a.output$summary[3,1]*x))
plot(x, Y)
lines(x,pred)

```

RC Since no information was available on the expected values for alpha, beta and gamma, we initially started with a uniform distribution from -10 to +10 for each of the three parameters. Since the estimated values for these parameters as well as the 95% credible intervals are far from -10 and +10, it was decided that the initial values and the distributions were appropriate.   

```{r GuR, include=F, echo=F, message=F, error=F}
set.seed(1)
Model4a.GuR.output = nimbleMCMC(Model4, niter=6000, nburnin=1000, thin=1, nchains=5)#RC changed this to 1 to ease the interpretation of the Gelman and Rubin plots. Link: plot.https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/
#RC changed niter=5000 to niter:6000 to be consistent with the model above

Model4a.GuR.List = mcmc.list(list(
  as.mcmc(Model4a.GuR.output$chain1),
  as.mcmc(Model4a.GuR.output$chain2),
  as.mcmc(Model4a.GuR.output$chain3),
  as.mcmc(Model4a.GuR.output$chain4),
  as.mcmc(Model4a.GuR.output$chain5)
))

```

Traceplots and autocorrelation plots for the parameters in the model can be found in figures \@ref(fig:traceplot4) and \@ref(fig:acf4), respectively. Additionally, the Gelman and Rubin convergence diagnostic for each of the parameters can be found in table \@ref(tab:Model4GuRTable), while the mulivariate potential scale reduction factor was found to be `r round(gelman.diag(Model4a.GuR.List)[[2]],3)`; the corresponding plots for the Gelman and Rubin diagnostic each of the parameters are displayed in figure \@ref(fig:Model4GuRPlot). RC Gelman and Rubin convergence diagnostic plots showed that a burn-in of 4000 is a reasonable choice, since scale-reduction is sufficient before 4000 iterations were reached. For that reason, the total number of iterations was increased to 9000, with a burn-in of 4000 iterations.  






```{r , include=F, echo=F, message=F, error=F}
#RC
set.seed(1)
Model4.GuR.output = nimbleMCMC(Model4, niter=9000, nburnin=4000, thin=1, nchains=5)

Model4.GuR.List = mcmc.list(list(
  as.mcmc(Model4.GuR.output$chain1),
  as.mcmc(Model4.GuR.output$chain2),
  as.mcmc(Model4.GuR.output$chain3),
  as.mcmc(Model4.GuR.output$chain4),
  as.mcmc(Model4.GuR.output$chain5)
))

```

RC After this adjustment, the multivariate potential scale reduction factor was found to be `r round(gelman.diag(Model4.GuR.List)[[2]],3)`.
The results reported in the remainder of this document will be based on 9000 iterations of which 4000 are burn-ins.

```{r GuR, include=F, echo=F, message=F, error=F}

set.seed(1)
Model4.output = nimbleMCMC(Model4, 
                           niter=9000,
                           nburnin=4000,
                           thin=10, 
                           summary=T)

Model4.output$summary

```

## 4.3 Produce summary statistics and kernel density plots for the posterior samples of the regression parameters. Check the Monte Carlo standard error for each of the parameters to assess the accuracy of your estimates.

Summary statistics of the parameters for the posterior distribution can be found in table \@ref(tab:postSummary4) while the plots for the density of the parameters of the posterior distribution are displayed in figure \@ref(fig:densplot4).

## 4.4 Are the parameters a posteriori correlated? Check by making pairwise scatterplots of the MCMC samples.

A visual exploration of correlation between the parameters of the posterior distribution can be found in figure \@ref(fig:corrPlots4). Spearman's correlations ($\rho$) and corresponding p-values for pairwise correlations were `r pairwise.correlation(Model4.output$samples[,1], Model4.output$samples[,2])` for $\alpha$ and $\beta$, `r pairwise.correlation(Model4.output$samples[,1], Model4.output$samples[,3])` for $\alpha$ and $\gamma$, `r pairwise.correlation(Model4.output$samples[,1], Model4.output$samples[,4])` for $\alpha$ and $\sigma$, `r pairwise.correlation(Model4.output$samples[,2], Model4.output$samples[,3])` for $\beta$ and $\gamma$, `r pairwise.correlation(Model4.output$samples[,2], Model4.output$samples[,4])` for $\beta$ and $\sigma$, and `r pairwise.correlation(Model4.output$samples[,3], Model4.output$samples[,4])` for $\gamma$ and $\sigma$.

```{r corMatrix4, include=F}

cor(as.data.frame(Model4.output$samples), method="spearman") %>%
  kable(align="cccc", linesep="", booktabs=T, caption = "Correlation matrix between the parameters of the posterior distribution.") %>%
  kable_styling(bootstrap_options = c("condensed")) %>% 
  column_spec(1, bold=T)

```


## 4.5 Check the sensitivity of the priors. Change the priors by increasing or decreasing the variance. Does this affect the results?

```{r , include=F}
set.seed(1)
Model4.lo.code <- nimbleCode(
  {
    # likelihood specification
    for (i in 1:N) {
      Y[i] ~  dnorm( alpha / (1 + beta*exp(gamma *x[i])),tau) 
    }
    # prior information
    alpha ~ dunif(-5,5)
    beta ~ dunif(-5,5)
    gamma ~ dunif(-5,5)
    sigma ~ dunif(0,100)
         tau <- 1/(sigma^2)
   #RC changed this because sigma ~ dnorm(0, 10) seems not a realistic distribution of sigma (mean zero and 50% probability of having a negative sigma) moreover growth curve plot looked awful. I used example code from the course (example1_nimble)
  })

Model4.lo = nimbleModel(Model4.lo.code,
                     constants=model4.constant,
                     data=model4.data,
                     inits=model4.inits)

Model4.lo.output = nimbleMCMC(Model4.lo, niter=9000, nburnin=4000, thin=10, summary=T)

#RC code to check observed vs predicted - can be deleted
pred <- Model4.lo.output$summary[1,1] / (1 + Model4.lo.output$summary[2,1] *exp(Model4.lo.output$summary[3,1]*x))
plot(x, Y)
lines(x,pred)

set.seed(1)
Model4.hi.code <- nimbleCode(
  {
    # likelihood specification
    for (i in 1:N) {
      Y[i] ~  dnorm( alpha / (1 + beta*exp(gamma *x[i])),tau) 
    }
    # prior information
    alpha ~ dunif(-100,100)
    beta ~ dunif(-100,100)
    gamma ~ dunif(-100,100)
    sigma ~ dunif(0,10000)
         tau <- 1/(sigma^2)
    #RC changed this because sigma ~ dnorm(0, 10) seems not a realistic distribution of sigma (mean zero and 50% probability of having a negative sigma) moreover growth curve plot looked awful
  })

Model4.hi = nimbleModel(Model4.hi.code,
                     constants=model4.constant,
                     data=model4.data,
                     inits=model4.inits)

Model4.hi.output = nimbleMCMC(Model4.hi, niter=9000, nburnin=4000, thin=10, summary=T)

#RC code to check observed vs predicted#RC code to check observed vs predicted - can be deleted
pred <- Model4.hi.output$summary[1,1] / (1 + Model4.hi.output$summary[2,1] *exp(Model4.hi.output$summary[3,1]*x))
plot(x, Y)
lines(x,pred)


```

A sensitivity analysis of changes of the variance for the prior assumptions on the parameter estimates can be found in table \@ref(tab:SensSigma). RC estimates did not change considerably with changing assumptions for the priors, which is desirable.    


```{r SensSigma, include=T, echo=FALSE, error=FALSE, message=FALSE, warning=FALSE}

rbind(Model4.lo.output$summary,
      Model4.output$summary,
      Model4.hi.output$summary) %>%
  "colnames<-"(c("mean", "median", "SD", " ", " ")) %>%
  kable(align="ccccc", linesep ="", booktabs=T, caption = "Summary statistics for the parameters of the posterior distribution.") %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  pack_rows("Reduced variance", 1, 4) %>%
  pack_rows("Original variance", 5, 8) %>%
  pack_rows("Increased variance", 9, 12) %>%
  column_spec(1, bold=T) %>%
  add_header_above(c(" "=4, "95% credible interval"=2))


```


## 4.6 Make a plot of the posterior growth curve. Visualize also the uncertainty of the plotted curve. How does this compare to the observed data?  

```{r Model4fit}

alpha = Model4.output$summary[[1,1]]
beta = Model4.output$summary[[2,1]]
gamma = Model4.output$summary[[3,1]]

f = function(x) alpha/(1 + beta*exp(gamma*x))

ObsPred = as.data.frame(
  cbind(x, Y)
)

ObsPred = mutate(ObsPred, Yobs = 0)
for(i in 1:nrow(ObsPred)){
  ObsPred[i,"Yobs"] = alpha/(1+beta*exp(gamma * ObsPred[i,"x"]))
}

ggplot() +
  geom_point(aes(x=x, y=Y)) +
  geom_function(fun = f)

```

## 4.7 Give the posterior predictive distribution of the length of a new animal at age 25.  

```{r }

x = c( 1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0,
       7.0, 8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
       13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5)
N = 27

Y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
      2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43, 2.47,
      2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57)

model4.data <- list('x' = x, 'Y'=Y)
model4.constant <- list('N' = N)

model4.inits <- list(sigma=1, alpha=0, beta=0, gamma=0 )

# MODEL SPECIFICATION 
Model4.code <- nimbleCode(
  {
    # likelihood specification
    for (i in 1:N) {
      Y[i] ~  dnorm( alpha / (1 + beta*exp(gamma *x[i])),tau)
      Y25 ~ dnorm(alpha/(1+beta*exp(gamma*25)),tau)
    }
    # prior information
    alpha ~ dunif(-10,10)
    beta ~ dunif(-10,10)
    gamma ~ dunif(-10,10)
    sigma ~ dflat()
    tau <- 1/(sigma^2)
  })

# SET UP MODEL OBJECT 

Model4 <- nimbleModel(Model4.code, 
                     constants=model4.constant,
                     data=model4.data,
                     inits=model4.inits)
Model4.compiled <- compileNimble(Model4)


# Generate MCMC samples and save output for specified variables
out <- nimbleMCMC(Model4,
                  niter=6000,nburnin=1000,thin=10,
                  summary=TRUE)

# Posterior summary statistics
out$summary

#POSTERIOR DISTRIBUTION
library(coda)
densplot(as.mcmc(out$samples))


```

\newpage
# Exercise No. 5

## Explain the model that is given in the code. What is the likelihood (write it out)? What prior is assumed for the parameters?

This model assumes that the number of people demonstrating HIV-seroconversion within this study ($r_i$) follows a binomial distribution with parameters $\theta_i$ for people randomized to either the vaccine ($i=1$) or placebo ($i=2$), with a prior distribution of $\theta$ described by a $\beta$-distribution with parameters $\alpha$=$\beta$=0.5.  
The likelihood function is described by:
$$ lik(\theta_i|n) = \binom{n_i}{r_i} r_i^{\theta_i} \cdot (1-\theta_i)^{n_i-r_i} $$


## Use 10,000 MCMC iterations with a burn-in of 1,000 MCMC iterations for estimation of the parameters. What is the MCMC error? Check convergence of the MCMC chain for all parameters. Discuss.

```{r, include=F}

n = c(8197, 8198)
N = sum(n)
r = c(51, 74)
R = sum(r)

alpha = 0.5
beta = 0.5
theta = c(r[1]/n[1], r[2]/n[2])

Model5.constants = list("n" = n)
Model5.data = list("r" = r)
Model5.inits = list("theta" = theta)

Model5.Code = nimbleCode({
  for (i in 1:2)
    {
      r[i] ~ dbin(theta[i],n[i])
      theta[i] ~ dbeta(0.5,0.5)
  }
})


Model5 = nimbleModel(Model5.Code, 
                     constants=Model5.constants,
                     data=Model5.data,
                     inits=Model5.inits)

Model5.compiled = compileNimble(Model5)

Model5.out = nimbleMCMC(Model5,
                  niter=10000,
                  nburnin=1000,
                  thin=1,
                  summary=T)

rm(Model5.compiled)
```


## Give the posterior mean and median for the proportions. Give also 95 credible intervals.

Posterior mean and median together with a 95% credible interval can be found in table.

###  Obtain a 95 credible interval for VE. How can we interpret this interval?

```{r, include=F}

VE.matrix = matrix(nrow=nrow(Model5.out$samples))

for(i in 1:nrow(VE.matrix)){
  VE.matrix[i] = as.numeric(round(1 - Model5.out$samples[i,1]/Model5.out$samples[i,2],3))
}

VE.matrix = VE.matrix %>% 
  as.data.frame() %>%
  "colnames<-"(c("V"))

CI_lo = round(quantile(VE.matrix$V, 0.025),3)
CI_hi = round(quantile(VE.matrix$V, 0.975),3)

VE.sample = mcmc(data=VE.matrix)
#HPDinterval(VE.sample, alpha=0.05)

```

A 95% equal-tail credible interval from the distribution of VE is [`r CI_lo`; `r CI_hi`]. This interval contains 95% of the a posteriori most probable values for VE. Being an equal-tail interval, there are, however, VE-values outside this credible interval, that are more likely than some parameters inside this interval and vice versa. The more narrow 95% highest-posterior-density interval is given by [`r HPDinterval(VE.sample, alpha=0.05)[[1]]`; `r HPDinterval(VE.sample, alpha=0.05)[[2]]`].

### Give a plot of the posterior density of VE

Plot for the posterior density of VE can be found in figure \@ref(fig:posteriorPlotVE).

```{r, include=F}

Gruppenweite = 0.01
minVE = as.numeric(min(VE.matrix$V))

VE.frequencies = as.data.frame(table(VE.matrix))
VE.frequencies = mutate(VE.frequencies, Prop = as.numeric(VE.frequencies$Freq/9000))
VE.frequencies = mutate(VE.frequencies, group = minVE - as.numeric(VE.frequencies$V))

VE.frequencies$V = as.numeric(VE.frequencies$V)


```

## What is the posterior probability that the proportion $\theta_1$ of infected amongst those vaccinated is actually smaller than the proportion $\theta_2$ of infected among those that take the placebo?

```{r, include=F}

VE.neg = VE.matrix %>% subset(V < 0) %>% table() %>% sum()

Prop.VE.neg = VE.neg/nrow(VE.matrix)

```

The posterior probability of $\theta_1 < \theta_2$ equals the probability of VE < 0, corresponding to `r round(Prop.VE.neg,3)`.

## Give a plot of the posterior predictive distribution for a future 100 subjects that are vaccinated, as well as the posterior predictive distribution for a future 100 subjects that take the placebo. Give also a 95% posterior predictive set.

## Do a sensitivity analysis of the prior. How do results change in case you use the following priors:

### Use a non-informative prior for both $\theta_1$ and $\theta_2$.

```{r, include=F, message=FALSE, warning=FALSE}

Model.Sens1.Code = nimbleCode({
  for (i in 1:2)
    {
      r[i] ~ dbin(theta[i],n[i])
      theta[i] ~ dbeta(1,1)
  }
})

Model.Sens1 = nimbleModel(Model.Sens1.Code, 
                     constants=Model5.constants,
                     data=Model5.data,
                     inits=Model5.inits)

Model.Sens1.compiled = compileNimble(Model.Sens1)

Model.Sens1.out = nimbleMCMC(Model.Sens1,
                  niter=10000,
                  nburnin=1000,
                  thin=1,
                  summary=T)



```

Results of the posterior summaries for the estimates for both parameters under initial assumptions as well as the assumption of weakly informative (flat) priors for the distribution or $\theta_1$ and $\theta_2$ can be found in table \@ref(tab:Model5Sens1). Of notice, the effects of the changes in the prior are negligible.  

### Use a prior that reflects a previous study ($r_1$=56, $n_1$=8,202 for people who have taken the vaccine; $r_2$=76, $n_2$=8,200 for people who have taken placebo).

```{r, include=F, message=FALSE, warning=FALSE}

n = c(8202, 8200)
N = sum(n)
r = c(56, 76)
R = sum(r)

model.constants = list("n" = n)
model.data = list("r" = r)
model.inits = list("theta" = theta)

Model.Sens2.Code = nimbleCode({
  for (i in 1:2)
    {
      r[i] ~ dbin(theta[i],n[i])
      theta[i] ~ dbeta(1,1)
  }
})

Model.Sens2 = nimbleModel(Model.Sens2.Code, 
                     constants=model.constants,
                     data=model.data,
                     inits=model.inits)

Model.Sens2.compiled = compileNimble(Model.Sens2)

Model.Sens2.out = nimbleMCMC(Model.Sens2,
                  niter=10000,
                  nburnin=1000,
                  thin=1,
                  summary=T)

rm(Model.Sens2)
rm(Model.Sens2.compiled)

```

Results of the posterior summaries for the estimates for both parameters under initial assumptions as well as the assumptions derived from a previous study for the distribution or $\theta_1$ and $\theta_2$ can be found in table \@ref(tab:Model5Sens2). Again, the effects of the changes in the prior are negligible. 

\newpage

# Figures

\newpage

```{r traceplot4, fig.cap="Traceplots of the parameters alpha, beta, gamma, and sigma after MCMC-sampling from the posterior distribution with 6,000 iteraions and a burn-in of 1,000."}

par(mfrow=c(2,2))

traceplot(as.mcmc(Model4a.output$samples))

```


```{r traceplot4, fig.cap="Traceplots of the parameters alpha, beta, gamma, and sigma after MCMC-sampling from the posterior distribution with 9,000 iteraions and a burn-in of 4,000."}

par(mfrow=c(2,2))

traceplot(as.mcmc(Model4.output$samples))

```

\newpage

```{r acf4, fig.cap="Plots of autocorrelation between the parameters in the model. Results for 6000 iterations of which 1000 burn-ins", fig.width=10}

acf(as.mcmc(Model4a.output$samples))

```

```{r acf4, fig.cap="Plots of autocorrelation between the parameters in the model with 9000 iterations and 4000 burn-ins.", fig.width=10}

acf(as.mcmc(Model4.output$samples))

```

\newpage

```{r Model4GuRPlot, fig.cap="Gelman and Rubin convergence diagnostic plots: 6000 iterations of which 1000 burn-ins."}

gelman.plot(Model4a.GuR.List)
```



```{r Model4v2GuRPlot, fig.cap="Gelman and Rubin convergence diagnostic plots with 9000 iterations in total, of which 4000 burn-ins."}
gelman.plot(Model4.GuR.List)
```

\newpage

```{r densplot4, fig.cap="Density plot of the posterior parameters.", fig.width=10}

par(mfrow=c(2,2))
densplot(as.mcmc(Model4.output$samples))

```

\newpage

```{r corrPlots4, include=T, fig.cap="Pairwise correlations between the parameters of the posterior distribution."}

par(mfrow=c(3,2))

p1 = ggplot(as.data.frame(Model4.output$samples)) + 
  geom_point(aes(x=alpha, y=beta)) +
  labs(x="alpha", y="beta") + theme_bw()
p2 = ggplot(as.data.frame(Model4.output$samples)) + 
  geom_point(aes(x=alpha, y=gamma))+
  labs(x="alpha", y="gamma") + theme_bw()
p3 = ggplot(as.data.frame(Model4.output$samples)) + 
  geom_point(aes(x=alpha, y=sigma))+
  labs(x="alpha", y="sigma") + theme_bw()
p4 = ggplot(as.data.frame(Model4.output$samples)) + 
  geom_point(aes(x=beta, y=gamma))+
  labs(x="beta", y="gamma") + theme_bw()
p5 = ggplot(as.data.frame(Model4.output$samples)) + 
  geom_point(aes(x=beta, y=sigma))+
  labs(x="beta", y="sigma") + theme_bw()
p6 = ggplot(as.data.frame(Model4.output$samples)) + 
  geom_point(aes(x=gamma, y=sigma))+
  labs(x="gamma", y="sigma") + theme_bw()

grid.arrange(p1,p2,p3,p4,p5,p6, nrow=3, ncol=2)

```




\newpage

```{r posteriorPlotVE, include=T, fig.cap="Caption to follow"}

ggplot(data=VE.matrix, aes(x=V)) +
  geom_histogram(aes(y = after_stat(count/ sum(count))), binwidth = 0.02) +
  theme_bw()

```

\newpage

# Tables

\newpage

```{r Model4GuRTable, includ=T, echo=F, message=F, error=F}
gelman.diag(Model4.GuR.List)[1] %>%
  kable(linesep="", booktabs = TRUE, align = "cc", digits = 3, 
        caption="Gelman and Rubin convergence diagnostics for each parameter in the model.") %>%
  kable_styling(bootstrap_options = c("condensed"), latex_options = "hold_position")
```

\newpage

```{r postSummary4}

Model4.output$summary %>%
  "colnames<-"(c("mean", "median", "SD", " ", " ")) %>%
  kable(align="ccccc", linesep ="", booktabs=T, caption = "Summary statistics for the parameters of the posterior distribution.") %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  column_spec(1, bold=T) %>%
  add_header_above(c(" "=4, "95% credible interval"=2))

```



\newpage

```{r Model5Summary}

Summary = as.matrix(Model5.out$summary[,c(1,2,4,5)]) %>%
  "colnames<-"(c("mean", "median", " ", " ")) 

kable(Summary, align="cccc", linesep ="", booktabs=T, caption="Means and medians with 95 credible intervals for the parameters of the posterior distribution") %>%
  kable_styling(bootstrap_options = c("condensed"), latex_options = c("hold_position")) %>%
  column_spec(1, bold=T) %>%
  add_header_above(c(" "=3, "95% credible interval"=2))

```

\newpage

```{r Model5Sens1}

rbind(
  Model5.out$summary,
  Model.Sens1.out$summary) %>% 
  kable(align="cccc", linesep="", booktabs=T, 
        caption="Results for sensitivity analysis on the priors, assuming no prior knowledge or deciding againgst including prior knowledge") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  pack_rows("Original assumptions", 1, 2) %>%
  pack_rows("Assumption of non-informative priors for theta", 3, 4)

```

\newpage

```{r Model5Sens2}

rbind(
  Model5.out$summary,
  Model.Sens1.out$summary) %>% 
  kable(align="cccc", linesep="", booktabs=T, 
        caption="Assuptions from results of previous studies") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  pack_rows("Original assumptions", 1, 2) %>%
  pack_rows("Assumption of information from previous study", 3, 4)

```
