## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ks)
library(fields) # for plotting purposes
library(mlOSP)
library(DiceKriging)
library(tgp)  # use lhs from there
library(randtoolbox)  # use sobol and halton QMC sequences
library(hetGP)
library(laGP)
library(ggplot2)
library(pander)
data("int_2d")

## -----------------------------------------------------------------------------
model2d <- list(x0 = rep(40,2),K=40,sigma=rep(0.2,2),r=0.06,div=0,T=1,dt=0.04,dim=2,sim.func=sim.gbm, payoff.func=put.payoff)
model2d$pilot.nsims <- 1000
model2d$look.ahead <- 1
model2d$cand.len <- 1000   # size of candidate set m_0 for acquisition function
model2d$max.lengthscale <- c(50,50)
model2d$min.lengthscale <- c(5,5)
model2d$tmse.eps <- 0

model2d$ucb.gamma <- 1.96

model2d$seq.design.size <- 100  # budget N = r_0 * k = 2500
model2d$batch.nrep <- 25  # initial replication r_0
model2d$total.budget <- 2500

model2d$init.size <- 20   # initial design size k_0
model2d$init.grid <- int_2d

model2d$tmse.eps <- 0
model2d$kernel.family <- "Matern5_2"  # kernel function for Gaussian Process
model2d$ucb.gamma <- 1.96
model2d$update.freq <- 10   # number of sequential design steps to update the GP surrogate
model2d$r.cand <- c(20, 30,40,50,60, 80, 120, 160) # r_L

## ----testing-homtp2d, message=FALSE, warning=FALSE, fig.width=6---------------
### GP + ADSA
set.seed(110)
model2d$batch.heuristic <- 'adsa'
model2d$ei.func <- "amcu"
oos.obj.adsa <- osp.seq.batch.design(model2d, method="hetgp")

### GP + ABSUR
set.seed(122)
model2d$batch.heuristic <- 'absur'
model2d$ei.func <- 'absur'
oos.obj.absur <- osp.seq.batch.design(model2d, method="hetgp")

### plot Figure 6
plt.2d.surf.batch(oos.obj.adsa$fit[[15]], oos.obj.absur$fit[[15]], oos.obj.adsa$batches[1:oos.obj.adsa$ndesigns[15] - 1, 15], oos.obj.absur$batches[1:oos.obj.absur$ndesigns[15] - 1, 15], "ADSA", "ABSUR", x=seq(25,50,len=201),y = seq(25, 50,len=201))

## ---- message=FALSE, warning=FALSE, fig.height=3.5, fig.width=6, fig.cap="Timing Value (background color) and Exercise Boundary (zero-contour) at $t=0.6$ using ADSA"----
oos.obj.adsa$ndesigns[15]  # number of unique designs
### plot Figure 6 right panel - ADSA
plt.2d.surf.with.batch(oos.obj.adsa$fit[[15]], 
                       oos.obj.adsa$batches[1:oos.obj.adsa$ndesigns[15] - 1, 15])

## ---- message=FALSE, warning=FALSE,fig.height=3.5, fig.width=6, fig.cap="Timing Value and Exercise Boundary using ABSUR"----
oos.obj.absur$ndesigns[15]  # number of unique designs
### plot Figure 6 left panel - ABSUR
plt.2d.surf.with.batch(oos.obj.absur$fit[[15]], 
                       oos.obj.absur$batches[1:oos.obj.absur$ndesigns[15] - 1, 15])

