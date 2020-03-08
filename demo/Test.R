load("~/Google Drive/GPstuff-4.6/int_2d.RData")

library(mlOSP)

model2d <- list(x0 = rep(40,2),K=40,sigma=rep(0.2,2),r=0.06,div=0,T=1,dt=0.04,dim=2,
                sim.func=sim.gbm, payoff.func=put.payoff)
MM <- (model2d$T/model2d$dt)
model2d$pilot.nsims <- 1000
model2d$look.ahead <- 1

library(randtoolbox)
sob30 <- sobol(55, d=2,scrambling=1)
sob30 <- sob30[ which( sob30[,1] + sob30[,2] <= 1) ,]  # a lot are on the diagonal
model2d$cand.len <- 1000
model2d$max.lengthscale <- c(20,20)
model2d$min.lengthscale <- c(0.2,0.2)
model2d$init.size <- 10   # initial design size

model2d$init.grid <- int_2d
model2d$tmse.eps <- 0
model2d$kernel.family <- "Gaussian"
model2d$ucb.gamma <- 1.96

model2d$seq.design.size <- 100
model2d$batch.nrep <- 20
model2d$ei.func <- "amcu"
library(mlOSP)
library(laGP)
set.seed(101); 
oos.obj <- osp.seq.design(model2d,method="homtp")
oos.obj$fit[[15]]$X0
library(ggplot2)
plt.2d.surf(oos.obj$fit[[15]])

xt <- seq(25, 50, len = 50)
xt1 <- rep(xt, 50)
xt2 <- rep(xt, each = 50)
xt <- data.frame(xt = xt1, x2 = xt2)

pred <- predict(x=cbind(xt1,xt2), object=oos.obj$fit[[15]])
pred_dat <- data.frame(xt1 = xt1, xt2 = xt2, est = pred$mean, var = pred$sd2)
designs <- data.frame(x = oos.obj$fit[[15]]$x, y = oos.obj$fit[[15]]$y)
# write.csv(pred_dat, file = "~/Dropbox/GP1_Plots/tp_bermudan_pred.csv", row.names = FALSE)
# write.csv(designs, file = "~/Dropbox/GP1_Plots/tp_bermudan_designs.csv", row.names = FALSE)

############################
load("~/Dropbox/int_2d.RData")
model2d$init.size <- 20   # initial design size

model2d$kernel.family <- "gauss"
model2d$ucb.gamma <- 1.96

model2d$seq.design.size <- 100
model2d$batch.nrep <- 20
model2d$ei.func <- "amcu"

set.seed(101)
oos.obj <- osp.seq.design(model2d,method="trainkm")

plt.2d.surf(oos.obj$fit[[15]])
