 knitr::opts_chunk$set(echo = TRUE, fig.path='figures/', dev=c('png', 'pdf'), dpi=100, fig.width=6.4, fig.height=3.6,cache=TRUE,warning=FALSE,message=FALSE)

library(ks)
library(fields) # for plotting purposes, use quilt.plot in 2D
library(mlOSP)
library(DiceKriging)
library(tgp)  # use lhs from there
library(randtoolbox)  # use sobol and halton QMC sequences
library(randomForest)
library(earth)
library(hetGP)
library(kernlab)
library(ggplot2)
library(kableExtra)
library(gridExtra)
library(mvtnorm)
library(nnet)
library(np)
library(laGP)  # for distance function

set.seed(262)
lsModel <- list(dim=1,r=0.05,sigma=0.2,div=0,dt=0.1)
rmc <- list(); paths <- list(); payoff <- list()

x0 <- 40;  Np <- 20; Nt <- 5; 
paths[[1]] <- sim.gbm( matrix(rep(x0, Np), nrow=Np,byrow=T), lsModel, lsModel$dt)
for (j in 2:(Nt+1))
    paths[[j]] <- sim.gbm( paths[[j-1]], lsModel, lsModel$dt)

payoff[[4]] <- pmax(40-paths[[5]],0)  # outputs Y(x). Inputs are paths[[4]] (X(4))
rmc[[4]] <- lm(y ~ poly(x,2,raw=TRUE), data.frame(y= payoff[[4]], x = paths[[4]]))


testx <- seq(25,58,len=50); 
predy <- predict.lm(rmc[[4]], new=data.frame(x=testx))

g1.base <- ggplot(data.frame(y= payoff[[4]], x = paths[[4]]), aes(x, y)) + geom_point(size=2.5) +
   theme_light() +  geom_smooth(method="lm", formula=y~poly(x, 2), size=1.5, fill=NaN) +
   labs(subtitle=paste("k=", 4, ": quadratic fit is ", round(coef(rmc[[4]])[1],3),
                       round(coef(rmc[[4]])[2],3),"x+", 
                       round(coef(rmc[[4]])[3],3), "x^2")) +
  xlab("X(4)=x") + ylab("Y(x)") + theme(plot.subtitle=element_text(size=8)) +
  scale_x_continuous(expand=c(0.01,0.01))  + scale_y_continuous(limits=c(-1,13)) 

# Show the immediate payoff and predicted cont
to.plot <- data.frame(contValue=pmax(40-testx, 0), payoff=exp(-lsModel$r*lsModel$dt)*predy,
                      valueFunc=pmax(exp(-lsModel$r*lsModel$dt)*predy,0 ),x=testx)
#plot(testx, ,lwd=3, col="red", xlab="S", ylab="Cont Value", cex.lab=1.3,type="l")
#lines(testx, , lwd=3, lty=2, col="blue")
#lines(testx, pmax(exp(-lsModel$r*lsModel$dt)*predy,0 ), lwd=3, col="blue")
find_stopBnd <- function(guessx, obj, disc) {
  predy <- predict.lm(obj, new=data.frame(x=guessx))*disc - (40-guessx)
  return(predy)
}
stopBnd <- array(0, dim=c(4,1))
stopBnd[4] <- uniroot(find_stopBnd, c(33,40), rmc[[4]], exp(-lsModel$r*lsModel$dt))$root
#lines(, c(-0.4, -0.4), col="green", lwd=5)

g2.base <- ggplot(to.plot, aes(x=x)) + 
  geom_line(aes(y=contValue, color="q"), size=1.5) +
  geom_line(aes(y=payoff, color="h"), size=1.5, linetype="twodash") +
  geom_line(data=data.frame(x=c(25, stopBnd[4]), y=c(-0.5, -0.5)), aes(x,y,color="s"),size=2.5) +
  geom_line(aes(y=valueFunc,color="h"),size=1.2) +
   theme_light() + scale_colour_manual("", 
                      breaks = c("q", "h", "s"),
                      values = c("red", "blue", "green"),
                      labels=c("q(4,x)", "h(4,x)", "Stopping Region" ))+
   labs(x="X(4)", y="q(4,x)") + theme(legend.position=c(.7, .75),legend.text = element_text(size = 8),
                                      axis.title=element_text(size=9)) +
  scale_x_continuous(expand=c(0.02,0.02))  + scale_y_continuous(limits=c(-0.5,13)) 
grid.arrange(g1.base,g2.base,nrow=1)

put1d.model <- c(K=40, payoff.func=put.payoff,  # payoff function
            x0=40,sigma=0.2,r=0.06,div=0,T=1,dt=0.04,dim=1, sim.func=sim.gbm,
            km.cov=4,km.var=1,kernel.family="matern5_2",  # GP emulator params
            look.ahead=1,pilot.nsims=0,batch.nrep=200,N=25)

train.grid.1d <- seq(16, 40, len=25)
km.fit <- osp.fixed.design(put1d.model,input.domain=train.grid.1d, method="km")

put1d.model$nk=20  # number of knots for the smoothing spline
spl.fit <- osp.prob.design(30000,put1d.model,method="spline")

check.x <- seq(25, 40, len=500)   # predictive sites
km.pred <- predict(km.fit$fit[[10]],data.frame(x=check.x), type="UK") 
to.plot.2 <- data.frame(km.fit=km.pred$mean,x=check.x,km.up=km.pred$upper95,km.down=km.pred$lower95,
                        sp.fit=predict(spl.fit$fit[[10]],check.x)$y)
ggplot(to.plot.2, aes(x=x)) + 
  geom_line(aes(y=km.fit), size=1.25,color="black") +
  geom_line(aes(y=sp.fit), size=1.25, color="purple") +
  geom_line(aes(y=km.up),size=0.5,color="black",linetype="twodash") +
  geom_line(aes(y=km.down),size=0.5,color="black",linetype="twodash") +
   theme_light() + geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(expand=c(0,0),limits=c(26,40))  + scale_y_continuous(expand=c(0,0), limits=c(-0.4,1.1)) +
   labs(x="X[t]", y="Timing Value") 
  #geom_rug(aes(x=x), data= data.frame(x=train.grid.1d),sides="b")

model2d <- list(look.ahead=1, K=40,x0=rep(40,2),sigma=rep(0.2,2),r=0.06,div=0,
                      T=1,dt=0.04,dim=2,sim.func=sim.gbm, payoff.func=put.payoff)

bas22 <- function(x) return(cbind(x[,1],x[,1]^2,x[,2],x[,2]^2,x[,1]*x[,2]))
model2d$bases <- bas22
prob.lm <- osp.prob.design(15000,model2d, method="lm")

model2d$N <- 150  # N_unique 
model2d$kernel.family <- "gauss" # squared-exponential kernel
model2d$batch.nrep <- 100
model2d$pilot.nsims <- 0

sob150 <- sobol(276, d=2) # Sobol space-filling sequence
# triangular approximation domain
sob150 <- sob150[ which( sob150[,1] + sob150[,2] <= 1) ,]  
sob150 <- 25+30*sob150  # Lower-left triangle in [25,55]x[25,55], see Fig 

sob.km <- osp.fixed.design(model2d,input.domain=sob150, method="km")

g1.lm <- plt.2d.surf( prob.lm$fit[[15]], x=seq(25,50, len=101), y=seq(25,50,len=101),  bases=model2d$bases) + guides(fill = guide_colourbar(barwidth = 0.4, barheight = 7))
g2.km <- plt.2d.surf( sob.km$fit[[15]], x=seq(25,50, len=101), y=seq(25,50,len=101))+ylab("") +
  guides(fill = guide_colourbar(barwidth = 0.4, barheight = 7))
g2.km$layers[[3]]$aes_params$size <- 1.4
grid.arrange(g1.lm,g2.km,nrow=1,widths=c(3,2.7))

nSims.2d <- 40000
nSteps.2d <- 25
set.seed(102)
test.2d <- list()
test.2d [[1]] <- model2d$sim.func( matrix(rep(model2d$x0, nSims.2d), nrow=nSims.2d, byrow=T), 
                               model2d, model2d$dt)
for (i in 2:(nSteps.2d+1))
   test.2d [[i]] <- model2d$sim.func( test.2d [[i-1]], model2d, model2d$dt)
oos.lm <- forward.sim.policy( test.2d, nSteps.2d, prob.lm$fit, model2d)
oos.km <- forward.sim.policy( test.2d, nSteps.2d, sob.km$fit,  model2d)
print( c(mean(oos.lm$payoff), mean(oos.km$payoff)) )  # estimates of check{V}(0,X(0))
# sanity check: estimated European option value
print(mean( exp(-model2d$r*model2d$T)*model2d$payoff.func(test.2d [[nSteps.2d]], model2d)))


ls.lm <- osp.prob.design(30000,model2d,method="lm",subset=1:15000)

model2d$batch.nrep <- 100
model2d$kernel.family="matern5_2"
model2d$pilot.nsims <- 1000
model2d$min.lengthscale <- c(1,1)
model2d$max.lengthscale <- c(20,20)
model2d$N <- 500

lattice136 <- as.matrix(expand.grid( seq(25,55,len=16), seq(25,55,len=16)))
lattice136 <- lattice136[ which( lattice136[,1] + lattice136[,2] <= 80) ,]

put2d.lattice <- osp.fixed.design(model2d,input.dom=lattice136, method="km")

model2d$qmc.method <- NULL
model2d$N <- 400  # space-filling inputs to generate. Only in-the-money ones are kept
put2d.lhsAdaptive<- osp.fixed.design(model2d,input.dom=0.04, method="km")

model2d$qmc.method <- randtoolbox::halton
model2d$N <- c(rep(300,8), rep(500,8), rep(800,8)) # design size across time-steps  
put2d.haltonRange <- osp.fixed.design(model2d,input.dom=-1, method="km")

g1.lattice <- plt.2d.surf(put2d.lattice$fit[[10]], x=seq(24,52,len=101),y=seq(24,52,len=101))
g1.lhs <- plt.2d.surf(put2d.lhsAdaptive$fit[[10]],x=seq(24,52,len=101),y=seq(24,52,len=101))
g2.lhs <- plt.2d.surf(put2d.lhsAdaptive$fit[[20]],x=seq(24,52,len=101),y=seq(24,52,len=101))
g1.halton <- plt.2d.surf(put2d.haltonRange$fit[[10]],x=seq(26,48,len=101),y=seq(26,48,len=101))
g2.halton <- plt.2d.surf(put2d.haltonRange$fit[[20]],x=seq(26,48,len=101),y=seq(26,48,len=101))
#g3.km <- plt.2d.surf(put2d.probRep.km$fit[[22]],x=seq(26,48,len=101),y=seq(26,48,len=101))
margin = theme(plot.margin = unit(c(1,0,0,0), "cm"))
#grid.arrange(g1.halton+margin, g2.halton+margin, nrow=1)
tmp <- ggplot_gtable(ggplot_build(g1.lattice))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]

g1.lhs <- g1.lhs + theme(legend.position = "none") + xlab("") + ylab("") +
  scale_fill_gradientn(limits = c(-0.5, 2), oob = scales::squish, colors=tim.colors(64))
g1.lattice <- g1.lattice + theme(legend.position = "none") +  xlab("")  +
   scale_fill_gradientn(limits = c(-0.5, 2), oob = scales::squish, colors=tim.colors(64))
g1.halton <- g1.halton + theme(legend.position = "none")  +
  scale_fill_gradientn(limits = c(-0.5, 2), oob = scales::squish, colors=tim.colors(64))
g2.halton <- g2.halton + theme(legend.position = "none") + ylab("") +
  scale_fill_gradientn(limits = c(-0.5, 2), oob = scales::squish, colors=tim.colors(64))
g1.halton$layers[[3]]$aes_params$size <- 1.4
g1.lhs$layers[[3]]$aes_params$size <- 1.4
g2.halton$layers[[3]]$aes_params$size <- 1.4
g1.lattice$layers[[3]]$aes_params$size <- 1.4
grid.arrange(g1.lattice, g1.halton, g1.lhs, g2.halton, legend, nrow=2,ncol=3,widths=c(4,4,1),
             layout_matrix=cbind(c(1,2),c(3,4),c(5,5)))

oos.1 <- forward.sim.policy(test.2d, nSteps.2d, put2d.lattice$fit,  model2d)
oos.2 <- forward.sim.policy(test.2d, nSteps.2d, put2d.lhsAdaptive$fit,  model2d)
oos.3 <- forward.sim.policy(test.2d, nSteps.2d, put2d.haltonRange$fit,  model2d)
print( c(mean(oos.1$payoff), mean(oos.2$payoff), mean(oos.3$payoff)) )

model2d$init.size <- 30   # initial design size
sob30 <- randtoolbox::sobol(55, d=2)  # build a Sobol space-filling design to initialize
sob30 <- sob30[ which( sob30[,1] + sob30[,2] <= 1) ,]  
sob30 <- 25+30*sob30
model2d$init.grid <- sob30

model2d$batch.nrep <- 25  # N_rep
model2d$seq.design.size <- 120  # final design size -- a total of 3000 simulations
model2d$ei.func <- "sur"  # Stepwise Uncertainty Reduction acquisition function
model2d$kernel.family <- "matern5_2"
model2d$km.cov <- c(15,15); model2d$km.var <- 1
put2d.sur.km <- osp.seq.design(model2d, method="trainkm")
oos.sur.km <- forward.sim.policy( test.2d, nSteps.2d, put2d.sur.km$fit, model2d)

model2d$ei.func <- "tmse"  # targeted mean-squared-error I_n(x)
model2d$tmse.eps <- 0.06   # tMSE parameter
model2d$kernel.family <- "Matern5_2"
put2d.tmse.hetgp <- osp.seq.design(model2d, method="hetgp")
oos.tmse.hetgp <- forward.sim.policy( test.2d, nSteps.2d, put2d.tmse.hetgp$fit, model2d)

model2d$ei.func <- "smcu"  # straddle maximum contour uncertainty
model2d$ucb.gamma <- 1  # sMCU parameter
model2d$kernel.family <- "Gaussian"
put2d.mcu.tp <- osp.seq.design(model2d, method="homtp") # homoskedastic TP 
oos.smcu.tp <- forward.sim.policy( test.2d, nSteps.2d, put2d.mcu.tp$fit, model2d)
print( c(mean(oos.sur.km$payoff), mean(oos.tmse.hetgp$payoff), mean(oos.smcu.tp$payoff)) )

g.tmse.hetgp <- plt.2d.surf(put2d.tmse.hetgp$fit[[10]],x=seq(26,48,len=101),y=seq(26,48,len=101))
g.mcu.tp <- plt.2d.surf(put2d.mcu.tp$fit[[10]],x=seq(26,48,len=101),y=seq(26,48,len=101))
g.sob30 <- plt.2d.surf(put2d.sur.km$fit[[10]],x=seq(26,48,len=101),y=seq(26,48,len=101))

g.sob30 <- g.sob30 + scale_fill_gradientn(limits = c(-0.5, 2), oob = scales::squish, colors=tim.colors(64)) + margin + theme(legend.position = "none")
g.tmse.hetgp <- g.tmse.hetgp+ scale_fill_gradientn(limits = c(-0.5, 2), oob = scales::squish, colors=tim.colors(64)) +ylab("") + theme(legend.position = "none") + margin
g.mcu.tp <- g.mcu.tp+ scale_fill_gradientn(limits = c(-0.5, 2), oob = scales::squish, colors=tim.colors(64)) + guides(fill = guide_colourbar(barwidth = 0.4, barheight = 7)) + ylab("")

g.mcu.tp$layers[[3]]$aes_params$size <- 1.4
g.sob30$layers[[3]]$aes_params$size <- 1.4
g.tmse.hetgp$layers[[3]]$aes_params$size <- 1.4
g.mcu.tp$theme$axis.title$size <- 10
g.sob30$theme$axis.title$size <- 10
g.tmse.hetgp$theme$axis.title$size <- 10

grid.arrange(g.sob30,g.tmse.hetgp,g.mcu.tp+margin,nrow=1,widths=c(2.2, 2.1, 2.9))

model2d$seq.design.size <- 100  # N_unique
model2d$batch.nrep <- 20      # Batch size N_rep
model2d$total.budget <- 2500  # total simulation budget N

model2d$kernel.family <- "gauss"  # GP kernel function: squared-exponential
model2d$update.freq <- 5      # frequency of re-fitting GP hyperparameters
model2d$r.cand <- c(20, 30, 40, 50, 60, 80, 120, 160) # parameter for ABSUR

set.seed(110)
model2d$batch.heuristic <- 'adsa'  # ADSA with AMCU acquisition function
model2d$ei.func <- 'amcu'
oos.obj.adsa <- osp.seq.batch.design(model2d, method="trainkm")

model2d$batch.heuristic <- 'absur'  # Adaptively Batched SUR
model2d$ei.func <- 'absur'
oos.obj.absur <- osp.seq.batch.design(model2d, method="trainkm")

library(RColorBrewer)
library(scales)
g1.adsa <- plt.2d.surf.with.batch(oos.obj.adsa$fit[[15]], 
                       oos.obj.adsa$batches[1:oos.obj.adsa$ndesigns[15] - 1, 15])
#oos.obj.absur$ndesigns[15]  # number of unique designs
### plot Figure 6 left panel - ABSUR
g1.absur <-plt.2d.surf.with.batch(oos.obj.absur$fit[[15]], 
                       oos.obj.absur$batches[1:oos.obj.absur$ndesigns[15] - 1, 15])
g1.absur$layers[[2]]$aes_params$size <- 1.45
g1.adsa$layers[[2]]$aes_params$size <- 1.45
grid.arrange(g1.adsa, g1.absur, nrow=1)

modelBrGl3d <- list(K=100, r=0.05, div=0.1, sigma=rep(0.2,3),T=3, dt=1/3,
  x0=rep(90,3),dim=3, sim.func=sim.gbm,payoff.func=maxi.call.payoff)

# Generate out-of-sample test set
set.seed(44)
nSims.3d <- 20000
nSteps.3d <- 9
test.3d <- list()
test.3d[[1]] <- sim.gbm( matrix(rep(modelBrGl3d$x0, nSims.3d), nrow=nSims.3d, byrow=T), modelBrGl3d)
for (i in 2:nSteps.3d)
   test.3d[[i]] <- sim.gbm( test.3d[[i-1]], modelBrGl3d)
# European option price
#mean( exp(-modelBrGl3d$r*modelBrGl3d$T)*modelBrGl3d$payoff.func(test.3d[[nSteps.3d]],modelBrGl3d))  

modelBrGl3d$rf.ntree = 200  # random forest parameters
modelBrGl3d$rf.maxnode=200
call3d.rf <- osp.prob.design(100000,modelBrGl3d,method="randomforest")
oos.rf <- forward.sim.policy(test.3d,nSteps.3d,call3d.rf$fit,modelBrGl3d,compact=TRUE)

modelBrGl3d$nn.nodes <- 50
call3d.nnet <- osp.prob.design(N=100000,modelBrGl3d, method="nnet")
oos.nnet <- forward.sim.policy( test.3d,nSteps.3d,call3d.nnet$fit,modelBrGl3d)

modelBrGl3d$N <- 800  # N_unique
modelBrGl3d$batch.nrep <- 25 # N_rep
lhs.rect <- matrix(0, nrow=3, ncol=2) # domain of approximation
lhs.rect[1,] <- lhs.rect[2,] <- lhs.rect[3,] <- c(50,150)
modelBrGl3d$qmc.method <- randtoolbox::sobol  # space-filling using QMC sequence

call3d.lhsFixed.rvm <- osp.fixed.design(modelBrGl3d,input.domain=lhs.rect, method="rvm")
oos.rvm <- forward.sim.policy(test.3d, nSteps.3d, call3d.lhsFixed.rvm$fit, modelBrGl3d)

modelBrGl3d$np.kertype <- "gaussian"
modelBrGl3d$np.kerorder <- 2
modelBrGl3d$np.regtype <- "lc"

call3d.sobFixed.np <- osp.fixed.design(modelBrGl3d,input.domain=lhs.rect, method="npreg")
oos.np <- forward.sim.policy( test.3d,nSteps.3d,call3d.sobFixed.np$fit,modelBrGl3d)

modelBrGl3d$kernel.family <- "matern5_2"
# covariance function hyperparameters: process variance and lengthscales
modelBrGl3d$km.var=20; modelBrGl3d$km.cov=c(15,15,15)  
set.seed(1)
call3d.lhsFixed.km <- osp.fixed.design(modelBrGl3d,input.domain=lhs.rect, method="km")  
oos.km <- forward.sim.policy( test.3d,nSteps.3d,call3d.lhsFixed.km$fit, modelBrGl3d)

load("call3d-regr-stats2.RData")
den.pr <- c( mean(stats.rf[,1]), mean(stats.nn[,1]), mean(stats.np[,1]), mean(stats.rvm[,1]), mean(stats.km[,1]))
den.sd <- c( sd(stats.rf[,1]), sd(stats.nn[,1]), sd(stats.np[,1]), sd(stats.rvm[,1]), sd(stats.km[,1]))/5
den.runtime <- c( mean(stats.rf[,3]), mean(stats.nn[,3]), 60*mean(stats.np[,3]), mean(stats.rvm[,3]), 60*mean(stats.km[,3]))
#units(call3d.lhsFixed.rvm$timeElapsed) <- "secs"
#units(call3d.lhsFixed.km$timeElapsed) <- "secs"
#units(call3d.rf$timeElapsed) <- "secs"
#units(call3d.sobFixed.np$timeElapsed) <- "secs"
#units(call3d.nnet$timeElapsed) <- "secs"
 kable(data.frame(method=c("RF", "nnet", "npreg", "RVM", "GP-km"), 
                  #price=c( mean(oos.rf$payoff), mean(oos.nnet$payoff),mean(oos.np$payoff),
                  #         mean(oos.rvm$payoff),  mean(oos.km$payoff)),
                  #  runtime=c(as.numeric(call3d.rf$timeElapsed), as.numeric(call3d.nnet$timeElapsed), as.numeric(call3d.sobFixed.np$timeElapsed), as.numeric(call3d.lhsFixed.rvm$timeElapsed), as.numeric(call3d.lhsFixed.km$timeElapsed)),
                  price=den.pr,sd=den.sd,timedenali=den.runtime), format="markdown",
       caption="\\label{tbl:2dput}Comparison of 3D max-Call solvers with different regression emulators. Running times based on a Surface Book laptop with 8GB memory and Inter Core i-5 2.6GHZ processor.",
       col.names=c("Emulator", "Mean Price", "Std. Error", "Time (secs)"), digits=c(3,3,4,1))

earthParams <- c(earth.deg=2,earth.nk=200,earth.thresh=1E-8)  # passed to {earth}
call3d.tvr <- osp.tvr(N=80000, c(modelBrGl3d,earthParams), method="earth")
oos.tvr <- forward.sim.policy( test.3d,nSteps.3d, call3d.tvr$fit, modelBrGl3d )
print( mean(oos.tvr$payoff))

# polynomials of degree <= 2
bas2 <- function(x) return(cbind(x[,1],x[,1]^2,x[,2],x[,2]^2,x[,1]*x[,2],x[,3],x[,3]^2,
                                 x[,3]*x[,2],x[,1]*x[,3]))
# polynomials up to degree 3 + the payoff
bas3 <- function(x) return(cbind(x[,1],x[,1]^2,x[,2],x[,2]^2,x[,1]*x[,2],x[,3],x[,3]^2,
                                 x[,3]*x[,2],x[,1]*x[,3], x[,1]^3,x[,2]^3,x[,3]^3,
                             x[,1]^2*x[,2],x[,1]^2*x[,3],x[,2]^2*x[,1],x[,2]^2*x[,3],
                             x[,3]^2*x[,1],x[,3]^2*x[,2],x[,1]*x[,2]*x[,3],
                             maxi.call.payoff(x,modelBrGl3d)) )  # include the payoff

modelBrGl3d$bases <- bas2  # 10 coefficients to fit
lm.run2 <- osp.prob.design(300000,modelBrGl3d,method="lm")
oos.lm2 <- forward.sim.policy(test.3d, nSteps.3d, lm.run2$fit, modelBrGl3d)
modelBrGl3d$bases <- bas3  # 21 coefficients to fit
lm.run3 <- osp.prob.design(300000,modelBrGl3d,method="lm")
oos.lm3 <- forward.sim.policy(test.3d, nSteps.3d, lm.run3$fit, modelBrGl3d)

modelBrGl3d$bases <-function(x) {
  sortedx <- t(apply(x, 1, sort, decreasing = TRUE))  # sort coordinates in decreasing order
  return(cbind(sortedx[,1],sortedx[,1]^2, sortedx[,1]^3,sortedx[,1]^4, sortedx[,2], 
               sortedx[,2]^2, sortedx[,3], sortedx[,1]*sortedx[,2], sortedx[,1]*sortedx[,3] ))
}
lm.run4 <- osp.prob.design(300000,modelBrGl3d,method="lm")

oos.lm.sorted <- forward.sim.policy(test.3d, nSteps.3d, lm.run4$fit, modelBrGl3d)
print(mean(oos.lm.sorted$payoff))

modelBrGl3d$nChildren <- 5
bw.run <- osp.probDesign.piecewisebw(300000,modelBrGl3d,test=test.3d)

modelBrGl3d$pilot.nsims <- 1000
modelBrGl3d$batch.nrep <- 60
modelBrGl3d$kernel.family <- "Matern5_2" # different naming compared to km 

modelBrGl3d$N <- 500   
modelBrGl3d$qmc.method <- randtoolbox::halton
put3d.hetgp <- osp.fixed.design(modelBrGl3d,input.dom=0.02, method="hetgp")
oos.hetgp <- forward.sim.policy(test.3d, nSteps.3d, put3d.hetgp$fit, modelBrGl3d)
print(round(mean(oos.hetgp$payoff),digits=4))

## g1 <-plt.2d.surf(put2d.haltonAdaptive.hetgp$fit[[6]],x=seq(26,48,len=101),y=seq(26,48,len=101)) +
##   theme(legend.position = "none") + xlab("") +
##   scale_fill_gradientn(limits = c(-0.5, 4), oob = scales::squish, colors=tim.colors(64))
## g2 <-plt.2d.surf(put2d.haltonAdaptive.hetgp$fit[[14]],x=seq(26,48,len=101),y=seq(26,48,len=101)) +
##  theme(legend.position = "none") + xlab("") + ylab("") +
##   scale_fill_gradientn(limits = c(-0.5, 4), oob = scales::squish, colors=tim.colors(64))
## g3 <-plt.2d.surf(put2d.haltonAdaptive.hetgp$fit[[22]],x=seq(26,48,len=101),y=seq(26,48,len=101)) +
##   xlab("") + ylab("") + scale_fill_gradientn(limits = c(-0.5, 4), oob = scales::squish, colors=tim.colors(64))
## g1$layers[[3]]$aes_params$size <- 1.4
## g2$layers[[3]]$aes_params$size <- 1.4
## g3$layers[[3]]$aes_params$size <- 1.4
## g1$theme$axis.title$size=10
## g2$theme$axis.title$size=10
## g3$theme$axis.title$size=10
## grid.arrange(g1, g2, g3, nrow=1, widths=c(2.5,2.35,3))
## #summary(put2d.haltonAdaptive.hetgp$fit[[14]])

modelBrGl3d$lagp.type="alcray"
modelBrGl3d$lagp.end=40
put3d.lagp <- osp.prob.design(7500,modelBrGl3d,method="lagp",subset=1:2500)

models=c('M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9')
modelsDim <- c(1,1,2,2,2,3,5,5,5)
modelPayoffs <-c("Put", "Put", "Basket Put", "Max Call", "SV Put", "Max Call", 
                 "Max Call", "Max Call", "Basket Put")
modelSimulators <- c("GBM", "GBM", "GBM", "GBM", "SV Heston", "GBM", "GBM", "GBM", "GBM Cor")
modelRefs <- c("[@LS]", "[@LS]", "[@LS]", "[@LS]", "[@Rambharat11]", "[@BroadieCao08]", "[@BroadieCao08]", "[@BroadieCao08]", "[@Lelong19]" ) 
modelSteps <- c(25,25,25,9,50,9,9,9,20)
modelNotes <- c("Classic 1D Put", "Same as M1 but out-of-the-money", "2D symmetric Put", "2D Call with dividends", "Stochastic volatility model with a Put payoff ", "3D Symmetric max-Call", "5D Symmetric Call", "5D Asymmetric max-Call out-of-the-money", "5D asymmetric correlated Put")
modelInfo <- data.frame(Model=models,dim=modelsDim,steps=modelSteps,Payoff=modelPayoffs,Dynamics=modelSimulators,Notes=modelNotes,Ref=modelRefs)

kable(modelInfo, format="markdown", caption="\\label{tbl:osp-instances}Benchmarked OSP Instances",
      booktabs=T, digits = 2) %>% 
  kable_styling(bootstrap_options="striped", full_width=F, latex_options = c("striped"),
                position="center") %>%
column_spec(6, width = "11em")  %>% column_spec(4, width = "6em")

load(file="solver-summary-denali2.RData" )
kable(all.solvers[1:9,], "latex", caption="\\label{tbl:bench-prices}Benchmarked Bermudan option prices for a given run. All methods share a common test set for each instance.",
       col.names=c('Model','S1-LM', 'S2-RF', 'S3-MARS', 'S4-TvR', 'S5-NNet', 'S6-BW', 'S7-GP', 'S8-ADSA', 'S9-Seq', 'S10-hetGP'),
      booktabs=T, digits = 2) %>% 
  kable_styling(bootstrap_options="striped", full_width=F, latex_options = c("striped", "scale_down"),
                position="center")


kable(solvers.time[1:9,], "latex", caption="\\label{tbl:bench-times}Benchmarked Algorithm Running Times (secs)", col.names=c('Model','S1-LM', 'S2-RF', 'S3-MARS', 'S4-TvR', 'S5-NNet', 'S6-BW', 'S7-GP', 'S8-ADSA', 'S9-Seq', 'S10-hetGP'),
      booktabs=T, digits = 1) %>% 
  kable_styling(bootstrap_options="striped", full_width=F, latex_options = c("striped", "scale_down"),
                position="center")

load("stats3d-call.RData")
colnames(stats.3dcall) <- c("base","N","run","inSample","outOfSample", "runTime")
s3 <- data.frame(stats.3dcall)
s3$N <- factor( s3$N, labels =c("40K", "80K", "160K", "320K"))
s3$base <- factor( s3$base, labels =c("quadr", "cubic"))
g.lm <- ggplot(s3,aes(y=outOfSample,x=N,fill=base)) + geom_boxplot()  + #facet_wrap(~base)
  geom_point(position=position_jitterdodge(jitter.width=0.25),size=0.8) + theme_light() +
  coord_cartesian(ylim=c(11.1,11.21))  +
  xlab("# of Paths") + ylab("Out of Sample Mean Payoff") +
  theme(axis.title = element_text(size=8), legend.title=element_text(size=8),legend.text=element_text(size=8)) +  
   scale_fill_brewer(palette="BuPu")
load("nn-tune.Rdata")
colnames(nn.3d.tune) <- c("Nodes","N","run","inSample","outOfSample", "runTime")
n3 <- data.frame(nn.3d.tune)
n3$N <- factor( n3$N, labels =c("40K", "80K", "160K"))
n3$Nodes <- factor( n3$Nodes, labels =c("25", "50", "100"))
g.nn <- ggplot(n3,aes(y=outOfSample,x=N,fill=Nodes)) + geom_boxplot()  + #facet_wrap(~base)
  geom_point(position=position_jitterdodge(jitter.width=0.25),size=0.8) + theme_light() +
  xlab("# of Paths") + ylab("Out of Sample Mean Payoff") +
  coord_cartesian(ylim=c(10.93,11.21)) + 
  theme(axis.title = element_text(size=8),legend.title=element_text(size=8),legend.text=element_text(size=8)) +  #+ stat_boxplot(geom = "errorbar",width=0.25) 
   scale_fill_brewer(palette="BuPu") + guides(color=guide_legend(title="Bases")) 
   
grid.arrange(g.lm, g.nn, nrow=1)

## load("rvm-tune.Rdata")
## colnames(rvm.3d.tune) <- c("Kernel","N","run","outOfSample", "runTime")
## r3 <- data.frame(rvm.3d.tune)
## r3$N <- factor( r3$N, labels =c("500", "800", "1200"))
## r3$Kernel <- factor( r3$Kernel, labels =c("rbfdot","laplacedot","anovadot"))
## g.r <- ggplot(r3,aes(y=runTime,x=N,fill=Kernel)) + geom_boxplot()  + #facet_wrap(~base)
##   geom_point(position=position_jitterdodge(jitter.width=0.25),size=0.85) + theme_light() +
##   xlab("# of Paths") + ylab("Out of Sample Payoff") +
##   coord_cartesian(ylim=c(10.95,11.21)) +
##   theme(axis.title = element_text(size=9),legend.title=element_text(size=8)) +  #+ stat_boxplot(geom = "errorbar",width=0.25)
##    scale_fill_brewer(palette="BuPu") + guides(color=guide_legend(title="Bases"))
## g.r

set.seed(10)
swingModel <- list(dim=1, sim.func=sim.gbm, x0=100,
            swing.payoff=put.payoff, n.swing=3,K=100, 
            sigma=0.3, r=0.05, div=0,
            T=1,dt=0.02,refract=0.1,
            N=800,pilot.nsims=1000,batch.nrep=25)
swingModel$nk=16  # number of knots for the smoothing spline
spl.swing <- swing.fixed.design(swingModel,input.domain=0.03, method ="spline")

set.seed(10); test.swing <- list()  # 25000 forward scenarios
test.swing[[1]] <- sim.gbm( matrix(rep(swingModel$x0, 25000),nrow=25000), swingModel)
for (i in 2:50)
  test.swing[[i]] <- swingModel$sim.func( test.swing[[i-1]], swingModel)

oos.spl3 <- swing.policy(test.swing,50,spl.swing$fit,swingModel,offset=1,n.swing=3)
mean(oos.spl3$totPayoff)

oos.spl2 <- swing.policy(test.swing,50,spl.swing$fit,swingModel,offset=1,n.swing=2)
oos.spl1 <- swing.policy(test.swing,50,spl.swing$fit,swingModel,offset=1,n.swing=1)

sim.corGBM <- function( x0, model, dt)
{   # build a matrix of rho*sigma_i*sigma_j, plus correct the diagonal to be sigma_i^2
    sigm <- model$rho*kronecker(model$sigma, t(model$sigma)) +  
            (1-model$rho)*diag(model$sigma^2)  
    
    newX <- x0*exp( rmvnorm(nrow(x0), sig=sigm*dt, 
                            mean= (model$r- model$div- model$sigma^2/2)*dt) )
    return (newX)
}

modelBecker <- list(dim=5,sigma=0.08*(1:5), r= 0.05, div=0.1, rho=0, 
                    x0 = rep(90,5), T=3, K=100, dt=1/3, 
                    sim.func=sim.corGBM, payoff.func=maxi.call.payoff)

hetGP.params <- list(max.lengthscale=rep(40,5),batch.nrep=100,
                     kernel.family="Matern5_2",pilot.nsims=1000,
                     look.ahead=1,N=500)
modelBecker <- c(modelBecker,hetGP.params)
beckerFit <- osp.fixed.design(modelBecker,input.dom=0.02, method="hetgp")

nSims.becker <- 100000
nSteps.becker <- 9
set.seed(102)
test.Becker <- list()
test.Becker[[1]] <- modelBecker$sim.func( matrix(rep(modelBecker$x0, nSims.becker),
                nrow=nSims.becker, byrow=T), modelBecker, modelBecker$dt)
for (i in 2:nSteps.becker)
   test.Becker[[i]] <- modelBecker$sim.func( test.Becker[[i-1]], modelBecker, modelBecker$dt)

oos.becker <- forward.sim.policy( test.Becker, nSteps.becker, beckerFit$fit, modelBecker)$payoff


#par(mfrow=c(3,2),mar=c(5,4,2,2))
g.base <- list(); counter <- 1
for (k in 4:2) {
  
  rmc[[k]] <- lm(y ~ poly(x,2,raw=TRUE), data.frame(y= payoff[[k]], x = paths[[k]]))
  #plot(paths[[k]], payoff[[k]],xlim=c(29,54), ylim=c(-2,10), xlab="S(t)", ylab="Continuation Value", 
  #   pch=19, cex=1.4, cex.lab=1.2, cex.axis=1.2, bty="n", main=paste("k=", k," with fit ", 
  #        round(coef(rmc[[k]])[1],3), round(coef(rmc[[k]])[2],3),"x+", round(coef(rmc[[k]])[3],3), "x^2"))
  g.base[[counter]] <- ggplot(data.frame(y= payoff[[k]], x = paths[[k]]), aes(x, y)) + geom_point(size=2.5) +
   theme_light() +  geom_smooth(method="lm", formula=y~poly(x, 2), size=1.5, fill=NaN) +
   labs(subtitle=paste("k=", k, ": quadratic fit is ", round(coef(rmc[[k]])[1],3),
                       round(coef(rmc[[k]])[2],3),"x+", 
                       round(coef(rmc[[k]])[3],3), "x^2")) +
  xlab("") + ylab("Y(x)") + theme(plot.subtitle=element_text(size=8),axis.title=element_text(size=8)) +
  scale_x_continuous(expand=c(0.01,0.01))  + scale_y_continuous(limits=c(-1,13)) 

  counter <- counter +1
  #coef(rmc[[k]])
  predy <- predict.lm(rmc[[k]], new=data.frame(x=testx))
  #lines(testx, predy, lwd=4, col="lightblue")
  #plot(testx, pmax(40-testx, 0),lwd=3, col="red", xlab="S", ylab="Cont Value", cex.lab=1.3,type="l")
  #lines(testx, pmax(exp(-lsModel$r*lsModel$dt)*predy,0 ), lwd=3, col="blue")
  #stopBnd[k] <- uniroot(find_stopBnd, c(33,40), rmc[[k]], exp(-lsModel$r*lsModel$dt))$root
  #lines(c(stopBnd[k], 55), c(-0.4, -0.4), col="green", lwd=5)
  to.plot <- data.frame(contValue=pmax(40-testx, 0), payoff=exp(-lsModel$r*lsModel$dt)*predy,
                      valueFunc=pmax(exp(-lsModel$r*lsModel$dt)*predy,0 ),x=testx)

  g.base[[counter]] <- ggplot(to.plot, aes(x=x)) + 
  geom_line(aes(y=contValue, color="q"), size=1.5) +
  geom_line(aes(y=payoff, color="h"), size=1.5, linetype="twodash") +
  #geom_line(data=data.frame(x=c(25, stopBnd[k]), y=c(-0.5, -0.5)), aes(x,y,color="s"),size=2.5) +
  geom_line(aes(y=valueFunc,color="h"),size=1.2) +
   theme_light() + scale_colour_manual("", 
                      breaks = c("q", "h"),
                      values = c("red", "blue"),
                      labels=c("q(k,x)", "h(k,x)"))+
   labs(x="", y="") + theme(legend.position=c(.7, .7),legend.text = element_text(size = 8),
                                      axis.title=element_text(size=8)) +
  scale_x_continuous(limits=c(30,55))  + scale_y_continuous(limits=c(-0.5,13)) 
  counter <- counter + 1
  
  payoff[[k-1]] <- pmax(exp(-lsModel$r*lsModel$dt)*rmc[[k]]$fitted.values, 40-paths[[k]], 0 )
}
grid.arrange(g.base[[1]],g.base[[2]],g.base[[3]],g.base[[4]],g.base[[5]],g.base[[6]],nrow=3)
mean(payoff[[1]])*exp(-lsModel$r*lsModel$dt)


BModel <- list()
BModel[[1]] <- list(dim=1,
            sim.func=sim.gbm,
            K=40, 
            payoff.func=put.payoff,
            x0=40,
            sigma=0.2,
            r=0.06,
            div=0,
            T=1,dt=0.04)


BModel[[2]] <- list(dim=1,
            sim.func=sim.gbm,
            K=40, 
            payoff.func=put.payoff,
            x0=44,
            sigma=0.2,
            r=0.06,
            div=0,
            T=1,dt=0.04)


BModel[[3]] <- list(dim=2,
                    K=40,
                    x0=rep(40,2),
                    sigma=rep(0.2,2),
                    r=0.06,div=0,
                    T=1,dt=0.04,
                    sim.func=sim.gbm, 
                    payoff.func=put.payoff)


BModel[[4]]<- list(dim=2,
                   K=100, 
                   r=0.05, 
                   div=0.1, 
                   sigma=rep(0.2,2),
                   T=3, dt=1/3,
                   x0=rep(110,2),
                   sim.func=sim.gbm,
                   payoff.func= maxi.call.payoff)


BModel[[5]] <- list(K=100,
                    x0=c(90, log(0.35)),
                    r=0.0225,div=0,sigma=1,
    T=50/252,dt=1/252,
    svAlpha=0.015,svEpsY=1,svVol=3,svRho=-0.03,svMean=2.95,
    eulerDt=1/2520, dim=2,
    sim.func=sim.expOU.sv,
    payoff.func =sv.put.payoff)

#putPr <- osp.probDesign.piecewisebw(40000,modelSV5)
 # get putPr$price= 16.81677


BModel[[6]]<- list(dim=3,
                   K=100, 
                   r=0.05, 
                   div=0.1, 
                   sigma=rep(0.2,3),
                   T=3, dt=1/3,
                   x0=rep(90,3),
                   sim.func=sim.gbm,
                   payoff.func= maxi.call.payoff)


BModel[[7]] <- list(dim=5, 
                sim.func=sim.gbm, 
                r=0.05, 
                div=0.1, 
                sigma=rep(0.2,5), 
                x0=rep(100, 5),  # also 70, 130
                payoff.func=maxi.call.payoff, 
                K=100, 
                T=3, 
                dt=1/3)

# paper_result=c(3.892,26.12,59.235) for S_0 = 70, 100, 130

BModel[[8]] <- list(dim=5, 
                sim.func=sim.gbm, 
                r=0.05, 
                div=0.1, 
                sigma=c(0.08,0.16,0.24,0.32,0.4), 
                x0=rep(70, 5),  
                payoff.func=maxi.call.payoff, 
                K=100, 
                T=3, 
                dt=1/3)

# paper_result=c(11.756,37.730,73.709) for S_0 = 70, 100, 130

BModel[[9]] <- list(dim=5, 
                sim.func=sim.gbm.cor, 
                r=0.05, 
                div=0, 
                sigma=0.2, 
                x0=rep(100, 5), 
                rho=0.2,
                K=100,
                payoff.func=put.payoff, 
                T=3, 
                dt=3/20) # 20 steps
#paper_result=4.254 


BModel[[10]] <- list(dim=10, 
                sim.func=sim.gbm.moving.ave, 
                r=0.05, 
                div=0, 
                sigma=0.2, 
                x0=c(100,0,0,0,0,0,0,0,0,0),  
                payoff.func=call.payoff, 
                K=100, 
                T=1, 
                dt=0.02) # 50 steps

# paper_result=c(5.214,11.403,40.440) for S_0 = 90, 100, 130
