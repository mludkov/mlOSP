## ----results="hide",warning=FALSE,message=FALSE,error=FALSE,echo=FALSE-------------------------------------
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

#' 
#' # Test Sets
#' 
#' Generate test sets for the 10 benchmarked models. Models are defined in the arxiv paper appendix.
## ----create-tests------------------------------------------------------------------------------------------
set.seed(44)
all.tests <- list()
nSims <- rep(10000,6)
nSteps <- rep(20,6)

for (j in 1:10) {
  nSims[j] <- 25000
  nSteps[j] <- BModel[[j]]$T/BModel[[j]]$dt
  test.j <- list()
  test.j[[1]] <- BModel[[j]]$sim.func( matrix(rep(BModel[[j]]$x0, nSims[j]), nrow=nSims[j], byrow=T), BModel[[j]])
  for (i in 2:nSteps[j])
    test.j[[i]] <- BModel[[j]]$sim.func( test.j[[i-1]], BModel[[j]])
  # European option price
  print(mean( exp(-BModel[[j]]$r*BModel[[j]]$T)*BModel[[j]]$payoff.func(test.j[[nSteps[j]]],BModel[[j]])) )
  all.tests[[j]] <- test.j
}
save(all.tests,file="allTests-mlosp.RData")

#' 
#' # Solvers
#' 
#' Iterate through the 10 solvers, one by one.
#' 
#' ## LM solver
#' 
#' Define the bases: 3rd degree in 1d and 2d, else quadratic only
## ----lm-bases----------------------------------------------------------------------------------------------
poly_base<-function(x,dim){
  if (dim==1){
    base=cbind(x[,1],x[,1]^2, x[,1]^3)
  } 
  if (dim==2){
    base=cbind(x[,1],x[,1]^2,x[,2],x[,2]^2,x[,1]*x[,2],x[,1]^3, x[,2]^3, x[,1]^2*x[,2], x[,1]*x[,2]^2)
  }
  if (dim==3){
    base=cbind(x[,1],x[,1]^2,x[,2],x[,2]^2,x[,1]*x[,2],x[,3],x[,3]^2,x[,1]*x[,3], x[,2]*x[,3])
  }
  if (dim==5){
    base=cbind(x[,1],x[,1]^2,x[,2],x[,2]^2,x[,1]*x[,2],x[,3],x[,3]^2,x[,1]*x[,3], x[,2]*x[,3], 
               x[,4],x[,4]^2,x[,5],x[,5]^2,x[,1]*x[,4],x[,1]*x[,5], x[,2]*x[,4], x[,2]*x[,5],
               x[,3]*x[,4], x[,3]*x[,5], x[,4]*x[,5])
  }
  
  return(base)
}

#' 
#' 
#' Run the Linear model solver
## ----lm-solver---------------------------------------------------------------------------------------------

lmPrice <-array(0,dim=c(10,3))
for (j in 1:9){
  set.seed(1)
  lmModel <- BModel[[j]]
  lmModel$bases <- function(x) return(cbind( poly_base(x, lmModel$dim), lmModel$payoff.func(x,lmModel)))
  if (lmModel$dim < 3) {
    nPaths <- 41000
  } else {
    nPaths <- 101000
    
  }
  lmSolve <- osp.prob.design(N=nPaths,lmModel, subset=1:1000,method="lm")
  oos.lm <- forward.sim.policy( all.tests[[j]], lmModel$T/lmModel$dt, lmSolve$fit, lmModel)
  units(lmSolve$timeElapsed) <- "secs"
  lmPrice[j,] <- c(lmSolve$p[1],mean(oos.lm$payoff), lmSolve$timeElapsed )
  
}


#' 
#' ## S2: MARS: Multivariate Adaptive Regression Splines solver
## ----mars-solver-------------------------------------------------------------------------------------------
print('Now working on MARS')
library(earth)

mrPrice <-array(0,dim=c(10,3))
for (j in 1:10){
  set.seed(2)
  mrModel <- BModel[[j]]
  mrModel$earth.deg = 2;  # earth parameters
  mrModel$earth.nk = 100;
  mrModel$earth.thresh = 1e-8
  
  if (mrModel$dim < 3) {
    nPaths <- 41000
  } else {
    nPaths <- 101000
    
  }
  mrSolve <- osp.prob.design(N=nPaths,mrModel, subset=1:1000,method="earth")
  oos.mr <- forward.sim.policy( all.tests[[j]], mrModel$T/mrModel$dt, mrSolve$fit, mrModel)
  units(mrSolve$timeElapsed) <- "secs"
  mrPrice[j,] <- c(mrSolve$p[1],mean(oos.mr$payoff), mrSolve$timeElapsed )
  
}

#' 
#' 
#' ## S3: Random Forest Solver
## ----random-forest-solver----------------------------------------------------------------------------------
print('Now working on Random Forest')
library(randomForest)
rfPrice <-array(0,dim=c(10,3))
for (j in 1:9){
  set.seed(3)
  rfModel <- BModel[[j]]
  rfModel$rf.ntree = 200  # random forest parameters
 
  if (rfModel$dim < 3) {
    nPaths <- 41000
    rfModel$rf.maxnode=100 # number of nodes per tree
  
  }
  if (rfModel$dim == 3) {
    nPaths <- 101000
    rfModel$rf.maxnode=200
  }
  if (rfModel$dim == 5) {
    nPaths <- 201000
    rfModel$rf.maxnode=200
  }

  rfSolve <- osp.prob.design(N=nPaths,rfModel, subset=1:1000,method="randomforest")
  oos.rf <- forward.sim.policy( all.tests[[j]], rfModel$T/rfModel$dt, rfSolve$fit, rfModel)
  units(rfSolve$timeElapsed) <- "secs"
  rfPrice[j,] <- c(rfSolve$p[1],mean(oos.rf$payoff), rfSolve$timeElapsed )
  
}

#' 
#' ## S4: Neural Net (single-layer) solver
## ----nnet-solver-------------------------------------------------------------------------------------------
print('Now working on nnet')
library(nnet)
nnPrice <-array(0,dim=c(10,3))

nnNodes <- c(20,20,40,40,40,50,50,50,50)
nn.paths <- c(rep(41000,5),101000,rep(101000,3))
for (j in 1:9){
  set.seed(11)
  nnModel <- BModel[[j]]
  nnModel$nn.nodes = nnNodes[j]  
 
  nnSolve <- osp.prob.design(N=nn.paths[j],nnModel, subset=1:1000,method="nnet")
  oos.nn <- forward.sim.policy( all.tests[[j]], nnModel$T/nnModel$dt, nnSolve$fit, nnModel)
  units(nnSolve$timeElapsed) <- "secs"
  nnPrice[j,] <- c(nnSolve$p[1],mean(oos.nn$payoff), nnSolve$timeElapsed )
  
}

#' ## S5: TvR scheme with MARS emulator
## ----tvr-solver--------------------------------------------------------------------------------------------
library(earth)
print('Now working on TvR')

tvrPrice <-array(0,dim=c(10,3))
for (j in 1:9){
  set.seed(10)
  tvrModel <- BModel[[j]]
  tvrModel$earth.deg = 2;  # earth parameters
  tvrModel$earth.nk = 100;
  tvrModel$earth.thresh = 1e-8
  
  if (tvrModel$dim < 3) {
    nPaths <- 41000
  } else {
    nPaths <- 101000
    
  }
  tvrSolve <- osp.tvr(N=nPaths,tvrModel, subset=1:1000,method="earth")
  oos.tvr <- forward.sim.policy( all.tests[[j]], tvrModel$T/tvrModel$dt, tvrSolve$fit, tvrModel)
  units(tvrSolve$timeElapsed) <- "secs"
  tvrPrice[j,] <- c(tvrSolve$p[1],mean(oos.tvr$payoff), tvrSolve$timeElapsed )
  print(j)
}

#' 
#' 
#' 
#' ## S6: Hierarchical adaptive partitioning solver
## ----bouchard-warin-solver---------------------------------------------------------------------------------
print('Now working on Bouchard-Warin')
bwPrice <-array(0,dim=c(10,3))
for (j in 1:9){
  set.seed(5)
  bwModel <- BModel[[j]]
  if (bwModel$dim < 3) {
    bwModel$nChildren <- 8
    nPaths <- 40000
  } 
  if (bwModel$dim == 3) {
    bwModel$nChildren <- 5
    nPaths <- 100000
  }
  if (bwModel$dim == 5) {
    bwModel$nChildren= 4
    nPaths = 204800
  }  
  bwSolve <- osp.probDesign.piecewisebw(nPaths,bwModel,test=all.tests[[j]])
  units(bwSolve$timeElapsed) <- "secs"
  bwPrice[j,] <- c(bwSolve$price,mean(bwSolve$test),bwSolve$timeElapsed)
  
}


#' ## S7: GP train-km emulator with LHS space-filling design
## ----lhs-trainkm-solver------------------------------------------------------------------------------------
print('Now working on LHS km')
library(DiceKriging)
library(tgp)


lhsPrice <-array(0,dim=c(10,2))
lhsDesSize <- list()
lhs.inputs <- c(rep(400,5),800,1000,1000,2500)
for (j in 1:9){
  set.seed(4)
  lhsModel <- BModel[[j]]
  lhsModel$pilot.nsims <- 1000
  lhsModel$batch.nrep <- 100
  lhsModel$kernel.family="matern5_2"
  lhsModel$qmc.method <- NULL
  lhsModel$N <- lhs.inputs[j]

  lhsSolve <- osp.fixed.design(lhsModel, input.dom=0.02,method="trainkm")
  oos.lhs <- forward.sim.policy( all.tests[[j]], lhsModel$T/lhsModel$dt, lhsSolve$fit, lhsModel)
  units(lhsSolve$timeElapsed) <- "secs"
  lhsPrice[j,] <- c(mean(oos.lhs$payoff), lhsSolve$timeElapsed )
  lhsDesSize[[j]] <- rep(0, length(lhsSolve$fit))  # save the design size at each step
  for (jj in 1:length(lhsSolve$fit))
    lhsDesSize[[j]][jj] <- dim(lhsSolve$fit[[jj]]@X)[1]
  
}


#' 
#' ## S8: GP-km emulator with ADSA sequential batched design
#' 
## ----adaptive-batch-solver---------------------------------------------------------------------------------
print('Now Working on ADSA')
adsaPrice <- array(0,dim=c(10,2))

adsaParams <- list( ucb.gamma=1.96, tmse.eps=0, kernel.family="matern5_2",
                    look.ahead=1, cand.len=1000,pilot.nsims=1000, update.freq=5)

gbmFlag <-c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE,TRUE,TRUE,FALSE)
adsaC.batch <- c(20,20,10,10,10,10,10,10,10)
adsa.len.max <- c(10, 10,20,40,40,40,40,40,40)
adsa.len.min <- c(1,1,2,3,3,3,3,3,3)
adsa.desSize <- c(200,200,200, 200,200,250,500,500,500)
adsa.batch<- c(25,25,40,40,40,50,50,50,50)
adsaInitSize <- c(20,20,20,20,20,40,80,80,80)
adsaInitGrid <- list()
adsaInitGrid[[1]] <- array(seq(25,42,len=20),dim=c(20,1))
adsaInitGrid[[2]] <- array(seq(25,42,len=20),dim=c(20,1))
adsaInitGrid[[3]] <- 25+30*sobol(20,2)
adsaInitGrid[[4]] <- 80+70*sobol(20,2)
adsaInitGrid[[5]] <- sobol(20,2)
adsaInitGrid[[5]][,1]  <- 75 + 27*sobol(20,2)[,1]
adsaInitGrid[[5]][,2] <- -1.5 + 1*sobol(20,2)[,2]
adsaInitGrid[[6]] <- 80+60*sobol(40,3)
adsaInitGrid[[7]] <- 80+70*sobol(80,d=5)
adsaInitGrid[[8]] <- 60+t(c(40,50,60,70,90)*t(sobol(80,d=5)))
adsaInitGrid[[9]] <- 70+60*sobol(80,d=5)


adsa.total.budget <- c(10000, 10000, 10000, 10000,10000, 20000,40000,40000,40000)
adsaDesSize <- list()

library(ks)

for (j in 1:9){
  set.seed(7)
  adsaModel <- c(BModel[[j]],adsaParams)
  adsaModel$batch.nrep <- adsa.batch[j]
  adsaModel$max.lengthscale <- rep( adsa.len.max[j], adsaModel$dim)
  adsaModel$min.lengthscale <- rep( adsa.len.min[j], adsaModel$dim)
  if (j == 5) {
    adsaModel$max.lengthscale <- c(40,2)
    adsaModel$min.lengthscale <- c(3,0.1)
  }
  
  adsaModel$batch.heuristic <- 'adsa'
  adsaModel$ei.func <- 'amcu'
  adsaModel$total.budget <- adsa.total.budget[j]
  adsaModel$seq.design.size <- adsa.desSize[j]
  adsaModel$init.size <- adsaInitSize[j]
  adsaModel$init.grid <- adsaInitGrid[[j]]
  adsaModel$c.batch <- adsaC.batch[j]

  adsaModel$r.cand <- c(20, 30,40,50,60, 80, 120, 160, 200, 250)

  adsaSolve <- osp.seq.batch.design(adsaModel, method="trainkm", is.gbm=gbmFlag[j])
  oos.adsa <- forward.sim.policy( all.tests[[j]], adsaModel$T/adsaModel$dt, adsaSolve$fit, adsaModel)
  units(adsaSolve$timeElapsed) <- "secs"
  adsaDesSize[[j]] <- adsaSolve$ndesigns
  adsaPrice[j,] <- c(mean(oos.adsa$payoff), adsaSolve$timeElapsed)
  
}
#plt.2d.surf( adsaSolve$fit[[7]], x=seq(90,130, len=101), y=seq(90,130,len=101), ub=10)

#' 
#' ## S10:  Heteroskedastic GP model with a fixed simulation design (pseudo-regression):
#' 
## ----hetgp-fixed-sobol-solver------------------------------------------------------------------------------
print('Now working on hetGP')
hetGP.params <- list(kernel.family="Gaussian",pilot.nsims=0)
batch.sizes <- c(20,20,40,40,40,50,50,50,50)
hgp.len.min <- c(1,1,2,2,2,3,3,3,3)
hgp.len.max <- c(20,20,40,40,40,40,40,40,40)

require(randtoolbox)
fixed.designs <- list()
fixed.designs[[1]] <- 25+16*sobol(100)
fixed.designs[[2]] <- 25+16*sobol(100)
fixed.designs[[3]] <- 25+30*sobol(276,d=2)
fixed.designs[[3]] <- fixed.designs[[3]][ which( fixed.designs[[3]][,1] + fixed.designs[[3]][,2] <= 80) ,]  # a lot are on the diagonal
fixed.designs[[4]] <- 80+70*sobol(250,d=2)
fixed.designs[[5]] <- 75+27*sobol(250,d=2)
fixed.designs[[5]][,1] <- 75 + 27*sobol(250,d=2)[,1]
fixed.designs[[5]][,2] <- -1.5 + 1*sobol(250,d=2)[,2]
fixed.designs[[6]] <- 80+70*sobol(500,d=3)
fixed.designs[[7]] <- 80+70*sobol(800,d=5)
fixed.designs[[8]] <- 60+t(c(40,50,60,70,90)*t(sobol(800,d=5)))
fixed.designs[[9]] <- 70+60*sobol(1000,d=5)  # bigger because Put

hgpDesSize <- list()
hgpPrice <-array(0,dim=c(10,2))
for (j in 1:9){
  set.seed(7)
  hgpModel <- c(BModel[[j]],hetGP.params)
  hgpModel$batch.nrep <- batch.sizes[j]
  hgpModel$max.lengthscale <- rep( hgp.len.max[j], hgpModel$dim)
  hgpModel$min.lengthscale <- rep( hgp.len.min[j], hgpModel$dim)
  if (j == 5) {
    hgpModel$max.lengthscale <- c(40,1)
    hgpModel$min.lengthscale <- c(3,0.1)
  }
  if (hgpModel$dim == 1) {
    hgpModel$N <- length(fixed.designs[[j]])
  } else {
    hgpModel$N <- dim(fixed.designs[[j]])[1]
  }
  
  hgpSolve <- osp.fixed.design(hgpModel, input.dom=fixed.designs[[j]],method="hetgp")
  oos.hgp <- forward.sim.policy( all.tests[[j]], hgpModel$T/hgpModel$dt, hgpSolve$fit, hgpModel)
  units(hgpSolve$timeElapsed) <- "secs"
  hgpPrice[j,] <- c(mean(oos.hgp$payoff),hgpSolve$timeElapsed )
  hgpDesSize[[j]] <- rep(0, length(hgpSolve$fit))  # save the design size at each step
  for (jj in 1:length(hgpSolve$fit))
    hgpDesSize[[j]][jj] <-  dim(hgpSolve$fit[[j]]$X0)[1]
  
}

#' 
#' ## S9: GP-km emulator with sequential SUR design
#' 
## ----seq-csur-km-------------------------------------------------------------------------------------------
print('Now working on seq')
des.size <- c(100, 100, 150, 150, 150, 250,300,320,200)
init.n <- c(20, 20,  30, 30, 40, 50,100,100,100)
batch.n <- c( 20, 20, 40, 40, 40, 100,125,125,125)
seq.len.min <- c(1,1,2,2,2,3,3,3,3)
seq.len.max <- c(20,20,40,40,40,40,40,40,40)

init.grids <- list()
init.grids[[1]] <- matrix(25+16*sobol(20),ncol=1)
init.grids[[2]] <- matrix(25+16*sobol(20),ncol=1)
init.grids[[3]] <- 25+30*sobol(30,d=2)
init.grids[[4]] <- 80+70*sobol(30,d=2)
init.grids[[5]] <- 75+27*sobol(30,d=2)
init.grids[[5]][,1] <- 75 + 27*sobol(30,d=2)[,1]
init.grids[[5]][,2] <- -1.5 + 1*sobol(30,d=2)[,2]
init.grids[[6]] <- 80+70*sobol(50,d=3)
init.grids[[7]] <- 80+70*sobol(100,d=5)
init.grids[[8]] <- 40+70*sobol(100,d=5)
init.grids[[9]] <- 70+60*sobol(100,d=5)
seqPrice <- array(0, dim=c(10,2))
seq.candLen <- c(500,1000,500,500,1000,2000,2000,4000,2000)
lhs.size <- c(rep(0.02,6),0.05, 0.02, 0.05)

for (j in 1:9){
  set.seed(9)
  seqGP.params <- list(kernel.family = "Matern5_2", ei.func="sur", 
                       batch.nrep =batch.n[j], cand.len=seq.candLen[j], pilot.nsims=1000, lhs.rect=lhs.size[j],
                     init.size=init.n[j], init.grid = init.grids[[j]], seq.design.size=des.size[j])
  seqModel <- c(BModel[[j]],seqGP.params)
  seqModel$max.lengthscale <- rep( seq.len.max[j], seqModel$dim)
  seqModel$min.lengthscale <- rep( seq.len.min[j], seqModel$dim)
  if (j == 5) {
    seqModel$max.lengthscale <- c(40,1)
    seqModel$min.lengthscale <- c(3,0.1)
  }
  seqSolve <- osp.seq.design(seqModel, method="hetgp")
  oos.seq <- forward.sim.policy( all.tests[[j]], seqModel$T/seqModel$dt, seqSolve$fit, seqModel)
  units(seqSolve$timeElapsed) <- "secs"
  seqPrice[j,] <- c(mean(oos.seq$payoff),seqSolve$timeElapsed )
  print(j)
}

#' 
#' 
#' ## Aggregate all results
#' 
#' 
## ----solver-table------------------------------------------------------------------------------------------
all.solvers <- data.frame(models=c('M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'),
                          lm.poly=lmPrice[,2],rf=rfPrice[,2], mars=mrPrice[,2],
                          tvr.mars=tvrPrice[,2],nnet=nnPrice[,2],
                          bou.war=bwPrice[,2], lhs.trainkm = lhsPrice[,1],
                          adsa=adsaPrice[,1], seq.sur=seqPrice[,1], hetgp.fix=hgpPrice[,1])
solvers.time <- data.frame(models=c('M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'),
                            lm.time = lmPrice[,3], rf.time=rfPrice[,3], mr.time =mrPrice[,3],
                            tvr.time=tvrPrice[,3], nn.time=nnPrice[,3],
                           bw.time = bwPrice[,3], lhs.time=lhsPrice[,2],
                           adsa.time=adsaPrice[,2],seq.time=seqPrice[,2],het.time = hgpPrice[,2])
solvers.insample <- data.frame(models=c('M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'),
                               lm.in=lmPrice[,1],rf.in=rfPrice[,1], mars.in=mrPrice[,1],
                               tvr.in=tvrPrice[,1],nnet.in=nnPrice[,1],bw.in=bwPrice[,1]) 

#save(bwSolve,lhsSolve,lmSolve,mrSolve,rfSolve,all.solvers, file="solvers-mlosp.RData") -- 200MB put in Downloads
save(solvers.time, all.solvers,solvers.insample,file="solver-summary.RData" )
kable(all.solvers, caption="Benchmarked Bermudan Option Prices") %>% 
  kable_styling(bootstrap_options="striped", full_width=F)
kable(solvers.time, caption="Benchmarked Algorithm Running Times") %>% 
  kable_styling(bootstrap_options="striped", full_width=F)

#' 
#' 
#' # Swing Solver
#' 
## ---- swing-km,  eval=FALSE--------------------------------------------------------------------------------
## swingModel$kernel.family="matern5_2"
## swingModel$qmc.method <- NULL
## swingModel$N <- 120;
## km.swing2 <- swing.fixed.design(swingModel,input.domain=c(65,67,69,seq(70,101,len=117)), method ="trainkm")
## # next one is much smaller since only in-the-money is kept
## #km.swing <- swing.fixed.design(swingModel,input.domain=0.02, method ="trainkm")
## # Below are benchmarks for the corresponding Put
## #km.put <- osp.fixed.design(swingModel,input.domain=0.02,method="trainkm")
## #spl.put <- osp.prob.design(60000,swingModel,subset=1:20000,method="spline")
## #spl.tvr <- osp.tvr(60000,swingModel,subset=1:20000,method="spline")
## #earthParams <- c(earth.deg=2,earth.nk=200,earth.thresh=1E-8)

