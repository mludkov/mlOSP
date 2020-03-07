# load("~/Google Drive/GPstuff-4.6/int_2d.RData")
load("~/Dropbox/int_2d.RData")
library(mlOSP)
library(laGP)
library(ggplot2)

source("~/Dropbox/mlOSP-master/R/ospSeqBatchDesign.R")
source("~/Dropbox/mlOSP-master/R/batchDesign_utils.R")

model2d <- list(x0 = rep(40,2),K=40,sigma=rep(0.2,2),r=0.06,div=0,T=1,dt=0.04,dim=2,
                sim.func=sim.gbm, payoff.func=put.payoff)
MM <- (model2d$T/model2d$dt)
model2d$pilot.nsims <- 1000
model2d$look.ahead <- 1

# option price for test points 
NN <- 16000
MM <- 25
set.seed(101)
mygr <- list()
mygr[[1]] <- model2d$sim.func( matrix(rep(model2d$x0, NN), nrow=NN, byrow=T),
                               model2d, model2d$dt)
for (i in 2:(MM+1)) {
  mygr[[i]] <- model2d$sim.func( mygr[[i-1]], model2d, model2d$dt)
}

# sanity check: European option value
option.payoff <- put.payoff
print(mean(exp(-model2d$r*model2d$T)*option.payoff(mygr[[MM]], model2d)))

model2d$cand.len <- 1000
model2d$max.lengthscale <- c(20,20)
model2d$min.lengthscale <- c(0.2,0.2)
model2d$init.size <- 20   # initial design size
model2d$init.grid <- int_2d
model2d$tmse.eps <- 0

model2d$ucb.gamma <- 1.96

model2d$seq.design.size <- 100
model2d$batch.nrep <- 20
model2d$ei.func <- "amcu"
set.seed(101)

model2d$r.cand =c(20,30,40,50,60, 80, 120, 160)

# tp + fb
model2d$batch.heuristic <- 'fb'
model2d$ei.func <- 'amcu'
model2d$kernel.family <- 'Gaussian'

oos.obj.fb.tp <- osp.seq.batch.design(model2d, method="homtp")
oos.obj.fb.tp$fit[[15]]$X0

plt.2d.surf(oos.obj.fb.tp$fit[[15]])

# tp + rb
model2d$batch.heuristic <- 'rb'
model2d$kernel.family <- "Gaussian"

oos.obj.rb.tp <- osp.seq.batch.design(model2d, method="homtp")
oos.obj.rb.tp$fit[[15]]$X0

plt.2d.surf(oos.obj.rb.tp$fit[[15]])

# tp + mlb
model2d$batch.heuristic <- 'mlb'

oos.obj.mlb.tp <- osp.seq.batch.design(model2d, method="homtp")
oos.obj.mlb.tp$fit[[15]]$X0
plt.2d.surf(oos.obj.mlb.tp$fit[[15]])

# tp + absur
model2d$batch.heuristic <- 'absur'
model2d$ei.func <- 'absur'

oos.obj.absur.tp <- osp.seq.batch.design(model2d, method="homtp")
oos.obj.absur.tp$fit[[15]]$X0
plt.2d.surf(oos.obj.absur.tp$fit[[15]])

# tp + adsa
model2d$batch.heuristic <- 'adsa'
model2d$ei.func <- 'amcu'

oos.obj.adsa.tp <- osp.seq.batch.design(model2d, method="homtp")
oos.obj.adsa.tp$fit[[15]]$X0
plt.2d.surf(oos.obj.adsa.tp$fit[[15]])

# tp + ddsa
model2d$batch.heuristic <- 'ddsa'
model2d$ei.func <- 'amcu'

oos.obj.ddsa.tp <- osp.seq.batch.design(model2d, method="homtp")
oos.obj.ddsa.tp$fit[[15]]$X0
plt.2d.surf(oos.obj.ddsa.tp$fit[[15]])

# hetgp + fb
model2d$batch.heuristic <- 'fb'
model2d$ei.func <- 'amcu'
model2d$kernel.family <- 'Gaussian'

oos.obj.fb.hetgp <- osp.seq.batch.design(model2d, method="hetgp")
oos.obj.fb.hetgp$fit[[15]]$X0

plt.2d.surf(oos.obj.fb.hetgp$fit[[15]])

# hetgp + rb
model2d$batch.heuristic <- 'rb'

oos.obj.rb.hetgp <- osp.seq.batch.design(model2d, method="hetgp")
oos.obj.rb.hetgp$fit[[15]]$X0

plt.2d.surf(oos.obj.rb.hetgp$fit[[15]])

# hetgp + mlb
model2d$batch.heuristic <- 'mlb'

oos.obj.mlb.hetgp <- osp.seq.batch.design(model2d, method="hetgp")
oos.obj.mlb.hetgp$fit[[15]]$X0
plt.2d.surf(oos.obj.mlb.hetgp$fit[[15]])

# hetgp + absur
model2d$batch.heuristic <- 'absur'
model2d$ei.func <- 'absur'

oos.obj.absur.hetgp <- osp.seq.batch.design(model2d, method="hetgp")
oos.obj.absur.hetgp$fit[[15]]$X0
plt.2d.surf(oos.obj.absur.hetgp$fit[[15]])

# hetgp + adsa
model2d$batch.heuristic <- 'adsa'
model2d$ei.func <- 'amcu'

oos.obj.adsa.hetgp <- osp.seq.batch.design(model2d, method="hetgp")
oos.obj.adsa.hetgp$fit[[15]]$X0
plt.2d.surf(oos.obj.adsa.hetgp$fit[[15]])

# hetgp + ddsa
model2d$batch.heuristic <- 'ddsa'
model2d$ei.func <- 'amcu'

oos.obj.ddsa.hetgp <- osp.seq.batch.design(model2d, method="hetgp")
oos.obj.ddsa.hetgp$fit[[15]]$X0
plt.2d.surf(oos.obj.ddsa.hetgp$fit[[15]])

# trainkm + fb
model2d$batch.heuristic <- 'fb'
model2d$ei.func <- 'amcu'
model2d$kernel.family <- 'gauss'

oos.obj.fb.km <- osp.seq.batch.design(model2d, method="trainkm")
oos.obj.fb.km$fit[[15]]@X

plt.2d.surf(oos.obj.fb.km$fit[[15]])

# trainkm + rb
model2d$batch.heuristic <- 'rb'
model2d$kernel.family <- 'gauss'

oos.obj.rb.km <- osp.seq.batch.design(model2d, method="trainkm")
oos.obj.rb.km$fit[[15]]@X

plt.2d.surf(oos.obj.rb.km$fit[[15]])

# trainkm + mlb
model2d$batch.heuristic <- 'mlb'
model2d$kernel.family <- 'gauss'

oos.obj.mlb.km <- osp.seq.batch.design(model2d, method="trainkm")
oos.obj.mlb.km$fit[[15]]@X

plt.2d.surf(oos.obj.mlb.km$fit[[15]])

# trainkm + absur
model2d$batch.heuristic <- 'absur'
model2d$ei.func <- 'absur'
model2d$kernel.family <- 'gauss'

oos.obj.absur.km <- osp.seq.batch.design(model2d, method="trainkm")
oos.obj.absur.km$fit[[15]]@X

plt.2d.surf(oos.obj.absur.km$fit[[15]])

# trainkm + adsa
model2d$batch.heuristic <- 'adsa'
model2d$ei.func <- 'amcu'
model2d$kernel.family <- 'gauss'

oos.obj.adsa.km <- osp.seq.batch.design(model2d, method="trainkm")
oos.obj.adsa.km$fit[[15]]@X

plt.2d.surf(oos.obj.adsa.km$fit[[15]])

# trainkm + ddsa
model2d$batch.heuristic <- 'ddsa'
model2d$ei.func <- 'amcu'
model2d$kernel.family <- 'gauss'

oos.obj.ddsa.km <- osp.seq.batch.design(model2d, method="trainkm")
oos.obj.ddsa.km$fit[[15]]@X

plt.2d.surf(oos.obj.ddsa.km$fit[[15]])

# batch design
oos.fb.hetgp <- forward.sim.policy( mygr, MM, oos.obj.fb.hetgp$fit, model2d)
oos.rb.hetgp <- forward.sim.policy( mygr, MM, oos.obj.rb.hetgp$fit, model2d) 
oos.mlb.hetgp <- forward.sim.policy( mygr, MM, oos.obj.mlb.hetgp$fit, model2d) 
oos.absur.hetgp <- forward.sim.policy( mygr, MM, oos.obj.absur.hetgp$fit, model2d) 
oos.adsa.hetgp <- forward.sim.policy( mygr, MM, oos.obj.adsa.hetgp$fit, model2d) 
oos.ddsa.hetgp <- forward.sim.policy( mygr, MM, oos.obj.ddsa.hetgp$fit, model2d)

oos.fb.tp <- forward.sim.policy( mygr, MM, oos.obj.fb.tp$fit, model2d) 
oos.rb.tp <- forward.sim.policy( mygr, MM, oos.obj.rb.tp$fit, model2d) 
oos.mlb.tp <- forward.sim.policy( mygr, MM, oos.obj.mlb.tp$fit, model2d) 
oos.absur.tp <- forward.sim.policy( mygr, MM, oos.obj.absur.tp$fit, model2d) 
oos.adsa.tp <- forward.sim.policy( mygr, MM, oos.obj.adsa.tp$fit, model2d) 
oos.ddsa.tp <- forward.sim.policy( mygr, MM, oos.obj.ddsa.tp$fit, model2d) 

oos.fb.km <- forward.sim.policy( mygr, MM, oos.obj.fb.km$fit, model2d) 
oos.rb.km <- forward.sim.policy( mygr, MM, oos.obj.rb.km$fit, model2d) 
oos.mlb.km <- forward.sim.policy( mygr, MM, oos.obj.mlb.km$fit, model2d) 
oos.absur.km <- forward.sim.policy( mygr, MM, oos.obj.absur.km$fit, model2d) 
oos.adsa.km <- forward.sim.policy( mygr, MM, oos.obj.adsa.km$fit, model2d) 
oos.ddsa.km <- forward.sim.policy( mygr, MM, oos.obj.ddsa.km$fit, model2d) 

print(c(mean(oos.fb.tp$payoff), mean(oos.rb.tp$payoff), mean(oos.mlb.tp$payoff), 
        mean(oos.absur.tp$payoff), mean(oos.adsa.tp$payoff), mean(oos.ddsa.tp$payoff)))

print(c(mean(oos.fb.hetgp$payoff), mean(oos.rb.hetgp$payoff), mean(oos.mlb.hetgp$payoff), mean(oos.absur.hetgp$payoff),
        mean(oos.adsa.hetgp$payoff), mean(oos.ddsa.hetgp$payoff)))

print(c(mean(oos.fb.km$payoff), mean(oos.rb.km$payoff), mean(oos.mlb.km$payoff), 
        mean(oos.absur.km$payoff), mean(oos.adsa.km$payoff), mean(oos.ddsa.km$payoff)))
