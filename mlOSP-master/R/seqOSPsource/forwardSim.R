######################################
#' Simulate forward to sample h(X_{tau}) using GLOBAL LINEAR REGRESSION
#' @param M     is the number of forward steps
#' @param X     is a vector containing the initial states
#' @param coefs is a list of lm objects specifying the classifier based on a regression
#' @param model is a list containing all the model parameters
#' @export
forward.sim.regr <- function( x,M,coefs,model)
{
    curX <- model$sim.func( x,model,model$dt)
    payoff <- option.payoff(x, model$K)

    contNdx <- 1:length(x)
    i <- 1

    while (i <=M & length(contNdx) > 0) {
        rule <- boundary.policy( curX[contNdx], coefs[[i]])
        payoff[contNdx]  <- exp(-(i-1)*model$dt*model$r)*option.payoff( curX[contNdx], model$K)
        contNdx <- contNdx[which( rule > 0 | payoff[contNdx] == 0) ]
        curX[contNdx] <- model$sim.func( curX[contNdx],model,model$dt)
        i <- i+1
    }
    return( payoff)
}


################################################
#' Simulate h(X_tau) using MIKE'S HOMEGROWN TREES
#'
#' @param x     is a vector of starting values
#' @param M     is number of time steps to forward simulate
#' @param fit   is a list of fits containing the classifiers to be used
#' @param model is the list of all model parameters
#' @export
forward.sim.tree <- function( x,M,fit,model)
{
    curX <- model$sim.func( x,model,model$dt)
    payoff <- option.payoff(x, model$K)
    tau <- rep(0, length(x))
    sims <- vector(mode="list",len=model$look.ahead+1)
    save.ndx <- vector(mode="list",len=model$look.ahead+1)
    fvalue <- vector(mode="list",len=model$look.ahead+1)

    contNdx <- 1:length(x)
    i <- 2

    while (i <=(M+1) & length(contNdx) > 0) {
        indices <- get.tree.index( curX[contNdx],fit[i,c('lower','upper')])
        rule <- fit[[i,'mean']][indices]

        payoff[contNdx]  <- exp(-(i-1)*model$dt*model$r)*option.payoff( curX[contNdx], model$K)
        tau[contNdx] <- (i-1)*model$dt

        # continue if the expected gain is positive or payoff is zero
        contNdx <- contNdx[which( rule > 0 | payoff[contNdx] == 0) ]
        sims[[min(i,model$look.ahead+1)]] <- curX[contNdx]
        save.ndx[[min(i,model$look.ahead+1)]] <- contNdx

        curX[contNdx] <- model$sim.func( curX[contNdx],model,model$dt)
        i <- i+1
    }
    for (i in 2:(model$look.ahead+1))   # payoff for a trajectory starting at x^n_{t+i} which was still alive 
    then
       fvalue[[i]] <- payoff[ save.ndx[[i]] ]*exp((i-1)*model$dt*model$r)

    return( list(payoff=payoff,sims=sims,fvalue=fvalue,tau=tau))
}
