#### Election forecast

rm(list=ls())

# Get data from csv
polls <- read.csv("all_polls.csv")
lvpolls <- polls[polls[,5]=='100% likely' | polls[,5]=='Extremely likely',]
dlvpolls <- lvpolls[lvpolls[,6]=='Clinton'|lvpolls[,6]=='Trump',]

# Get unique dates
dates=unique(polls[,1])
npolls=length(dates) # number of distinct polls

# Get unique states, alphabetical
states = sort(unique(polls[,2]))
nstates = length(states) # number of states + DC

# Make vector with number of electoral votes in each state
evotes = c(3,9,6,11,55,9,7,3,3,29,16,4,6,4,20,11,6,8,8,
           11,10,4,16,10,10,6,3,15,3,5,4,14,5,6,29,18,7,
           7,20,4,9,3,11,38,6,13,3,12,10,5,3)

# Get outcome of polls
nn = matrix(-99,nstates,npolls) # Record total decided likely voters here
outcomes = matrix(-99,nstates,npolls) # Record successes here
for ( ii in 1:length(states) ) {
  for ( jj in 1:length(dates) ) {
    n = sum(dlvpolls$Date == dates[jj] & dlvpolls$Geography == states[ii])
    y = sum(dlvpolls$Date == dates[jj] & dlvpolls$Geography == states[ii] &
              dlvpolls$Question..2.Answer == 'Clinton')
    nn[ii,jj] = n
    outcomes[ii,jj] = y
  }
}

# Define prior on p_i
p.hypers <- matrix(2,nstates,2) # column 1 is shape, col 2 is rate for beta dist
prev_elections_table <- 
  read.csv("modern_results_by_state.csv")
npe <- dim(prev_elections_table)[2] -1 # Number of previous elections
prev_elections <- prev_elections_table[,2:(npe+1)]
# dem_rep_vics will hold number of dem and rep victories in previous elections
dem_rep_vics <- matrix(-99,nstates,2)
for ( i in 1:nstates ) { 
  dem_rep_vics[i,] = c(sum(prev_elections[i,]), npe-sum(prev_elections[i,]) )
  }
#repvics = rep(10,nstates) - demvics
p.hypers = p.hypers + dem_rep_vics

# MCMC settings
M = 1e5 # Total number of samples (including burn-in)
burn_in = 5e3 # floor(M/3)

# MCMC loop function
forecast <- function(a0, 
                     M,
                     burnIn = 0,
                     p=rep(0.5,51), 
                     ns = nn, 
                     poll.outcomes = outcomes, 
                     e.votes = evotes, 
                     p.hyper = p.hypers, 
                     sigma = 1) {
  
  # Get some initial settings
  alpha0  = p.hyper[,1]
  beta0   = p.hyper[,2]
  nstates = dim(poll.outcomes)[1]
  a0.vec  = rep(a0,npolls)^seq(npolls-1,0) # Get vector of decreasing powers of a0
  g0      = log(a0/(1-a0)) # logit transform to eliminate boundary constraints
  # Get parameters for initial beta draw
  alpha   = alpha0 + poll.outcomes %*% a0.vec
  beta    = beta0 + (ns - poll.outcomes) %*% a0.vec
  
  # Set up vehicles for records:
  a0.rec = rep(-99,M)
  accept.rec = rep(-99,M)
  r.rec = rep(-99,M)
  p.rec = matrix(-99,M,nstates)
  state.results.rec = matrix(-99,M,nstates)
  forecast.rec = rep(-99,M)
  
  for (ii in 1:M) {
    
    ### MH step to get a0
    g0.s = rnorm(1,g0,sigma) # Draw new g0 candidate
    a0.s = exp(g0.s)/(1+exp(g0.s)) # Reverse logit trans
    a0.s.vec = rep(a0.s,npolls)^seq(npolls-1,0)
    alpha.s = alpha0 + poll.outcomes %*% a0.s.vec
    beta.s  = beta0 + (ns - poll.outcomes) %*% a0.s.vec
    
    mh.lnum <- sum( dbeta(p, alpha.s, beta.s, log = T) + log(a0.s*(1-a0.s)) ) # log of numerator of MH acceptance ratio
    mh.lden <- sum( dbeta(p, alpha, beta, log = T) + log(a0*(1-a0)) ) # log of denominator of MH acceptance ratio
    r = exp(mh.lnum - mh.lden) # acceptance ratio
    
    accept = 0 # logical to tell us whether new candidate accepted
    if (runif(1) < r) {
      a0 = a0.s
      alpha = alpha.s
      beta = beta.s
      g0 = g0.s
      a0.vec = a0.s.vec
      accept = 1
    }
    
    ### Draw p_i from conditional dist for each state
    p = rbeta(nstates, alpha, beta)
    
    ### Generate state outcomes
    state.results = rbinom(nstates,1,p)
    
    ### Generate election outcome
    forecast=0
    hrc.evotes = state.results %*% e.votes
    if (hrc.evotes >= 270) {forecast = 1} # Clinton victory
    
    ### Record things
    a0.rec[ii] = a0
    accept.rec[ii] = accept
    r.rec[ii] = r
    p.rec[ii,] = p
    state.results.rec[ii,] = state.results
    forecast.rec[ii] = forecast
    
    ### Get periodic update
    if (ii %% 1000 == 0) {
      par(mfrow=c(2,2))
      plot(a0.rec[1:ii])
      plot(p.rec[1:ii,10])
      plot(state.results.rec[(ii-100):ii,40])
      plot(forecast.rec[(ii-100):ii])
    }
  }
  
  nonBurnIn = M - burnIn
  return(list("forecasts"       = tail(forecast.rec,nonBurnIn), 
              "state.forecasts" = tail(p.rec,nonBurnIn),
              "a0"              = tail(a0.rec,nonBurnIn), 
              "accept"          = tail(accept.rec,nonBurnIn), 
              "state.results"   = tail(state.results.rec,nonBurnIn)))
  
}

# Run MCMC and get prediction
res <- forecast(.9,M,burn_in)
prediction <- sum(res$forecasts)/length(res$forecasts)
prediction # This is the estimated probability of Clinton victory


library(coda)

# Get results for each state
state.predictions <- apply(res$state.forecasts,2,mean)
pred.mcmc = as.mcmc(res$state.forecasts) # Coerce the vector into a MCMC object
# pred.mcmc.hpd gives a HPD interval for the probability of each state
# going for Clinton.
pred.mcmc.hpd = round(HPDinterval(pred.mcmc, prob = 0.9),3)                  # Find 95% HPD interval for theta using the CODA function
# The following loop makes a table with one line per state, along with
# expected victor in that state, estimated prob. of Clinton victory, 
# and .9 HPD interval of Clinton victory. In addition, if a state's 
# .9 HPD interval inclues .5, then that state is labelled as a 
# swing state.
for (ii in 1:nstates) {
  safety = '  Swing  ' 
  if (0.5 > pred.mcmc.hpd[ii,2]) {safety = '         '}
  if (0.5 < pred.mcmc.hpd[ii,1]) {safety = '         '}
  winner = '  Clinton  '
  if (state.predictions[ii]<0.5) {winner = '  Trump    '}
  print( paste( states[ii], 
                winner,
                formatC(state.predictions[ii],format='f',digits=3),
                safety, 
                formatC(pred.mcmc.hpd[ii,1],format='f',digits=3), 
                formatC(pred.mcmc.hpd[ii,2],format='f',digits=3)),
         quote=F)
}

par(mfrow=c(1,1))

## Examine presidential prediction
pres.mcmc = as.mcmc(res$forecasts)
#Check autocorrelation of presidential prediction
autocorr.plot(pres.mcmc)
# Check effective sample size of presidential prediction
effectiveSize(pres.mcmc)

## Examine state predictions
plot(pred.mcmc)
# Check autocorrelation of state p_i
autocorr.plot(pred.mcmc)
# Check effective sample size of state p_i
effectiveSize(pred.mcmc)
# Now let's do that for some particular states:
plot(pred.mcmc[,41]) # SC
autocorr.plot(pred.mcmc[,5]) # California
effectiveSize(pred.mcmc[,10]) # Florida

## Examine a0
a0.mcmc = as.mcmc(res$a0)
# Check autocorrelation of a0
autocorr.plot(a0.mcmc)
# Check effective sample size of a0
effectiveSize(a0.mcmc)
mean(a0.mcmc)
HPDinterval(a0.mcmc,prob=0.95)

plot(res$a0, typ = 'l')
plot(res$state.forecasts[,10], typ = 'l')


