library(Rcpp)
library(coda)
library(ggplot2)
library(grid)
library(plyr)
library(gridExtra)
library(MASS)


sourceCpp("../src/model.cpp")
source("../R/model.R")
source("../R/mcmc.R")
source("../R/mcmc_aux.R")
source("../R/analyses.R")
source("../R/plots.R")

n <- 50
tmax <- 1000
n_samples <- 2
times <- seq(0,tmax,by=1)
strains <- c("A","B","C")

nstrains <- 3
cr12 <- 0.3
cr13 <- 0.1
cr23 <- 0.4
tp <- 21
m <- 0.002
a_propn <- 0.7
b_propn <- 0.1
c_propn <- 0.6

a_range <- c(200,600)
b_range <- c(100,400)
c_range <- c(300,800)

sample_start <- 0

mu <- 4
mu_sigma <- 2
error <- 5
S <- 0.79
EA <- 0.2


simulated <- simulate_and_save(n, sample_start, tmax, n_samples, strains, 
                               c(cr12,cr13,cr23),tp,m, c(a_propn,b_propn,c_propn),
                               list(a_range,b_range,c_range),mu,mu_sigma,
                               c("sigma"=error,"S"=S,"EA"=EA,"max_titre"=13),"test",TRUE)

all_data <- simulated[["all_data"]]
all_ind_pars <- simulated[["all_ind_pars"]]
listDat <- simulated[["listDat"]]
pop_pars <- simulated[["pop_pars"]]

p1 <- plot_serology(listDat[,c("t1","A_V1","t2","A_V2")],"A",TRUE)
p2 <- plot_serology(listDat[,c("t1","B_V1","t2","B_V2")],"B",TRUE)
p3 <- plot_serology(listDat[,c("t1","C_V1","t2","C_V2")],"C",TRUE)

for(i in 1:length(all_data)) all_data[[i]] <- unname(as.matrix(all_data[[i]]))

mcmcPars <- c("iterations"=50000,"popt"=0.1,"opt_freq"=500,"thin"=1,"burnin"=1000,"adaptive_period"=50000,"save_block"=500,
              "tuning"=40000,"N_pop_start"=0)

fixed_pop <- c(0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1)
fixed_ind <- c(0,0,0,0)

run_metropolis_MCMC(all_data, times,mcmcPars,"test",
                    fixed_pop,fixed_ind,pop_pars,all_ind_pars)
chains <- NULL
for(i in 1:n) chains[i] <- paste("test_",i,"_chain.csv",sep="")
results <- calculate_statistics(chains,all_ind_pars,all_data,sum(mcmcPars[c("burnin","adaptive_period")]))
write.csv(results,"test_results.csv")
y1 <- generate_counts(results)
y2 <- generate_measures(y1)

