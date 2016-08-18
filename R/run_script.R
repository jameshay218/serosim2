launch_sim_mcmc_job <- function(all_params, strains, mcmcPars, runLocal=TRUE, filename, wd){
    fixed_pop <- c(0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1)
    fixed_ind <- c(0,0,0,0)

    if(runLocal) setwd("~/net/home/serosim/outputs/")
    else setwd("Q:/serosim/outputs")

    if(!dir.exists(wd)) dir.create(wd)
    setwd(wd)

    addNoise <- all_params[["addNoise"]]
    
    n <- all_params[["n"]]
    tmax <- all_params[["tmax"]]
    tstart <- all_params[["tstart"]]
    
    n_samples <- all_params[["n_samples"]]
    times <- seq(0, tmax, by=1)

    nstrains <- length(strains)

    cr12 <- all_params[["cr12"]]
    cr13 <- all_params[["cr13"]]
    cr23 <- all_params[["cr23"]]

    tp <- all_params[["tp"]]
    m <- all_params[["m"]]

    a_propn <- all_params[["a_propn"]]
    b_propn <- all_params[["b_propn"]]
    c_propn <- all_params[["c_propn"]]

    a_range <- all_params[["a_range"]]
    b_range <- all_params[["b_range"]]
    c_range <- all_params[["c_range"]]

    mu <- all_params[["mu"]]
    mu_sigma <- all_params[["mu_sigma"]]
    error <- all_params[["sigma"]]
    S <- all_params[["S"]]
    EA <- all_params[["EA"]]

    max_titre <- all_params[["max_titre"]]

    simulated <- simulate_and_save(n, tstart,tmax, n_samples, strains, 
                               c(cr12,cr13,cr23),tp,m, c(a_propn,b_propn,c_propn),
                               list(a_range,b_range,c_range),mu,mu_sigma,
                               c("sigma"=error,"S"=S,"EA"=EA,"max_titre"=13),addNoise)

    all_data <- simulated[["all_data"]]
    all_ind_pars <- simulated[["all_ind_pars"]]
    listDat <- simulated[["listDat"]]
    pop_pars <- simulated[["pop_pars"]]
    

    plot_serology(listDat[,c("t1",paste(strains[1],"_V1",sep=""),"t2",paste(strains[1],"_V2",sep=""))],strains[1],TRUE)
    plot_serology(listDat[,c("t1",paste(strains[2],"_V1",sep=""),"t2",paste(strains[2],"_V2",sep=""))],strains[2],TRUE)
    plot_serology(listDat[,c("t1",paste(strains[3],"_V1",sep=""),"t2",paste(strains[3],"_V2",sep=""))],strains[3],TRUE)

    for(i in 1:length(all_data)) all_data[[i]] <- unname(as.matrix(all_data[[i]]))

    run_metropolis_MCMC(all_data, times,mcmcPars,"test",
                        fixed_pop,fixed_ind,pop_pars,all_ind_pars)

    chains <- NULL
    for(i in 1:n) chains[i] <- paste(filename,"_",i,"_chain.csv",sep="")
    
    results <- calculate_statistics(chains,all_ind_pars,all_data,sum(mcmcPars[c("burnin","adaptive_period")]))
    write.csv(results,paste(filename,"_results.csv",sep=""))
    y1 <- generate_counts(results)
    y2 <- generate_measures(y1)
    return(y2)
}
