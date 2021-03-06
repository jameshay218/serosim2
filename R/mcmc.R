#' @export
run_metropolis_MCMC <- function(all_data,
                              times,
                              mcmcPars=c("iterations"=1000,"popt"=0.44,"opt_freq"=50,"thin"=1,"burnin"=100,"adaptive_period"=100,"save_block"=500,"tuning"=100,"N_pop_start"=50),
                              filename,
                              pop_fixed,
                              ind_fixed,
                              all_pop_pars,
                              all_ind_pars,
                              findY0s=FALSE,
                              oddLikelihood=FALSE){
    iterations <- mcmcPars["iterations"]
    tuning <- mcmcPars["tuning"]
    popt <- mcmcPars["popt"]
    opt_freq<- mcmcPars["opt_freq"]
    thin <- mcmcPars["thin"]
    burnin <- mcmcPars["burnin"]
    adaptive_period<- mcmcPars["adaptive_period"]
    save_block <- mcmcPars["save_block"]
    N_pop_start <- mcmcPars["N_pop_start"]
    non_tuning_pop <- max(N_pop_start - burnin,0)
############################################################################
    ## MCMC parameter extraction
    ############################################################################ 
    TUNING_ERROR<- 0.1

    if(tuning > adaptive_period) tuning <- adaptive_period

    tmax <- max(times)
    print("MCMC paramters extracted")
    ############################################################################


    ############################################################################
    ## Also need to monitor population parameters
    ############################################################################
    pop_step <- 1 # population parameter proposal step size
    pop_pars <- c("mu_mu"=runif(1,1,8),"mu_sigma"=runif(1,1,5),
                  "error_sigma"=5,"S"=0.79,"EA"=0.2,
                  "cr12"=runif(1,0.1,0.5),"cr13"=runif(1,0.1,0.5),"cr23"=runif(1,0.1,0.5),
                  "cr12_sigma"=1,"cr13_sigma"=1,"cr23_sigma"=1,
                  "tp"=21,"tp12"=21,"tp13"=21,"tp23"=21,
                  "tp_sigma"=1,"tp12_sigma"=1,"tp13_sigma"=1,"tp23_sigma"=1,                  
                  "m"=runif(1,0.0005,0.005),"m12"=runif(1,0.0005,0.005),"m13"=runif(1,0.0005,0.005),"m23"=runif(1,0.0005,0.005),
                  "m_sigma"=1,"m12_sigma"=1,"m13_sigma"=1,"m23_sigma"=1
                  )
    
    
    ## Initial covariance matrix
    tmp_pars <- pop_pars
    all_pop_pars[which(pop_fixed==0)] <- pop_pars[which(pop_fixed==0)]
    ##tmp_pars[which(pop_fixed == 1)] <- 1
    covMat <- diag(0.1*tmp_pars^2)
    
    mcmc_pop_filename <- paste(filename,"_pop_chain.csv",sep="")

    ## Setup empty matrices for saving MCMC
    pop_chain_colnames <- c("sampno",names(pop_pars),"lnlike")
    empty_pop_chain <- matrix(nrow=adaptive_period,ncol=length(pop_chain_colnames))
    empty_save_pop_chain <- matrix(nrow=save_block,ncol=length(pop_chain_colnames))
    colnames(empty_pop_chain) <- colnames(empty_save_pop_chain) <- pop_chain_colnames
                              
    POP_PARAMETER_CONTROL <- list("pars"=all_pop_pars,"fixed_pars"=pop_fixed,"filename"=mcmc_pop_filename,"cov"=covMat,
                                  "step_scale"=pop_step,"save_chain"=empty_save_pop_chain,"opt_chain"=empty_pop_chain,
                                  "accepted"=0,"iter"=0,"curr_lik"=0,"chain_i"=1,"opt_i"=1)
    print("Pop parameters set up")
    ############################################################################

    
    likelihood_ind <- make_likelihood(NULL,NULL, all_pop_pars, all_data[[1]],times, oddLikelihood)
    
    ############################################################################
    ## Each individual needs to be monitored
    ############################################################################
    n_individuals <- length(all_data)

    ## Empty matrices to save MCMC steps for each individual
    ind_pars <- c("mu_i"=8,"ti1"=100,"ti2"=100,"ti3"=100)
    if(findY0s) ind_pars <- c(ind_pars, "y0_1"=0,"y0_2"=0,"y0_3"=0)
    chain_colnames <- c("sampno",names(ind_pars),"lnlike")

    ind_step <- 1 # initial MCMC step size

    ## Create list of lists for all individuals
    pop_lik <- 0 ## Need to calculate initial population likelihood
    ind_liks <- numeric(n_individuals)
    
    ALL_INDIVIDUALS_CONTROL <- list()
    print("Setting up individual parameters")
    for(i in 1:(n_individuals)){
        mcmc_chain_file <- paste(filename,"_",i,"_chain.csv",sep="") ## Each individual has a save file
        tmp_accepted <- tmp_iter <- 0
        
        ## Create empty chain to store every iteration for the adaptive period
        chain <- matrix(nrow=adaptive_period,ncol=length(chain_colnames))
        
        ## Create empty chain to store "save_block" iterations at a time
        save_chain <- matrix(nrow=save_block,ncol=length(chain_colnames))

        colnames(chain) <- colnames(save_chain) <- chain_colnames

        ## Starting parameter values and covariance matrix for mvrnorm proposals
        startPars <- all_ind_pars[[i]]$pars
        ## Starting positions for non-fixed parameters
        ##startPars[which(ind_fixed==0)] <- runif(length(startPars[which(ind_fixed==0)]),1,15)
        startPars[1] <- runif(1,1,12)
        startPars[2:4] <- runif(3,1,tmax)
        tmpPars <- startPars
        
        tmp_y0s <- all_ind_pars[[i]]$y0s
        if(findY0s){
            tmp_y0s <- rpois(3,0.1)
            tmpPars <- c(tmpPars, 0.1,0.1,0.1)
        }
        
        covMat <- diag(0.1*tmpPars^2)
        tmp_y0s <- all_ind_pars[[i]]$y0s

        ## Calculate starting likelihood
        ini_lik <- likelihood_ind(startPars,tmp_y0s,POP_PARAMETER_CONTROL$pars,all_data[[i]])
        pop_lik <- pop_lik + ini_lik
        ind_liks[i] <- ini_lik
        
        ## Save starting parameters and likelihood
        tmp_table <- as.data.frame(array(dim=c(1,length(chain_colnames))))
        tmp_table[1,] <- c(0,as.numeric(tmpPars),ini_lik)
        colnames(tmp_table) <- chain_colnames
        write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

        ALL_INDIVIDUALS_CONTROL[[i]] <- list("pars"=startPars,"fixed_pars"=ind_fixed,"y0s"=tmp_y0s,"filename"=mcmc_chain_file,
                                             "accepted"=tmp_accepted,"iter"=tmp_iter,"opt_chain"=chain,
                                             "save_chain"=save_chain,"cov"=covMat,"step_scale"=ind_step,
                                             "curr_lik"=ini_lik,"chain_i"=1,"opt_i"=1)
    }
    print("Individual parameters set up")
    ############################################################################

    ############################################################################
    #### Initial likelihood for pop parameters
    ############################################################################
    POP_PARAMETER_CONTROL$curr_lik <- pop_lik
    ## Set up initial csv file
    chain_colnames <- c("sampno",names(all_pop_pars),"lnlike")
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(0,as.numeric(all_pop_pars),pop_lik)
    colnames(tmp_table) <- chain_colnames
    ## Write starting conditions to file
    write.table(tmp_table,file=mcmc_pop_filename,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)
    print("Initial population parameters saved")
    ############################################################################













    
    ############################################################################
    ############################################################################
    ## ACTUAL CHAIN NOW
    ############################################################################
    ############################################################################
    ## Go through chain
    proposed_ind_liks <- numeric(n_individuals)
    print("Starting MCMC chain....")
    for (i in 1:(iterations+adaptive_period+burnin)){
        pop_lik <- 0
        
        ############################################################################
        ## INDIVIDUAL UPDATES
        ############################################################################
        ## Update for each individual
        for(j in 1:n_individuals){
            chain_i <- ALL_INDIVIDUALS_CONTROL[[j]]$chain_i
            opt_i <- ALL_INDIVIDUALS_CONTROL[[j]]$opt_i
            ALL_INDIVIDUALS_CONTROL[[j]]$iter <- ALL_INDIVIDUALS_CONTROL[[j]]$iter + 1

        
            ######################################################
            ## MULTIVARIATE NORMAL PROPOSAL
            ######################################################
            #proposal <- mvr_proposal(ALL_INDIVIDUALS_CONTROL[[j]]$pars, ALL_INDIVIDUALS_CONTROL[[j]]$cov,ALL_INDIVIDUALS_CONTROL[[j]]$step_scale,ind_fixed)
            if(findY0s){
                tmp <- mvr_proposal(ALL_INDIVIDUALS_CONTROL[[j]]$pars, ALL_INDIVIDUALS_CONTROL[[j]]$cov,ALL_INDIVIDUALS_CONTROL[[j]]$step_scale,ind_fixed, ALL_INDIVIDUALS_CONTROL[[j]]$y0s)
                proposal <- tmp[["proposed"]]
                prop_y0s <- tmp[["y0s"]]
            } else {
                proposal <- mvr_proposal(ALL_INDIVIDUALS_CONTROL[[j]]$pars, ALL_INDIVIDUALS_CONTROL[[j]]$cov,ALL_INDIVIDUALS_CONTROL[[j]]$step_scale,ind_fixed, NULL)
                prop_y0s <- ALL_INDIVIDUALS_CONTROL[[j]]$y0s
            }
            ######################################################
            ## Check that proposals are within allowable range
            if(ind_prior(proposal,tmax,prop_y0s) > -0.5){
                ######################################################
                ## LIKELIHOOD CALCULATION FOR THIS INDIVIDUAL
                ######################################################
                lnlike <- proposed_ind_liks[j] <- likelihood_ind(proposal, prop_y0s, POP_PARAMETER_CONTROL$pars,all_data[[j]])
                log_prob <- min(lnlike - ALL_INDIVIDUALS_CONTROL[[j]]$curr_lik,0)
                ######################################################

                ######################################################
                ## Check for acceptance. If better, or if worse and proportional to how worse accept the move
                if(log(runif(1)) < log_prob){
                    ALL_INDIVIDUALS_CONTROL[[j]]$pars <- proposal
                    if(findY0s) ALL_INDIVIDUALS_CONTROL[[j]]$y0s <- prop_y0s
                    ALL_INDIVIDUALS_CONTROL[[j]]$curr_lik <- ind_liks[j] <- lnlike  ## Record new likelihood
                    ALL_INDIVIDUALS_CONTROL[[j]]$accepted <- ALL_INDIVIDUALS_CONTROL[[j]]$accepted + 1
                }

                ## Change step size to reach better acceptance rate if within adpative part. Also only adapt if proposal was valid
                if(i > burnin && i < (burnin + adaptive_period)){
                    ALL_INDIVIDUALS_CONTROL[[j]]$step_scale <- rm_scale(ALL_INDIVIDUALS_CONTROL[[j]]$step_scale, i-burnin, popt, log_prob, adaptive_period)
                    pcur <- ALL_INDIVIDUALS_CONTROL[[j]]$accepted/ALL_INDIVIDUALS_CONTROL[[j]]$iter
                    #ALL_INDIVIDUALS_CONTROL[[j]]$step_scale <- scaletuning(ALL_INDIVIDUALS_CONTROL[[j]]$step_scale,popt,pcur)
                }
            }
            ##################################################
            ## ADAPTIVE PERIOD
            ##################################################
            ## If within adaptive period, record chain more regularly for adaptation
            if(i < (adaptive_period+burnin) && i > burnin){
                toSave <- ALL_INDIVIDUALS_CONTROL[[j]]$pars
                if(findY0s) toSave <- c(toSave, ALL_INDIVIDUALS_CONTROL[[j]]$y0s)
                ALL_INDIVIDUALS_CONTROL[[j]]$opt_chain[opt_i,] <- c(i,toSave, ALL_INDIVIDUALS_CONTROL[[j]]$curr_lik)
                ALL_INDIVIDUALS_CONTROL[[j]]$opt_i <- ALL_INDIVIDUALS_CONTROL[[j]]$opt_i + 1
                
                ## If this is an optimisation step, scale the step size and update the covariance matrix
                ## Note that we don't tune covariance matrix for first 10% of tuning steps to give time to build up some moves
                if(ALL_INDIVIDUALS_CONTROL[[j]]$opt_i%%opt_freq == 0 && ALL_INDIVIDUALS_CONTROL[[j]]$opt_i > (0.10*tuning+burnin) && i < (burnin + tuning)){
                    covMat <- cov(ALL_INDIVIDUALS_CONTROL[[j]]$opt_chain[1:(ALL_INDIVIDUALS_CONTROL[[j]]$opt_i-1),2:(ncol(chain)-1)]) ## Only take covariance for the actual parameters, and only up to current step
                    ALL_INDIVIDUALS_CONTROL[[j]]$cov <- covMat
                }
            }
            ##################################################

            ## Add to overall population likelihood
            pop_lik <- pop_lik + ALL_INDIVIDUALS_CONTROL[[j]]$curr_lik
            ######################################################
            

            ######################################################
            ## If this iteration is meant to be recorded, save it to the MCMC chain
            ######################################################
            if(i %% thin == 0){
                toSave <- ALL_INDIVIDUALS_CONTROL[[j]]$pars
                if(findY0s) toSave <- c(toSave, ALL_INDIVIDUALS_CONTROL[[j]]$y0s)
                ALL_INDIVIDUALS_CONTROL[[j]]$save_chain[chain_i,] <- c(i,toSave, ALL_INDIVIDUALS_CONTROL[[j]]$curr_lik)
                ALL_INDIVIDUALS_CONTROL[[j]]$chain_i <- ALL_INDIVIDUALS_CONTROL[[j]]$chain_i + 1
                ## If enough rows have been recorded, write these to csv and reset records
                if(ALL_INDIVIDUALS_CONTROL[[j]]$chain_i > save_block){
                    write.table(ALL_INDIVIDUALS_CONTROL[[j]]$save_chain,file=ALL_INDIVIDUALS_CONTROL[[j]]$filename,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
                    ALL_INDIVIDUALS_CONTROL[[j]]$chain_i <- 1
                }
            }
        }
        POP_PARAMETER_CONTROL$curr_lik <- pop_lik
        if(i%%save_block == 0) print(i)        

        ############################################################################
        ## POPULATION UPDATES
        ############################################################################
        ## If after the point of population updates
        if(i > N_pop_start){
            POP_PARAMETER_CONTROL$iter <- POP_PARAMETER_CONTROL$iter + 1
            chain_i <- POP_PARAMETER_CONTROL$chain_i
            opt_i <- POP_PARAMETER_CONTROL$opt_i
            
            ## Proposal
            pop_proposal <- mvr_proposal(POP_PARAMETER_CONTROL$pars,POP_PARAMETER_CONTROL$cov,POP_PARAMETER_CONTROL$step_scale,pop_fixed)

            ## If proposal is allowable
            if(pop_prior(pop_proposal) > -0.5){

                ## Sum up likelihoods for each individual
                for(j in 1:length(all_data)){
                    proposed_ind_liks[j] <- likelihood_ind(ALL_INDIVIDUALS_CONTROL[[j]]$pars, ALL_INDIVIDUALS_CONTROL[[j]]$y0s, pop_proposal, all_data[[j]])
                }
                lks <- sum(proposed_ind_liks)
                log_prob <- min(lks-POP_PARAMETER_CONTROL$curr_lik,0)
                
                if(log(runif(1)) < log_prob){
                    POP_PARAMETER_CONTROL$curr_lik <- lks
                    POP_PARAMETER_CONTROL$pars <- pop_proposal
                    for(j in seq_along(ALL_INDIVIDUALS_CONTROL)){
                        ALL_INDIVIDUALS_CONTROL[[j]]$curr_lik <- ind_liks[j] <- proposed_ind_liks[j]
                    }
                    POP_PARAMETER_CONTROL$accepted <- POP_PARAMETER_CONTROL$accepted + 1
                    
                }
                
                ## If this is an optimisation step, scale the step size and update the covariance matrix
                if(i < (adaptive_period+burnin) && i > burnin){
                    #pcur <- POP_PARAMETER_CONTROL$accepted/POP_PARAMETER_CONTROL$iter
                    #scale <- scaletuning(POP_PARAMETER_CONTROL$step_scale,popt,pcur)
                    scale <- rm_scale(POP_PARAMETER_CONTROL$step_scale, (i-N_pop_start), popt, log_prob, (adaptive_period-N_pop_start))
                    ##POP_PARAMETER_CONTROL$step_scale <- min(scale,0.1)
                    POP_PARAMETER_CONTROL$step_scale <- scale
                    #print(paste("Scale: ", scale,sep=""))
                    #print(paste("Pcur: " ,pcur,sep=""))
                    POP_PARAMETER_CONTROL$iter <- POP_PARAMETER_CONTROL$accepted <- 0
                }
            }
            
            if(i < (adaptive_period+burnin) && i > burnin){
                POP_PARAMETER_CONTROL$opt_chain[opt_i,] <- c(i,POP_PARAMETER_CONTROL$pars, POP_PARAMETER_CONTROL$curr_lik)
                POP_PARAMETER_CONTROL$opt_i <- POP_PARAMETER_CONTROL$opt_i + 1
                if(POP_PARAMETER_CONTROL$opt_i %% opt_freq == 0 && POP_PARAMETER_CONTROL$opt_i > (0.1*tuning + burnin) && i < (burnin + tuning)){
                                        #print(paste("Pcur: ",POP_PARAMETER_CONTROL$accepted/POP_PARAMETER_CONTROL$iter,sep=""))
                    covMat <- cov(POP_PARAMETER_CONTROL$opt_chain[1:(POP_PARAMETER_CONTROL$opt_i-1),2:(ncol(POP_PARAMETER_CONTROL$opt_chain)-1)])
                    #print(covMat[which(pop_fixed==0),which(pop_fixed==0)])
                    POP_PARAMETER_CONTROL$cov <- covMat
                }
            }
            ######################################################
            ## If this iteration is meant to be recorded, save it to the MCMC chain
            ######################################################
            if(i %% thin == 0){
                POP_PARAMETER_CONTROL$save_chain[chain_i,] <- c(i,POP_PARAMETER_CONTROL$pars, POP_PARAMETER_CONTROL$curr_lik)
                POP_PARAMETER_CONTROL$chain_i <- POP_PARAMETER_CONTROL$chain_i + 1
                ## If enough rows have been recorded, write these to csv and reset records
                if(POP_PARAMETER_CONTROL$chain_i > save_block){
                    write.table(POP_PARAMETER_CONTROL$save_chain,file=POP_PARAMETER_CONTROL$filename,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
                    POP_PARAMETER_CONTROL$chain_i <- 1
                }
            }
        }
    }
    return(TRUE)
}
