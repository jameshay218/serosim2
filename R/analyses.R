#' Save ind pars
#'
#' Formats a list of individual parameters to save as table
#' @param ind_pars the list of parameters
#' @param tmax max time corresponding to not infected
#' @return a table of pars
#' @export
save_ini_pars <- function(ind_pars, tmax){
    all_pars <- NULL
    for(i in seq_along(ind_pars)){
        all_pars <- rbind(all_pars, c(i, signif(ind_pars[[i]]$pars[1],3),ind_pars[[i]]$pars[2:4],ind_pars[[i]]$y0s))
    }
    colnames(all_pars) <- c("individual","mu_i","tis1","tis2","tis3","y0_1","y0_2","y0_3")
    all_pars[,c("tis1","tis2","tis3")][all_pars[,c("tis1","tis2","tis3")] == tmax] <- "Not Infected"
    
    return(all_pars)
}

#' Save data
#'
#' Formats list of data to single table
#' @param data the list of data frames
#' @return a single table of formatted data
#' @export
save_all_data <- function(data){
    allDat <- NULL
    for(i in seq_along(data)){
        tmp <- as.matrix(data[[i]])
        tmpDat <- NULL
        viruses <- colnames(tmp)
        for(j in 2:ncol(tmp)){
            tmpDat <- rbind(tmpDat, as.matrix(data[[i]])[,c(1,j)])
        }
        tmpDat <- cbind(i,rep(viruses[2:length(viruses)],each=nrow(tmp)),tmpDat)
        allDat <- rbind(allDat,tmpDat)
    }
    colnames(allDat) <- c("Individual","Strain","Time","Titre")
    return(allDat)    
}

#' Check for seroconversion
#'
#' Checks for seroconversion within an individual data frame
#' @param data data for each individual (ie. each individual has a matrix of sampling times and titres)
#' @param y0s vector of initial titre values (ie. at t = 0)
#' @return a vector of infected/not infected for each strain
#' @export
seroconversion_check <- function(data,y0s){
    tmpDat <- rbind(y0s,data[,2:ncol(data)])
    results <- character(ncol(tmpDat))
    for(i in 1:ncol(tmpDat)){
        results[i] <- "NOT INFECTED"
        for(j in 1:(nrow(tmpDat)-1)){
            if(tmpDat[j+1,i] >= tmpDat[j,i] + 2){
                results[i] <- "INFECTED"
            }
        }
    }
    return(results)
}

#' Get posterior propn
#'
#' Calculates proportion of posterior estimate from chain that indicates infection
#' @export
get_infection_ratios <- function(chain, burnin, threshold, tmin,tmax){
    infection_names <- c("ti1","ti2","ti3")
    chain <- chain[chain$sampno > burnin,infection_names]

    results <- data.frame(nrow=length(infection_names),ncol=2)
    for(i in seq_along(infection_names)){
        propn <- length(chain[(chain[,infection_names[i]] < tmax),infection_names[i]])/nrow(chain)
        label <- "NEGATIVE"
        if(propn > threshold) label <- "POSITIVE"
        results[i,] <- c(propn,label)
    }
    colnames(results) <- c("posterior","infection")
    return(results)
    
    
}


#' @export
straight_to_results <- function(chain_names, all_ind_pars, all_data, burnin,threshold=0.5){
    return(generate_measures(generate_counts(calculate_statistics(chain_names,all_ind_pars,all_data,burnin,threshold))))
}

#' Calculate infection statistics from MCMC chains
#' 
#' Make sure you're in the working directory with all of the chains
#' @export
calculate_statistics <- function(chain_names, all_ind_pars, all_data, burnin,threshold=0.5){
    #'###########################
    #'# ANALYSIS
    #'###########################
    allResults <- NULL
    
    ## Go through each individual
    for(i in seq_along(chain_names)){
        tmpPars <- all_ind_pars[[i]]$pars[2:4]
        tmax <- max(all_data[[i]][,1])
        tmin <- min(all_data[[i]][,1])
        actualResults <- !(tmpPars >= tmax)
        
        tmp <- get_infection_ratios(data.table::fread(chain_names[i],data.table=FALSE),burnin,threshold,tmin,tmax)
        seroconversion <- seroconversion_check(all_data[[i]],all_ind_pars[[i]]$y0s)
        allResults <- rbind(allResults, cbind("individual"=i,"virus"=c("A","B","C"),tmp,actualResults, seroconversion))
    }
    colnames(allResults) <- c("Individual","Virus","Posterior Proportion","Model Predicted","Actual Infection", "Seroconversion Predicted")


    allResults$model_result <- NULL
    allResults$seroconversion_result <- NULL

    for(i in 1:nrow(allResults)){
        if(allResults[i,"Actual Infection"] == FALSE){ ## If no infection
            if(allResults[i,"Model Predicted"]=="NEGATIVE"){ ## If no predicted infection
                allResults[i,"model_result"] <- "TRUE NEGATIVE"
            } else {
                allResults[i,"model_result"] <- "FALSE POSITIVE"
            }
            
            if(allResults[i,"Seroconversion Predicted"] == "NOT INFECTED"){
                allResults[i,"seroconversion_result"] <- "TRUE NEGATIVE"
            } else {
                allResults[i,"seroconversion_result"] <- "FALSE POSITIVE"
            }
        } else { ## If infection
            if(allResults[i,"Model Predicted"]=="POSITIVE"){
                allResults[i,"model_result"] <- "TRUE POSITIVE"
            } else {
                allResults[i,"model_result"] <- "FALSE NEGATIVE"
            }
            
            if(allResults[i,"Seroconversion Predicted"] == "INFECTED"){
                allResults[i,"seroconversion_result"] <- "TRUE POSITIVE"
            } else {
                allResults[i,"seroconversion_result"] <- "FALSE NEGATIVE"
            }
        }
    }
    return(allResults)
}

#' @export
generate_counts <- function(results){
    sero_results <- table(results[,"seroconversion_result"])
    model_results <- table(results[,"model_result"])
    final <- data.frame("FALSE NEGATIVE" = c(sero_results["FALSE NEGATIVE"], model_results["FALSE NEGATIVE"]),
                        "FALSE POSITIVE" = c(sero_results["FALSE POSITIVE"], model_results["FALSE POSITIVE"]),
                        "TRUE NEGATIVE" = c(sero_results["TRUE NEGATIVE"], model_results["TRUE NEGATIVE"]),
                        "TRUE POSITIVE" = c(sero_results["TRUE POSITIVE"], model_results["TRUE POSITIVE"]),
                        row.names=NULL
                        )
    final[is.na(final)] <- 0
    colnames(final) <- c("FALSE NEGATIVE","FALSE POSITIVE","TRUE NEGATIVE","TRUE POSITIVE")
    return(as.matrix(final))
}

#' @export
generate_measures <- function(counts){
    sero_results <- counts[1,]
    model_results <- counts[2,]
    sensitivity_sero <- unname(sero_results["TRUE POSITIVE"]/(sero_results["TRUE POSITIVE"]+sero_results["FALSE NEGATIVE"]))
    specificity_sero <- unname(sero_results["TRUE NEGATIVE"]/(sero_results["TRUE NEGATIVE"]+sero_results["FALSE POSITIVE"]))
    precision_sero <- unname(sero_results["TRUE POSITIVE"]/(sero_results["TRUE POSITIVE"]+sero_results["FALSE POSITIVE"]))

    sensitivity_model <- unname(model_results["TRUE POSITIVE"]/(model_results["TRUE POSITIVE"]+model_results["FALSE NEGATIVE"]))
    specificity_model <- unname(model_results["TRUE NEGATIVE"]/(model_results["TRUE NEGATIVE"]+model_results["FALSE POSITIVE"]))
    precision_model <- unname(model_results["TRUE POSITIVE"]/(model_results["TRUE POSITIVE"]+model_results["FALSE POSITIVE"]))

    results <- data.frame("Method"=c("Seroconversion","Model"),"Sensitivity"=c(sensitivity_sero,sensitivity_model),"Specificity"=c(specificity_sero,specificity_model),"Precision"=c(precision_sero,precision_model))
    return(results)
}

#' @export
simulate_and_save <- function(n, sample_start, tmax, samples, strains, crs, tp, m, infection_propns, infection_ranges, mu, mu_sigma, error_pars,filename="test",addNoise=FALSE, useY0s=FALSE){
    
    startTime <- sample_start
    endTime <- tmax
    nstrains <- length(strains)
    
    infection_risks <- NULL
    for(i in seq_along(infection_ranges)){
        range <- infection_ranges[[i]]
        propn <- infection_propns[i]
        infection_risks[[i]] <- generate_toy_inc(propn,tmax,range[1],range[2])
    }

    pop_pars <- c("mu_mu"=mu,"mu_sigma"=mu_sigma,
              "error_sigma"=unname(error_pars["sigma"]),"S"=unname(error_pars["S"]),"EA"=unname(error_pars["EA"]),
              "cr12"=crs[1],"cr13"=crs[2],"cr23"=crs[3],
              "cr12_sigma"=1,"cr13_sigma"=1,"cr23_sigma"=1,
              "tp"=tp,"tp12"=tp,"tp13"=tp,"tp23"=tp,
              "tp_sigma"=1,"tp12_sigma"=1,"tp13_sigma"=1,"tp23_sigma"=1,                  
              "m"=m,"m12"=m,"m13"=m,"m23"=m,
              "m_sigma"=1,"m12_sigma"=1,"m13_sigma"=1,"m23_sigma"=1
              )
    cr <- matrix(c(1,crs[1],crs[2],crs[1],1,crs[3],crs[2],crs[3],1),nrow=nstrains,ncol=nstrains)
    tp <- matrix(rep(tp,nstrains^2),nrow=nstrains,ncol=nstrains)
    m <- matrix(rep(m,nstrains^2),nrow=nstrains,ncol=nstrains)

    processParams <- list("mu"=mu,"mu_sigma"=mu_sigma,"tp"=tp,"cr"=cr,"m"=m)
    
    y0s <- matrix(0,ncol=3,nrow=n)
    
    samplingTimes <- matrix(0,ncol=samples,nrow=n)
    for(i in 1:nrow(samplingTimes)){
        samplingTimes[i,] <- sort(sample(seq(startTime,endTime,by=1),samples))
    }
    
    y <- overall_simulation(
        individuals = n,
        strains = strains,
        strainIncidences=infection_risks,
        samplingTimes=samplingTimes,
        processParams=processParams,
        y0s=y0s,
        startTime=0,
        endTime=endTime,
        addNoise=addNoise,
        noiseParams=error_pars[c("S","EA","max_titre")],
        discreteData=TRUE,
        multiple_strains,
        add_noise_func,
        TRUE
        )

    all_data <- y$all_data
    all_ind_pars <- y$all_ind_pars
    listDat <- y$listDat

    param_save <- save_ini_pars(all_ind_pars,tmax+1)
    data_save <- save_all_data(all_data)

    filename1 <- paste(filename,"_params.csv",sep="")
    filename2 <- paste(filename,"_data.csv",sep="")
    filename3 <- paste(filename,"_pop_params.csv",sep="")

    save_pop_pars <- c(pop_pars, infection_propns, do.call("c",infection_ranges))
    names(save_pop_pars) <- c(names(pop_pars),paste(strains,"_propn",sep=""),rep(c("start","end"),length(strains)))

    write.csv(param_save,filename1)
    write.csv(data_save,filename2)
    write.csv(t(save_pop_pars),filename3)
    
    return(list("all_data"=all_data,"all_ind_pars"=all_ind_pars,"listDat"=listDat,"pop_pars"=pop_pars))
}



