#' @export
rm_scale <- function(step_scale, mc, popt,log_prob, N_adapt)
{
	dd <- exp(log_prob)
	if( dd < -30 ){ dd <- 0 }
	dd <- min( dd, 1 )

	rm_temp <- ( dd - popt )/( (mc+1)/(0.01*N_adapt+1) )
	
	out <- step_scale*exp(rm_temp)
	
	out <- max( out, 0.02 )
	out <- min( out, 2)
	out
}

#' @export
mvr_proposal <- function(pars, cov, step,fixed,y0s=NULL){
    indices <- which(fixed==0)
    if(length(pars[indices]) == 0) return(pars)
    proposed <- pars

    if(!is.null(y0s)){
        proposed <- c(proposed, y0s)
        y0_indices <- c(length(proposed)-2,length(proposed)-1,length(proposed))
        indices <- c(indices, y0_indices)
    }
    #proposed[indices]  <- pars[indices] + step*mvrnorm(n=1,mu=rep(0,length(pars[indices])),cov[indices,indices])
    proposed[indices] <- mvrnorm(n=1,mu=proposed[indices],Sigma=step*cov[indices,indices])

    if(!is.null(y0s)){
        y0s <- proposed[y0_indices]
        proposed <- proposed[!(seq_along(proposed) %in% y0_indices)]
        return(list("proposed"=proposed,"y0s"=y0s))
    }
    
    return(proposed)
}

#' @export
y0_proposal <- function(y0s){
    newY0s <- rpois(3, y0s)
    return(y0s)
}

#' @export
ind_prior <- function(pars, tmax, y0s=NULL){
    if(any(pars < 0)) return(-1)
    if(pars[2] != tmax+2 && abs(pars[2] - pars[3]) <= 21) return(-1)
    if(pars[3] != tmax+2 && abs(pars[3] - pars[4]) <= 21) return(-1)
    if(pars[4] != tmax+2 && abs(pars[2] - pars[4]) <= 21) return(-1)
    if(any(pars[2:4] > 2000)) return(-1)
    if(pars[1] > 50) return(-1)
    if(!is.null(y0s)){
        if(any(y0s < 0)) return(-1)
    }
    return(1)
}

#' @export
pop_prior <- function(pars){
    if(any(pars < 0)) return(-1)
    if(any(pars > 100)) return(-1)
    if(any(pars[6:8] > 1)) return(-1)
    return(1)
}


#' @export
make_likelihood <- function(ind_pars, y0s, pop_pars, data, times, oddLikelihood=FALSE){
    if(!oddLikelihood){
        f <- function(ind_pars, y0s, pop_pars, data){
            ## Population parameters
            tp_pop <- matrix(pop_pars[c("tp","tp12","tp13","tp12","tp","tp23","tp13","tp23","tp")],nrow=3,ncol=3)
            cr_pop <- matrix(c(1,pop_pars["cr12"],pop_pars["cr13"],pop_pars["cr12"],1,pop_pars["cr23"],pop_pars["cr13"],pop_pars["cr23"],1),nrow=3,ncol=3)
            m_pop <- matrix(pop_pars[c("m","m12","m13","m12","m","m23","m13","m23","m")],nrow=3,ncol=3)

            ## Error pars
            error_pars <- pop_pars[c("error_sigma","S","EA")]
            
            ## Individual parameters
            tis <- ind_pars[c("tis1","tis2","tis3")]
            mu_ind <- ind_pars["mu_i"]


            ## Mixed effects parameters
            mu_pop <- pop_pars["mu_mu"]
            mu_pop_sigma <- pop_pars["mu_sigma"]
            
            ## Population sigmas (if we want mixed effects for these)
            cr_sigmas <- c(pop_pars["cr12_sigma"],pop_pars["cr13_sigma"],pop_pars["cr23_sigma"])
            tp_sigmas <- c(pop_pars["tp_sigma"],pop_pars["tp12_sigma"],pop_pars["tp13_sigma"],pop_pars["tp23_sigma"])
            m_sigmas <- c(pop_pars["m_sigma"],pop_pars["m12_sigma"],pop_pars["m13_sigma"],pop_pars["m23_sigma"])
            
            ## Call to an Rcpp function
            return(posterior_1(tis, y0s, mu_ind, cr_pop, tp_pop, m_pop, data, error_pars, mu_pop, mu_pop_sigma))
        }
        return(f)
    } else {
f <- function(ind_pars, y0s, pop_pars, data){
            ## Population parameters
            tp_pop <- matrix(pop_pars[c("tp","tp12","tp13","tp12","tp","tp23","tp13","tp23","tp")],nrow=3,ncol=3)
            cr_pop <- matrix(c(1,pop_pars["cr12"],pop_pars["cr13"],pop_pars["cr12"],1,pop_pars["cr23"],pop_pars["cr13"],pop_pars["cr23"],1),nrow=3,ncol=3)
            m_pop <- matrix(pop_pars[c("m","m12","m13","m12","m","m23","m13","m23","m")],nrow=3,ncol=3)

            ## Error pars
            error_pars <- pop_pars[c("error_sigma","S","EA")]
            
            ## Individual parameters
            tis <- ind_pars[c("tis1","tis2","tis3")]
            mu_ind <- ind_pars["mu_i"]


            ## Mixed effects parameters
            mu_pop <- pop_pars["mu_mu"]
            mu_pop_sigma <- pop_pars["mu_sigma"]
            
            ## Population sigmas (if we want mixed effects for these)
            cr_sigmas <- c(pop_pars["cr12_sigma"],pop_pars["cr13_sigma"],pop_pars["cr23_sigma"])
            tp_sigmas <- c(pop_pars["tp_sigma"],pop_pars["tp12_sigma"],pop_pars["tp13_sigma"],pop_pars["tp23_sigma"])
            m_sigmas <- c(pop_pars["m_sigma"],pop_pars["m12_sigma"],pop_pars["m13_sigma"],pop_pars["m23_sigma"])
            
            ## Call to an Rcpp function
            return(posterior_4(tis, y0s, mu_ind, cr_pop, tp_pop, m_pop, data, error_pars, mu_pop, mu_pop_sigma))
        }
        return(f)

    }
}

#' @export
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    return(step)
}


