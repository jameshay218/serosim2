
#' Boost/wane model
#' 
#' Data generating function for a multiple strain, boosting and waning serological model.
#' @param tis vector of times of infection eg. c(0, 45, 70)
#' @param y0s vector of baseline HAI values for each infecting strain
#' @param mu_pars matrix of boosting parameters
#' @param tp_pars matrix of time to peak parameters
#' @param m_pars matrix of waning parameters
#' @param times vector of times in days.
#' @return matrix of times and longitundinal simulation for each infecting strain
#' @seealso \code{\link{single_strain}}
#' @export
multiple_strains <- function(tis, y0s, mu_pars,tp_pars,m_pars, times){
    #' Create matrix with a column for each infection and a column for times
    dat <- matrix(ncol=length(tis)+1, nrow=length(times))
    dat[,1] <- times

    #' Sort infection times into ascending order and store indices of sorted order. Add this to list of process parameters if needed
    actual_ti <- sort(tis,index.return=TRUE)
    indices <- actual_ti$ix
    tis <- actual_ti$x

    #' For each strain, calculate the single trajectory
    for(i in 2:(length(tis)+1)){
        dat[,i] <- single_strain(mu_pars[i-1,indices], tp_pars[i-1,indices],m_pars[i-1,indices], tis, y0s[i-1],times)
    }
    return(dat)    
}


#' Single strain boost/wane model
#'
#' Simulates the time series antibody kinetics for a single strain from arbritary infection events
#' @param mu_pars vector of boosting parameters
#' @param tp_pars vector of time to peak parameters (days)
#' @param m_pars vector of waning rate parameters
#' @param ti_pars vector of infection times
#' @param times vector of times to run model over
#' @return vector of HAI values at each day of the simulation
#' @seealso \code{\link{multiple_strains}}
#' @export
single_strain <- function(mu_pars, tp_pars, m_pars, ti_pars, y0, times){
    ii <- 1
    y <- NULL
    lower_bound <- -10

    for(i in 1:length(ti_pars)){
        if(i == length(ti_pars)){
            final_t = times[length(times)]
        } else {
            final_t <- ti_pars[i+1]
        }
        ti <- ti_pars[i]
        mu <- mu_pars[i]
        tp <- tp_pars[i]
        m <- m_pars[i]
        while(ii < (length(times)+1) & times[ii] <= final_t){
            t <- times[ii]
            tmp <- 0
            if(t <= ti) tmp <- y0
            else if(t > ti & t <= (ti + tp)) tmp <- (mu/tp)*t - (mu/tp)*ti + y0
            else tmp <- -m*t + m*(ti+tp) + mu + y0        
            y[ii] <- tmp
            if(y[ii] < lower_bound) y[ii] <- lower_bound
            ii <- ii + 1
        }
        if(i < length(ti_pars)){
            if(ti_pars[i+1] <= ti + tp) y0 <- (mu/tp)*ti_pars[i+1] - (mu/tp)*ti + y0
            else y0 <- m*(ti + tp - ti_pars[i+1]) + mu + y0

            if(y0 < lower_bound) y0 <- lower_bound
        }
    }
    return(y)
}
