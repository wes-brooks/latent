#' Estimate parameters of a latent variable model
#' 
#' Estimate the parameters of a latent variable model by the method of maximum likelihood, using conjugate gradient descent to maximize likelihood and an EM algorithm to handle censored data.
#' 
#' @param data matrix of observed counts, number of rows equal to the number of observations and columns equal to the number of categories.
#' @param min.detect minimum limit of detection for the FIB assay - values below this threshold are censored
#' @param event vector of event assignment labels, one entry per row of data. Indicates the event to which each row of data belongs.
#' @param specific vector of TRUE/FALSE values, one entry per column of data. Each entry indicates whether the corresponding FIB species is human-specific.
#' @param verbose should the function provide verbose output about its progress? Defaults to TRUE.
#' 
#' @return A list of results
#' 
#' @export
latent = function(data, min.detect, event, specific=NULL, verbose=TRUE) {
    #Initial parameters:
    xx = c(rep(as.integer(!specific), length(unique(event))), rep(1, nrow(data) + ncol(data)))
    finished = FALSE
    
    f.new = log.lik(data, xx, event)
    f.old = -Inf
    tol = .Machine$double.eps %>% sqrt
    tol = 1e-5
    check=Inf
    
    while (!finished) {
        #These iterations restart conjugacy:
        converged = FALSE
        while (!converged) {
            
            #Prepare to iterate conjugate gradient descent:
            i=0
            f.outer = f.old
            f.old = -Inf
            t = 1
            while(f.new>f.old && !converged && i<length(xx)) {
                i = i+1
        
                dir.new = score(data, xx, event, specific=specific)
                dir.new = dir.new / sqrt(sum(dir.new^2))
                
                #First iteration, ignore conjugacy - thereafter, use it.
                #s.new is the vector of the new step (in parameter space)
                if (i==1) {  
                    s.new = dir.new
                } else {
                    conj = (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
                    s.new = dir.new + conj * s.old
                }
                
                #Find the optimal step size
                #Backtracking: stop when the loss function is majorized
                condition = (log.lik(data, xx + t*s.new, event) < f.new + sum((t*s.new)*dir.new) + 1/(2*t)*sum((t*s.new)^2))[1]
                while( condition ) {
                    condition = (log.lik(data, xx + t*s.new, event) < f.new + sum((t*s.new)*dir.new) + 1/(2*t)*sum((t*s.new)^2))[1]
                    
                    #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
                    if (is.na(condition)) {
                        converged = TRUE
                        condition = FALSE
                    }
                    
                    t = 0.8*t
                }
                
                #Find the optimal step
                step = s.new * t
                
                #Make t a little bigger so next iteration has option to make larger step:
                t = t / 0.8 / 0.8
                p = xx + step
                
                #save for next iteration:
                dir.old = dir.new
                s.old = s.new
                
                #Only save the new parameters if they've decreased the loss function
                f.proposed = log.lik(data, p, event)
                if (f.proposed > f.old)
                    xx = p
                f.old = f.new
                f.new = f.proposed
            }
            
            if (verbose) cat(paste("Likelihood objective: ", f.new, "\n", sep=""))
            
            
            if ((f.new - f.outer) < tol * f.outer) converged = TRUE
        }
        
        d = length(unique(event))
        p = ncol(data)
        
        alpha = xx[1:(d*p)]
        beta = xx[d*p + 1:p]
        gamma = xx[((d+1)*p+1):length(xx)]
        data.new = E.step(alpha, beta, gamma, data, min.detect, event)
        
        check.old = check
        check = sum((data.new - data)^2, na.rm=TRUE) / sum(data^2, na.rm=TRUE)
        data = data.new
        
        if (check.old - check < (abs(check.old) + tol)*tol) {
            if (tol<=sqrt(.Machine$double.eps))
                finished = TRUE
            tol = max(tol/2, sqrt(.Machine$double.eps))
            cat(paste("Iterating with tol=", tol, "\n", sep=""))
        }
    }
    
    #Compile the results and return
    result = list()
    result$data = data
    result$min.detect = min.detect
    result$event = event
    result$specific = specific
    result$alpha = matrix(alpha, d, p, byrow=TRUE)
    colnames(result$alpha) = names(result$beta)
    result$beta = beta
    result$gamma = gamma
    
    return(result)
}