#'Expectation step
#'
#'Compute the expectations of censored FIB counts, resitricted to be less than the censoring threshold.
#'
#'@param alpha vector of intercepts where each entry is the intercept for one species/event combination
#'@param beta slope vector with one entry per FIB species, which indicates the marginal increase in that species' log-mean for a unit increase in the contamination index
#'@param gamma contamination vector with one entry per row of data, indicating the contamination index for the corresponding row
#'@param data matrix of FIB counts with one row per sample and one column per FIB species
#'@param event vector with one entry per row of data where each entry indicates the event with which the data row is associated
#'
#'@return A matrix of FIB counts, where the censored counts are replaced by their expectations under the current parameters.
#'
E.step = function(alpha, beta, gamma, data, event) {
    #Basic constants:
    n = nrow(data)
    p = ncol(data)
    d = length(unique(event))
    
    #Event-level alphas:
    alpha.local = matrix(0, n, p)
    for (k in 1:d) {
        indx = which(event==unique(event)[k])
        alpha.local[indx,] = matrix(alpha[p*(k-1) + 1:p], length(indx), p, byrow=TRUE)
    }
    
    for (j in 1:ncol(data))
        for (t in 1:nrow(data))
            if (!is.na(data[t,j]))
                if (data[t,j] <= cens[j]) {
                    #Get the probability that the count is 0,1,...,MLD
                    #Where MLD is the minimum level of detection (censoring threshold).
                    cens.log.lik = dpois((0:cens[j]), exp(gamma[t]*beta[j] + alpha.local[t,j]), log=TRUE)
                    cens.log.lik = cens.log.lik - max(cens.log.lik, na.rm=TRUE)
                    
                    #Compute the expected count, conditional on count being no greater than censoring threshold
                    #If the total probability below the censoring threshold is indistinguishable from zero,
                    #then set the expectation to the censoring threshold.
                    if (sum(exp(cens.log.lik))==0) {
                        data[t,j] = cens[j]
                    } else data[t,j] = sum((0:cens[j]) * exp(cens.log.lik), na.rm=TRUE) / sum(exp(cens.log.lik), na.rm=TRUE)
                }
    y
}
