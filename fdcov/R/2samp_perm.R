# 2-sample permutation test for the equality of covariance operators
# 
# This method performs a test for the equality of the covariance operators of 2 data samples.
# 
# @param dat n X p data matrix of n samples of p long vectors.
# @param grp n long vector of group labels.
# @param iter Number of permutations. Defaults to 1000.
# @param dist Distance between covariance operators. Can be 'sq' (square-root distance), 'tr' (trace distance),'pr' (Procrustes distance), 'hs'(Hilbert-Schmidt distance) or 'op' (operator distance). Defaults to 'sq'.
# @param cent If FALSE, the mean functions of the groups are supposed to be different, therefore data are centred before performing the test. Defaults to FALSE.
# 
# @return The function returns the p-value associated to the test.
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
# 
# @references Pigoli, Davide, John A. D. Aston, Ian L. Dryden, and Piercesare Secchi (2014). "Distances and inference for covariance operators." Biometrika: asu008.
# 
# @examples 
# \dontrun{
# ## Phoneme data
# 
# library(fdcov)
# library(fds)
# 
# # Create data set 
# data(aa); data(sh)
# dat = cbind(aa$y[,1:20],sh$y[,1:20])
# dat = t(dat)
# grp = c(rep(1,20),rep(2,20))
# 
# # Test the equality of the covariance operators
# p = twosample.perm(dat,grp)
# p # p-value
# }
# 
# 
twosample.perm = function(dat, grp, iter = 1000, dist = 'sq', cent = FALSE, load = FALSE){
    
    table_groups = table(grp) # groups table
    C = length(table_groups)  # number of groups
    try(if(C!=2) stop('The number of groups must be equal to 2. If k>2, the  function ksample.perm can be used.'))
    N = length(grp)
    
    ### Step 1: data centring
    
    if(cent == FALSE){
        nonalign_dat = dat
        for(i in 1:length(grp)){
            dat[i,] = nonalign_dat[i,]-colMeans(nonalign_dat[grp==grp[i],],na.rm=TRUE)
        }
    }
    
    T = rep(5,iter+1) # test statistics vector initialisation
    T[1] = distCov(cov(dat[grp==1,],use='pairwise'),cov(dat[grp==2,],use='pairwise'),dist)
    
    if(load) pb = txtProgressBar(min = 0, max = iter, style = 3) # create progress bar
    for (bb in 1:iter){
        dat.perm = dat[sample(N),]
        T[bb+1] = distCov(cov(dat.perm[grp==1,],use='pairwise'),cov(dat.perm[grp==2,],use='pairwise'),dist)
        if(load) setTxtProgressBar(pb, bb-1) # update progress bar 
    }
    if(load) close(pb) # close progress bar
    
    return(mean(T[-1]>=T[1]))
}