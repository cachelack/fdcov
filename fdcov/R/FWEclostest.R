# Closed testing procedure for strong FWE control
# 
# @param T (iter+1) X k matrix of permutation test statistics
# @param comb Combining function
# 
# @return A vector of adjusted p-values
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it} 
# 
# @export

FWE.clos<-function(T, comb = 'dire', loading = FALSE){
  
    #library(matlab)
    
    iter = dim(T)[1]      # number of permutations
    k = dim(T)[2]         # number of partial tests
    rows = ncycles = 2^k  # total number of tests
 
    x = matrix(0,rows,k)
  
    if(loading) pb = txtProgressBar(min = 0, max = k, style = 3) # create progress bar
    for (i in 1:k){  
        nreps = rows/ncycles
        ncycles = ncycles/2 
        zo = matrix(c(0,1),nreps,2,byrow=TRUE)
        zoc = rbind(as.matrix(zo[,1]),as.matrix(zo[,2]))
        settings = matlab::repmat(zoc,c(1,ncycles)) 
        x[,k-i+1] = settings
        if(loading) setTxtProgressBar(pb, k) # update progress bar 
    }
    if(loading) close(pb) # close progress bar
    x = x[-1,]
  
    print(T[1:10,])
    
    T2 = matrix(0,iter,(rows-1))       
    for(j in 1:(rows-1)){ 
        # for (i in 1:(iter+1)){
            T2[,j] = switch(comb,
                                dire = apply(as.matrix(T[,x[j,]==1]),1,sum), # [i]
                                fish = apply(as.matrix(T[,x[j,]==1]),1,comb.fish), # [i]
                                lipt = apply(as.matrix(T[,x[j,]==1]),1,comb.lipt), # [i]
                                warning('The selected combining function is not available'))
        # }
    }  
    
    print(T2[1:5,1:5])
    
    rawP = apply(T2,2,perm.pval)
    adjP = rep(5,k)
    for(l in 1:k) adjP[l] = max(rawP[x[,l]==1])
  
    return(adjP)
}
