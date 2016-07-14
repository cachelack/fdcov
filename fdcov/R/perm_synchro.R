# Synchronised permutations of a data set
# 
# @param dat n x p data set matrix
#
# @return This function returns a matrix containing the partial test statistics of each pairwise comparison for each permutation applied to the data set
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
#

perm.sync = function(dat, grp, iter, dist, load = FALSE){

    table_groups = table(grp) # groups table
    C = length(table_groups)  # number of groups
    K = C*(C-1)/2             # number of partial tests
    p = dim(dat)[2]           # number of samplings per function
    
    T = array(5,dim = c((iter+1),K)) # test statistics vector initialisation
    
    ### Step 2: compute test statistic 
    
    cont = 1
    for(i in 1:(C-1)){ 
        for(j in (i+1):C){ # for each pair of groups
            # compute test statistic for initial data
            T[1,cont] = distCov(cov(dat[grp==i,],use='pairwise'),cov(dat[grp==j,],use='pairwise'),dist)
            cont = cont+1
        }
    }
    
    ### Steps 3 and 4: apply iter permutations and compute the test statistic for each permuted data set

    # Build pseudomatrix
    
    n = (dim(dat)[1])/C
    X = array(0,dim = c((2*n),p,K)) # matrix 2*n X p X number of groups
    cont = 1
    for(i in 1:(C-1)){ # for each pair of groups
        for(j in (i+1):C){
            X[,,cont] = rbind(dat[grp==i,],dat[grp==j,]) # fill in the matrix
            cont = cont+1
        }
    }
    
    
    # Apply 'iter' permutations to the pseudomatrix and compute test statistics
    
    if(load) pb = txtProgressBar(min = 0, max = iter, style = 3) # create progress bar
    for(bb in 2:(iter+1)){
        X.perm = X[sample(2*n),,] # select permutation      
        cont = 1
        for(i in 1:(C-1)){ # for each pair of groups
            for(j in (i+1):C){ # compute test statistic for the permuted dataset
                T[bb,cont] = distCov(cov(X.perm[c(1:n),,cont],use = 'pairwise'),cov(X.perm[-c(1:n),,cont],use = 'pairwise'),dist)
                cont = cont+1
            }
        }
        if(load) setTxtProgressBar(pb, bb-1) # update progress bar 
    } # end iter
    if(load) close(pb) # close progress bar
    
    return(T)
}