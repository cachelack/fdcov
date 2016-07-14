#' Multiple-sample permutation test for the equality of covariance operators of functional data
#' 
#' The method performs a test for the equality of the covariance operators of multiple data samples.
#' It can also perform all of the pairwise comparisons between the groups and compute a p-value for each of them.
#' This feature is useful when the global null hypothesis is rejected, so one may want to find out which samples have different covariances. 
#' 
#' @param dat n X p data matrix of n samples of p long vectors.
#' @param grp n long vector of group labels.
#' @param iter Number of permutations. Defaults to 1000.
#' @param perm Type of permutation, can be 'sync' (if all the data samples are of the same size) or 'pool'. Defaults to 'sync'
#' @param dist Distance between covariance operators. Can be 'sq' (square-root distance), 'tr' (trace distance),'pr' (Procrustes distance), 'hs'(Hilbert-Schmidt distance) or 'op' (operator distance). Defaults to 'sq'.
#' @param adj p-value adjustment. Defaults to TRUE.
#' @param comb Can be 'tipp' (for Tippett), 'maxT', 'dire' (direct), 'fish' (Fisher) or 'lipt' (Liptak). Defaults to 'tipp'.
#' @param part If FALSE, the function computes only the global p-value; otherwise it computes also all the p-values corresponding to the pairwise comparisons. Defaults to FALSE.
#' @param cent If FALSE, the mean functions of the groups are supposed to be different, therefore data are centred before performing the test. Defaults to FALSE.
#' @param load Boolean flag, which if TRUE, prints a loading bar.
#' 
#' @return If \code{part} is set to FALSE, the output is the p-value associated to the global test. If \code{part} is TRUE, the function returns also all the p-values of the pairwise comparisons.
#' 
#' @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
#' 
#' @references Pigoli, Davide, John A. D. Aston, Ian L. Dryden, and Piercesare Secchi (2014). "Distances and inference for covariance operators." Biometrika: 101(2):409–422.
#' 
#' @examples 
#' \dontrun{
#' ## Phoneme data
#' 
#' library(fdcov)
#' library(fds)
#' 
#' # Create data set 
#' data(aa); data(ao); data(dcl);data(iy);data(sh)
#' dat = cbind(aa$y[,1:20],ao$y[,1:20],dcl$y[,1:20],iy$y[,1:20],sh$y[,1:20])
#' dat = t(dat)
#' grp = c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))
#' 
#' # Test the equality of the covariance operators
#' p = ksample.perm(dat, grp, iter=100, part = TRUE)
#' p$global # global p-value
#' p$partial # partial p-values
#' }
#' \dontshow{
#' library(fdcov)
#' library(fds)
#' 
#' # Create data set 
#' data(aa); data(ao); data(dcl);data(iy);data(sh)
#' dat = cbind(aa$y[,1:20],ao$y[,1:20],dcl$y[,1:20],iy$y[,1:20],sh$y[,1:20])
#' dat = t(dat)
#' grp = c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))
#' 
#' # Test the equality of the covariance operators
#' p = ksample.perm(dat,grp,iter=2)
#' }
#' @export

ksample.perm = function(dat, grp, iter = 1000, perm = 'sync', dist = 'sq', adj = TRUE,
                        comb = 'tipp', part = FALSE, cent = FALSE, load = FALSE){

    table_groups = table(grp) # groups table
    C = length(table_groups)  # number of groups
    if(C < 2) stop('The number of groups must be at least 2.')
    if(C == 2) return(list(global = twosample.perm(dat, grp, iter, dist, load)))
       
    ### Step 1: data centring
  
    if(cent == FALSE){
        nonalign_dat = dat
        for(i in 1:length(grp)){
            dat[i,] = nonalign_dat[i,]-colMeans(nonalign_dat[grp == grp[i],],na.rm=TRUE)
        }
    }
     
    ### Steps 2, 3 and 4: apply iter permutations and compute the test statistic for each permuted data set
    
    T = switch (perm,
                    sync = perm.sync(dat, grp, iter, dist, load),
                    pool = perm.pool(dat, grp, iter, dist, load),
                    stop('The selected permutation strategy is not available'))  
  
    ### Global p-values
  
    if (comb == 'tipp' | comb == 'fish' | comb == 'lipt' | adj == FALSE) P = perm.t2p(T)
    
    P.glob = switch (comb,
                        tipp = mean(min(P[1,])       >= apply(P[-1,],1,min)),
                        maxT = mean(max(T[1,])       <= apply(T[-1,],1,max)),
                        dire = mean(sum(T[1,])       <= apply(T[-1,],1,sum)),
                        fish = mean(comb.fish(P[1,]) <= apply(P[-1,],1,comb.fish)),
                        lipt = mean(comb.lipt(P[1,]) <= apply(P[-1,],1,comb.lipt)),
                        stop('The selected combining function is not available'))
          
    ### Partial p-values
    
    if(part == TRUE){
        if(adj == TRUE){ # Corrected p-values
            P.part = switch (comb,
                                tipp = FWE.tipp(P),
                                maxT = FWE.maxT(T),
                                dire = FWE.clos(T, comb, load),
                                fish = FWE.clos(P, comb, load),
                                lipt = FWE.clos(P, comb, load))
        }else{P.part = matrix(P[1,],nrow = C*(C-1)/2,ncol = 1) # Raw p-values
        }
        
        # Group names definition
        name = character()
        cont = 1
        for(i in 1:(C-1)){
            for(j in (i+1):C){
                name[cont] = paste(i,j,sep="-")
                cont = cont+1
            }
        }
        
        P.part = data.frame(dist = T[1,],p_value = P.part, signif = perm.sig(P.part))
        rownames(P.part) = name
       
        return(list(global = P.glob,partial = P.part))
    }
    else{return(list(global = P.glob))} 
}
