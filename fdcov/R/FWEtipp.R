# Step-down Tippett procedure for strong FWE control
# 
# @param P (iter+1) X k matrix of permutation p-values
# 
# @return A vector of adjusted p-values
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it} 

FWE.tipp = function(P){

    ord = order(P[1,],decreasing=FALSE)             # sort the observed p-values in increasing order and store the order
    P.ord = P[,ord]                                 # put the columns of matrix P in the new order orde
    
    k = dim(P)[2]                                   # number of tests
    p.ris = array(5,dim=c(k,1))                     # create vector of adjusted p-values
  
    Pcomb = apply(P.ord,1,min)                      # combine vectors of p-values with Tippett's comb. fct.
    p.ris[1] = p.glob = mean(Pcomb[-1]<=Pcomb[1])   # first adjusted p-value corresponds with the global p-value
  
    if(k>2){                                        # apply tippett step-down algorithm for p-value adjustement
        for(j in 2:(k-1)){
            T = apply(P.ord[,j:k],1,min) 
            p.ris[j] = max(mean(Pcomb[-1]<=Pcomb[1]),p.ris[(j-1)])
        }
    }
    
    p.ris[k] = max(P.ord[1,k],p.ris[k-1])             # last adjusted p-value
    p.ris[ord] = p.ris                                # put the ajusted p-values in the right order
    rownames(p.ris) = colnames(P)
    
  return(p.ris)
}