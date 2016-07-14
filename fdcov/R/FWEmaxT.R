# Step-down Tippett procedure for strong FWE control
# 
# @param T (iter+1) X k matrix of permutation test statistics
# 
# @return The global p-value and a vector of adjusted p-values
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it} 

FWE.maxT = function(T){

  ord = order(T[1,],decreasing=TRUE)        # put the vector of observed test staiistics in decreasing order and store the order
  T.ord = T[,ord]                           # put the columns of matrix T in the new order 'ord'
  
  k = dim(T)[2]                             # number of tests
  p.ris = array(5,dim=c(k,1))               # create vector of adjusted p-values
  
  # Compute smallest p-value
  Tcomb = apply(T.ord,1,max)                  # combine vectors of p-values with max comb. fct.
  p.ris[1] = p.glob=mean(Tcomb[-1] >= Tcomb[1]) # the first adjusted p-value corresponds with the global p-value
  
  # Compute the other p-values
  if(k>2){                                  # apply general step-down algorithm for p-value adjustement
    for(j in 2:(k-1)){
        Tcomb = apply(T.ord[,j:k],1,max) 
        p.ris[j] = max(mean(Tcomb[-1] >= Tcomb[1]),p.ris[(j-1)])
      }
  }
  
  # Compute greatest p-value
  Tcomb = T.ord[,k] 
  p.ris[k] = max(mean(Tcomb[-1] >= Tcomb[1]),p.ris[k-1])    # last adjusted p-value
  
  # Put the ajusted p-values in the correct order
  p.ris[ord] = p.ris                                     
  
  rownames(p.ris) = colnames(T)
  return(p.ris)
}