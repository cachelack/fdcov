# Test statistics to p-values
# 
# @param T Matrix of test statistics
# 
# @return Matrix containing p-values
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
# 

perm.t2p<-function(T){

  #if(is.null(dim(T))){T<-array(T,dim=c(length(T),1))}
  oth<-2#seq(1:length(dim(T)))[-1]

  B<-dim(T)[1]-1							
  p<-dim(T)[2]	
  #if(length(dim(T))==3){C<-dim(T)[3]}

  rango<-function(x){
    r=1-rank(x[-1],ties.method="min")/B+1/B
    return(c(mean(x[-1]>=x[1]),r))
  }

  P=apply(T,oth,rango)
  return(P)
}