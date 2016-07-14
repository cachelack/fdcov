# Liptak combining function
# 
# @param v Vector of p-values
# 
# @return This function returns a global p-value computed according to the Fisher combining function
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
# 
# @references Pesarin, Fortunato, and Luigi Salmaso (2010). Permutation tests for complex data: theory, applications and software. John Wiley & Sons.
# 
# @export
# 
comb.lipt = function(v){
    
    iter = 1000 # trovare un modo per metterlo nei parametri
    q = rep(0,length(v))
    for (i in 1:length(v)){
        p =  (v[i] + 1/(2*iter))/(1 + 1/iter)
        q[i] = qnorm(p)
    }
    
    return(-sum(q))
}