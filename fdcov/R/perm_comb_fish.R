# Fisher combining function
# 
# @param v Vector of p-values
# 
# @return This function returns a global p-value computed according to the Fisher combining function
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
# 
# @references Pesarin, Fortunato, and Luigi Salmaso (2010). Permutation tests for complex data: theory, applications and software. John Wiley & Sons.
# 

comb.fish = function(v){
    
    return(-2*sum(log(v)))
    
}