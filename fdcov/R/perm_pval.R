# Permutation p-value
# 
# Given a vector of test statistics, where the first one is considered to be the one corresponding to the observed data set, computes the p-value of the test
# 
# @param t Vector of test statistics
# @param extr To set which values are more extreme than the one observed. Can be \code{greater} or \code{lesser}.
# 
# @return This function returns the p-value of the test, that is the percentages of test statistics in the vector that are more extreme than the one observed
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
# 
perm.pval = function(t, extr = 'greater'){
     if(extr == 'greater') p = mean(t >= t[1],na.rm=TRUE) 
     else                  p = mean(t <= t[1],na.rm=TRUE)
     return(p)
    }
