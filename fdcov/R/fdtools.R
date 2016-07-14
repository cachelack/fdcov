
#
# set max entry in each row to 1, rest to 0
#

setMax21 <- function( arr )
{
  bst= apply(arr,1,max);
  arr= arr/bst;
  return(floor(arr));
}

# Subtract off mean vector
fdtool_centerData <- function(funct)
{
  n      = nrow(funct);
  mu     = colSums(funct);
  return( t(t(funct) - mu/n) );
}

# Sum a covset of n matrices
# covset : array (n X p X p) 
# return : (p X p) matrix 
fdtool_sumArray <- function( covset )
{
  return(
     apply( covset, c(2,3), sum )
  );
}
