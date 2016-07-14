
########################
#  Rademacher Averages #
########################


#  Generate a random Rademacher vector
#  Input:
#    n, number of iid Rademacher's to generate
#  Output: 
#    n-long vector of iid Rademacher's
rrademacher <- function( n ){
  return( 2*rbinom( n, 1, 0.5)-1 );
}

# Same as above only geneates 1 and i
# instead of 1 and -1
rrademacherImg <-function( n ){
  return( rbinom( n, 1, 0.5)*(-1+1i)+1 );
}


fdtool_radMeanImg <- function( funct, cov1, iter, p )
{
  n   = nrow(funct);
  res = 0;
  for( j in 1:iter )
  {
    r   = rrademacherImg( n );
    temp= fdtool_centerData(funct)*r;
    rad = (t(temp)%*%temp)/n - (sum(r^2)*cov1)/n;
    res = res + pschnorm(Re(rad),p);
  }
  return( res/iter );
}

fdtool_radMeanListImg <- function( functList, covList, iter, p, tCov=0 )
{
  cnt = length( functList );
  n   = lapply( functList,nrow );
  res = 0;
  for( j in 1:iter )
  { 
    rad = 0;
    radp= 0;
    radm= 0;
    for( i in 1:cnt )
    {
      r    = rrademacherImg( n[[i]] );
      temp = fdtool_centerData(functList[[i]])*r;
      if(p==1/2){
        tmp2  = temp*(r==1);
        radp  = radp + (t(tmp2)%*%tmp2 - 
               sum(r^2)*covList[[i]])/n[[i]];
        tmp2  = temp*(r==1i)*(-1i);
        radm  = radm + (t(tmp2)%*%(tmp2) - 
               sum(r^2)*covList[[i]])/n[[i]];
      } else {
        rad  = rad + pschnorm( (t(temp)%*%temp - 
               sum(r^2)*tCov)/n[[i]],p);
      }
    }
    if(p==1/2){
      res  = res + pschnorm(sqrtMat(radp)-sqrtMat(radm),2);
    } else {
      res  = res + (rad);
    }
  }
  return( res/iter );
}

