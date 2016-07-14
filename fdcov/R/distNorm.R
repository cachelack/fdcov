###########################
# Distance Functions for  #
# covariance operators    #
###########################

# General distances function
# 
# @param mat1 First covariance matrix
# @param mat2 Second covariance matrix
# @param dist Distance between covariance operators. Can be 'sq' (square-root), 'tr' (trace),'pr' (Procrustes), 'hs'(Hilbert-Schmidt) or 'op' (operator).
# @return Distance.
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# 
distCov <- function( mat1, mat2, type )
{
   switch( type,
     sq = distSqrt(mat1,mat2),
     tr = distTrac(mat1,mat2),
     pr = distProc(mat1,mat2),
     hs = distHsno(mat1,mat2),
     op = distOper(mat1,mat2)
   );
}

# Trace Class distance
# 
# @param mat1 First covariance matrix
# @param mat2 Second covariance matrix
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
distTrac <- function( mat1, mat2 )
{
  return( pschnorm( mat1-mat2,1 ) );
}

# Hilbert-Schmidt distance
# 
# @param mat1 First covariance matrix
# @param mat2 Second covariance matrix
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
distHsno <- function( mat1, mat2 )
{
  return( pschnorm( mat1-mat2, 2 ) );
}

# Operator norm distance
# 
# @param mat1 First covariance matrix
# @param mat2 Second covariance matrix
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
distOper <- function( mat1, mat2 )
{
  return( pschnorm( mat1-mat2, -1 ) );
}

# Square Root distance
# 
# @param mat1 First covariance matrix
# @param mat2 Second covariance matrix
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
distSqrt <- function( mat1, mat2 )
{
  smat1 = sqrtMat(mat1);
  smat2 = sqrtMat(mat2);
  return( pschnorm( smat1-smat2, 2 ) );
}

# Procrustes distance
# 
# @param mat1 First covariance matrix
# @param mat2 Second covariance matrix
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
distProc <- function( mat1, mat2 )
{
  smat1 = sqrtMat(mat1);
  smat2 = sqrtMat(mat2);
  matC  = t(smat2)%*%smat1;
  svdC  = svd(matC);
  matR  = svdC$u%*%t(svdC$v);
  return( pschnorm( smat1-smat2%*%matR, 2 ) );
}

###########
#  Norms  #
###########

# Hilbert-Schmidt (Frobenius) Norm
# 
# @param sig covariance matrix (i.e. symmetric positive definite)
# 
# @param HS Norm of sig
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
hsnorm <- function( sig )
{
  return(sqrt(sum(abs(sig)^2)));
}

# p-Schatten Norm
# 
# @param sig covariance matrix (i.e. symmetric positive definite)
# @param p   [1,Inf] or 1/2
# 
# @return p-Schatten Norm of sig
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
pschnorm <- function( sig, p )
{
  if( p==2 )
    return( hsnorm(sig) );
  if(p==1/2)
    return( sqrt(pschnorm(sig,1)) );
  eigval = eigen( sig, symmetric=TRUE, only.values=TRUE );
  if( p==-1||is.infinite(p) )
    return( max( abs(eigval$values) ) );
  return( sum(abs(eigval$values)^p)^(1/p) );
}

# Computes Square Root of matrix A
# 
# @param A matrix
# 
# @return Square root of A
# 
# @author Adam B Kashlak \email{ak852@cam.ac.uk}
# 
# @export
# 
sqrtMat <-function(A)
{
  eig = eigen( A );
  val = eig$values;
  #val = pmax( eig$values, rep(0,nrow(A)) );
  d   = sqrt( as.complex(val) );
  D   = diag( d );
  V   = eig$vectors;
  return( V%*%D%*%t(V) );
}