
#' k-sample test for equality of covariance operators
#'
#' \code{ksample.com} performs a k-sample test for equality
#' of covariance operators using concentration
#' inequalities.
#'
#' This function tests for the equality of k covariance 
#' operators given k sets of functional data.  It makes 
#' use of Talagrand's concentration inequality in the Banach
#' space setting. The argument p specifies the p-Schatten 
#' norm used in the test. As detailed in Kashlak et al (2016),
#' the most power is achieved using the trace class norm (p=1),
#' which is the default value.
#'
#' This test is inherently conservative as it constructed by
#' concatenating many concentration inequalities together.
#' Consequently, the method may be tuned by adjusting the 
#' arguments scl1 and scl2 to achieve the desired empirical
#' size for the users specific data set.  Otherwise, it can
#' be used as a quick first pass before a more powerful but
#' more computational test, such as specifically \code{ksample.perm},
#' is run.  More information on tuning this method can be 
#' found in the reference.
#' 
#' @param  dat (n X m) data matrix of n samples
#'         of m long vectors.
#' @param  grp n long vector of group labels.
#' @param  p p-Schatten norm in [1,Inf], Default is 1. (see Details)
#' @param  alpha the desired size of the test, Default is 0.05.
#' @param  scl1 scales the deviation part of the
#'         concentration inequality. (see Details)
#' @param  scl2 scales the Rademacher part of the
#'         concentration inequality. (see Details)
#' @return Boolean value for whether or not the test
#'         believes the alternative hypothesis is true.
#'         ( i.e. Does there exist at least two categories of the k 
#'         whose covariance operators are not equal? )
#' @references
#'   Kashlak, Adam B, John AD Aston, and Richard Nickl (2016).
#'   "Inference on covariance operators via concentration
#'   inequalities: k-sample tests, classification, and clustering via
#'   Rademacher complexities", (in review)
#' @author Adam B Kashlak \email{ak852@cam.ac.uk}
#' @examples
#' # Load in phoneme data
#' library(fds)
#' # Setup data arrays
#' dat1 = rbind( t(aa$y)[1:20,], t(sh$y)[1:20,] );
#' dat2 = rbind( t(aa$y)[1:20,], t(ao$y)[1:20,] );
#' dat3 = rbind( dat1, t(ao$y)[1:20,] );
#' # Setup group labels
#' grp1 = gl(2,20);
#' grp2 = gl(2,20);
#' grp3 = gl(3,20);
#' # Compare two disimilar phonemes (should return TRUE)
#' ksample.com(dat1,grp1);
#' # Compare two similar phonemes (should return FALSE)
#' ksample.com(dat2,grp2);
#' # Compare three phonemes (should return TRUE)
#' ksample.com(dat3,grp3);
#' @export

ksample.com <- function( 
  dat, grp, p=1, alpha=0.05, scl1=1,scl2=1  
){
  # Check p
  if( p<1 ){
    print("Error: p must be in [1,Inf]");
    return();
  }
  # Distance from category covs to pooled cov
  lhs = 0;
  # pooled weak variance
  pvar= 0;
  # pooled sample covariance
  pCov= cov(dat);
  # total sample size
  n   = nrow(dat);
  # group sample sizes
  ngrp= c();
  # create lists for Rademacher computation
  indx    = 1;
  datList = list();
  covList = list();
  # iterate over groups
  for( i in unique(grp) ){
    sDat= subset(dat,grp==i);
    sCov= cov(sDat);
    ngrp= c(ngrp,nrow(sDat));
    lhs = lhs + pschnorm( sCov-pCov, p );
    pvar= pvar+ pschnorm( sCov, p )^2*nrow(sDat);
    datList[[indx]] = sDat;
    covList[[indx]] = sCov;
    indx = indx+1;
  }
  # Compute pooled weak variance
  pvar = sqrt(2*pvar/n)
  # Compute Rademacher average
  # TODO: create lists above
  rad  = fdtool_radMeanListImg( datList, covList, 1, p, pCov )*scl2;
  # Compute sum 1/n_i
  invS = sum( 1/ngrp );
  # Compute rhs of inequality
  rhs  = rad + scl1*( 
     pvar*sqrt(-2*log(2*alpha)*invS) +
     pvar*log(2*alpha)*invS/3
  );
  return( lhs>rhs );
}


