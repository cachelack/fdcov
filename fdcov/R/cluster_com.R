

#################################
# Covariance Clustering via     #
# Concentration Inequalities    #
#################################

#' @title Functional data clustering via concentration inequalities
#'  
#' @description
#' \code{cluster.com} clusters sets of functional data via their 
#' covariance operators making use of an EM style algorithm with
#' concentration inequalities.  
#' 
#' @details
#' This function clusters individual curves or sets of curves
#' by considering the distance between their covariance operator
#' and each estimated category covariance operator.  The implemented
#' algorithm reworks the concentration inequality based 
#' classification method \code{classif.com} into an EM style algorithm.
#' This method iteratively updates the probability of a given observation
#' belonging to each of the k categories.  These probabilities are 
#' in turn used to update the category means.  This process continues
#' until either the total number of iterations is reached or a
#' computed likelihood begins to decrease signaling the arrival of 
#' a local optimum. 
#' 
#' If the argument labl is NULL, then every curve is clustered
#' separately.  If labl contains factors used to group the curves,
#' then each set of curves is classified as one group.  For example,
#' if you have multiple speakers and multiple speech samples from 
#' each speaker, you can group the data from each speaker together
#' in order to cluster based on each speakers' covariance operator
#' rather than based on each speech sample individually.
#' 
#' If the flag SOFT is set to TRUE, then soft clustering occurs.
#' In this case, given k different labels, a k-long probability vector
#' is returned for each observation whose entries correspond to the 
#' probability that the observed function belongs to a specific label.
#' 
#' @param  dat (n X m) data matrix of n samples
#'         of m long vectors.
#' @param  labl An optional vector of n labels to group curves. (see Details)
#' @param  grpCnt Number of clusters into which to split the data.
#' @param  iter Number of iterations for EM algorithm.
#' @param  SOFT Boolean flag for whether or not category probabilities
#'         should be returned.
#' @param  PRINTLK Boolean flag, which if TRUE, prints 
#'         likelihood values for each iteration.
#' @param  LOADING Boolean flag, which if TRUE, prints
#'         a loading bar.
#' @param  IGNORESTOP Boolean flag, which if TRUE, will ignore 
#'         early stopping conditions and cause the EM algorithm to
#'         run for the total amount of desired iterations.
#' @return \code{cluster.com} returns a vector a labels with one entry for 
#'         each row of data corresponding to one of the k categories
#'         ( or an array of probability vectors if SOFT=TRUE ).
#' @references
#'   Kashlak, Adam B, John A D Aston, and Richard Nickl (2016).
#'   "Inference on covariance operators via concentration
#'   inequalities: k-sample tests, classification, and clustering via
#'   Rademacher complexities", in review
#' @author Adam B Kashlak \email{ak852@cam.ac.uk}
#' @examples
#' \dontrun{
#'  # Load phoneme data 
#'  library(fds);
#'  # Setup data to be clustered
#'  dat  = rbind( t(aa$y[,1:20]),t(iy$y[,1:20]),t(sh$y[,1:20]) );
#'  # Cluster data into three groups
#'  clst = cluster.com(dat,grpCnt=3);
#'  matrix(clst,3,20,byrow=TRUE);
#'  
#'  # cluster groups of curves
#'  dat  = rbind( t(aa$y[,1:40]),t(iy$y[,1:40]),t(sh$y[,1:40]) );
#'  lab  = gl(30,4);
#'  # Cluster data into three groups
#'  clst = cluster.com(dat,labl=lab,grpCnt=3);
#'  matrix(clst,3,10,byrow=TRUE);
#' }
#' \dontshow{
#'  # Load phoneme data 
#'  library(fds);
#'  # Setup data to be clustered
#'  dat  = rbind( t(aa$y[,1:5]),t(iy$y[,1:5]) )
#'  clst = cluster.com(dat);
#' }
#' @export 

cluster.com <- function( 
  dat, labl=NULL, grpCnt=2, iter=30, SOFT=FALSE, 
  PRINTLK=TRUE, LOADING=FALSE, IGNORESTOP=FALSE
){
  # Group observations if desired
  if(is.null(labl))
    labl = gl(nrow(dat),1)
  ulab = unique(labl);
  n    = length(ulab);
  # generate initial probabilities
  prb = rdirichlet( n, rep(1/2,grpCnt) );
  bst = prb; # best probability vectors
  # calculate inital proportions
  prp = colSums(prb)/n;
  ## calculate inital cov ops and norms
  sig = fdem_processLabels( dat,grpCnt,labl,prb );
  
  # Calculate Pooled Covariance Operator and 
  # hence, the pooled Weak Variance
  pool= matrix(0,ncol(dat),ncol(dat));
  for( ii in 1:grpCnt ){
    pool = pool + sig[[ii]]$ave*prp[ii];
  }
  pool = sqrt(2)*pschnorm(pool,1);
  # Set Weak Variances to the Pooled WVar
  for( ii in 1:grpCnt ){ 
    sig[[ii]]$wvar = pool;
  }

  # Vector of Likelihoods
  lk = rep(0,iter);
  # Boolean Vector?
  fxd= rep(FALSE,n);
  # iterate
  if(LOADING)
    pb  = txtProgressBar(min = 1, max = iter, style = 3);
  for( it in 2:iter ){
    indx=1;
    for( lab in ulab ){  # For each curve or group
      # Calc prob of the group
      cpOut   = fdcls_calcProb( 
        subset(dat,labl==lab), sig, "oper",RETRAW=TRUE 
      );
      # Update Probability
      prb[indx,] = cpOut$prb;
      indx       = indx+1;
      # Update log likelihood
      lk[it]  = lk[it] + cpOut$like;
    }
    # If lk begins to increase, we found a local optimum. Stop!
    if( !IGNORESTOP && it>2 && lk[it]>lk[it-1] )
      break();
    # Update choice of bst
    bst = prb;
    # If lk is zero, no more changes will occur. Stop!
    if( !IGNORESTOP && lk[it]==0 )
      break();
    if(PRINTLK)
      print(lk[it]);
    # Update category proportions
    prp = rbind(prp, colSums(prb)/n);
    # Process labels again
    sig = fdem_processLabels( dat,grpCnt,labl,prb );
    # Computed pooled Weak Variance
    pool= matrix(0,ncol(dat),ncol(dat));
    for( ii in 1:grpCnt ){
      pool = pool + sig[[ii]]$ave*prp[it,ii];
    }
    pool = sqrt(2)*pschnorm(pool,1);#*(decay^(it-1));
    ## Print Results
    #print(paste(lk[it],pool,radLine));
    # Set all labels to pooled WVar
    for( ii in 1:grpCnt ){ 
      sig[[ii]]$wvar = pool;
    }
    if(LOADING)
      setTxtProgressBar(pb,it);
  }
  if(LOADING)
    close(pb);
  
  if(SOFT){
    return(bst);
  }
  grp = apply(bst,1,which.max);
  return( grp );
}

#
#  Turn weighted functional data into
#  cov oper, norm, and rad norm
#

fdem_processLabels <- function( dat,grpCnt,labl,prb )
{
  res = list();
  for( i in 1:grpCnt ){
    # functional mean and sd
    mn  = colSums( dat*prb[labl,i] )/sum(prb[labl,i]);
    std = nrow(dat)*apply( dat*prb[labl,i],2,sd )/sum(prb[labl,i]);
    dati= t((t(dat)-mn));
    # Average of data in category i
    ave = fdem_calcCov( dati,grpCnt,labl,prb[,i] );
    # Weak variance of data in category i
    wvr = 2*pschnorm( ave, 1 );
    # Rademacher Average of data in category i
    rad = fdem_radMeanImg( dati, ave, grpCnt, labl, prb[,i], 10, 1 );
    # Save results for category i
    res[[i]] = list( ave=ave,wvr=wvr,rad=rad,mn=mn,std=std );
    #print(paste(wvr,rad))
  }
  return(res);
}


#
#  Random Dirichlet Draw
#

rdirichlet <- function(n, a)
{
  len = length(a);
  res = rgamma( len*n, a, 1 );
  res = t(matrix(res,len,n));
  res = res/rowSums(res);
  return(res);
}

#
# Calc weighted Rad average
#

fdem_radMeanImg <- function( dat, cov1, grpCnt, labl, prb, iter, p )
{
  prb = prb[labl];
  n   = length(unique(labl));
  res = 0;
  for( j in 1:iter )
  {
    r   = rrademacherImg( n );
    r   = r[labl];
    temp= fdem_wcenter(dat,grpCnt,prb)*sqrt(prb)*r;
    rad = (t(temp)%*%temp - sum(r^2*prb)*cov1)/sum(prb);
    res = res + pschnorm(Re(rad),p);
  }
  return( res/iter );
}

#
# Calc weighted covariance
#

fdem_calcCov <- function(dat,grpCnt,labl,prb)
{
  prb = prb[labl];
  dat = fdem_wcenter(dat,grpCnt,prb)*sqrt(prb);
  res = (t(dat)%*%dat)/sum(prb);
  return(res);
}

#
#  Center data with weights
#

fdem_wcenter <- function( dat, mm, prb )
{
  return( (t(t(dat)-colSums(dat*prb)/sum(prb))) );
}

