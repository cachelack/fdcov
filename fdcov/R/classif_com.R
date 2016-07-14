#################################
# Covariance Classification via #
# Concentration Inequalities    #
#################################

#' @name classifier-com
#' @rdname classify_com
#'  
#' @title Functional data classifier via concentration inequalities
#'  
#' @description
#' \code{classif.com} trains a covariance operator based 
#' functional data classifier that makes use of concentration inequalities.  
#' \code{predict.classif.com} uses the previously trained classifier
#' to classify new observations.
#' 
#' @details
#' These functions are used to train a functional data classifier
#' and to predict the labels for a new set of observations.
#' This method classifies based on the distances between each 
#' groups' sample covariance operator.  A simplified version of
#' Talagrand's concentration inequality is used to achieve this.
#' 
#' If the flag SOFT is set to TRUE, then soft classification occurs.
#' In this case, given k different labels, a k-long probability vector
#' is returned for each observation whose entries correspond to the 
#' probabilities that the observed function belongs to each specific label.
#' 
#' @param  dat (n X m) data matrix of n samples
#'         of m long vectors.
#' @param  datGrp A vector of group labels.
#' @param  object A concentration-of-measure classifier object
#'         of class inheriting from \code{classif.com}.
#' @param  SOFT Boolean flag, which if TRUE, returns soft
#'         classification for each observation.
#' @param  LOADING Boolean flag, which if TRUE, prints
#'         a loading bar.
#' @param  ... additional arguments affecting the predictions produced.
#' @return \code{classif.com} returns a functional data classifier
#'         object.  \code{predict.classif.com} returns a vector
#'         of n labels ( or an array of n probability vectors if  
#'         SOFT=TRUE )
#' @references
#'   Kashlak, Adam B, John AD Aston, and Richard Nickl (2016).
#'   "Inference on covariance operators via concentration
#'   inequalities: k-sample tests, classification, and clustering via
#'   Rademacher complexities", (in review)
#' @author Adam B Kashlak \email{ak852@cam.ac.uk}
#' @examples
#' \dontrun{
#   # Load phoneme data
#'  library(fds);
#'  # Setup training data
#'  dat1 = rbind( 
#'    t(aa$y[,1:100]), t(ao$y[,1:100]), t(dcl$y[,1:100]), 
#'    t(iy$y[,1:100]), t(sh$y[,1:100]) 
#'  );
#'  # Setup testing data
#'  dat2 = rbind( 
#'    t(aa$y[,101:400]), t(ao$y[,101:400]), t(dcl$y[,101:400]), 
#'    t(iy$y[,101:400]), t(sh$y[,101:400]) 
#'  );
#'  
#'  datgrp = gl(5,100);
#'  clCom = classif.com( datgrp, dat1 );
#'  grp = predict( clCom, dat2, LOADING=TRUE );
#'  acc = c(
#'    sum( grp[1:300]==1 ), sum( grp[301:600]==2 ), sum( grp[601:900]==3 ), 
#'    sum( grp[901:1200]==4 ), sum( grp[1201:1500]==5 )
#'  )/300;
#'  print(rbind(gl(5,1),signif(acc,3)));
#' }
#' \dontshow{
#   # Load phoneme data
#'  library(fds);
#'  # Setup training data
#'  dat1 = rbind( 
#'    t(aa$y[,1:5]), t(ao$y[,1:5]), t(dcl$y[,1:5]) 
#'  );
#'  # Setup testing data
#'  dat2 = rbind( 
#'    t(aa$y[,11:15]), t(ao$y[,11:15]), t(dcl$y[,11:15])
#'  );
#'  
#'  datgrp = gl(3,5);
#'  clCom = classif.com( datgrp, dat1 );
#'  grp = predict( clCom, dat2 );
#' }
NULL

#' @rdname classify_com
#' @export 
classif.com <- function( datGrp, dat )
{
  grps    = unique(datGrp);
  len     = length(grps);
  datList = list();
  for( i in 1:len ){
    datList[[i]] = subset(dat,datGrp==grps[i]);
  }
  res = fdcls_processLabels( datList, grps, "oper" );
  class(res) <- "classif.com";
  return(res);
}

#' @rdname classify_com
#' @export 
predict.classif.com <- function( 
  object, dat, SOFT=FALSE, LOADING=FALSE, ...
){
  if(!inherits(object, "classif.com")) 
    warning("calling predict.classif.com() on incorrect object!");
  clCOM = object;
  len = nrow(dat);
  res = c();
  if(LOADING)
    pb  = txtProgressBar(min = 0, max = len, style = 3);
  for( i in 1:len ){
    res = rbind(res,fdcls_calcProb( dat[i,], clCOM, "oper" ));
    if(LOADING)
      setTxtProgressBar(pb,i);
  }
  if(LOADING)
    close(pb);
  if(SOFT) 
    return(res);
  #res = setMax21(res);
  grp = apply(res,1,which.max);
  return( grp );
}

#
#  Given unlabled data, construct label probabilities from data
#  Input:
#    dat,     Functional Data
#    labList, list containing label information

fdcls_calcProb <- function( dat, labList, type, RETRAW=FALSE )
{
  n   = nrow(dat);
  if(is.null(n)) n=1;
  len = length(labList);
  phi = rep(0,len);
  # Compute unnormalized log label probabilities
  for( i in 1:len ){
    phi[i] = fdcls_computeLogPhi( n, dat, labList[[i]], type );
  }
  # Normalize phi
  mat = phi-t(matrix(phi,len,len));
  res = colSums(exp(mat));
  if(RETRAW){
    lk = 2*phi-t(matrix(phi,len,len));
    # For numerical stability
    lmn= apply(lk,2,max);
    lk = t(t(lk)-lmn);
    lk = exp(lk);
    lk = colSums(lk)*exp(lmn);
    lk = log((sum(subset(1/lk,is.finite(1/lk)))));
    if(!is.finite(lk)) lk=0;
    return(list(prb=1/res,like=lk));
  }
  return( 1/res );
}

#
# Compute log unnormalized label probability, 
# i.e. log P( trueLabel=lab | data, lab_info )
# Input:
#   n,   
#   dat, Data
#   lab, Label

fdcls_computeLogPhi <- function( n, dat, lab, type )
{
  if( n==1 ){
    tAve = outer(as.double(dat-lab$mn),as.double(dat-lab$mn))
  } else {
    tAve = t(t(dat)-lab$mn);
    tAve = t(tAve)%*%tAve;
    #tAve = fdcls_computeMean( t(t(dat)-lab$mn), type );
  }

  # r = || sum f_i - g || - R_n
  r = fdcls_computeNorm(tAve-lab$ave,type)-lab$rad;

  return(
    -n*r^2/( 2*( lab$wvr ) )
  );
}

#
# Create Label objects
#

#
#  Turn Labelled functional data in to 
#  cov oper, norm, and rad norm
#  Input:
#    datList, list of labeled functional data
#

fdcls_processLabels <- function( datList, grps, type )
{
  len = length(datList);
  res = list();
  for( i in 1:len ){
    # functional mean and sd
    mn  = apply( datList[[i]],2,mean );
    std = apply( datList[[i]],2,sd );
    datList[[i]] = t((t(datList[[i]])-mn));
    # Average of data in category i
    ave = fdcls_computeMean( datList[[i]], type );
    # Weak variance of data in category i
    wvr = fdcls_computeWVar( datList[[i]], ave, type );
    # Rademacher Average of data in category i
    rad = fdcls_computeRade( datList[[i]], ave, type );
    # Save results for category i
    res[[i]] = list( ave=ave,wvr=wvr,rad=rad,mn=mn,std=std,grp=grps[i] );
    #print(paste(wvr,rad))
  }
  return(res);
}

# Compute mean ( mean of operator is covariance )
fdcls_computeMean <- function( data, type ){
  switch( type,
    fnct = colSums(data)/nrow(data),
    oper = cov(data)
  );
}

# Compute weak variance
fdcls_computeWVar <- function( data, ave, type ){
  switch( type,
    fnct = pschnorm( cov(data) - ave%*%t(ave), 1 ),
    oper = 2*pschnorm( ave, 1 )  
  );
}

# Compute Rademacher average
fdcls_computeRade <- function( data, ave, type ){
  switch( type,
    fnct = sum(abs(colSums( 
      fdtool_centerData(data)*rrademacher(ncol(data)) 
    )/nrow(data)))/ncol(data),
    oper = fdtool_radMeanImg( data, ave, 100, 1 )
  );
}

# Use trace-class norm only
fdcls_computeNorm <- function( ave, type ){
  switch( type,
    fnct = sum(abs(ave))/length(ave),
    oper = pschnorm( ave, 1 )
  );
}

