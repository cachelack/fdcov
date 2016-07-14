#' Plot partial p-values
#'
#' \code{perm.plot} plots all of the partial comparison p-values in a matrix.
#'
#' @param p Output of function perm.test, if part = TRUE.
#' @param k Number of groups, must be greater than 2.
#' @param lab Group labels. Defaults to 1, 2, ..., k.
#' @param save Boolean variable that indicates if the plot must be saved as an .eps. Defaults to FALSE.
#' @param name If \code{save} is TRUE, this is the filename of the plot. Defaults to \code{pvalues.eps}.
#'
#' @return \code{perm.plot} plots the partial p-values in a matrix.
#' 
#' @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
#' 
#' @references Pigoli, Davide, John A. D. Aston, Ian L. Dryden, and Piercesare Secchi (2014). "Distances and inference for covariance operators." Biometrika: 101(2):409–422.
#'
#' @examples 
#' \dontrun{
#' ## Phoneme data
#' 
#' library(fdcov)
#' library(fds)
#' 
#' # Create data set
#' data(aa); data(ao); data(dcl);data(iy);data(sh)
#' dat=cbind(aa$y[,1:20],ao$y[,1:20],dcl$y[,1:20],iy$y[,1:20],sh$y[,1:20])
#' dat=t(dat)
#' grp=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))
#' 
#' # Test the equality of the covariance operators
#' p=ksample.perm(dat,grp,iter=100,only.glob=FALSE)
#' 
#' # Plot partial p-values
#' perm.plot(p,5, lab=c('aa','ao','dcl','iy','sh'))
#' }
#' 
#' @export

perm.plot = function(p, k, lab = NULL, save = FALSE, name = 'pvalues.eps'){
  
    if(k<3) stop('The number of groups must be at least 3.')
    
    pmatrix = matrix(1,k,k)
    cont = 1
    for(i in 1:(k-1)){
        for(j in (i+1):k){
            pmatrix[i,j] = p$partial$p_value[cont]
            cont = cont+1
        }
    }
    
    if(!is.null(lab)){
        colnames(pmatrix) = rownames(pmatrix) = lab
    }
    col = c(0,0,0,0,'#b30000','#e34a33','#fc8d59','#fdcc8a','#fef0d9')

    corrplot::corrplot(t(pmatrix),method = "color",type = "lower",tl.col='black',addCoef.col = "black",
             is.corr = FALSE, cl.lim = c(0,1), col = col, tl.pos = 'ld')
    if(save){
        setEPS()
        postscript(name)
        corrplot::corrplot(t(pmatrix),method = "color",type = "lower",tl.col = 'black',addCoef.col = "black",
           is.corr = FALSE, cl.lim = c(0,1), col=col, tl.pos = 'ld')
        dev.off()
    }

}
