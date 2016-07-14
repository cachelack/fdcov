# P-value significance
# 
# @param p P-value
# 
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
# 


perm.sig<-function(p){

  ### Significance function
    
  s<-character(length=dim(p)[1])
  
  for(i in 1:dim(p)[1]){
    if(p[i]<0.1){
      s[i]="."
      if(p[i]<0.05){
        s[i]="*"
        if(p[i]<0.01){
          s[i]="**"
          if(p[i]<0.001){
            s[i]="***"
            
          }
        }
      }
    }
  }
  return(s)
}# end sig.