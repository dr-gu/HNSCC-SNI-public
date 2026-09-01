   
                         
   
                
                                    
                                                                                    
                                                                                             
   
            
                                                          
                                                     
                                                          
                                                                 
                                                               
  
               
                         
   
            
                                
   
            
                              

doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
                
    
                   
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
                        
    yr <- (yr - mean(yr)) / sd(yr)
    
                                 
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
                      
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}
