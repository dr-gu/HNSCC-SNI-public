         
                                                        
    
                       
     
               
                                                                                       
                                                                                            
                                                         
   
                                                   
                
                                                                                            
                                  
                                     
                                           
                                                                                                     
                                                         
                                       
                                                                                                
  
            
                                                                                                   
                                                                                                    
                                         
   
  
            
                  
                                     
                         
                            
                                  
         
                                      
                     
                          
                    
                    
                                   
                            
   
   
                   
                                                                    
   
                    
                                                                                                              
   
                    
                                                                                                  
                                                                 
   
                         
                                              
   
                   
                                     
   
               
                                       
   
                       
                                                                               
   
                  
                                                                                                        
   
                    
                                                                                       
   

CIBERSORT.server <- function(codePath, filePath, fileName= "",signaturePath, saveDir, perm = 100,CHIPorRNASEQ="RNA-seq", sample="multiple",dataType="TPM"){

                                               
	
  source(sprintf("%s/2.CIBERSORT_main.R",codePath));
  source(sprintf("%s/3.CIBERSORT_func.R",codePath));
  source(sprintf("%s/4.CoreAlg.R",codePath));
  source(sprintf("%s/5.doPerm.R",codePath));
                                                   
  
  print("The Infiltration source codes is import!");

	
	
                                                                                                       
	               
	
  mixture_file = get(load(paste0(filePath,"/expression.from.users.tpm.RData")))
  print("The user data is import!");
	
                                                                                                                   
	
  if(CHIPorRNASEQ=="RNA-seq"){
    
    signature_name<-"LM14_name3.0.txt"
    
  }else if(CHIPorRNASEQ=="Microarray"){
    
    signature_name<-"LM22_name.txt"
    
    }
  
                               
  sig_matrix <- read.table(sprintf("%s/%s",signaturePath, signature_name), sep = "\t", stringsAsFactors = FALSE,header=T,row.names=1,check.names=F)
  
  print("The signature matrix is import!");
  
                                             
  Result <- CIBERSORT_server(mixture_file= mixture_file,
                             sig_matrix = sig_matrix ,
                             saveDir=saveDir,
                             perm=perm, 
                              CHIPorRNASEQ=CHIPorRNASEQ,
                             sample=sample)
  
  print("CIBERSORT.server is end");
	
  write("CIBERSORT.server", file=sprintf("%s/processResult.txt",saveDir),append=T);
}
