                 
                                          
                                                                                                     
                                                                                          
                                                                               
                                                                    
                                                                              
                                                                                   
                                                                                                 
                                                                                                    
                                                                                    

TIP_process <- function(
  codePath,
  filePath,
  fileName,
  saveDir,
  sampleNumber,
  permTimes,
  type.of.data,
  format.of.file,
  sample)
  {
  print("start to source function!")
  print(paste(codePath,"3.ImmuneInfiltration/1.CIBERSORT.server.R",sep="/"))
  print(paste(codePath,"2.ImmuneActivityScore/immunityScore.server.R",sep="/"))
  source(paste(codePath,"3.ImmuneInfiltration/1.CIBERSORT.server.R",sep="/"))
  source(paste(codePath,"2.ImmuneActivityScore/immunityScore.server.R",sep="/"))
  
  print("start to run immunityScore.server!")
  immunityScore.server(codePath =paste(codePath,"2.ImmuneActivityScore",sep="/"), 
                       filePath = filePath,
                       fileName = fileName,
                       saveDir = filePath,
                       sampleNumber = sampleNumber,
                       permTimes = permTimes,
                       type.of.data = type.of.data,
                       format.of.file = format.of.file);
  
  print("start to run CIBERSORT.server!")
  CIBERSORT.server(codePath =paste(codePath,"3.ImmuneInfiltration",sep="/"),
                   filePath = filePath,
                   fileName = fileName,
                   signaturePath = paste(codePath,"3.ImmuneInfiltration",sep="/"),
                   saveDir = filePath,
                   perm = permTimes,
                   CHIPorRNASEQ = type.of.data,
                   sample= sample,
                   dataType = format.of.file)
  
  print("the TIP_integration is end !")
}
