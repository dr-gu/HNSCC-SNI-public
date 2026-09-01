                     
                                                                             
                                                                                                     
                                                                                          
                                                                               
                                                                    
                                                                              
                                                                                   
                                                                                                 
                                                                                                    
                                                                                     
                                                              
                                                                             
                                                   

TIP_integration <- function(
codePath,
filePath,
fileName,
saveDir,
sampleNumber,
permTimes,
type.of.data,
format.of.file,
sample,
CancerType, 
Samples, 
email){
  print("The filePath...")
  print(filePath)
  sink(file = sprintf("%s/%s",saveDir, "allPrint.txt"), append = TRUE)
  
  print("The filePath...")
  print(filePath)
  print("The saveDir...")
  print(saveDir)
  
               
  print("start to run checkError!")
  process.check1 <- tryCatch({
    source(paste(codePath,"1.MainFunction/ErrorProcess.R",sep="/"))
    checkError(filePath = filePath, fileName = fileName, format.of.file = format.of.file, 
               type.of.data = type.of.data, codePath = codePath, saveDir = saveDir)
   
  }, error=function(e) {
    cat("checkError processing wrong...\n")
    print(e$message)
    write(e$message, file = sprintf("%s/ErrorString.txt",saveDir), append = TRUE)
    return("error")  
  })
                   
  if(!is.null(process.check1)){
    sink();
    stop();
  }else{
    cat("checkError processing over...\n")
    }
  
  
                         
  process.check2 <- tryCatch({
    source(paste(codePath,"1.MainFunction/processResult.server.R",sep="/"))
    processResult.server(filePath, CancerType, type.of.data, format.of.file, sample,sampleNumber,Samples, email)
    
  }, error=function(e) {
    cat("processResult.server processing wrong...\n")
    print(e$message)
    write(e$message, file = sprintf("%s/ErrorString.txt", saveDir), append = TRUE)
    return("error")  
  })
  
  if(process.check2 == "error"){
    sink();
    stop();
  }else{
    cat("processResult.server over...\n")
  }
  
                
  process.check3 <- tryCatch({
    source(paste(codePath,"1.MainFunction/TIP_process.R",sep="/"))
    TIP_process(codePath, filePath, fileName, saveDir, sampleNumber, permTimes,
                type.of.data, format.of.file, sample)
    
  }, error=function(e) {
    cat("TIP_process processing wrong...\n")
    print(e$message)
    write(e$message, file = sprintf("%s/ErrorString.txt",saveDir),append = TRUE)
    return("error")  
  })
  if(process.check3 == "error"){
    sink();
    stop();
  }else{
    cat("TIP process over...\n")
  }
  
  sink()
}
  