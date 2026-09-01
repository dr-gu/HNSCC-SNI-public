                          
   
                                                                                                     
                                                       
   
                                                                                                    
                                                                                         
                                                                              
                                                                   
                                                                             
                                                                                  
                                                                                                 
                                                                                                  
                                          
               
               
   
                   

immunityScore.server <- function(codePath, filePath, fileName, saveDir, sampleNumber, permTimes, type.of.data, format.of.file){
	
	
	
                                            
	load(sprintf("%s/signature.GeneSymbol.list.RData",codePath));
	load(sprintf("%s/annotation.RData",codePath));
	load(sprintf("%s/genename_length3.0.Rdata",codePath));                                                          
	source(sprintf("%s/Count2TPM3.0.R",codePath));
  
	source(sprintf("%s/ProcessMultipleSample.R",codePath));                                                                        
	source(sprintf("%s/ProcessSingleSample.R",codePath));                                                                     
	
	source(sprintf("%s/ssgsea.core.R",codePath)); 
	source(sprintf("%s/ssGSEAPermutation.R",codePath));
	
  source(sprintf("%s/makeHeatmapData.R",codePath));                                                
	source(sprintf("%s/PCAScatterPlot.R",codePath)); 
	source(sprintf("%s/score_boxplot.R",codePath));
	source(sprintf("%s/changeStepNames.R",codePath));
	
  source(sprintf("%s/filterOutliers.R",codePath)); 
  
	print("source is end");
	print("set.seed");
	set.seed(1:100);
	
	library(pheatmap);
	
  if(sampleNumber > 1){
	  example <- get(load(sprintf("%s/%s",saveDir, "expression.afterIDConvert.RData")))
    print("read down!")
	  process.check <- tryCatch({
	    example.result <- ProcessMultipleSample(expression.from.users = example, 
	                                            save.dir = saveDir, 
	                                            signatureList = signature.GeneSymbol.list, 
	                                            perm.times = permTimes, 
	                                            signature.annotation = annotation,
	                                            type.of.data=type.of.data,
	                                            format.of.file=format.of.file,
	                                            gene.length.path = codePath)
	  }, error=function(e) {
	    errorString <- "Other"
	    cat("Processing wrong...\n")
	    write(errorString, file = sprintf("%s/ErrorString.txt",saveDir))
	                                                               
	    write(e$message, file = sprintf("%s/ErrorInssGSEA.txt", saveDir))   
	  })
	  						
	}else{
	  example <- get(load(sprintf("%s/%s",saveDir, "expression.afterIDConvert.RData")))
	  print("read down!")
	  process.check <- tryCatch({
	    example.result <- ProcessSingleSample(expression.from.users = example, 
	                                          save.dir = saveDir, 
	                                          signatureList = signature.GeneSymbol.list, 
	                                          perm.times = permTimes, 
	                                          signature.annotation = annotation,
	                                          type.of.data=type.of.data,
	                                          format.of.file=format.of.file,
	                                          gene.length.path=codePath)
	   
	  }, error=function(e) {
	    errorString <- "Other"
	    cat("ssGSEA processing wrong...\n")
	    write(errorString, file = sprintf("%s/ErrorString.txt",saveDir))
	    write(e$message, file = sprintf("%s/ErrorInssGSEA.txt", saveDir))
	       
	  })
	}
  print("immunityScore.server is end");
  write("immunityScore.server", file=sprintf("%s/processResult.txt",saveDir),append=T);
  write(format(Sys.time(), "%Y.%m.%d.%H.%M.%S"), file=sprintf("%s/processResult.txt",saveDir),append=T);
    
}