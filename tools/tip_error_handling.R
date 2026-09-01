                   
                                                                                              
                                                                                            
                                   
   
                                                                                          
                                                                               
                                                                                                     
                                                                                                   
                                                                                                 
                                                                                                                     
                            

DealWithError <- function(filePath, fileName, codePath, format.of.file, type.of.data, saveDir){
                                                                                           
  input.check <- tryCatch({
    expression.profile <- read.table(sprintf("%s/%s",filePath, fileName), sep = "\t", stringsAsFactors = FALSE, header =TRUE, 
                                     row.names = 1, check.names=F, na.strings = NULL)
                      
  }, error=function(e){
    print(e$message)
    
    return(NULL)    
  })
  
  
               
  if(is.null(input.check)){
               
    input.check2 <- tryCatch({
      expression.profile <- read.table(sprintf("%s/%s",filePath, fileName), sep = "\t", stringsAsFactors = FALSE, header =TRUE, 
                                       check.names=F, na.strings = NULL)
      if(sum(duplicated(expression.profile[,1])) > 0){errorString <- "DuplicateNameError"; return(errorString);}
    }, error=function(e){
      print(e$message)
      
      return(NULL)    
    })
    
    if(is.null(input.check2)){errorString <- "InputError"; return(errorString);}
    
    
  }else{
       SampleNumber.check <- ncol(expression.profile)
       
                                                                                                           
       if(SampleNumber.check == 0){errorString <- "SeparatorError"; return(errorString);}
                                                                                                           
                                
       gene.info <- get(load(sprintf("%s/2.ImmuneActivityScore/human_gene2ensembl2symbol_list.RData", codePath)));
       id <- as.character(rownames(expression.profile))
       names(id) <- id
       
       inter_id <- intersect(id, gene.info[[1]][,1])
       if(length(inter_id) > 0){
         id[inter_id] <- gene.info[[1]][inter_id, 2]
       }
       
       inter_id<-intersect(id, gene.info[[2]][,1])
       if(length(inter_id) > 0){
         id[inter_id] <- gene.info[[2]][inter_id, 2]
       }
       
       rownames(expression.profile) <- id
       print("Rownames to geneSymbol done!")
       
       load(sprintf("%s/2.ImmuneActivityScore/annotation.RData",codePath));
       overlap.num <- length(intersect(rownames(expression.profile), annotation[,3]))
       if(overlap.num < 20){errorString <- "GeneNameError"; return(errorString);}
                                                                                                           
      if(all(!is.na(as.numeric(colnames(expression.profile)[-1])))){
        errorString <- "HeaderNameError"; return(errorString);
                                                                                                           
      }else if(sum(complete.cases(expression.profile))!=nrow(expression.profile) | sum(is.na(expression.profile))!=0 | 
         length(which(expression.profile==""))!=0 | length(which(expression.profile==" "))!=0){
        errorString <- "DefaultValueError"; return(errorString);
      }
  }
  if(!dir.exists(saveDir)){ dir.create(saveDir) }
  save(expression.profile, file = paste0(saveDir, "/expression.afterIDConvert.RData"))
  return(NULL)
}



                
                                                  
                                                                                          
                                                                               
                                                                                                     
                                                                                                   
                                                                                                 
                                                                                                                             

checkError <- function(filePath, fileName, codePath, format.of.file, type.of.data, saveDir){
  result <- DealWithError(filePath = filePath, fileName = fileName, codePath = codePath, 
                          format.of.file = format.of.file, type.of.data = type.of.data, saveDir=saveDir)
  print(result)
  if(!is.null(result)){
    if(!dir.exists(saveDir)){ dir.create(saveDir) }
    write(result, file = sprintf("%s/ErrorString.txt",saveDir))
  }
}
