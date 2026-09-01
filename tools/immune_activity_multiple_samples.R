                           
                                                                                              
                                                                                           
                                                                                                                 
                                                              
   
                                                                                 
                                                                     
                                                                                        
                                                                                     
                                                                                                                      
                           
                                                                                                  
                                                                                                    
                                                                                           

ProcessMultipleSample <- function(expression.from.users, save.dir, 
                                  signatureList = signature.GeneSymbol.list, 
                                  perm.times = 100, 
                                  signature.annotation = annotation,
                                  type.of.data = c("Microarray", "RNA-seq"), 
                                  format.of.file = c("TPM","Count"),
                                  gene.length.path){
  
  if(!dir.exists(save.dir)){ dir.create(save.dir) }
            
                                                                                                                
            
                                             
  separate.remain <- sapply(signatureList, function(i){
    num <- length(intersect(i, rownames(expression.from.users)))
    return(num)
  })
  
  show.num <- paste0("(", separate.remain, "/", sapply(signatureList, length), ")")
  names(show.num) <- names(signatureList)
  
                                                  
  positive.l <- grep("positive", names(signatureList))
                                    
  mix.show.num <- data.frame(separate.remain = separate.remain, signature.num = sapply(signatureList, length),
                             step = c(1,2,2,3,3,4:20,21,21,22,22,23,23))
  
  
  l0 <- which(separate.remain == 0)
                                                                                           
  a <- which(l0 %in% positive.l) 
  b <- which(l0 %in% (positive.l-1)) 
  if(length(a > 0) | length(b > 0)){
    delete.signature <- c(l0[-c(a,b)],l0[a],l0[a]-1,l0[b],l0[b]+1)
  }else{
    delete.signature <- l0
  }
  
  if(length(delete.signature) > 0){
    signatureList <- signatureList[-delete.signature]
    
    mix.show.num <- mix.show.num[-delete.signature,]
    
    get.separate <- unlist(strsplit(names(delete.signature),"\\."))
    delete.step <- intersect(get.separate[grep("Step", get.separate)], c("Step2","Step3","Step5","Step6","Step7"))
    delete.step <- gsub("Step","",delete.step)
                                     
    if(length(delete.step) > 0){
      signature.annotation <- signature.annotation[-which(signature.annotation[,4] %in% delete.step), ]
    }
    
  } 
  print(delete.signature)
  print(show.num)
  
  mix.show.num$step <- as.factor(mix.show.num$step)
  mix.show.num <- paste0("(", tapply(mix.show.num[,1], mix.show.num$step, sum), "/", 
                         tapply(mix.show.num[,2], mix.show.num$step, sum), ")")
  positive.l <- grep("positive", names(signatureList))
  if(length(positive.l) > 0){
    names(mix.show.num) <- names(signatureList)[-positive.l]
    names(mix.show.num) <- gsub(".negative", "", names(mix.show.num))
  }else{
    names(mix.show.num) <- names(signatureList)
  }
  
  mix.show.num <- paste(names(mix.show.num), mix.show.num, sep = " ")
  
  
  
            
                                                                                                               
            
                         
                                                                                                
  if(type.of.data == "RNA-seq" & format.of.file == "TPM"){
    save(expression.from.users, file = paste0(save.dir, "/expression.from.users.tpm.RData"))
    normalized.expression.from.users <- as.matrix(log2(expression.from.users+1)+1)
    
  }else if(type.of.data == "RNA-seq" &  format.of.file == "Count"){
                           
    normalized.expression.from.users <- count2Tpm(countName = expression.from.users, filePath = gene.length.path, sample="multiple")
    print(paste("count to TPM down---dim(normalized.expression.from.users)", dim(normalized.expression.from.users)))
    save(normalized.expression.from.users, file = paste0(save.dir, "/expression.from.users.tpm.RData"))
    normalized.expression.from.users <- as.matrix(log2(normalized.expression.from.users+1)+1)
    
  }else{                                             
    save(expression.from.users, file = paste0(save.dir, "/expression.from.users.tpm.RData"))
    if(max(expression.from.users) < 50){
      normalized.expression.from.users <- as.matrix(expression.from.users+1)
    }else{
      normalized.expression.from.users <- as.matrix(log2(expression.from.users+1)+1)
    }
   
  }
  
                                                                                              
  sample.number <- ncol(expression.from.users)
  permutation.exp <- matrix(rep(normalized.expression.from.users, times = (perm.times+1)), 
                            nrow = nrow(normalized.expression.from.users))
                                                  
  for(i in 1:perm.times){
    a <- ssGSEA.MultipleSamples.Permutation(exp.profile = normalized.expression.from.users)
    permutation.exp[, i*sample.number + 1:sample.number] <- a
  }
  rownames(permutation.exp) <- rownames(normalized.expression.from.users)
  
  exp.n <- ncol(permutation.exp)
  exp.m <- nrow(permutation.exp)
                                 
  superposition.rank.matrix <- superpositionRank(raw.matrix = permutation.exp, m = exp.m, n = exp.n)
                                                    
  permutation.exp.rank <- rankMatrixByCol(superposition.rank = superposition.rank.matrix, m = exp.m, n = exp.n)
  
                                              
  permutation.score <- t(sapply(signatureList, function(sig){
    gSetIdx <- which(rownames(permutation.exp) %in% sig)
                                                                                      
    t1 <- numeric(exp.m)
    t1[gSetIdx] <- 1
    geneSet.position <- matrix(rep(t1, times = exp.n), ncol = exp.n)
                                           
    score.for.oneset <- oneSetssGSEA.ES(exp.rank = permutation.exp.rank, superposition.rank.es = superposition.rank.matrix, 
                                        gSet.position.matrix = geneSet.position, alpha=0.25, m.num = exp.m, n.num = exp.n)
    
    return(score.for.oneset)
  }))
  print("permutation.score over!")
  
  
  
            
                                                                                                                           
            
                                                                                                                     
                                                             
  for(t in 1:perm.times){
    permutation.oneExp.score <- permutation.score[, t*sample.number + 1:sample.number]
    r <- which(permutation.score[,1:sample.number] * permutation.oneExp.score < 0) 
    permutation.oneExp.score[r] <- NA
    permutation.score[, t*sample.number + 1:sample.number] <- permutation.oneExp.score
  }
  
                           
  zScore.matrix <- sapply(1:sample.number, function(s){
                                                                                       
    perm.matrix.forOneSample <- permutation.score[,((0:perm.times) * sample.number + s)]
                          
    zScore.forOneSample <- apply(perm.matrix.forOneSample, 1, function(x) {
      if(length(which(!is.na(x))) < 2){
        zScore <- x[1]
      }else{
        
        zScore <- scale(x)[1]
      }
      return(zScore)
    })
    return(zScore.forOneSample)
  })
  print("zScore over!")
  
                                                
  positive.l <- grep("positive", rownames(zScore.matrix)) 
  if(length(positive.l) > 0){
                                                
    zScore.matrix[positive.l-1, ] <- zScore.matrix[positive.l-1, ]*(-1)
    direction <- 1:length(signatureList)
    direction[positive.l] <- direction[positive.l-1] 
  }else{
    direction <- 1:length(signatureList)
  }
  print("direction over!")
  
                                                                  
  ssGSEA.normalized.score <- rowsum(zScore.matrix, direction)
  if(length(positive.l) > 0){
    rownames(ssGSEA.normalized.score) <- rownames(zScore.matrix)[-positive.l]
    rownames(ssGSEA.normalized.score) <- gsub(".negative", "", rownames(ssGSEA.normalized.score))
  }else{
    rownames(ssGSEA.normalized.score) <- rownames(zScore.matrix)
  }
  colnames(ssGSEA.normalized.score) <- colnames(expression.from.users)
  
  ssGSEA.normalized.score <- round(ssGSEA.normalized.score, 3)
  save(ssGSEA.normalized.score, file = paste0(save.dir, "/ssGSEA.normalized.score.RData"))
  
  
  
            
                                                                                                                                     
            
                                                                                                                          
                
  final.mix.show.num <- paste('["',paste(rev(mix.show.num), collapse='","'), '"]', sep="")
  write(final.mix.show.num, file = paste0(save.dir, "/final.mix.show.num.txt"))
  print("final.mix.show.num over!")
                   
  samplePaste <- paste('["',paste(colnames(expression.from.users), collapse='","'), '"]', sep="")
  write(samplePaste, paste0(save.dir,"/samplePaste.txt"))
                              
  print("ssGSEA.normalized.score summary: ")
  print(summary(as.vector(ssGSEA.normalized.score)))
  zScore.ssGSEA.normalized.score <- t(apply(ssGSEA.normalized.score, 1, scale))
  colnames(zScore.ssGSEA.normalized.score) <- colnames(expression.from.users)
  print("zScore.ssGSEA.normalized.score summary: ")
  print(summary(as.vector(zScore.ssGSEA.normalized.score)))
  zScore.ssGSEA.normalized.score <- filterOutliers(zScore.ssGSEA.normalized.score)
  ssGSEA.normalized.score.coordinate <- makeHeatmapData(zScore.ssGSEA.normalized.score)
  write(ssGSEA.normalized.score.coordinate, file = paste0(save.dir, "/ssGSEA.normalized.score.coordinate.txt"))
                                 
  ssGSEA <- data.frame(Steps = rownames(ssGSEA.normalized.score), ssGSEA.normalized.score)
  colnames(ssGSEA)[-1] <- colnames(expression.from.users)
  write.table(ssGSEA, file = paste0(save.dir, "/ssGSEA.normalized.score.txt"), quote = FALSE, 
              row.names = FALSE, col.names = TRUE, sep = "\t")
  print("ssGSEA over!")
                             
  overall.score <- colSums(ssGSEA.normalized.score)
  overall.score <- round(overall.score, 3)
  
  a <- data.frame(Samples = colnames(ssGSEA.normalized.score), Overall.activity = overall.score)
  write.table(a, file = paste0(save.dir, "/overall.score.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  write.table(range(overall.score), file = paste0(save.dir, "/overall.score.MinMax.txt"), quote = FALSE,
              row.names = FALSE, col.names = FALSE, sep = "\n")
  print("overall score over!")
 
                                                                                                                                 
                
                
  step4.lo <- grep("Step4", rownames(ssGSEA.normalized.score))
  if(length(step4.lo) > 0){
    indicator <- 1:nrow(ssGSEA.normalized.score)
    indicator[step4.lo] <- indicator[step4.lo[1]]
    step4.com.score <- apply(ssGSEA.normalized.score, 2, function(x){tapply(x, indicator, mean)})
    step4.com.score <- round(step4.com.score, 3)
    rownames(step4.com.score) <- sort(c(rownames(ssGSEA.normalized.score)[-step4.lo], "Step4"))
  }else{
    step4.com.score <- ssGSEA.normalized.score
  }
 
  rownames(step4.com.score) <- changeStepNames(rownames(step4.com.score))
  colnames(step4.com.score) <- colnames(ssGSEA.normalized.score)
  nu <- BarPlotFormat(step4.com.score, save.dir = save.dir)
  print("BarPlotFormat over!")
  
                                        
  a <- sapply(setdiff(1:nrow(ssGSEA.normalized.score), step4.lo), function(step){ 
    b <- sapply(1:ncol(ssGSEA.normalized.score), function(sample){ 
      write.table(ssGSEA.normalized.score[step, sample], 
                  file = paste0(save.dir, "/", colnames(ssGSEA.normalized.score)[sample], "_", 
                                rownames(ssGSEA.normalized.score)[step], ".txt"), 
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    })
  })
  print("step123567 over!")
                 
  if(length(step4.lo) > 0){
    b <- sapply(1:ncol(ssGSEA.normalized.score), function(sample){ 
      radar <- data.frame(name = gsub("Step4.", "", rownames(ssGSEA.normalized.score)[step4.lo]), 
                          sample = ssGSEA.normalized.score[step4.lo, sample])
      write.table(t(radar), file = paste0(save.dir, "/", colnames(ssGSEA.normalized.score)[sample], "_Step4.txt"), 
                  quote = FALSE, row.names = TRUE, col.names = FALSE, sep = ",")
      })
  }
  print("step123567 over!")
  
                                                                                                                                     
             
  anno <- unique(signature.annotation[,3:4])
  inter <- intersect(signature.annotation[,3], rownames(expression.from.users))
  print(paste0("intersect number: ", length(inter)))
  anno <- anno[which(anno[,1] %in% inter),] 
  print(paste("dim(anno)---", dim(anno)))
  
  annotation_row_length <- tapply(anno[,2], anno[,2], length)
  print(paste("sum of annotation_row_length: ", sum(annotation_row_length)))
  
  expression.profile <- as.matrix(expression.from.users[anno[,1], ]) 
  print(paste("dim(expression.profile): ",dim(expression.profile)))
  storage.mode(expression.profile) <- "numeric"
  save(expression.profile, file = paste0(save.dir, "/expression.profile.RData"))
  
  
  heatmap.expression <- expression.profile
  rownames(heatmap.expression) = 1:nrow(heatmap.expression)
  
  annotation_row = data.frame(Steps = factor(rep(paste0("Step", names(annotation_row_length)), times = annotation_row_length)))
  rownames(annotation_row) = rownames(heatmap.expression)
  ann_colors <- list(Steps = c(Step1 = "#66C2A5", Step2="#FC8D62", Step3="#8DA0CB", Step4="#E78AC3", 
                               Step5 = "#A6D854", Step6="#FFD92F", Step7="#E5C494")) 
  
                  
  heatmap.expression <- log2(heatmap.expression+1)
  print("heatmap.expression summary before filter outliers: ")
  print(summary(as.vector(heatmap.expression)))
  print(which(apply(heatmap.expression, 1, sd) == 0))
  scale.expression <- t(apply(heatmap.expression,1,scale))
  colnames(scale.expression) <- colnames(heatmap.expression)
  print("scale.expression summary before filter outliers: ")
  print(summary(as.vector(scale.expression)))
  scale.expression <- filterOutliers(scale.expression)
  
  png(file = paste0(save.dir, "/SignatureGenes.Expression.png"))
  b <- pheatmap(scale.expression, cluster_row = FALSE, cluster_col = TRUE, 
                col=colorRampPalette(c("navy", "white", "firebrick3"))(100),
                annotation_row = annotation_row, 
                main = "Signature Genes Expression", labels_row = "", labels_col = NULL, scale = "none", 
                annotation_colors = ann_colors, border_color = NA)
  dev.off()
  
                             
  heatmap.expression <- expression.profile
  rownames(heatmap.expression) = 1:nrow(heatmap.expression)
  
  signature.expression <- data.frame(GeneSymbol = anno[,1], Steps = anno[,2], round(heatmap.expression, 3))
  signature.expression[,1] <- as.character(signature.expression[,1])
  colnames(signature.expression)[-(1:2)] <- colnames(expression.from.users)
  write.table(signature.expression, file = paste0(save.dir, "/SignatureGenes.Expression.txt"), quote = FALSE,
              row.names = FALSE, col.names = TRUE, sep = "\t")
  print("signature.expression over!")
  
                                
  heatmap.expression <- data.frame(GeneSymbol = rownames(expression.profile), heatmap.expression)
  heatmap.expression <- unique(heatmap.expression)
  print(paste("dim(heatmap.expression)----", dim(heatmap.expression)))
  data1 <- t(log2(heatmap.expression[,-1]+1))
  pca_data<- prcomp(data1, center = TRUE, scale. = FALSE) 
  
  
  cumsum.pca <- cumsum(pca_data$sdev^2)/sum(pca_data$sdev^2)
  print(paste("cumsum(pca_data$sdev^2)-----", cumsum.pca[1:3]))
  rownames(pca_data$x) <- colnames(expression.from.users)
  save(pca_data, file = paste0(save.dir, "/pca.result.RData"))
  
                           
  PCASatterPlot(pca.matrix = round(pca_data$x[,1:2],3), save.file = paste0(save.dir, "/PCA.SatterPlot.txt"))
  print("PCA over!")
  
                                                         
  box.expression <- signature.expression
  box.expression[,-(1:2)] <- log2(box.expression[,-(1:2)]+1)
  Score_boxplot(Exp.class.matrix = box.expression, saveDir = save.dir)
  print("boxplot over!")
  
  
            
                                                                                                                         
            
  result <- list(ssGSEA.normalized.score = ssGSEA.normalized.score, signature.expression = expression.profile, 
                 separate.remain = separate.remain, 
                 show.num = show.num, delete.signature = delete.signature, mix.show.num = mix.show.num)
  save(result, file = paste0(save.dir, "/multiSample.result.RData"))
  
  return(list(separate.remain = separate.remain, show.num = show.num, delete.signature = delete.signature,
              ssGSEA.normalized.score = ssGSEA.normalized.score, signature.expression = expression.profile))
}
