                          
                                                                     
                                                                                          
                                                              
                                                                                                 
                                                                                                   
                                                                                    
                                                                              
                                                                             
                                                   

processResult.server <- function(filePath, CancerType, type.of.data, format.of.file, sample, sampleNumber, Samples, email){

	if(email==""){email = "no Email";}
	inputInfor <- c(CancerType, type.of.data, format.of.file, sample,sampleNumber,Samples, email, format(Sys.time(), "%Y.%m.%d.%H.%M.%S"));
	
	write.table(inputInfor, file=sprintf("%s/processResult.txt",filePath),quote = FALSE, row.names = FALSE, col.names= FALSE);
	print("All R compute finishes successfully.");
}
