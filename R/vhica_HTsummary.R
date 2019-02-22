vhica_HTsummary <- function(TEName,imageObj,vhicaObj,DivRate){
  #Creating the Summary Final Comparison file
  Finalcomparison<-1
  Finalcomparison<-as.null(Finalcomparison)
  #Getting the stats matrix from the imajeObj 
  upper_matriz<-imageObj$stats
  #Keeping only the upper matrix
  upper_matriz[lower.tri(upper_matriz)]<-NA 
  #Filter matrix by logbase10(0.05) and transforming it in logical matrix
  LogicalFilter<-upper_matriz<(-1.3) 
  #Extracting the comparison which presented significant p-value 
  significantfilterNumber<-subset(as.data.frame(as.table(upper_matriz)),Freq<(-1.3))
  if(dim(significantfilterNumber)[1]==0){
    return(NULL)
  }
  #Extracting the comparison which presented significant p-value to get species pairs as dataframe
  significantfilter<-subset(as.data.frame(as.table(LogicalFilter)),Freq=="TRUE") 
  
  #Extracting from transposon matrix the ones given as arguments in the function
  transposons<-subset(vhicaObj$div,grepl(TEName,seq) & dS!="NA") 
  for(i in 1:dim(transposons)[1]){
    elemento<-unlist(strsplit(as.character(transposons[i,1]),"[.]"))
    Size<-length(elemento)
    if(Size==2){
      transposons[i,3]<-paste(transposons[i,3],".",elemento[2],sep="")
      transposons[i,4]<-paste(transposons[i,4],".",elemento[2],sep="")
    }else if(Size==3){
      transposons[i,3]<-paste(transposons[i,3],".",elemento[2],sep="")
      transposons[i,4]<-paste(transposons[i,4],".",elemento[3],sep="")
    }
  }
  
  #Walking through the TRUE elements of the logical matrix
  for(i in 1:dim(significantfilter)[1]){
    
    #Extraction species names 1 and 2
    member1<-as.character(significantfilter[i,1]) 
    member2<-as.character(significantfilter[i,2])
    
    #Extracting specific p-values from comparison between species 1 and 2 for the TE given as argument in the function
    comparisonNumero<-subset(significantfilterNumber,((Var1==member1 | Var1==member2) & (Var2==member1 | Var2==member2) & Var1!=Var2))
    
    comparison<-subset(transposons,((sp1==member1 | sp1==member2) & (sp2==member1 | sp2==member2) & sp1!=sp2))
    
    #If no DivRate given no estimates in Mya
    if(DivRate==-1){
      pvalor<-10^(comparisonNumero$Freq)
      comparison<-cbind(comparison,pvalor)
      #Adding information in the final object
      Finalcomparison<-rbind(Finalcomparison,comparison) 
    }else{
      #Calculating DivRate
      time<-comparison[1,][2]/(2*DivRate)
      colnames(time) <- c("Time(Mya)")
      comparison2<-cbind(comparison[1,],time)
      pvalor<-10^(comparisonNumero$Freq)
      comparison2<-cbind(comparison2,pvalor)
      #Adding information in the final object
      Finalcomparison<-rbind(Finalcomparison,comparison2) 
      
    }
      
  }
  return(Finalcomparison)
}
