vhica_HTsummary <- function(TEName,imageObj,vhicaObj,DivRate=-1){
  #Creating the Summary Final Comparison file
  Finalcomparison<-NULL
  #Getting the stats matrix from the imajeObj 
  upper_matriz<-imageObj$stats
  #Keeping only the upper matrix
  upper_matriz[lower.tri(upper_matriz)]<-NA 
  #Extracting the comparison which presented significant p-value 
  significantfilterNumber<-subset(as.data.frame(as.table(upper_matriz)),Freq<(-1.3))
  if(dim(significantfilterNumber)[1]==0){
    return(NULL)
  }
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
  for(i in 1:dim(significantfilterNumber)[1]){
    
    #Extraction species names 1 and 2
    member1<-as.character(significantfilterNumber[i,1]) 
    member2<-as.character(significantfilterNumber[i,2])
    
    #Extracting specific p-values from comparison between species 1 and 2 for the TE given as argument in the function
    comparisonNumero<-subset(significantfilterNumber,((Var1==member1 | Var1==member2) & (Var2==member1 | Var2==member2) & Var1!=Var2))
    comparison<-subset(transposons,((sp1==member1 | sp1==member2) & (sp2==member1 | sp2==member2) & sp1!=sp2))
    pvalor<-10^(comparisonNumero$Freq)
    #If no DivRate given no estimates in Mya
    if(DivRate==-1){
      comparison<-cbind(comparison,pvalor)
    }else{
      #Calculating DivRate
      time<-comparison[2]/(2*DivRate)
      colnames(time) <- c("Time(Mya)")
      comparison<-cbind(comparison,time)
      comparison<-cbind(comparison,pvalor)
    }
    #Adding information in the final object
    Finalcomparison<-rbind(Finalcomparison,comparison) 
    
  }
  return(Finalcomparison)
}
