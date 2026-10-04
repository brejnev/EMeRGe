#### midaUtil.R : Contains that are specific to MIDA pipeline ##### ----  

#### Section 1: ====== Function definition ====== 

#'Converts unlabeled data into vector
#'@param unlabeled_vector vector of average isotologue M+0, M+1, ..., M+n, for unlabeled samples.  MIDs are averaged over unlabeled samples 
#'
make.corrMat<-function(unlabeled_vector){
    l=length(unlabeled_vector)
    #Convert vector into matrix
    unlabeled<-matrix(1:l,nrow=l,ncol=l,byrow=TRUE)
    unlabeled
    for(i in 1:dim(unlabeled)[1])
    {
        for (j in 1:dim(unlabeled)[2])
        {
          unlabeled[i,j]<-unlabeled_vector[abs(j-i)+1]
        }
    }
    #Convert into triangular matrix
    unlabeled[lower.tri(unlabeled)] <- 0
    return(unlabeled)
}


#'Gets isotopologue names eg. "Citrate_M.0" "Citrate_M.1" "Citrate_M.2" ... (0, 1, 2 ...) from compound names
#'@param inpt a character vector containing isotopologue names
#'
getMs<-function(inpt){  
      ## Testing lines
      # inpt=curCols;        
      Ms=unlist(lapply(as.character(inpt),function(n){
            m=strsplit(n,"_M\\.")[[1]][2]
            return(m)}))
      return(Ms) 
}


#'Gets the group name from the sample names (biological replicates prefix) 
#'@param rname names or rownames=sample names
#'
getRepSetNames<-function(rnames){
    unique(unlist(lapply(rnames,function(n){
        v=strsplit(n,"_")[[1]]
      r=paste(v[1:(length(v)-1)],collapse="_")
      return(r)})))
}

#'Takes a matrix and generate a barplot 
#'@param mtx  the input matrix
#'@param colz  list of colors generate for each isotopologue
#'@param main_ main title of the figure 
#'@param ylab_ y-axis lable (i.e Total Enrichment)
#'
ggBar<-function(mtx,colz,main_="", ylab_=""){ 
          ## Testing lines         
          # mtx=mt; colz=colorz; main_="Total Enrichment"; ylab_="Total Enrichment"
        
          df<-melt(mtx)
          
          compd_n<- gsub("_M\\.\\d{1,3}" ,"",as.character(df$X2)[[1]])
          
          df$Isotope= sapply(as.character(df$X2),function(n){v=strsplit(n,"_")[[1]]; m=v[length(v)]; r=gsub("\\.", "+", m); return(r) } )
          
          df$Isotope=factor(df$Isotope,levels=rev(unique(df$Isotope)))
          df$X1=factor(df$X1,levels=rownames(mtx))
          options(scipen=10000)
          
          colz1=colz[unique(as.character(df$Isotope))]
          g=ggplot(df, aes(fill=Isotope, y=value, x=X1)) + 
                    geom_bar(position="stack", stat="identity",colour="black",size=.3) +
                    #scale_fill_manual(name="Isotope",values=colz)+
                    scale_fill_manual(name="Isotope",values=colz1)+
                    
                    #scale_fill_manual(name="Isotope",values=colz[unique(as.character(df$Isotope))])
                    ggtitle(paste0(compd_n," \n(",main_,")")) +
                    theme(plot.title=element_text(face="bold",hjust=0.5),axis.text.x = element_text(angle = 90, vjust = .5,colour = "black"),axis.text.y = element_text(colour = "black"),)+theme(plot.title=element_text(face="bold",hjust=0.5),panel.border = element_blank(),panel.background =element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))+
                    xlab("Samples") +
                    ylab(ylab_) 
          return(g)
}

#'Takes a matrix and generate a barplot with standard error
#'@param mtx  the input matrix
#'@param sderr standard error vector
#'@param colz  list of colors generate for each isotopologue
#'@param main_ main title of the figure 
#'@param ylab_ y-axis lable (i.e Total Enrichment)
#'
ggBar.sde<-function(mtx,sderr,colz,main_="", ylab_=""){ 
        
          df<-melt(mtx)
          
          #compd_n<-strsplit(as.character(df$X2)[[1]],".Results")[[1]]
          #df$Isotope=getMs(df$X2,compd_n)
          
          compd_n<- gsub("_M\\.\\d{1,3}" ,"",as.character(df$X2)[[1]])
          df$Isotope= sapply(as.character(df$X2),function(n){v=strsplit(n,"_")[[1]]; m=v[length(v)]; r=gsub("\\.", "+", m); return(r) } )
          
          
          df$Isotope=factor(df$Isotope,levels=rev(unique(df$Isotope)))
          df$X1=factor(df$X1,levels=rownames(mtx))
          options(scipen=10000)
          
          colz1=colz[unique(as.character(df$Isotope))]
          
          g=ggplot(df, aes(fill=Isotope, y=value, x=X1)) + 
                    geom_bar(position="stack", stat="identity",colour="black",size=.3) +
                    scale_fill_manual(name="Isotope",values=colz1)+
                    ggtitle(paste0(compd_n," \n(",main_,")")) +
                    theme(plot.title=element_text(face="bold",hjust=0.5),axis.text.x = element_text(angle = 90, vjust = .5,colour = "black"),axis.text.y = element_text(colour = "black"),)+theme(plot.title=element_text(face="bold",hjust=0.5),panel.border = element_blank(),panel.background =element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))+
                    xlab("Samples") +
                    ylab(ylab_) 
          
          for (i in 1:dim(mtx)[1]){
                    nm=rownames(mtx)[i]
                    yval=sum(df$value[which(df$X1==nm)])
                    g <-g +geom_segment(aes_string(x = i, y =yval , xend = i, yend = yval+sderr[i]),size=.3)
                    g <-g +geom_segment(aes_string(x = i-.1, y =yval+sderr[i] , xend = i+.1, yend = yval+sderr[i]),size=.3)
          }
          
          return(g)
}

#'Computes the standard error from Total FE  (sde colculated by sums of M+0, M+1, ..., M+n, all biological replicates)
#'@param df data frame of Total FE 
#'@param mt data frame of sample group name vs compound MID (mt is for one compound)
#' 
get.sde<-function(df,mt){
        
      # test lines
           # df=DF_MID_FE;  mt=mt.rep
        
      rnames.st=rownames(mt);cnames=colnames(mt)
      sub.df=df[,cnames]
      sde<-numeric()
      for (rnm in rnames.st) {
            ix=which(startsWith(rownames(sub.df),rnm))
            v=rowSums(data.matrix_bm(sub.df[ix,]))
            
            if (length(v)==1){  # if we have one biological replicate
                  sde=c(sde,0)   
            }else{  # if we have several biological replicates
                  sde=c(sde, sd(v)/sqrt(length(v))) 
            }
        
        
      } 
      return(sde)
}

#'Filters  noisy compounds after natural abundance correction by: removing compound with FE <95 or >105  in the M+0 (12C)  in 2/3 of the unlabeled samples
#'@param df input data frame (corrected data)
#'@param compoundz list of compound
#'@param L_ row indices for unlabeled / labeled samples
#'@param FEpercThresh_ FE percentage threshold (default 5)
#'
M0_FE_QC<-function(df,compoundz, L_, FEpercThresh_){
      ## Testing lines
      #  df=OutputData; compoundz=compNotesDF$Compound; L_=L; FEpercThresh_=5
      
      colnames_ = CompNameWOmid_Mida(dfIn=df)
      
      for (i in 1:length(compoundz)){
              #cat(i,"\n")
              cix=grep(compoundz[i], colnames_, ignore.case=T) 
              m=as.numeric(df[L_[["unlabeled"]],cix[1]])
              
              if (length(which(m < (100-FEpercThresh_) | m > (100+FEpercThresh_) ))/length(m)>=2/3){
                    df[,cix]<-NA
              }
      }
      
      return(df)
}



#'Imputes of missing values using the minimum detection values per compound per sample group
#'Exclude compounds below baseline detection (1000) i.e. if () 2/3 of unlabeled samples M+0 are below 1000) or all labeled sample entries are < 1000, here 1000 is adjustable threshold "baselineD"
#'@param df input data frame
#'@param compoundz list of compound names
#'@param baseLineD base-line/threshold for compound detection (default = 1000)
#'
imputMin<-function(df,compoundz,baseLineD){    ### This routine expects the order of isotopes for a particular compounds as M+0, M+1, M+2, ... M+n
  ## Testing lines
  # df=InputData1; compoundz=compNotesDF$Compound; baseLineD=baseLineDet_
  # df_bak <- InputData1
  
  snames=unique(df$Labeling)
  
  # J colnames_ <- colnames(df) %>% str_remove_all("_M\\.\\d{1,3}", "")
  colnames_=CompNameWOmid_Mida(dfIn=df)
  # print("compoundz")
  # print(compoundz)
  # print('colnames_')
  # print(colnames_)  

  # print(df%>%colnames())
  
  for (i in 1:length(compoundz)){
        #cat(i,"\n") ;cat(compoundz[i],"\n")
        
        cix=grep( paste0("^",compoundz[i],"$"), colnames_ , ignore.case = TRUE) # J LLeucine, Lleucine | NAcetylglutamate | Nacetylglutamate inconsistency in compound nameing
        print("compoundz")
        print(compoundz[i])
        print("cix - get here. line 197")
        print(cix)
        print("colnames")
        print(df[,cix] %>% colnames)
        m= data.matrix_bm(df[,cix]) # error! df[, cix] -> vector instead of data.frame.
        print('m is calculated.')
        uns=grep( "unlabeled",df$Labeling ,ignore.case = T)
        
        #Check if compound is bellow basline detection (M+0)
        
        if (length(which(m[uns,1]<baseLineD))/length(uns) > 2/3  ){ 
                message(paste0("\t\t*",compoundz[i]," : falls bellow detection baseline ",baseLineD," for M+0 unlabeled sample. It will be removed"))
                df[,cix]=NA
                next
        
        }else if ( length(which(m[-uns,] >= baseLineD))==0){
                message(paste0("\t\t*",compoundz[i]," : falls bellow baseline detection ",baseLineD," in all labeled samples. It will be removed"))
                df[,cix]=NA
                next
        
        }
        else{
              for (j in 1:length(snames)){ 
                  rix=which(startsWith(tolower(df$Labeling),tolower(snames[j])))
                  rix=grep(tolower(snames[j]), tolower(df$Labeling))
                  
                  m1=m[rix,]
                  if(sum(m1)==0){
                          cat(compoundz[i]," x ",snames[j]," : all entries are zero thus are skipped\n")
                          next
                  }else{
                  
                    #min accross sample for individual isotopologue
                    
                        if (is.matrix(m1)){
                            m2=apply(m1,2,function(cl){
                                  ret=c()
                                  if (all(cl==0)){
                                        ret=cl
                                  }else{ 
                                        min_=min(cl[which(cl>0)])
                                        cl[which(cl==0)]<- min_
                                        ret=cl
                                  }
                                  return(ret)
                            })
                            df[rix,cix]<-m2
                        }
                    }
              }
          
        }
      }
      return(df)
}


#'Removes compound with low detection
#'@param df input dataframe
#'
removeCompNA<-function(df){
      compdNA=which(colSums(is.na(df)) > 0)
      if (length(compdNA)>0){
            df <- df[ ,-compdNA]
            message(paste0("\t\t*",length(compdNA)," column(s) removed"))  
      }else{
            message("\t\t* No column removed")
      }
      return(df)
}


##  Adjusts the compound List after filtration step
## Adjusts/removes compounds from the compound note dataframe to match what remain in the Input or processed Dataframe after filtration
#'@param df input dataframe from raw/pre-processed data
#'@param compsDF compound data frame from compound notes
#'@metaCols column index for metadata in the input/processed data table (dataframe)  : required argument since all table has a metadata column eg. sample name group name
#'
adj.CompDF<-function(df,compsDF,metaCols){
        ## Testing lines
        # df= DF_TotLabel_FE; compsDF=compNotesDF; metaCols=c(1,2)
        # df=InputData
        # df=InputData1
        # df=OutputData
        
  
        colnames_<-CompNameWOmid_Mida(dfIn=df)
        if (length(metaCols)>0){colnames_=colnames_[-metaCols]} 
        
        retained_=intersect(compsDF$Compound, unique(colnames_))
        removed_=setdiff(compsDF$Compound, unique(colnames_))
        
        compsDF=compsDF[retained_,]
        
        message("\t\t", length(removed_)," compound(s) removed: ", paste0(removed_,collapse=", "))
        return(compsDF)
}


## Adjusts/removes compounds from Input Dataframe to match what remains in the compound note data frame
#'@param df_ input dataframe from raw/pre-processed data
#'@param compdsDF_ compound dataFrame from the compound notes input file
#'@metaCols_ column index for metadata, these must not be removed, when removing filtred compounds
#'
adj.DataDF<-function(df_, compdsDF_, metaCols_=c()){
       ## Test lines
        df_=InputData; compdsDF_=compNotesDF; metaCols_=c(1)
        # df_=df; compdsDF_=compoundzDF; metaCols_=c(1)
        
       compds_=compdsDF_$Compound                        
       dm1=dim(df_)[2]        
       colnames_=CompNameWOmid_Mida(dfIn=df_)
       
         
       idx=unlist(lapply(compds_,function(cmp){ r= grep( pattern=paste0("^",cmp,"$"),colnames_,ignore.case = T); return(r)}))

       df_=df_[,c(metaCols_,idx)]
       
       message("\t\t",dm1-dim(df_)[2]," column(s) removed.")
       return(df_)
}

## Substitutes NAs with Zeros in a data frame
#'@param df_ input data frame
#'
sub_Na_Zero<-function(df_){
      for (i in 1: nrow(df_)){
            for (j in 2:ncol(df_)){
                  if (is.na(df_[i,j])){df_[i,j]="0"}
                  if (df_[i,j]==""){df_[i,j]="0"}
            }
      }
      return(df_)
}

#'Substitutes negative values with Zeros
#'@param df input data frame
#'
sub_Ng_Zero<-function(df){
      for (i in 1:nrow(df)){
            for (j in 2:ncol(df)){
                  if (df[i,j]< 0){df[i,j]=0}
            }
      }
      return(df)
}

#'Calculates coefficient of variation give a numeric vector or matrix
#'@param X input numeric vector or matrix 
#'
cv=function(X){r=(sd(X,na.rm=TRUE)/mean(X,na.rm=TRUE))*100;return(r)}

#' Takes a data frame table and generate the coefficient of variation (sum all isotopologues and computing cv)
#' @param df the data frame to compute coefficient of variation for
#'@param compoundzDF compound data frame from compound notes
#'@param metaCols indices of the columns containing the metadata  in the df dta fram
#'
coev_df<-function(df,compoundzDF,metaCols=c()){
        
        ## Testing lines
        #  df=OutDf; compoundzDF=compNotesDF; metaCols=c(1,2)
        #  df=InputData1; compoundzDF=compNotesDF; metaCols=c(1)
        
        colnames_=CompNameWOmid_Mida(dfIn=df)
        if (length(metaCols)>0){colnames_=colnames_[-metaCols]}
        uniqueColnames=unique(colnames_)
        
        if (length(compoundzDF$Compound)>length(uniqueColnames)){
                compDF = adj.CompDF(df, compoundzDF,metaCols=metaCols)
                dataDF=df
        } else if (length(compoundzDF$Compound)<length(uniqueColnames)){
                dataDF = adj.DataDF(df, compoundzDF,metaCols=metaCols)                
                compDF=compoundzDF
        }else{
                dataDF=df
                compDF=compoundzDF
        }
        
        compoundz=compDF$Compound
        
        sNames=getRepSetNames(rownames(dataDF))
        
        M.cv=matrix(NA,nrow=length(sNames), ncol=length(compoundz), dimnames=list(sNames, compoundz))

        for (i_ in 1:length(sNames)){
                r.ix=grep(sNames[i_],rownames(dataDF))
                for (j_ in 1:length(compoundz)){
                        compound=compoundz[j_]
                        c.ix=grep(compound, colnames(dataDF))
                        M.cv[i_,compound]=cv(rowSums(data.matrix_bm(dataDF[r.ix,c.ix])))
                }
        } 
        
        return(M.cv)
}




#'Return the coordinates of where arrow should start and stop between two successive barplots on the TCA figure
#'@param LisFrom
#'@param LisTo
#'@param sideFrom the side of the barplot where the arrow should originate from
#'@param sideTo the side of the barplot where the arrow should point
#'
arrEnds<-function(LisFrom,LisTo,sideFrom=0,sideTo=0){
      #LisFrom=LFig[["Suc"]]; LisTo=LFig[["AKG"]]; sideFrom=4;sideTo=1
      
      if (sideFrom==0 |sideTo==0){stop("Please specify the sides")}
      vFrom=c()
      if (sideFrom==1){ vFrom=c(mean(c(LisFrom[c(1,3)])),LisFrom[2])}
      if (sideFrom==2){ vFrom=c(LisFrom[1],mean(LisFrom[c(2,4)]))}
      if (sideFrom==3){ vFrom=c( mean(LisFrom[c(1,3)]),LisFrom[4])}
      if (sideFrom==4){ vFrom=c( LisFrom[3],  mean(LisFrom[c(2,4)]))}
      vTo=c()
      if (sideTo==1){ vTo=c(mean(c(LisTo[c(1,3)])),LisTo[2])}
      if (sideTo==2){ vTo=c(LisTo[1],mean(LisTo[c(2,4)]))}
      if (sideTo==3){ vTo=c( mean(LisTo[c(1,3)]),LisTo[4])}
      if (sideTo==4){ vTo=c( LisTo[3],  mean(LisTo[c(2,4)]))}
      v=c(vFrom,vTo)
      return(v)
}

#'Makes QC plot based on CV
#'@param df input data frame
#'@param main title of the figure
#'
cvPlt<-function(df,main=""){

  #colorz=colorRamps::primary.colors()[c(1,15,4,5,12,10,13,6,7,8,9,11,14,16,17,18,19,20,21,22,23)]
  
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colorz = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  
  
  d.m<-melt(df);names(d.m)=c("SampleName","Compound","cv")
  d.m$SampleName=factor(d.m$SampleName,levels=rownames(df))
  
  g<-ggplot(d.m, aes(x=Compound, y=cv, fill=SampleName)) +
           geom_bar(stat="identity", colour="black", position="dodge", size=0.25, width=0.8, alpha=0.8) +
           
           #scale_fill_manual(values=s.colorz[1:dim(df)[1]])+
           scale_fill_manual(values=colorz)+
           #scale_fill_manual(values=c("dodgerblue4","green","darkorchid4","magenta","brown","gray"))+
           labs(fill="Groups" )+
            ggtitle(main)+
           theme(plot.title=element_text(face="bold",hjust=0.5),axis.text.x = element_text(angle = 90, vjust = .5,colour = "black"),axis.text.y = element_text(colour = "black"),)+theme(plot.title=element_text(face="bold",hjust=0.5),panel.border = element_blank(),panel.background =element_blank(),panel.grid.major = element_line(color="gray",size=0.15),panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))+
           coord_flip()
  return(g)
}

#'Makes a figure for the TCA cycle using generated bargraphs in pdf
#'@param inpath1 path to where the pdf figure of total enrichment are located
#'
tcaPlt<-function(inpath1){
        #test line
                # inpath1="/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/MIDA_test/MC00097-B16-NewMtd_out/extra"
        
        mrg=1
        arrowL=1.7 #arrow length
        spaceH=1   #Horizontal spacing
        spaceV=1   #vertical spacing
        fig.W=11   #figure width
        fig.H=10   #figure height
        
        height_=fig.H*7 + mrg*4 + arrowL*7  + spaceV*7;  height_
        width_=fig.W*5 + mrg*2 + arrowL*5  + spaceH*5;  width_
        
        #row1
        stY=0
        stX=(width_-fig.W)/2 ;                                 Suc=c(stX,stY,stX+fig.W,stY+fig.H);Suc
        stX=(width_-fig.W)/2 + fig.W*2 + spaceH*2 + arrowL*2 ; Gln=c(stX,stY,stX+fig.W,stY+fig.H);Gln
        #row2
        stY= fig.H + spaceH + arrowL
        stX=(width_-fig.W)/2 - (fig.W + spaceH + arrowL)     ; Fum=c(stX,stY,stX+fig.W,stY+fig.H);Fum
        stX=(width_-fig.W)/2 + fig.W + spaceH + arrowL       ; AKG=c(stX,stY,stX+fig.W,stY+fig.H);AKG
        stX=(width_-fig.W)/2 + fig.W*2 + spaceH*2 + arrowL*2 ; Glt=c(stX,stY,stX+fig.W,stY+fig.H);Glt
        #row3
        stY= fig.H*2 + spaceH*2 + arrowL*2
        stX=(width_-fig.W)/2 - (fig.W + spaceH*2 + arrowL)   ;  Mal=c(stX,stY,stX+fig.W,stY+fig.H);Mal
        stX=(width_-fig.W)/2 + fig.W + spaceH*2 + arrowL     ;  IsoC=c(stX,stY,stX+fig.W,stY+fig.H);IsoC
        #row4
        stY= fig.H*3 + spaceH*3 + arrowL*3 ;
        stX=(width_- (2*fig.W+ spaceH + arrowL ))/2        ;  OAA=c(stX,stY,stX+fig.W,stY+fig.H);OAA
        stX=(width_+ (spaceH + arrowL ))/2                 ;  Cit=c(stX,stY,stX+fig.W,stY+fig.H);Cit
        #row5
        stY= fig.H*4 + spaceH*2+  arrowL*3
        stX=(width_-fig.W)/2 - (fig.W*2 + spaceH + arrowL*2)    ; Asn=c(stX,stY,stX+fig.W,stY+fig.H);Asn
        stX=(width_-fig.W)/2 - (fig.W + spaceH + arrowL)    ; Asp=c(stX,stY,stX+fig.W,stY+fig.H);Asp
        stX=(width_-fig.W)/2 + fig.W*2 + spaceH + arrowL*2  ; Ita=c(stX,stY,stX+fig.W,stY+fig.H);Ita
        
        #row6
        stY= fig.H*5 + spaceH*5+  arrowL*4 
        stX=(width_-fig.W)/2 - (fig.W + spaceH + arrowL)  ; Lac=c(stX,stY,stX+fig.W,stY+fig.H);Lac
        stX=(width_-fig.W)/2                              ; Pyr=c(stX,stY,stX+fig.W,stY+fig.H);Pyr
        stX=(width_-fig.W)/2 + fig.W + spaceH + arrowL    ; Ala=c(stX,stY,stX+fig.W,stY+fig.H);Ala
        #row7
        stY= fig.H*6 + spaceH*6+  arrowL*5 
        stX=(width_-fig.W)/2                              ; Glu6P=c(stX,stY,stX+fig.W,stY+fig.H);Glu6P
        #row8
        stY= fig.H*7 + spaceH*7+  arrowL*6 
        stX=(width_-fig.W)/2                              ; Glu=c(stX,stY,stX+fig.W,stY+fig.H);Glu
        
        
        # 
        LFig=list(Glu,Glu6P,Pyr,Lac,Ala,Cit,Ita,IsoC,AKG,Glt,Gln,Suc,Fum,Mal,OAA,Asp,Asn)
        LTca<-list("Glucose","Glucose-6-P","Pyruvate","Lactate","L.Alanine","Citrate","Itaconate","Isocitrate","alpha.Ketoglutarate","L.Glutamate","L.Glutamine","Succinate","Fumarate","Malate",
                   "Oxalacetate","Aspartic.Acid","L.Asparagine")
        names(LFig)<-names(LTca)<-c("Glu","Glu6P","Pyr","Lac","Ala","Cit","Ita","IsoC","AKG","Glt","Gln","Suc","Fum","Mal","OAA","Asp","Asn")
        
        par(mar=rep(mrg,4))
        plot(0,0, type="n",xaxt="n", yaxt="n", xlab="", ylab="",main="TCA cycle",
             xlim=c(0,width_), ylim=c(0,height_))
        
        for (i in 1:length(LFig)){
            AbrComp=names(LFig)[i]
            comp=LTca[[AbrComp]]
            
            if (file.exists(paste0(inpath1,"/",comp,"_TE.repset.pdf"))){
                im=image_read_pdf(paste0(inpath1,"/",comp,"_TE.repset.pdf"),density=300)
                rasterImage(im, LFig[[i]][1],LFig[[i]][2], LFig[[i]][3],LFig[[i]][4])  
            }else{
                text(LFig[[i]][1]+fig.W/2,LFig[[i]][2]+fig.H/2,labels =comp,cex=.5)
            }
        }
        
        arr_col='blue'
        LArr<-list()
        LArr[[1]]<-c(arrEnds(LFig[["AKG"]],LFig[["Suc"]],1,4),.4)
        LArr[[2]]<-c(arrEnds(LFig[["Suc"]],LFig[["Fum"]],2,1),.4)
        LArr[[3]]<-c(arrEnds(LFig[["AKG"]],LFig[["Glt"]],4,2),0)
        LArr[[4]]<-c(arrEnds(LFig[["Glt"]],LFig[["Gln"]],1,3),0)
        LArr[[5]]<-c(arrEnds(LFig[["Fum"]],LFig[["Mal"]],3,1),.2)
        LArr[[6]]<-c(arrEnds(LFig[["Mal"]],LFig[["OAA"]],3,2),.4)
        LArr[[7]]<-c(arrEnds(LFig[["Cit"]],LFig[["IsoC"]],4,3),.4)
        LArr[[8]]<-c(arrEnds(LFig[["OAA"]],LFig[["Asp"]],2,1),.4)
        LArr[[9]]<-c(arrEnds(LFig[["OAA"]],LFig[["Cit"]],4,2),.2)
        LArr[[10]]<-c(arrEnds(LFig[["Cit"]],LFig[["Ita"]],4,2),.2)
        LArr[[11]]<-c(arrEnds(LFig[["Pyr"]],LFig[["Cit"]],1,3),-.4)
        LArr[[12]]<-c(arrEnds(LFig[["Pyr"]],LFig[["Ala"]],4,2),0)
        LArr[[13]]<-c(arrEnds(LFig[["Pyr"]],LFig[["Lac"]],2,4),0)
        LArr[[14]]<-c(arrEnds(LFig[["IsoC"]],LFig[["AKG"]],1,3),.1)
        LArr[[15]]<-c(arrEnds(LFig[["Glu6P"]],LFig[["Pyr"]],1,3),0)
        LArr[[16]]<-c(arrEnds(LFig[["Glu"]],LFig[["Glu6P"]],1,3),0)
        LArr[[17]]<-c(arrEnds(LFig[["Asp"]],LFig[["Asn"]],2,4),0)
        
        for (i in 1:length(LArr)){
          iArrows( LArr[[i]][1],LArr[[i]][2],LArr[[i]][3],LArr[[i]][4], h.lwd=1, sh.lwd=1, curve=LArr[[i]][5], sh.col=arr_col,size=.2,width=1.2)
        }
}       

#'Generates my own progress bar
#'@param prg_ current execution step
#'@param LEN_ number of execution steps
#'
get.prog<-function(prg_,LEN_){
        n=ceiling(prg_*100/LEN_)
        bar=paste0(rep("|",floor(n/10)),collapse = "")
        
        strng=paste0("|",bar,paste0(rep(".",10-nchar(bar)),collapse = "")," ",n,"%|")
        assign("prg", prg_+1,envir=parent.frame(n=1))
        return(strng)
}

#'Generates sample names using the Study factors from the metadata file
#'@param fileP full path to the metadata file
#'
sFactors2Names_Mida<-function(fileP){
        
        ## Test line
        # fileP=paste0(inPath_,"/",mtdFile)
        
        metaData=read_excel(path=fileP)
        metaData=as.data.frame(metaData)
        
        sF.newix=grep("factor",tolower(colnames(metaData)))
        
        ## Remove the study factor column of all NA's
        sF.ix=c();  for (ix in sF.newix){if (!all(is.na (metaData[,ix]))){sF.ix=c(sF.ix,ix)}}
        
        Replicate=c()
        Counter=list()
        LabelingS=c()
        
        for (i in 1:dim(metaData)[1]){
                
                if (length(which(is.na(metaData[i,sF.ix])))==0){    ## get sample names
                        curName=paste0(metaData[i,sF.ix], collapse = "_")
                        if (!curName %in% names(Counter)){Counter[[curName]]<-1}else{Counter[[curName]]=Counter[[curName]]+1}
                        Replicate=c(Replicate, paste0(curName,"_",Counter[[curName]]))
                        LabelingS=c(LabelingS,metaData[i,"Labeling"])
                }
                else if  (length(which(is.na(metaData[i,sF.ix])))==length(sF.ix)){ 
                          Replicate=c(Replicate, metaData$CoreID[i] )           ## get QC and Blank
                          LabelingS=c(LabelingS,metaData[i,"Labeling"])
                          
                }else{
                        message("Input error: missing study factors. Please check your metadata file, in study factor colunms\n","column ",i," : ", paste0(metaData[i,sF.ix],collabse=" "))
                        stop(" --- MOLAR has fail to completed the analysis ---")
                }
        }
        
        Df=cbind(metaData[1:2], Replicate,"Labeling"=LabelingS,metaData[3:dim(metaData)[2]])
}

#'Reads the input files for MIDA
#'@param inPath_ full path of the input directory containing all three inputfiles of MIDA
#'@param prefix_ prefix of study name. The input file name without "_MIDA_RawData.xlsx", "_MIDA_CmpdNotes.xlsx" or "_MIDA_Metadata.xlsx" suffix
#'
read_INPUT_files_Mida<-function(inPath_=inPath,  prefix_=prefix){
  
        # test lines
        # inPath_=inPath_; prefix_=prefix_
        
        inputFile=paste0( prefix_,"_MIDA_RawData.xlsx")
        cmpdFile=paste0( prefix_,"_MIDA_CmpdNotes.xlsx")
        mtdFile=paste0( prefix_,"_MIDA_Metadata.xlsx")
        
        prefix_List_=list(gsub("\\.xlsx$","",inputFile), gsub("\\.xlsx$","",cmpdFile),gsub("\\.xlsx$","",mtdFile))
        
        ## Check files : metadata and compound table
        if(!file.exists(paste0(inPath_,"/",inputFile))) {ret=paste0("\t\tError :", basename(inputDataFile_), " has no raw data file associated with it!"); return(ret)}
        if(!file.exists(paste0(inPath_,"/",mtdFile))) {ret=paste0("\t\tError :", basename(inputDataFile_), " has no metadata file associated with it!"); return(ret)}
        if(!file.exists(paste0(inPath_,"/",cmpdFile))) {ret=paste0("\t\tError :", basename(inputDataFile_), " has no compound file associated with it!"); return(ret)}

        # Read raw data 
        InputData_=read_excel(paste0(inPath_,"/",inputFile),col_names = TRUE)
        InputData_=as.data.frame(InputData_)
        colnames(InputData_)[1:5]=c("Sample","Data File","Type","Level","Acq. Date-Time")
        InputData_=InputData_[-1,]
        InputData_=as.data.frame(InputData_)
        rownames(InputData_)=InputData_[,1]
        
        # Read compound data
        compNotesDF_=read_excel(  paste0(inPath_,"/",cmpdFile),skip =13) #,col_names = TRUE)
        compNotesDF_=as.data.frame(compNotesDF_)
        
        midsIx<-grep("M\\+\\d{1,2}$",compNotesDF_$Compound,ignore.case = T)
        if (length(midsIx)>0){
                compNotesDF_=compNotesDF_[-midsIx ,]
        }
        
        rownames(compNotesDF_)=compNotesDF_$Compound
        
        # Read metadata
        ## Generates sample names using the Study factors from the metadata file
        metaData_=sFactors2Names_Mida(paste0(inPath_,"/",mtdFile))
        
        inMetaFile_<<-paste0(inPath_,"/", prefix_,"_MIDA_MetadataE.csv")
        write.csv(metaData_,inMetaFile_,row.names = F,quote=F)
        metaData_=as.data.frame(metaData_)
        rownames(metaData_)=metaData_[ ,"Replicate"]
        
        return(list(InputData=InputData_, compNotesDF=compNotesDF_, metaData=metaData_, prefix=prefix_, prefix_List=prefix_List_))
}


#'Get the compound names list without the MID letters M+0, M+1, M+2 ...
#'@param dfIn input vector of compound isotopologues
#'
CompNameWOmid_Mida<-function(dfIn){
      ## Test lines        
      # dfIn=df_   
      tColnames_=gsub(pattern="_M.\\d{1,3}$","",colnames(dfIn) )
      return(tColnames_)
}


#'Given an MID table, it sum up the values of all isotopologues (MIDA format)
#'@param df_ input dataframe
#'@param comoundz vector containing compound names
#'
sum_MIDs_Mida<-function(df_,compoundz){
    ## Test line
    #  df_=InputData_sub; compoundz=compounds
    
    df1=df_[,c(1,2)]
    df2=df_[,-c(1,2)]
    
 
    tColnames_<-CompNameWOmid_Mida(dfIn=df2)
    
    list_out= lapply(compoundz, function(cmp){
        ##cix=which(tColnames_s==cmp)
        cix=grep( paste0("^",cmp,"$"),tColnames_,ignore.case = T)
        m=data.matrix_bm(df2[,cix])
        rowSums(m)
    })
    
    DF_out=do.call(cbind, list_out)
    colnames(DF_out)=compoundz
    DF_out=cbind(df1,DF_out)
    
    return(DF_out)
}


#'Gets the name of compounds to remove (those with no group with at >5 TotalLabel_FE least 2/3 of samples)
#'@param df input dataframe
#'@param  MinLabelThresh_ minimum threshold for total labeling (sum 13C or sum of M+1,M+2, ..., M+n)
#'
compToFilter_sum_13C<-function(df, MinLabelThresh_){
      ## Test lines
      # df=DF_TotLabel_FE;  MinLabelThresh_=5
      
      groups=unique(df$group)
      comps<-c()
      for (i in 3:dim(df)[2]){
        
            cpName = colnames(df)[i]
            
            cpCol=sapply(groups,function(grp){
                  ix=which(df$group==grp)    
                  col_=df[ix,i]
                  
                  ret=0
                  if ((length(which(col_ > MinLabelThresh_ ))/length(col_)) >=2/3){   #getting group with 2/3 samples with >5%
                        ret=1
                  }
                  return(ret)
            })
            if (max(cpCol)==0){
                  comps=c(comps,cpName)
            }
      }
      return(comps)
}


#'Renames the start of a compound name that starts with a digit.  ie. 5 -> X.5
#'@param vect the input vector containg compound names
#'
RenameStart<-function(vect){
        for (ij in 1:length(vect)){
                if (length(grep("^\\d",vect[ij]))>0){
                        vect[ij]=paste0("X.",vect[ij])
                }
        }
        return(vect)
}

#### Section 3: ======   Main execution function of the MIDA program    ==========

#'Main execution of the MIDA pipeline
#'@param inputDataFile_ full path of the input file
#'@param prefix_  prefix of the current output file  e.g "MC00028"
#'@param inPath_ the path of the input directory
#'@param baseLineDet_ base-line/threshold for compound detection (default = 1000)
#'@param FEpercThresh_ FE percentage threshold (default 5). It is used to filter compounds with 2/3 of unlabeled samples having FE <95 or >105, after correction for natural abundance
#'@param MinLabelThresh_ minimum threshold for total labeling (sum 13C or sum of M+1,M+2, ..., M+n), if 2/3 for sample are below, the compound is removed
#'@param topMetaboFDR__ FDR threshold to determining metabolites with significant difference between groups 
#'@param normMethod__ normalization method ("none", "quantile normalization", "normalization by reference feature" etc.),default set to "none"
#'@param transMethod__ data transformation method ("log normalization", "cubic-root normalization")
#'@param scalingMethod__ data scaling method (mean-centering, auto-scaling, pareto-scaling or range-scaling)
#'@param KeggPathwayAnalysisCheck_ determines whether or not do perform pathway analysis
#'@param KeggPathwayEnrichMethod_  Kegg pathway enrichment method (MSEA,QEA or ORA)
#'@param KeggPathwayOrganism_ the organism being studied
#'@param KeggPathwayFDR_ FDR threshold for significantly up-regulated pathways between two groups
#'@param run.Type_ indicates whether it is a tracer or profiling experiment
#'
mainFunMida<-function(
                  inputDataFile_,
                  prefix_,
                  inPath_,
                  
                  baseLineDet_,
                  FEpercThresh_,
                  MinLabelThresh_,
                  
                  topMetaboFDR_,
                  normMethod_,
                  transMethod_,
                  scalingMethod_,
                  
                  KeggPathwayAnalysisCheck_,
                  KeggPathwayEnrichMethod_,
                  KeggPathwayOrganism_,
                  KeggPathwayFDR_,
                  

                  run.Type_="tracer"
                 
){
     
     
        # *** Testing lines ----
        # ****** Start executing this block below to systematically run each line/section of the code in order to debug the program.
        testing="NO"
        if (testing=="YES") {
          
                  app.Path<<-"/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/"
                  app.Path<<-"/Volumes/projects_secondary/mnp/tool.dev/EMeRGe1.0/"
                  app.Path<<- '/Users/Jay.Kim/Rproject/EMeRGe1.0'
                  setwd(app.Path)
                  
                  
                  source("idxUtils.R")
                  source("molarUtils.R")
                  source("midaUtils.R")
                  
                  .on.public.web<<-FALSE
                  api.base <<- "api.metaboanalyst.ca"
                  api.base <<- "http://api.xialab.ca"
                  
                  methods_regex <<- "HILIC_[AB]_|LIPID_|T3pos_|PHpos_|PHneg_|PHTCA_"
                  methods_regex2<<-gsub("_\\|","\\|",methods_regex); methods_regex2=gsub("_$","",methods_regex2)
          
                  idx_PB_Blanks_pattern<<-"Blank|pb|qc|pos.id|neg.id|sample_"
                  
                  inputDataFile_="E:/RHome/Programs/EMeRGe1.0/testData/MIDA_test/MC00058_MultiTraceTCells-NewMtd_MIDA_CmpdNotes.xlsx"
                  inputDataFile_= "/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/MIDA_test/MC00097-B16-NewMtd_MIDA_CmpdNotes.xlsx"
                  inputDataFile_= "/Users/Brejnev.Muhire/Downloads/GC_FAMES_MC00131_MIDA_Metadata.xlsx"
                  inputDataFile_= "/Users/Brejnev.Muhire/Downloads/GC_FAMES_MC00139_Cells_MIDA_RawData.xlsx"
                  inputDataFile_= "/Users/Brejnev.Muhire/Downloads/MC00115_GC_MIDA_RawData.xlsx"
                  # good example file
                  inputDataFile_ <- "/Users/Jay.Kim/Rproject/EMeRGe1.0/test_data/MIDA/MC00140_Cells/GC_MC00140_MIDA_RawData.xlsx"
                  inputDataFile_ <- "/Users/Jay.Kim/Rproject/EMeRGe1.0/promblematic_data_mida/GC_MC00175_TBDMS_MIDA_RawData.xlsx"
                  
                 
                  inPath_<<-dirname(inputDataFile_)
                  prefix_<<-gsub(pattern="_MIDA_RawData\\.xlsx$|_MIDA_Metadata\\.xlsx$|_MIDA_CmpdNotes\\.xlsx$" ,"", basename(inputDataFile_))
                  
                  baseLineDet_=1000
                  FEpercThresh_=5
                  MinLabelThresh_=5
                  
                  topMetaboFDR_=0.001; normMethod_="None"; transMethod_="LogNorm"; scalingMethod_="MeanCenter"
                  preproSignalCorr.Flg_=FALSE; signalCorr.Flg_=FALSE; qc.blank.Flg_=TRUE; NA.thr.Flg_=TRUE; imputation.Flg_=TRUE
                  KeggPathwayEnrichMethod_="QEA"; KeggPathwayAnalysisCheck_=TRUE
                  KeggPathwayOrganism_="mmu"; KeggPathwayFDR_=0.1
                  
                  run.Type_="tracer"
                  currModule<<-"MIDA"
                  
                  cD<<-read_INPUT_files_Mida(inPath_=inPath_,  prefix_=prefix_)
                  
        }
      
        ## Checking of program settings ----
        options(stringsAsFactors = F, scipen=999)
        

        message("#---------------------  Current settings  ----------------------------------------#")
        message("\tWorking directory: ", inPath_)
        message(paste0("\tInput file:  ",inputDataFile_)) 
        message("\tBaseline detection: ", baseLineDet_)
        message("\tUnlabeled FE % threshold: ", FEpercThresh_)
        message("\tMinimum % labeling threshold: ", MinLabelThresh_)
        message("#---------------------------------------------------------------------------------#\n")
        
        ## Checking that the studySettings are valid
        load(file=paste0(inPath_,"/",prefix_,"_studySettings.RData"))
        
        
        ## Create the output directory
        subDir=paste0(inPath_,"/",prefix_,"_out/")
        subDir2=paste0(subDir,"extra")
        
        
        if (!dir.exists(subDir)){dir.create(subDir)} else {unlink(subDir, recursive = TRUE); dir.create(subDir)}
        setwd(subDir); dir.create(subDir2,recursive=TRUE)
        
        
        message(paste0("   * Current output directory : ", paste0(subDir)))
        
        
        LEN=11
        prg=1
        
        
        message(paste0("\n\n*-------------  Correction for ", prefix_ ," ...  --------------"))
                
        outF_MID_FE=paste0(prefix_,"_MID_FE.csv")    # Isotope natural abundance corrected values
        outF_PoolSize=paste0(prefix_,"_PoolSize.csv")
        outF_TotalLabel_FE=paste0(prefix_,"_TotalLabel_FE.csv")
                
        outHtmp="Heatmap.pdf"
                
        
        ### Step 1: Pre-processing of input file ----
        # Reading the tab- or comma-separatated values data and metadata
        message(paste0("\n",get.prog(prg,LEN)," * STEP 1 --- Reading the compound abundance table"))
        
        # Input files for uncorrected data, metadata and compound note from cD list of tables made in the server function
        InputData <<- cD$InputData;       
        metaData <<- cD$metaData     
        compNotesDF<<-cD$compNotesDF
        
        
        #Reformat the isotopologue names
        #reformat MIDA and cmpdNote compound Isotopologue names
        #isotopologues
        colnames(InputData)=gsub(" Results$|Results$","",colnames(InputData))
        colnames(InputData)=gsub(" M\\+|M\\+","_M.",colnames(InputData))
        for (jj in 6:dim(InputData)[2]){
                nn=colnames(InputData)[jj]
                if (length(grep("_M\\.\\d{1,2}$",nn))==0){
                        colnames(InputData)[jj]=paste0(nn,"_M.0")  # looks like replce space in colnames with _Mxxx - do InputData %>% colnames
                }      
        }
        
        #All compound names in input and  cmpdNote
        colnames(InputData)=gsub("\\s|'",".",colnames(InputData))
        colnames(InputData)=gsub("-","_",colnames(InputData))
        colnames(InputData)=RenameStart(colnames(InputData))
        
        compNotesDF$Compound=gsub("\\s|'",".",compNotesDF$Compound)
        compNotesDF$Compound=gsub("-","_",compNotesDF$Compound)
        compNotesDF$Compound=RenameStart(compNotesDF$Compound) # J any starting digit get a "X." in front of it (5-Glucose -> X.5-Glucose)
        rownames(compNotesDF)=compNotesDF$Compound

        
        ## Initial check for compound data ---- 
        # see whether the compounds in the compound note file list match those in the data table
        # And select only useful columns { Sample Name 'Name', and all compounds isotopogue columns }
        ex=c(); 
        ex.i=c()
                
        for (i in 1:length(compNotesDF$Compound)){
                cmpd=compNotesDF$Compound[i]
                        
                L=grep(pattern=paste0("^",cmpd),colnames(InputData),ignore.case = T)
                        
                ex.i<-c(ex.i,L)
                if (length(L)==0){ex=c(ex,0)}else{ex=c(ex,1)}
        } # J looks like looking for colnames from Inputdata match the compNotesDF$Compund (compNotesDF$Compond - reference)
                
        if (sum(ex)<length(compNotesDF$Compound)){
                message("\tThe following compounds are not found in the uncorrected input table,\nPlease check your files:")
                message(compNotesDF$Compound[which(ex==0)])
                message(paste (inputDataFile_," is not processed \n------------------"))
                return()
                
        }else{   #keep relevant columns
                rmv.idx=which(!(1:length(colnames(InputData))%in%ex.i))
                InputData=InputData[,-rmv.idx]
        }
                
                
        # Check the ordering of Isotopes M+0 M+1 M+3 (to make sure correction i.e matrix multiplication does not go wrong)
        allCols=colnames(InputData)
        ixs=c()
        for (namCPD in compNotesDF$Compound){
                
                ix=grep(pattern=paste0("^",namCPD),allCols,ignore.case = T)
                curCols=allCols[ix]  
                Ms=getMs(curCols)
                MsN=as.numeric(Ms)
                if (!identical(sort(MsN),MsN)){
                       message("\tWarning: Isotopologue in wrong order found:")
                       ixs=c(ixs,ix)
                       retF="Error"
                }
        }
        if (length(ixs)>0){
                message("\tThe wrong order of Isotopologues detected!\nPlease check your files and try again.")
                for (n in allCols[ixs]){cat(n,"\n")}
                message(paste0(inputDataFile_," : is not processed \n-----------------------"))
                return()
        }
                
                
        # Make selection of compounds or use all compounds  [those indicated in the compound file]
        compNotesDF_bak <- compNotesDF # J back up this dataframe
        compNotesDF=compNotesDF[- grep("x",compNotesDF$Response,ignore.case=T),] # J filter out Response == x
        InputData1 <- adj.DataDF(InputData,compNotesDF, metaCols=c())  
                
        #============ starting here, allow to be adj.sting  compounds or compNotesDF ----
        
        message(paste0("\t\tList of included compounds (",length(compNotesDF$Compound),") :",paste0(compNotesDF$Compound,collapse = ", ")))
   
       
        ## Read/add the metada ----
        message(paste0("\n",get.prog(prg,LEN)," * STEP 2 --- Reading the associated metadata"))
        
        
        if (!identical(sort(metaData$CoreID), sort(rownames(InputData1)))){
              
              meta=setdiff(metaData$CoreID,rownames(InputData1))
              raw=setdiff(rownames(InputData),metaData$CoreID)
              
              if (length(meta) > 0){
                      message (paste("   - Metedata unmatched names : ",paste0(meta,collapse = ", ")))
              }
              if (length(raw)>0){
                      message ("   - InputData unmatched names: ",paste0(raw,collapse = ", "))
              }
              message("   - The sample names in the InputData table do not match those in the metaData table.\nPlease check the files and start again.")
              setwd(getCurrentFileLocation()) 
              return()
              
        }
        
        
        Sample.Names=unlist(lapply(rownames(InputData1),function(rn){ metaData$Replicate[which(metaData$CoreID==rn)]}))
        Labeling.Info=unlist(lapply(rownames(InputData1),function(rn){ metaData$Labeling[which(metaData$CoreID==rn)]}))
        
        rownames(InputData1)=Sample.Names
        # J 
        InputData1_bak <- InputData1
        InputData1=cbind(Labeling=Labeling.Info,InputData1)
        
        
        # Also select relevant rows for  unlabeled and labeled - don't select "blank"
        rix<-grep(pattern="^unlabeled|^labeled",InputData1$Labeling,ignore.case = T)
        InputData1=InputData1[rix,]
        
        message(paste0("\n",get.prog(prg,LEN)," * STEP 3 --- Substitute missing values with zeros"))
        InputData1_sub_na_zero_bak <- InputData1
        # J replace <NA> with 0
        InputData1=sub_Na_Zero(InputData1)
       
        message(paste0("\n",get.prog(prg,LEN)," * STEP 4 --- Imputation of missing values using the minimum for isotope in each group"))
        message(paste0("\n",get.prog(prg,LEN),"      --- Removing all low detection compound(s)"))
        
        print("compoundz")
        print(compNotesDF$Compound)
        
        print("colnames_")
        print(InputData1 %>% colnames)
        ## the below line causes an error. 'x' in selecting a method for function 'as.matrix' : argument of length 0
        InputData1=imputMin(df=InputData1, compoundz=compNotesDF$Compound, baseLineD=baseLineDet_) 
        InputData1=removeCompNA(InputData1)   ;   if ( dim(InputData1)[2]<2){ message("\n      All compounds were filtered out. Sorry the analysis cannot proceed!");  return()}
        compNotesDF<-adj.CompDF(InputData1,compNotesDF,metaCols = c(1))  ### removes compounds that were removed in the OutputData table
        
        
        message(paste0("\n",get.prog(prg,LEN)," * STEP 5 --- Correction for remaining isotopes (",length(compNotesDF$Compound ) ,") natural abundance"))
        
        
        ## Get the indeces  of labeled and unlabeled
        sNames=unique(InputData1$Labeling)
        L=list()
        counter=1
        for (i in 1:length(sNames)){
                L[[counter]]<-which(InputData1$Labeling==sNames[i])
                counter=counter+1
        }
        names(L)<-tolower(sNames)
        
        OutputData<-InputData1
        OutputData[1:dim(OutputData)[1],2:dim(OutputData)[2]]=0
        
        
        ## Isotope natural abundance correction ----
        lowDet=c()
        for (compound in compNotesDF$Compound){
                
                message(paste0("\t\t Correcting for  ",compound))
          
                # Making the correction matrix based on unlabeled samples
                colIx=grep(pattern=paste0("^",compound),colnames(InputData1),ignore.case=T )
                
                if (length(colIx)==0){next}
                
                #unlabeledV= colMeans(data.matrix_bm(InputData1[L[[1]],colIx]))
                
                tmMat=InputData1[L[["unlabeled"]],colIx]
                tmMat <- data.matrix_bm(tmMat[,1:dim(tmMat)[2]])
                unlabeledV= colMeans(tmMat)
                
                LowDetection=FALSE
                if (sum(unlabeledV)==0 ){      #(length(which(unlabeledV==0))>0 ){  
                        LowDetection=TRUE; message(paste0("\t\t * Warning * : Low detection for ",compound)) 
                }else{
                        corrMat=make.corrMat(unlabeledV)
                        if  (rankMatrix(corrMat)[1]< dim(corrMat)[1]){
                                LowDetection=TRUE; message(paste0("\t\t * Warning * : The rank natural abundance correction matrix [CM] is < the no. of its rows or columns. CM cannot be inverted for ",compound))
                        }else{
                                Inv<-solve(corrMat)     
                        }
                }
                
                # make correction of labeled samples
                for (i in 1:length(L)){
                        #Take inverse of unlabled sample matrix
                        for (ri in L[[i]]){
                                if (LowDetection==TRUE){
                                        OutputData[ri,colIx]<- 0
                                        lowDet=c(lowDet,compounds)
                                }else{
                                        labeledV=as.numeric(InputData1[ri,colIx])
                                        
                                        if (sum(labeledV)==0){
                                                OutputData[ri,colIx]<- 0
                                                message(paste0("\t\t * Warning * : Low detection for ",compound));
                                                lowDet=c(lowDet,compound)
                                        }else{
                                                #Matrix multiplication
                                                outVec<-labeledV%*%Inv
                                                #Normalize
                                                outVecN<-(outVec*100)/sum(outVec)
                                                OutputData[ri,colIx]<- round(outVecN,10)
                                        }
                                }
                        }
                }
        }
        
        
        
        ## Post processing of corrected values ----
        message(paste0("\n",get.prog(prg,LEN)," * STEP 6 --- Generating tables  ..."))
        # Put back the negative values with zeros
        OutputData<-sub_Ng_Zero(OutputData)    
        
        # Filtration of compounds with eg. FE <95 or >105 (FEpercThresh_=5%) in the M+0 (12C)  in 2/3 of the unlabeled samples  ----
        message (paste0("\t\tQC : Unlabeled FE % threshold =",FEpercThresh_," | removing the compound(s) with fractional enrichment <",(100-FEpercThresh_), " or >",((100+FEpercThresh_))," in unlabeled samples"))
        
        OutputData<-M0_FE_QC( OutputData, compNotesDF$Compound, L, FEpercThresh_)
        
        OutputData=removeCompNA(OutputData); if ( dim(InputData1)[2]<2){ message("\n      All compounds were filtered out. Sorry the analysis cannot proceed!");  return()}
        compNotesDF<-adj.CompDF(OutputData,compNotesDF,metaCols = c(1)) ### remove compounds that were removed in the OutputData table
        message("\t\tremaining compounds (",length(compNotesDF$Compound),")")
        
        # Remove the unlabeled samples
        OutputData=OutputData[-L[["unlabeled"]],]
        
        ## regenerate samples for groups with less than 3 replicates   ** IMPORTANT ** ----
        
  
        ## Make table: MID_FE Total_Label and PoolSize  ----
        
        # Table1 MID FE
        sample=rownames(OutputData)
        group=gsub("_\\d{1,3}$","",sample)
        DF_MID_FE=cbind(sample,group,OutputData[,2:dim(OutputData)[2]])
       
        # Table2 PoolSize 
        InputData_sub= cbind( DF_MID_FE[,c("sample","group")],    InputData1[rownames(DF_MID_FE),colnames(DF_MID_FE)[3:dim(DF_MID_FE)[2]] ])
        DF_PoolSize<-sum_MIDs_Mida(df=InputData_sub,compoundz=compNotesDF$Compound)          ## Summing up with M+0,M+1 ... M+n or uncorrected
       
        # Table3 Total Label FE
        #ix=grep( "M\\.\\d{1,3}\\.Results$|M+d{1,3} Results$|",colnames(DF_MID_FE))
        ix=grep("M\\.0$",colnames(DF_MID_FE))
        
        DF_MID_FE_C13=DF_MID_FE[,-ix]
        DF_TotLabel_FE=sum_MIDs_Mida(df=DF_MID_FE_C13, compoundz=compNotesDF$Compound)
        
        
        # Filtration by sum(M+1, M+2 … M+4) , sum(13C) >  5 at least in one biological group ----
        
        message(paste0("\t\tQC: filtration by  Total labeling Fractional Enrichment, sum(13C) >  ", MinLabelThresh_," at least in one biological group"))
        compToRemove<-compToFilter_sum_13C(df=DF_TotLabel_FE, MinLabelThresh_)
        
        if (length(compToRemove)>0){
                
              ie=sapply(compToRemove,function(comp){ix=which(colnames(DF_TotLabel_FE)==comp); return(ix)})
              DF_TotLabel_FE=DF_TotLabel_FE[,-ie]
              
              iz=sapply(compToRemove,function(comp){ix=which(colnames(DF_PoolSize)==comp); return(ix)})
              DF_PoolSize=DF_PoolSize[,-iz]
            
              im=as.numeric(unlist(sapply(compToRemove,function(comp){ix=grep(comp,colnames(DF_MID_FE)); return(ix)})))
              DF_MID_FE=DF_MID_FE[,-im]
                      
              compNotesDF<-adj.CompDF(DF_TotLabel_FE,compNotesDF,metaCols = c(1,2) ) ### removes compounds that were removed in the OutputData table

              message("\t\tremaining compounds (",length(compNotesDF$Compound),")")
              
              
              if ( dim(InputData1)[2]==2 | dim(DF_PoolSize)[2]==2 |  dim(DF_MID_FE)[2]==2){ message("\n      All compounds were filtered out. Sorry the analysis cannot proceed!");  return()}
              
        }
        
        
        
        # Write the three first tables
        message( "\n        - Generating table :    ",outF_MID_FE)
        write.csv(DF_MID_FE,file=paste0(subDir, outF_MID_FE),row.names = F)
        message( "        - Generating table :    ",outF_PoolSize)
        write.csv(DF_PoolSize, file=paste0(subDir, outF_PoolSize),row.names = F)
        message( "        - Generating table :    ",outF_TotalLabel_FE,"\n")
        write.csv(DF_TotLabel_FE,file=paste0(subDir, outF_TotalLabel_FE),row.names = F)
        
        
        # Re-formating the MID_FE, TotalLabel_FE and  PoolSize tables ----
        #  Use the study settings for (sample removal, group exclusion, or study factor scaling)
        #  Adding a third replicates for the sample with Two group
     
        r=format_Tb_all(subDir_=subDir, prefix1__=prefix_, suffix_="_MID_FE.csv", studySettings_=studySettings, mod_="MIDA"); if (!is.data.frame(r)){message(r); return()}
        MID_FE_global<<-r
        
        
        ## Stop here if we only have one compound in the data set.
        if(length(unique(CompNameWOmid_Mida(MID_FE_global[,-c(1,2)])))<2){ message("\n      Only one compound is remaining in the dataset. The analysis ends here.");  return()}
        
        
        r=format_Tb_all(subDir_=subDir, prefix1__=prefix_, suffix_="_TotalLabel_FE.csv", studySettings_=studySettings, mod_="MIDA"); if (!is.data.frame(r)){message(r); return()}
        TotalLabel_FE_global<<-r
        
        r=format_Tb_all(subDir_=subDir, prefix1__=prefix_, suffix_="_PoolSize.csv", studySettings_=studySettings, mod_="MIDA"); if (!is.data.frame(r)){message(r); return()}
        PoolSize_global<<-r
        
       
        ## Make tables:  MID_PoolSize, Total_Labaling_PoolSize ---
        # Table4 MID Poolsize
        message("\n")
        r=make_MID_PoolSize_tb(subDir_=subDir, currPrefix_=prefix_, suffix_MID_="_MID_FE.csv", suffix_PS_="_PoolSize.csv",suffix_Out_="_MID_PoolSize.csv", mod_="MIDA" )
        MID_PoolSize_global<<-r
        
        
        # Table5 Total labeling Poolsize 
        r=make_TotalLabel_PoolSize_tb(subDir_=subDir, currPrefix_=prefix_, suffix_PS_="_PoolSize.csv", suffix_TlFE_="_TotalLabel_FE.csv",suffix_Out_="_TotalLabel_PoolSize.csv" );  if (!is.data.frame(r)){message(r); return()}       
        TotalLabel_PoolSize_global<<-r

      
        
        ## Make PRISM tables ----
        message("\n\t\tMaking table for PRISM")
        p_MID=paste0(subDir, prefix_, "_MID_FE.csv")  
        p_TPOOL=paste0(subDir, prefix_, "_PoolSize.csv") 
        
        dummy=make_PRISM_Tables(p_MID_=p_MID,  p_TPOOL_=p_TPOOL, subDir_=subDir , studySettings_=studySettings)
        
        
        # Generating first sets figures ----
        message(paste0("\n",get.prog(prg,LEN)," * STEP 7 --- Generating figures  ..."))
        # heatmap
        htmp.TL=data.matrix_bm(DF_TotLabel_FE[,3:dim(DF_TotLabel_FE)[2]])
        
        
        cix=which(apply(htmp.TL,2,max)>5)
        
        if (length(cix)>1){
                htmp.TL.ZS=apply(htmp.TL[,cix],2,function(cl){zscore(cl)}) #a
                pdf(file=paste0(subDir, outHtmp),onefile = T,height=10,width=5)
                pheatmap(t(htmp.TL.ZS),main="Total Fractional Enrichment\nsum of M+1, M+2, M+3, ...",border_color = "NA",cluster_rows = FALSE,cluster_cols = FALSE)
                dev.off()
        }
        
        
        ### Within replicates groups coefficient of variation
        #OutDf= sub_Ng_Zero(DF_MID_FE) 
        
        # Re-assigning the 'DF_MID_FE' dataframe to the final reformatted data set 'MID_FE_global'
        DF_MID_FE=MID_FE_global; rownames(DF_MID_FE)=MID_FE_global$sample
        DF_MID_PoolSize=MID_PoolSize_global; rownames(DF_MID_PoolSize)=MID_PoolSize_global$sample
        
        
        message("\t\tComputing coefficient of variation for the uncorrected data")
        InputCV=coev_df(InputData1, compNotesDF, metaCols=c(1)); write.csv(InputCV,file=paste0(subDir2,"/", "CV_uncorr.csv"))  # Uncorrected data 
        message("\t\tComputing coefficient of variation for the corrected MID FE data")
        #OutCV=coev_df(OutDf, compNotesDF, metaCols=c(1,2)); write.csv(OutCV,file=paste0(subDir2,"/", "CV_MID_FE.csv"))         # Corrected data
        OutCV=coev_df(DF_MID_FE, compNotesDF, metaCols=c(1)); write.csv(OutCV,file=paste0(subDir2,"/", "CV_MID_FE.csv")) 
       
        #make cv plot
        if (length(which(is.na(OutCV)))< (dim(OutCV)[1]*dim(OutCV)[2])){
              g1=cvPlt(InputCV,main="Raw pool size\nCoefficient of variation")
              ggsave(paste0(subDir,"CV_Raw_PoolSize.pdf"),  g1 , height=11, width=6)
              
              g1=cvPlt(OutCV,main="Corrected MID FE\nCoefficient of variation")
              ggsave(paste0(subDir,"CV_corrected_MID_FE.pdf"),  g1 , height=11, width=6)
        }
        
        
        colorz=c("white",colorRamps::primary.colors(),"grey30","grey45")[c(1,15,4,5,12,10,13,6,7,8,9,11,14,16,17,18,19,20,21,22,23,24,25,26,27,28,2,29,3)] #[sample(1:30,30)]#plot(1:30,1:30,pch=19,cex=1,col=colz)
        #plot(1:29,col=c("white",colorRamps::primary.colors(),"grey30")[c(1,16,5,6,13,11,14,  4,7,8,9,10,12,15,17,18,19,20,21,22,23,24,25,26,27,28,2,29,3)],pch=19,cex=2.5)
        names(colorz)= unlist(lapply(1:length(colorz),function(i){paste0("M+",as.character(i-1))}))
        
        #Total enrichment figures
        gL1<-gL2<-gL3<-gL4<-c()
        for (compound in compNotesDF$Compound){
                # cat(compound,"\n")
                # Total Enrichment matrix for each compound 
                cIx=which(startsWith(tolower(colnames(DF_MID_FE)),tolower(compound)))
                if (length(cIx)==0){next}
                
                mt<-data.matrix_bm(DF_MID_FE[,cIx])
                mt[which(mt< 0)]=0
                g1=ggBar(mt,colorz,main_="Total Enrichment",ylab_="Total Enrichment")
                gL1[[compound]]<-g1
                
                
                cIxPs=which(startsWith(tolower(colnames(DF_MID_PoolSize)),tolower(compound)))
                if (length(cIxPs)==0){next}
                mtPs<-data.matrix_bm(DF_MID_PoolSize[,cIxPs])
                mtPs[which(mtPs< 0)]=0
                g2=ggBar(mtPs,colorz,main_="MID PoolSize",ylab_="Total Enrichment")
                gL2[[compound]]<-g2
                
                
                # --   Averaging biological replicates   -- #
                #  Total Enrichment matrix for each compound (ComBined replicates)
                sNames.rep=getRepSetNames(rownames(mt))
                mt.rep=matrix(0,nrow=length(sNames.rep),ncol=ncol(mt),dimnames = list(sNames.rep,colnames(mt)))
                for (i in 1:length(sNames.rep)){
                        rix=which(startsWith(rownames(mt),sNames.rep[i]))
                        if (length(rix)==1){mt.rep[i,]<-mt[rix,]} else { mt.rep[i,]=colMeans(mt[rix,])}   
                }
                #Standard Error calculation
                sde=get.sde(DF_MID_FE, mt.rep)
                g3=ggBar.sde(mt.rep,sde,colorz,main_="Total Enrichment",ylab_="Total Enrichment")
                gL3[[compound]]<-g3
                
                
                #if (nrow(mt.rep)<=12){ H=4;W=4}; if (nrow(mt.rep)>12 &nrow(mt.rep)<15){ H=4;W=8}; if (nrow(mt.rep)>=15 &nrow(mt.rep)<20){ H=4;W=20}; if (nrow(mt.rep)>=20 &nrow(mt.rep)<30){ H=4;W=25};if (nrow(mt.rep)>=30){ H=4;W=30}
                ggsave(file=paste0(subDir2,"/",compound,"_TE.repset.pdf"),height=3,width=3)
                
                sNames.repPs=getRepSetNames(rownames(mtPs))
                mtPs.rep=matrix(0,nrow=length(sNames.repPs),ncol=ncol(mtPs),dimnames = list(sNames.repPs,colnames(mtPs)))
                for (i in 1:length(sNames.repPs)){
                        rix=which(startsWith(rownames(mtPs),sNames.repPs[i]))
                        if (length(rix)==1){mtPs.rep[i,]<-mtPs[rix,]} else { mtPs.rep[i,]=colMeans(mtPs[rix,])}   
                }
                #Standard Error calculation
                sde=get.sde(DF_MID_PoolSize, mtPs.rep)
                g4=ggBar.sde(mtPs.rep,sde,colorz,main_="MID PoolSize",ylab_="Total Enrichment")
                gL4[[compound]]<-g4
               
                # # Percentage enrichment matrix for each compund 
                # mtPs.prcnt.rep= t(apply(mt.rep,1,function(row){round(100*(row/sum(row)),2)})) 
                # g4=ggBar(mt.prcnt.rep,colorz,main_="Percentage Enrichment",ylab_="Percentage Enrichnmet")
                # gL4[[compound]]<-g4
                

                
        }
       
        
        message("\t\tGenerating figures for MID Fractional Enrichment  (all replicates)")
        grobs = lapply(gL1, "+", theme(plot.margin = unit(c(.8,.9,.8,.9), "cm")))
        grobs_ <- lapply(grobs, ggplotGrob)
        ggsave(paste0(subDir,"Total.Enrichment.replicates.pdf"), marrangeGrob(grobs = grobs_ , nrow=3, ncol=2), height=11, width=9)
       
        
        message("\t\tGenerating figures for MID PoolSize (all replicates)")
        grobs = lapply(gL2, "+", theme(plot.margin = unit(c(.8,.9,.8,.9), "cm")))
        grobs_ <- lapply(grobs, ggplotGrob)
        ggsave(paste0(subDir,"MID_PoolSize.replicates.pdf"), marrangeGrob(grobs = grobs_ , nrow=3, ncol=2), height=11, width=9)
        
        
        message("\t\tGenerating figures for total enrichment (combining replicates)") 
        grobs = lapply(gL3, "+", theme(plot.margin = unit(c(.8,.9,.8,.9), "cm")))
        grobs_ <- lapply(grobs, ggplotGrob)
        ggsave(paste0(subDir,"Total.Enrichment.pdf"), marrangeGrob(grobs = grobs_ , nrow=3, ncol=2), height=11, width=9)
        
        message("\t\tGenerating figures for MID PoolSize (combining replicates)") 
        grobs = lapply(gL4, "+", theme(plot.margin = unit(c(.8,.9,.8,.9), "cm")))
        grobs_ <- lapply(grobs, ggplotGrob)
        ggsave(paste0(subDir,"MID_PoolSize.pdf"), marrangeGrob(grobs = grobs_ , nrow=3, ncol=2), height=11, width=9)
        
        message("\t\tGenerating TCA cycle figure") 
        ##TCA cycle figure
        pdf(file=paste0(subDir,"TCA.pdf"),width=9,height=11)
        tcaPlt(subDir2) 
        dev.off()
        
        
        ## Statistical analysis (with MetaboAnalyst) ----
        message(paste0("\n",get.prog(prg,LEN)," * STEP 8 --- Statistical analysis (with MetaboAnalyst) ..."))
       
        
        outPathDirs <-statistical_Analyses_all (  subDirCurr=subDir, 
                                                  
                                                  prefixCurr=prefix_,
                                                  run.TypeCurr=run.Type_,
                                                  
                                                  compNotesDF=compNotesDF,
                                                  studySettingz=studySettings,  
                                                  
                                                  topMetaboFDR_=topMetaboFDR_,
                                                  normMethod_=normMethod_,
                                                  transMethod_=transMethod_,
                                                  scalingMethod_=scalingMethod_,
                                                  
                                                  KeggPathwayAnalysisCheck_=KeggPathwayAnalysisCheck_,
                                                  KeggPathwayEnrichMethod_=KeggPathwayEnrichMethod_,
                                                  KeggPathwayOrganism_=KeggPathwayOrganism_,
                                                  KeggPathwayFDR_=KeggPathwayFDR_,
                                                  
                                                  lbl_="8"
                                                )
        
        
        
        ## Step 9 : Generating the final reports ----
        
        message(paste0("\n",get.prog(prg,LEN)," * STEP 9 --- Generating pdf report  ..."))
        
      
        generateFinalReport(outPathDirz=outPathDirs, subDir_=subDir)
        
        setwd(app.Path)
        cat(app.Path)
        
        
        message(paste0("\n",get.prog(prg,LEN),"\t\t         MIDA run for ", prefix_, " is completed successfully!                "))
        
        message(paste0("\n\t\t\t\t*---      ", format(Sys.time(), "%a, %b %Y %d, %X"),"      ---*        \n\n"))
       
}

