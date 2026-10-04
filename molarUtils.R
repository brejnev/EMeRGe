#### midaUtil.R : Contains function that are specific to MOLAR and those that are share between MOLAR and MIDA ##### ----  


#### Section 1: ====== Definition Functions ======

#*********
## This code we have used some function of StatTarget for QC-base signal correction.
## Functions with names ending with "_bm" are those which Brejnev Muhire changed 
## from their original versions in StatTarget R-package.
#*********


#'Imputation using the half min value per compound per group
#'@param x input data frame 
#'@param group list of group names 
#'
minHalfValue <- function(x, group) {
        group = as.factor(as.numeric(group))
        for (i in 1:dim(x)[1]) {
                for (j in 1:dim(x)[2]) {
                        if (is.na(x[i, j]) == TRUE) {
                                x[i, j] <- tapply(as.numeric(x[, j]), group,
                                                  min, na.rm = TRUE)[group[i]]/2
                        }
                }
        }
        return(x)
}

#'Imputation using the min value per compound per group,  remove compound with all NAs in a particular group
#'@param x input data frame 
#'@param group list of group names 
#'
minHalfValue_bm<-function(x,group){
        
        message("After correction, imputation using the half min value per compound per group")
        
        gNs=unique(group)
        for (j in 1:dim(x)[2]){
                for (gN in gNs){
                        rxi=which(group==gN)
                        v=x[rxi,j]
                        
                        if (sum(is.na(v))>0) { #In missing values are present, impute
                                if (sum(is.na(v))<length(v)){
                                        halfMin=min(v,na.rm=T)/2
                                        v[is.na(v)]=halfMin
                                        x[rxi,j]=v
                                } else{
                                        message(paste0("All missing values!:\n compound: ",
                                                       colnames(x)[j]," \n samples: ",
                                                       paste(as.character(rownames(x)[rxi]),collapse=", "),
                                                       sep=""))
                                }
                        }
                }
        }
        na.id=which(colSums(is.na(x)) > 0)
        
        message(paste0("Removed compounds = ", length(na.id),
                       " \n samples: ",paste(as.character(colnames(x)[na.id]),collapse="\n"),
                       sep=""))
        if (length(na.id)>0){
                x=x[,-na.id]
        }
        return(x)
}

#'Imputation using the min value per compound per group
#'@param x input data frame 
#'@param group list of group names 
#'
minValue <- function(x, group) {
        group = as.factor(as.numeric(group))
        for (i in 1:dim(x)[1]) {
                for (j in 1:dim(x)[2]) {
                        if (is.na(x[i, j]) == TRUE) {
                                x[i, j] <- tapply(as.numeric(x[, j]), group,
                                                  min, na.rm = TRUE)[group[i]]
                        }
                }
        }
        return(x)
}

#'Imputation using the median value per compound per group
#'@param x input data frame 
#'@param group list of group names 
#'
medianvalue <- function(x, group) {
        group = as.factor(as.numeric(group))
        for (i in 1:dim(x)[1]) {
                for (j in 1:dim(x)[2]) {
                        if (is.na(x[i, j]) == TRUE) {
                                x[i, j] <- tapply(as.numeric(x[, j]), group,
                                                  median, na.rm = TRUE)[group[i]]
                        }
                }
        }
        return(x)
}

#'Takes the compound table and filter out all compounds with more than Frule (eg. 0.2) of
#'Missing data in any of classes (sample group)
#'@param m input data frame 
#'@param degree factor or list of sample group names
#'@param Frule_ filtration rule threshold (maximum proportional of missing values in a particular group)
#'
FilterMV = function(m, degree,Frule_=0.2) {
     
        #--- This is for testing ---
        # m=imsamFP; degree=classF
        
        dx <- c()
        for (i in 1:ncol(m)) {
                freq <- as.vector(tapply(m[, i], degree, function(x) {
                        sum(is.na(x))/length(x)
                }))
                if (sum(freq > 1 - Frule_) > 0)
                        dx <- c(dx, i)
        }
        
        if (length(dx) > 0)
                m <- m[, -dx]
        return(m)
}

#'This filters compounds based on comparison of QCs vs Blanks (rule: QCs must me 5x > than blanks)
#'@param df input data frame
#'@param cuttoff Threshold for QCs vs blanks i.e. if cuttoff=5 it mean QCs 5x > blanks must be cuttoff X
#'
Filter_qc_blnk<-function(df,cutoff=5 ){
       
        mt=data.matrix(df[,3:dim(df)[2]])
        ncmp=dim(mt)[2]
        
        qc.i=grep("qc",tolower(rownames(mt)))
        bl.i=grep("mc-pb",tolower(rownames(mt)))
        
        if (length(qc.i)==0){message("No QC samples found, Filtration based on QC and Blank cannot proceed!") }
        if (length(bl.i)==0){message("No blank samples found, Filtration based on QC and Blank cannot proceed!")}
        
        
        if (length(qc.i)==1){qc.means=mt[qc.i,]} else { qc.means=colMeans(mt[qc.i,],na.rm=T)}
        if (length(bl.i)==1){bl.means=mt[bl.i,]}else{bl.means=colMeans(mt[bl.i,],na.rm=T)}
        
        
        ixc=which((qc.means/bl.means) > cutoff)
        mt1=mt[,ixc]
        if (ncmp > dim(mt1)[2]){
                message("     QC/Blank-based-filtration (threshold =",cutoff,"): ",ncmp-dim(mt1)[2], " compounds removed")
                message("     ",paste(colnames(mt[,-ixc]),collabse="\n     "))
        }else{
                message("     QC/Blank-based filtration (threshold =",cutoff,"): no compounds removed")
        }
        rt=cbind(df[,1:2],mt1); colnames(rt)[1:2]=colnames(df)[1:2]
        return(rt)
}

#'Generating the figures : The plot of RSD distribution and data table.
#'Loplot, a visible figure of signal correction loplot provide the visible figure of QC-based
#'correction.  x The file before QC-based correction.  z The file after QC-based correction.  i
#'The name of variable.  The plot of RSD distribution and data table.
#'@param x input data frame before QC-based correction
#'@param z  the data after QC-based correction
#'@param ii the name of variable

#loplot(datmp, loessDatT, i)

loplot<-function (x, z, ii) {
        
        #Line for testing
             #x=datmp; z=loessDatT; ii=i
        
        cn <- colnames(x)
        qcid <- grep("qc", tolower(cn))
        meantem <- mean(z[ii, qcid])
        sdtem <- sd(z[ii, qcid])
        rsdcutoffU <- 0.3
        rsdcutoffL <- 0.15
        sdLineU <- rsdcutoffU * meantem
        sdLineL <- rsdcutoffL * meantem
        cuttemU <- meantem + sdLineU
        cuttemD <- meantem - sdLineU
        cuttemUL <- meantem + sdLineL
        cuttemDL <- meantem - sdLineL
        meantemR <- mean(x[ii, qcid])
        sdtemR <- sd(x[ii, qcid])
        rsdcutoffRU <- 0.3
        rsdcutoffRL <- 0.15
        sdLineRU <- rsdcutoffRU * meantemR
        sdLineRL <- rsdcutoffRL * meantemR
        cuttemRU <- meantemR + sdLineRU
        cuttemRD <- meantemR - sdLineRU
        cuttemRUL <- meantemR + sdLineRL
        cuttemRDL <- meantemR - sdLineRL
        RSD30_CV = paste0(rownames(x)[ii], "_", ii, ".pdf")
        dirout.loplot <- paste0(getwd(), "/loplot")
        if (!file.exists(dirout.loplot)) {dir.create(dirout.loplot)}#, showWarnings = FALSE)}
        
        pdf(paste(dirout.loplot, RSD30_CV, sep = "/"), width = 6,height = 7)
        
        #graphics::layout(matrix(1:2, nrow = 2))
        par(mfrow=c(2,1),mar=c(4,4,4,5))
        numY <- 1:dim(x)[2]
        # graphics::plot(numY, x[ii, ], pch = 19, col = "green",
        plot(numY, x[ii, ], pch = 19, col = "green",
             ylab = c("Intensity"), xlab = c("Injection Order"),
             main = paste0(rownames(x)[ii]," - Raw peaks"))
        points(qcid, x[ii, qcid], pch = 19, col = "red")
        abline(h = cuttemRU, col = "#D55E00", lwd = 1.2, lty = 1)
        abline(h = cuttemRD, col = "#D55E00", lwd = 1.2, lty = 1)
        abline(h = cuttemRUL, col = rgb(0, 0, 0, 0.3), lwd = 1.2,lty = 2)
        abline(h = cuttemRDL, col = rgb(0, 0, 0, 0.3), lwd = 1.2,lty = 2)
        mtext("+SD@30%CV", side = 4, line = 0.3, at = cuttemRU,col = "#D55E00", cex = 0.5, las = 1)
        mtext("-SD@30%CV", side = 4, line = 0.3, at = cuttemRD,col = "#D55E00", cex = 0.5, las = 1)
        mtext("+SD@15%CV", side = 4, line = 0.3, at = cuttemRUL,col = "#D55E00", cex = 0.5, las = 1)
        mtext("-SD@15%CV", side = 4, line = 0.3, at = cuttemRDL,col = "#D55E00", cex = 0.5, las = 1)
        legend("top", c("Sample", "QC"), col = c("green","red"), lty = 1, pch = 19, bty = "n", cex = 0.75, horiz = TRUE)
        
        #graphics::plot(numY, z[ii, ], pch = 19, col = "green",
        plot(numY, z[ii, ], pch = 19, col = "green",
             ylab = c("Intensity"), xlab = c("Injection Order"),
             main = paste0(rownames(x)[ii]," - Corrected peaks"))
        points(qcid, z[ii, qcid], pch = 19, col = "red")
        abline(h = cuttemUL, col = rgb(0, 0, 0, 0.3), lwd = 1.2,lty = 2)
        abline(h = cuttemDL, col = rgb(0, 0, 0, 0.3), lwd = 1.2, lty = 2)
        abline(h = cuttemU, col = "#D55E00", lwd = 1.2, lty = 1)
        abline(h = cuttemD, col = "#D55E00", lwd = 1.2, lty = 1)
        mtext("+SD@30%CV", side = 4, line = 0.3, at = cuttemU,col = "#D55E00", cex = 0.5, las = 1)
        mtext("-SD@30%CV", side = 4, line = 0.3, at = cuttemD,col = "#D55E00", cex = 0.5, las = 1)
        mtext("+SD@15%CV", side = 4, line = 0.3, at = cuttemUL,col = "#D55E00", cex = 0.5, las = 1)
        mtext("-SD@15%CV", side = 4, line = 0.3, at = cuttemDL,col = "#D55E00", cex = 0.5, las = 1)
        legend("top", c("Sample", "QC"), col = c("green","red"), lty = 1, pch = 19, bty = "n", cex = 0.75,horiz = TRUE)
        dev.off()
}


## Generat a PCA plot, given a data frame and return a ggplot2 object
#'@param DF input data frame
#'@param fig.tittle figure main tittle
#'
pltPCA<-function(DF,fig.title){
        
        DF=DF[-grep("qc",tolower(DF$Replicate)),]
        
        M=DF[,9:dim(DF)[2]];rownames(DF)=DF$Replicate
        
        M=cbind(M,Group=DF$GroupName)
        
        PCA.res=PCA(M[,-dim(M)[2]], graph = FALSE)
        
        colorz=c("black","red","blue","green","darkcyan","chocolate4","brown3","chartreuse4","darkmagenta","darkslateblue","dodgerblue4",
                 "gray30","deeppink","goldenrod4","darkblue","brown4","gray55","darkorange","deepskyblue","gold3","gray35","gray65","gray65","gray80")[1:dim(M)[1]]
        g<-fviz_pca_ind(PCA.res,
                        geom.ind = "point", # show points only (nbut not "text")
                        col.ind = M$Group, # color by groups
                        palette = colorz, #c("#00AFBB", "#E7B800", "#FC4E07"),
                        addEllipses = FALSE, # Concentration ellipses
                        legend.title = "Groups",
                        title=fig.title
        )
        return(g)
}


#'Getting the RSD, CV and other basic statistics
#' it returns : a two-column data frame, where the first column takes the value of the statistics, and
# the second its name, e.g.  'nobs', 'NAs', 'Minimum', 'Maximum', '1. Quartile', '3.  Quartile',
# 'Mean', 'Median', 'Sum', 'SE Mean', 'LCL Mean', 'UCL Mean', 'Variance', 'Stdev', 'Skewness',
# 'Kurtosis')
#'@param x input data frame
#'@param ci  a numeric value setting the confidence interval.
#'
bStatCor_bm<-function (x, ci = 0.95)
{
        ### Getting the RSD, CV and other statistics
        
        ## For testing
          #x=I
        
        y = as.matrix(x)
        
        ## If there no column names (compound names), use x1, x2, ..., xn to name the columns
        if (is.null(colnames(y))) {
                Dim = dim(y)[2]
                if (Dim == 1) {
                        colnames(y) = paste(substitute(x), collapse = ".")
                }
                else if (Dim > 1) {
                        colnames(y) = paste(paste(substitute(x), collapse = ""),
                                            1:Dim, sep = "")
                }
        }
        
        ## Computing the confidence interval (around the mean)
        cl.vals = function(x, ci) {
                x = x[!is.na(x)]
                n = length(x)
                if (n <= 1) {return(c(NA, NA))}    ## if you only have one compound, return c(NA, NA)
                se.mean = sd(x)/n#sqrt(stats::var(x)/n)    ## standard error of the mean; computer the covariance matrix (betwen colum i,e compound)  -> sqrt(cov(x)/n),   cov= 1/n-1 [sum( (x-mean(x))*y-mean(y))]
                t.val = qt((1 - ci)/2, n - 1)
                mn = mean(x)
                lcl = mn + se.mean * t.val
                ucl = mn - se.mean * t.val
                c(lcl, ucl)
        }
        
        nColumns = dim(y)[2]
        ans = NULL
        for (i in 1:nColumns) {
                X = as.numeric(y[, i])
                X.length = length(X)
                X = X[!is.na(X)]
                X.na = X.length - length(X)
                ## Computing different statistics for each comppound
                z = c(X.length, X.na, min(X), max(X),
                      as.numeric(quantile(X,prob = 0.25, na.rm = TRUE)),
                      as.numeric(quantile(X,prob = 0.75, na.rm = TRUE)),
                      mean(X),
                      median(X),
                      sum(X),
                      sd(X)/length(X),#sqrt(stats::var(X)/length(X)),
                      cl.vals(X,ci)[1], cl.vals(X, ci)[2],
                      var(X),#stats::var(X),
                      sd(X),#sqrt(stats::var(X)),
                      sd(X)/mean(X))#sqrt(stats::var(X))/mean(X))
                
                znames = c("nobs", "NAs", "Minimum", "Maximum", "1. Quartile",
                           "3. Quartile", "Mean", "Median", "Sum", "SE Mean",
                           "LCL Mean", "UCL Mean", "Variance", "Stdev", "RSD")
                
                result = matrix(z, ncol = 1)
                row.names(result) = znames
                ans = cbind(ans, result)
        }
        colnames(ans) = colnames(y)
        data.frame(round(ans, digits = 6))
}

#Calculate and returns an RSD values for each compound accross all samples
RsdCal_bm<-function (tb, group = TRUE, DistPattern = TRUE, output = TRUE,smpl=NULL,corr=NULL)
{
        ## Calculate and returns an RSD values for each compound accross all samples
        ## If group=TRUE it produces the cumulative distribution of compound RSDs for each group
        ## If DistPattern=It print the cumulative distribution matrix
        ## If output = TRUE, it save the cumulative distribution data, of Samples RSD
        ## smpl = specifies whether computing RSD for samples or QC
        ## corr = specifies whether before or after correction
        ## if(is.null(smpl)|is.null(corr)){stop("Error: RsdCal_bm -> Please specify the smpl and/or corr parameters ")}
        
        # ----- lines for testing
        # tb=raw_temp; group = TRUE; DistPattern = TRUE; output = TRUE; smpl="AllGroups"; corr="before.corr"
        # tb=lo_temp; group = TRUE; DistPattern = TRUE; output = TRUE; smpl="AllGroups"; corr="before.corr"
        
        
        #xr = tb
        #xs = tb[, 9:ncol(tb)] ""
        xs = tb[, 5:ncol(tb)] ##getting the data without GroupNumber or GroupName
        
        ## RSD per groupes
        if (group) {
                #x = cbind(xr[, 4],xr[, 3], xr[, 1], xs)      # swap the first two columns -> (group, sample names)
                #x.nn =  cbind(tb[,4],tb[, 3], xs)
                
                x.nn =  cbind(tb[,4],tb[,3], xs)
                
                colnames(x.nn)[1:2]=c("GroupNumber","GroupName")
                
                ## Give numerical names to the QCs, if QCs are presents (and represetend by NAs)
                #ixNA=which(is.na(x.nn[,"GroupNumber"]))
                #if (length(ixNA)>0){x.nn[ixNA,"GroupNumber"]=max(x.nn[,"GroupNumber"],na.rm=T)+1}
                #sorted = x.nn[order(x.nn[, 1]), ]    # sort by group numbers
                
                ## write each group intestinty values into csv, save them into tmp directory
                dirout.g = paste(getwd(), "/tmp", sep = "")
                if (!file.exists(dirout.g)){dir.create(dirout.g)}
                #g=unique(sorted$GroupNumber)
                gs=unique(x.nn$GroupName) ## gs : group names
                for (g in gs){
                        #rix=which(sorted$GroupNumber==ij)
                        rix=which(x.nn$GroupName==g)
                        #outputfileg = paste("r.", ij, ".csv", sep = "")
                        outputfileg = paste("rsd.", g, ".csv", sep = "")
                        #write.csv( sorted[rix,3:ncol(sorted)], paste(dirout.g, outputfileg, sep = "/"), row.names = FALSE)
                        write.csv( x.nn[rix,3:ncol(x.nn)], paste(dirout.g, outputfileg, sep = "/"), row.names = FALSE)
                }
                
                ### Compute RSD accros samples or QCs (in each group), for every compound
                dirout.w = paste(getwd(), "/RSDresult",sep = "")
                if (!file.exists(dirout.w)){dir.create(dirout.w)}
                
                NoF = length(gs)
                #nrow = dim(sorted)[2] - 2
                nrow = dim(x.nn)[2] - 2
                
                groupRSD <- matrix(NA, ncol = NoF, nrow = nrow)
                colnames(groupRSD)<-gs
                
                for (i in 1:NoF) {
                        #ni = paste("r.", i, ".csv", sep = "")
                        ni = paste("rsd.", gs[i], ".csv", sep = "")
                        pwdi = paste(getwd(), "/tmp/", ni, sep = "")
                        I = read.csv(pwdi, header = TRUE)
                        #I = I[, -1]
                        bS = t(bStatCor_bm(I)[nrow(bStatCor_bm(I)), ])  ### Get only the Statistics of interest, (Only RSD is retained at this point)
                        bS = as.matrix(bS)                              ### Return dim = no. of compounds
                        groupRSD[, gs[i]] = bS
                }                                                 ### Return dim = no. of groupes x no. compounds
                
                ### bind columns of  the RSD for all groupes and RSD computed with all samples
                bSall = t(bStatCor_bm(xs)[nrow(bStatCor_bm(xs)), ])
                groupRSD <- cbind(groupRSD, bSall)
                
                #Name the columns of
                #colnames(groupRSD)=c(unlist(lapply(g,function(ix){unique(sorted$GroupName[which(sorted$GroupNumber==ix)])})),"All/Total")
                colnames(groupRSD)[dim(groupRSD)[2]]="All/Total"
                
                RSD_all = paste0("RSD_all_",corr, ".csv")
                write.csv(round(t(groupRSD),2), paste(dirout.w, RSD_all, sep = "/"))
                
                #rownames(groupRSD) <- NULL
                #colnames(groupRSD) <- NULL
                
                # Computes the cumulative distribution of compound RSD values for each group
                #distout <- matrix(NA, nrow = NoF + 1, ncol = 20)
                distout <- matrix(NA, nrow = dim(groupRSD)[2], ncol = 20)
                rownames(distout)<-colnames(groupRSD)
                
                #for (i in 1:(NoF + 1)) {
                for (cl in colnames(groupRSD)) {
                        #valueB <- groupRSD[, i]
                        valueB <- groupRSD[, cl]
                        groupdist <- sapply(seq(0, 95, 5), function(x) {                  #propotion of values
                                sum(valueB * 100 <= x + 5, na.rm = TRUE) * 100/length(valueB)
                        })
                        #plot(seq(0, 95, 5), groupdist,ylab="Freq (x100)",xlab="RSD*100",main=paste0Cumulative Distribution of RSD in group ",i))
                        #distout[i, ] <- matrix(groupdist, 1)
                        distout[cl, ] <- matrix(groupdist, 1)
                }
                
                #View & save data of cumulative distribution of samples or QCs RSD
                colnames(distout) <- paste("CV", "<", seq(5, 100, 5),"%", sep = "")
                #rownames(distout) <-c(unlist(lapply(g,function(ix){unique(sorted$GroupName[which(sorted$GroupNumber==ix)])})),"All/Total")
                
                if (DistPattern) {
                        print(distout)
                }
                
                if (output) {
                        #table
                        RSDdist_CV=paste0("RSDdist_",smpl,"_",corr,".csv");
                        write.csv(round(distout,2), paste(dirout.w, RSDdist_CV, sep = "/"))
                        
                        #figure
                        fig=paste0("RSDdist_",smpl,"_",corr,".pdf");
                        pdf(paste(getwd(),fig,sep="/"),width=8,height=5)
                        colz=colorRamps::primary.colors()
                        plot(seq(0, 95, 5),distout[1,],pch=19,ylim=c(0,100),xlim=c(0,100),col=colz[1],ylab="Freq x100 (No. compounds % < RSD )",xlab="RSD*100",main=paste0("Cumulative Distribution of RSD\n",smpl,"  (" ,corr,")"))
                        lines(seq(0, 95, 5),distout[1,],col=colz[1])
                        for (f in 2:dim(distout)[1]){
                                points(seq(0, 95, 5),distout[f,],pch=19,col=colz[f])
                                lines(seq(0, 95, 5),distout[f,],col=colz[f])
                        }
                        AUCs=sapply(1:nrow(distout), function(vi){round(sum(distout[vi,])/ (length(distout[vi,])*100),3)})
                        legend("bottomright",legend=paste0(rownames(distout),", AUC:",AUCs), pch=19,col=colz[1:nrow(distout)])
                        dev.off()
                }
                tmpfile = paste(getwd(), "/tmp/", sep = "")
                unlink(tmpfile, recursive = TRUE)
        }
        ## Just compute the RSD values for all samples
        else {
                stop("Please specify whether set group=TRUE")
                #bSall = t(bStatCor_bm(xs)[nrow(bStatCor_bm(xs)), ])
        }
        # Save RSD values for for all samples
        # RSD_all = paste0("RSD_all_",corr, ".csv")
        # write.csv(round(bSall,2), paste(dirout.w, RSD_all, sep = "/"))
        
        #return(bSall)   # return dim = no. compounds
}

#'QC - based correction
#'@param metaDF input metadata frame
#'@param samDF input data frame
#'@param Frule=  filtration rule threshold (maximum proportional of missing values in a particular group) 
#'@param imputeM imputation method ("half of the minimum","minimum", "median" or  "KNN")
#'@param MLmethod= machine learning method for QC-base signal correction ("Loess Regression","random forest classifier)
#'@param ntree number of tree for a random-forest classifier (default = 500) 
#'@param QCspan  smoothing parameter (default=0) 
#'@param loess.degree polynomial degree (default=2)
#'@param plot determines whether to generate the plots or not
#'@param signalCorr.Flg determines whether to perform or skip the QC-based signal correction
#'@param qc.blank.Flg determines whether to perform or skip qc/blank-based filtration
#'@param NA.thr.Flg threshold of proportion of missing value at which compound should be filtered out
#'@param imputation.Flg determines whether to perform or skip imputation
#'

shiftCor_bm<-function(metaDF,
                      samDF,
                      Frule=0.8, 
                      imputeM = "KNN",
                      MLmethod="QCRFSC", 
                      ntree = 500, 
                      QCspan = 0, 
                      loess.degree = 2,
                      plot = TRUE,
                      
                      signalCorr.Flg,
                      NA.thr.Flg,
                      imputation.Flg
)
{
        ## QC - based correction
        cat("\n")
        
        
        message("\n** Step 1: Data files checking starts...")
        
        ## --------  This section is for testing and must be commented afterwards --------##
        # metaDF=metaData2; samDF=InputData2; Frule = 1-NA.thr_; imputeM = imput.method_ ; MLmethod = corr.method_; ntree = 500;  QCspan = QC.span_; loess.degree = 2; plot = TRUE
        # signalCorr.Flg=signalCorr.Flg_; NA.thr.Flg=NA.thr.Flg_; imputation.Flg=imputation.Flg
        
        if (!identical(metaDF[,"Replicate"],samDF[,"Sample"])){stop("Error: Sample names in data and metadata do not match.")}
        
        ### Check the number of samples in both two input file (data vs metadata)
        message("    - ", dim(metaDF)[1], " Metadata samples  and ", dim(samDF)[1], " data samples", sep = "")
        message("    - The Meta samples list (*NA, missing data from the Profile File)")
        
        
        ### Match the sample names in both data and metadata dataframe
        #mcdat <- samDF[, 1][match(metaDF[, 2], samDF[, 1])]
        mcdat <- samDF[, "Sample"][match(metaDF[, "Replicate"], samDF[, "Sample"])]
        print(as.vector(mcdat))
        
        
        if (any(is.na(mcdat))) {
                stop("\n", "   - Missing data from the Profile File! Check your data please!!")
        }
        if (dim(samDF)[1] - dim(metaDF)[1] > 0) {
                cat("\n", "   - Warning: The sample size in Profile File is larger than Meta File! ")
                cat("\n")
        }
        if (dim(samDF)[1] - dim(metaDF)[1] < 0) {
                stop("\n", "   - The sample size in Profile File should be no less than Meta File!",
                     " Check your data please!!")
        }
        
        ### Check whether we have QC
        #samFP <- samDF[metaDF$sample, ]
        
        samFP <- samDF[metaDF$Replicate, ]
        
        if (sum(is.na(metaDF$GroupNum)) <= 0) {
                stop("\n", "   - There were no QC samples in your data!")
        }
        
        ### Check the classes (sample group) and batches
        metaDF_stat <- metaDF
        #rownames(metaDF_stat) <- metaDF[, 1] #  *** ! **
        qc_seq <- rownames(metaDF_stat)
        qc_seq_tmp <- grep("qc", tolower(qc_seq))
        
        sam_seq_tmp <- metaDF_stat[-c(qc_seq_tmp), ]
        
        # if (sum(is.na(metaDF_stat$batch)) > 0) { *** ! **
        #         stop("\n", "There were missing values in batch! Check your data please!\n")
        # }
        #if (sum(is.na(sam_seq_tmp$class)) > 0) {
        
        if (sum(is.na(sam_seq_tmp$GroupNum)) > 0) {
                stop("\n", "   - There were missing values (Unclassified data) in sample Groups Num! \n         Check your data please!\n")
        }
        
        ### Print class(sample group) and batch summaries
        #metaDF_stat$class[is.na(metaDF_stat$class)] <- "QC" *** ! **
        metaDF_stat$GroupNum[is.na(metaDF_stat$GroupNum)] <- "QC"
        
        #data_stat = aggregate(metaDF_stat$class, by = list(Category = metaDF_stat$class),FUN = length)  *** ! **
        #data_stat = aggregate(metaDF_stat$GroupNum, by = list(Category = metaDF_stat$GroupNum),FUN = length)
        data_stat = aggregate(metaDF_stat$GroupName, by = list(Category = metaDF_stat$GroupName),FUN = length)
        
        #batch_stat = aggregate(metaDF_stat$class, by = list(Category = metaDF_stat$batch),FUN = length)
        
        colnames(data_stat) <- c("Group", "No.")
        cat("\n", "Meta-information:", "\n\n")
        print(data_stat)
        #colnames(batch_stat) <- c("Batch", "No.")
        #print(batch_stat)
        
        
        ### Print number of samples, QC and Metabolites
        num_sam <- dim(samDF[, 3:ncol(samDF)])
        num_sam <- data.frame(num_sam)
        colnames(num_sam) <- c("no.")
        rownames(num_sam) <- c("QC and samples", "Metabolites")
        cat("\n", "   - Metabolic profile information:", "\n\n")
        print(num_sam)
        
        
        ### Create the matrix without metadata for processing (row="samples" and colnms="metabolites")
        samFP_temp <- as.matrix(samFP[, 3:ncol(samFP)])
        
        ### Replace negative by zeros, then zeros by NA
        if (length(samFP_temp[samFP_temp < 0L]) > 0) {
                samFP_temp[samFP_temp < 0L] <- 0L
        }
        samFP_temp[samFP_temp == 0L] <- NA
        
        ### Imputation of Missing
        message("\n** Step 2: Evaluation of missing values ...")
        message("   - The number of missing values before QC based signal correction: ", sum(is.na(samFP_temp)))
        imsamFP <- samFP_temp
        
        ## Filtering compounds based on Frule
        #classF <- as.factor(metaDF$class)   ### class = sample group    *** ! **
        classF <- as.factor(metaDF$GroupNum)   ### class = sample group
        
        classF = addNA(classF)
        
        
        if (NA.thr.Flg==FALSE){
                message("   - Per group missing value filtration skipped ")
                imsamFPF=imsamFP
        }else if(NA.thr.Flg==TRUE){
                message("   - Per group missing value filtration starts ...")
                imsamFPF = FilterMV(imsamFP, degree=classF,Frule_=Frule)  
                
                Frule_warning = paste("   - The number of filtered compounds based on the  ", (1-Frule) * 100, "% rule :", sep = " ")
                message(Frule_warning, "     ", dim(imsamFP)[2] - dim(imsamFPF)[2])
        }
        
        imsamFP = as.matrix(imsamFPF)
        
        ### Imputation using KNN, minimum, minimum/2, median
        
        ### KNN (k, # of nearest neigbor default is usually 10) therefore is not suitable for dataset with few
        ### samples per group
        
        
        if (imputation.Flg==FALSE){
                message ( "\n** Step 3: Imputation is skipped")
                inputedData=imsamFP
        }else{
                message("\n** Step 3: Inputation starts ...")
                if (imputeM == "KNN") {
                        message("   - The imputation method was set at 'KNN'")
                        mvd <- impute::impute.knn(imsamFP, k = 3, rowmax = 0.99, colmax = 0.99,maxp = 15000)
                        inputedData <- mvd$data
                }
                
                ## Imputation using the min value per compound per group
                else if (imputeM == "min") {
                        message("   - The imputation method was set at 'min'")
                        inputedData = minValue(imsamFP, classF)
                }
                
                ## Imputation using the half min value per compound per group
                else if (imputeM == "minHalf") {
                        message("   - The imputation method was set at 'minHalf'")
                        inputedData = minHalfValue(imsamFP, classF)
                }
                
                ## Imputation using the median value per compound per group
                else if (imputeM == "median") {
                        message("  - The imputation method was set at 'median'")
                        inputedData = medianvalue(imsamFP, classF)
                }
                
                message("   - The number of missing value after imputation: ",sum(is.na(inputedData)))
                if (sum(is.na(inputedData)) > 0) {
                        inputedData = minHalfValue(inputedData, classF)
                        message("   The number of missing value after the second imputation: ",
                                sum(is.na(inputedData)))
                }
                message("   - Imputation completed!")
        }
        
        
        if (signalCorr.Flg==FALSE){
                message("\n** Step 4: QC - based signal correction skipped ..." )
                raw_temp<- cbind(metaDF, inputedData)
                return(raw_temp)
        }
        else if (signalCorr.Flg==TRUE){
                ### QC-based correction
                message("\n** Step 4: QC-based signal correction ... (with StatTarget)")
                message("   - Signal correction method : ",MLmethod)
                
                dat <- as.matrix(t(inputedData))
                numX <- 1:dim(dat)[2]
                
                
                ### QC-based signal correction - random forest (QC-RFSC)
                ### The correction factor of each feature (compound) is obtained by learning an ensemble of regression trees
                ### each grown from randomly sub-sampled QCs samples
                
                
                if (MLmethod == "QCRFSC") {
                        
                        REGfit = function(x, y, ntree = ntree) {
                                ## ---- For testing ------
                                # x = dat; y = numX; ntree = 500
                                
                                ## Check whether the first and the last samples are the QC samples.
                                cn <- colnames(x)
                                x <- as.matrix(x)
                                st_QC <- grep("qc", tolower(cn[1]))
                                ed_QC <- grep("qc", tolower(cn[length(cn)]))
                                if (length(st_QC) == 0) {
                                        stop("\n   Wrong: the first sample must be QC sample; please check ......")
                                }
                                if (length(ed_QC) == 0) {
                                        stop("\n   Wrong: the sample at the end of sequence must be QC sample; \n          please check ......")
                                }
                                
                                qcid <- grep("qc", tolower(cn))
                                pb <- txtProgressBar(min = 1, max = dim(x)[1], style = 3)
                                
                                ## For each compound  (rows=compounds, cols=samples), call randomForest only on QCs samples
                                for (i in 1:dim(x)[1]) {
                                        temp <- randomForest(data.frame(qcid), as.numeric(x[i,qcid]), ntree = ntree)
                                        y <- data.frame(y)
                                        colnames(y) <- "qcid"
                                        rfP <- predict(temp, y)
                                        
                                        debugLINES=FALSE
                                        if (debugLINES){
                                                #Line for testing
                                                par(mfrow=c(2,2));
                                                plot(qcid,x[i, qcid],main=paste0(rownames(x)[i]," - QC peaks"),ylab="Signal",col="red");
                                                cols=rep("black",length(rfP)); cols[qcid]="red" ;
                                                plot(x[i,],main="Raw values",col=cols ,ylab="Signal",xlab="Samples" )
                                                plot(rfP,main="Correction factors",col=cols,ylab="Value",xlab="Samples" )
                                                plot(as.numeric(x[i, ])/rfP,main="Corrected values",col=cols ,ylab="Signal",xlab="Samples" )
                                        }
                                        
                                        ### Infer the correction parameter for the compound at every sample (i)
                                        x[i, ] <- as.numeric(x[i, ])/rfP  ### Using the vector correctional factors to correct the intensities (i.e normalization)
                                        setTxtProgressBar(pb, i)
                                }
                                close(pb)
                                loessDat = x
                                return(loessDat)
                        }
                        loessDat <- REGfit(x = dat, y = numX, ntree = ntree)
                }
                
                
                if (MLmethod == "QCRLSC") {
                        
                        #NOTE: QC span default is .75, to used the AutoFit function (Optimitation of QCspan parameter) atlead 10 QC sample are requrired
                        
                        if (QCspan > 0) {
                                
                                ###  QC-base Robust Loess Signal Correction (QC-RLSC)
                                ### QC-based signal correction - nonlinear local polynomial regression LOESS (QC-RLSC)
                                # locally estimated scatterplot smoothing (LOESS)
                                
                                loessFit = function(x, y, QCspan, degree) {
                                        #Line for test
                                        #x = dat; y = numX; QCspan=0.75; degree=2
                                        
                                        cn <- colnames(x)
                                        st_QC <- grep("qc", tolower(cn[1]))
                                        ed_QC <- grep("qc", tolower(cn[length(cn)]))
                                        if (length(st_QC) == 0) {
                                                stop("   - Input error: the first sample must be QC sample; please check ......")
                                        }
                                        if (length(ed_QC) == 0) {
                                                stop("   - Input error: the sample at the end of sequence must be QC sample; \n          please check ......")
                                        }
                                        qcid <- grep("qc", tolower(cn))
                                        pb <- txtProgressBar(min = 1, max = dim(x)[1],
                                                             style = 3)
                                        for (i in 1:dim(x)[1]) {
                                                loe <- stats::loess(x[i, qcid] ~ qcid, span = QCspan,
                                                                    degree =degree)
                                                yf <- stats::predict(loe, y)
                                                
                                                #Lines for testing
                                                debugLINES=FALSE
                                                if (debugLINES){
                                                        
                                                        par(mfrow=c(2,2));
                                                        plot(qcid,x[i, qcid],main=paste0(rownames(x)[i]," - QC peaks"),col="red",ylab="Signal");
                                                        cols=rep("black",length(yf)); cols[qcid]="red" ;
                                                        plot(x[i,],main="Raw values",col=cols ,ylab="Signal",xlab="Samples" )
                                                        plot(yf,main="Correction factors",col=cols,ylab="Value",xlab="Samples" )
                                                        plot(as.numeric(x[i, ])/yf,main="Corrected values",col=cols ,ylab="Signal",xlab="Samples" )
                                                }
                                                x[i, ] <- as.numeric(x[i, ])/yf
                                                
                                                setTxtProgressBar(pb, i)
                                        }
                                        close(pb)
                                        loessDat = x
                                        return(loessDat)
                                }
                                loessDat <- loessFit(x = dat, y = numX, QCspan = QCspan,
                                                     degree = loess.degree)
                        }
                        else if (QCspan <= 0) {
                                cat("\n", "Warning: The QCspan was set at '0'.")
                                message("\n   - The GCV was used to avoid overfitting the observed data")
                                
                                ### QC-based signal correction - Generalized Cross-Validation
                                # This uses LOESS approach described above and finds optimal value
                                # of the smotthing parameter (span) between (0.5 and 0.75, step 0.01 )
                                autoFit <- function(xl, y) {
                                        #Line for test
                                        #xl = dat; y = numX
                                        
                                        cn <- colnames(xl)
                                        st_QC <- grep("QC", tolower(cn[1]))
                                        ed_QC <- grep("QC", tolower(cn[length(cn)]))
                                        if (length(st_QC) == 0) {
                                                stop("   - the first sample must be QC sample; please check ......")
                                        }
                                        if (length(ed_QC) == 0) {
                                                stop("   - the sample at the end of sequence must be QC sample; \n          please check ......")
                                        }
                                        qcid <- grep("qc", tolower(cn))
                                        pb <- txtProgressBar(min = 1, max = dim(xl)[1],
                                                             style = 3)
                                        for (i in 1:dim(xl)[1]) {
                                                Sys.sleep(1e-06)
                                                loe1 <- loess(xl[i, qcid] ~ qcid)
                                                env <- environment()
                                                sploe <- function(sp) {
                                                        loe2 <- get("loe1", envir = env)
                                                        mod <- stats::update(loe2, span = sp)
                                                        CVspan = statTarget:::loessGCV(mod)[["gcv"]]
                                                }
                                                sp <- c(seq(0.5, 0.75, 0.01))
                                                CVspan = as.matrix(lapply(sp, sploe))
                                                CVspan[!is.finite(as.numeric(CVspan))] <- NA
                                                minG <- data.frame(sp, CVspan)
                                                minspan <- minG[which.min(minG[, 2]), 1]
                                                minspan
                                                loeN <- stats::update(loe1, span = minspan)
                                                yf <- predict(loeN, y)
                                                xl[i, ] <- as.numeric(xl[i, ])/yf
                                                setTxtProgressBar(pb, i)
                                        }
                                        close(pb)
                                        loessDat = xl
                                        return(loessDat)
                                }
                                loessDat <- autoFit(xl = dat, y = numX)
                        }
                }
                
                ## Convert the corrected matrix data to numeric type and miltiply by 1000
                ## Convert the matrix of raw data to numerics
                loessDatmp <- apply(loessDat, 2, function(x) as.numeric(as.character(x)))
                loessDatT <- loessDatmp * 1000
                rownames(loessDatT) <- rownames(dat)
                datmp <- apply(dat, 2, function(x) as.numeric(as.character(x)))
                rownames(datmp) <- rownames(loessDatT)
                
                ## Imputation for all zeros and negative values:
                ## Set all negative to zeros: negative values may arise from due low dection of some compounds.
                if (length(loessDatT[loessDatT < 0L]) > 0) {
                        loessDatT[loessDatT < 0L] <- 0
                }
                
                ## Set all zeros to NAs
                loessDatT[loessDatT == 0L] <- NA
                
                ## Impute the missing values (NAs) per compound per group using the half of the mimum values
                if (sum(is.na(loessDatT)) > 0) {
                        loessDatT <- minHalfValue_bm(t(loessDatT), classF)
                        loessDatT <- t(loessDatT)
                        
                        #remove the columns from InputData
                        ix=which(rownames(datmp)%in%rownames(loessDatT))
                        datmp=datmp[ix,]
                        
                        ic=which(colnames(inputedData)%in%rownames(loessDatT))
                        inputedData=inputedData[,ic]
                        
                        message("   - The number of missing value after QC-based signal correction: ",
                                sum(is.na(loessDatT)))
                        
                }
                
                dirsc.ID =dirout.Bs=dirout.As= getwd()
                ## Generating the scater plots showing uncorrect and corrected signal for sample and QCs
                if (plot) {
                        message("   - High-resolution images output...")
                        pbloplot <- txtProgressBar(min = 1, max = dim(datmp)[1], style = 3)
                        for (i in 1:dim(datmp)[1]) {
                                #statTarget:::loplot(datmp, loessDatT, i)
                                loplot(datmp, loessDatT, i)
                                setTxtProgressBar(pbloplot, i)
                        }
                        close(pbloplot)
                        ## Save a combined pdf
                        pdfs=list.files(pattern = "\\.pdf$",path = paste0(getwd(), "/loplot"))
                        pdf_combine(input=paste0(getwd(),"/loplot/",pdfs),output=paste0(getwd(),"/QC.corr.pdf"))
                }
                
                
                ## Compute the Relative Standard Deviation (CV) for the all sample groups (samples + QCs) {before correction}
                
                raw_temp<- cbind(metaDF, inputedData)
                pca1<-pltPCA(raw_temp,"Sample PCA - before correction")
                
                write.csv(raw_temp, paste(dirout.Bs, "AllGroups_raw.csv", sep = "/"))
                
                message("   - Calculation of CV distribution of raw peaks (samples  and QC groups)...") ## This is usefull for batch comparison before correction
                RsdCal_bm(raw_temp, group = TRUE, DistPattern = TRUE,output = TRUE,smpl="AllGroups",corr="before.cor")
                
                # Compute the Relative Standard Deviation (CV) for the all sample groups (samples + QCs) {after correction}
                lo_temp<- cbind(metaDF, t(loessDatT))
                pca2<-pltPCA(lo_temp,"Sample PCA - after correction")
                
                write.csv(lo_temp, paste(dirout.Bs, "AllGroups_corr.csv", sep = "/"))
                
                message("   - Calculation of CV distribution of raw peaks (samples  and QC groups) - after correction...") ## This is usefull for batch comparison before correction
                RsdCal_bm(lo_temp, group = TRUE, DistPattern = TRUE,output = TRUE,smpl="AllGroups",corr="after.cor")
                
                
                grobs = lapply(list(pca1,pca2), "+", theme(plot.margin = unit(c(.8,1,.8,1), "cm")))
                grobs_ <- lapply(grobs, ggplotGrob)
                ggsave(paste(dirout.Bs, "AllGroups_pca.pdf", sep = "/"), marrangeGrob(grobs = grobs_ , nrow=2, ncol=1), height=10, width=6)
                
                # Combined the RSD.pdfs
                
                
                # Write all corrected values in a CSV DF
                
                cat("\n", "Output Link:", getwd(), "\n")
                cat("\n", "Correction Finished! Time: ", date(),"\n")
                cat("\n #####################################\n")
                cat(" #  Software Version: statTarget 2.0  #\n")
                cat(" ######################################\n")
                setwd(dirsc.ID)
                
                ## Write the log DF containing the used parameters
                scPam1 <- c("Frule", "MLmethod", "ntree","QCspan", "degree", "imputeM")
                scPam2 <- c(Frule, MLmethod, ntree, QCspan, loess.degree, imputeM) #@@@@@@
                scpam <- data.frame(scPam1, scPam2)
                colnames(scpam) <- c("parameter", "value")
                write.table(scpam, paste0(getwd(), "/ParameterShiftCor.log"),row.names = FALSE)
                tmpfilesc = paste(getwd(), "/tmp", sep = "")
                unlink(tmpfilesc, recursive = TRUE)
                
                message("   - Signal correction completed.")
        }
        
        
        return(lo_temp)
}

#'Shows the order of injection
#'@param tb input data frame
#'
showOrder<-function(tb){
        qc=grep("qc",tolower(colnames(tb)))
        colz=rep("black",dim(tb)[2])
        colz[qc]="red"
        plot(as.numeric(tb[5,]),col=colz,main=paste0("Run oder - ", rownames(tb)[5]),lwd=2,xlab="Samples names",ylab="Signal")
        abline(v=qc,col="magenta",lwd=.5,lty=3)
        legend("topleft",pch=1,col=c("black","red"),legend=c("Samples","QCs"))
}


#'Computes Delta RSD before and after correction
#'@param RSD1path path to the RSD table before correction
#'@param RSD2path path to the RSD table after correction
#'@param pathDRSD output file
#'
deltaRSD<-function(RSD1path,RSD2path,pathDRSD){
        #test lines
        #RSD1path=paste0(subDir,"/RSDresult/RSD_all_before.cor.csv")
        #RSD2path=paste0(subDir,"/RSDresult/RSD_all_after.cor.csv")
        RSD1<-read.csv( RSD1path,stringsAsFactors = F)
        RSD2<-read.csv( RSD2path,stringsAsFactors = F)
        
        ix=grep("All/Total",RSD1[,1])
        RSD1=RSD1[-ix,]
        ix=grep("All/Total",RSD2[,1])
        RSD2=RSD2[-ix,]
        
        mRSD1=as.matrix(RSD1[,2:dim(RSD1)[2]]);rownames(mRSD1)=RSD1[,1]
        mRSD2=as.matrix(RSD2[,2:dim(RSD2)[2]]);rownames(mRSD2)=RSD2[,1] #2:dim(RSD2)[1]]
        
        dRSD=mRSD2-mRSD1
        
        
        colorz=c("black","red","blue","green","darkcyan","chocolate4","brown3","chartreuse4","darkmagenta","darkslateblue","dodgerblue4",
                 "gray30","deeppink","goldenrod4","darkblue","brown4","gray55","darkorange","deepskyblue","gold3","gray35","gray65","gray65","gray80")
        
        
        pdf(pathDRSD,width=20,height = 6)
        par(mar=c(10,2,2,2))#,mfrow=c(2,1))
        cix=1:dim(dRSD)[2]
        ylim=range(dRSD[,cix])
        
        plot(dRSD[1,cix],xaxt='n',col=colorz[1],cex=.8,ylim=ylim,ylab="dRSD",xlab="",main="Change in RSD after correction")
        for (i in 2:dim(dRSD)[1]){
                points(dRSD[i,cix],col=colorz[i],cex=.8)
        }
        
        abline(h=0,col="magenta",lty=2)
        legend("topright",pch=1,legend=rownames(dRSD),col=colorz[cix],cex=.5)
        axis(1,at=cix,labels = colnames(dRSD)[cix],las=2,cex.axis=.7)
        dev.off()
}


## Generates sample names using the Study factors from the metadata file
#'@param fileP full path to the meta data file
#'
sFactors2Names_Molar<-function(fileP){
        
        ## Test line
        # fileP=paste0(inPath_,"/",mtdFile)
        
        metaData=read_excel(path=fileP)
        metaData=as.data.frame(metaData)
        
        sF.newix=grep("factor",tolower(colnames(metaData)))
        
        ## Remove the study factor column of all NA's
        sF.ix=c();  for (ix in sF.newix){if (!all(is.na (metaData[,ix]))){sF.ix=c(sF.ix,ix)}}
        
        Replicate=c()
        Counter=list()
        
        for (i in 1:dim(metaData)[1]){
                
                if (length(which(is.na(metaData[i,sF.ix])))==0){    ## get sample names
                        curName=paste0(metaData[i,sF.ix], collapse = "_")
                        if (!curName %in% names(Counter)){Counter[[curName]]<-1}else{Counter[[curName]]=Counter[[curName]]+1}
                        Replicate=c(Replicate, paste0(curName,"_",Counter[[curName]]))
                }
                else if  (length(which(is.na(metaData[i,sF.ix])))==length(sF.ix)){ 
                        Replicate=c(Replicate, metaData$CoreID[i] )           ## get QC and Blank
                }else{
                        message("Input error: missing study factors. Please check your metadata file, in study factor colunms\n","column ",i," : ", paste0(metaData[i,sF.ix],collabse=" "))
                        stop(" --- MOLAR has fail to completed the analysis ---")
                }
        }
        
        #Df=cbind(metaData[1:6], Replicate, metaData[7:dim(metaData)[2]])
        Df=cbind(metaData[1:2], Replicate,metaData[3:dim(metaData)[2]])
        return(Df)
}


#'Reads the input files (Raw Data, Compound note, and metadata): it generates sample names based on study factors.
#'@param inPath_ path to the MOLAR input file
#'@param prefix_ the prefix on the inputfile (name without "_MOLAR_RawData.xlsx", "_MOLAR_CmpdNotes.xlsx" or "_MOLAR_Metadata.xlsx" )
#'
read_INPUT_files_Molar<-function(inPath_,  prefix_){
        
        # test lines
        # inPath_=inPath; prefix_=prefix
        
        inputFile=paste0( prefix_,"_MOLAR_RawData.xlsx")
        cmpdFile=paste0( prefix_,"_MOLAR_CmpdNotes.xlsx")
        mtdFile=paste0( prefix_,"_MOLAR_Metadata.xlsx")
        
        prefix_List_=list(gsub("\\.xlsx$","",inputFile), gsub("\\.xlsx$","",cmpdFile),gsub("\\.xlsx$","",mtdFile))
        
        ## Check files : metadata and compound table
        if(!file.exists(paste0(inPath_,"/",mtdFile))) {ret=paste0("\t\tError :", basename(inputDataFile_), " has no metadata file associated with it!"); return(ret)}
        if(!file.exists(paste0(inPath_,"/",cmpdFile))) {ret=paste0("\t\tError :", basename(inputDataFile_), " has no compound note file associated with it!"); return(ret)}
        if(!file.exists(paste0(inPath_,"/",inputFile))) {ret=paste0("\t\tError :", basename(inputDataFile_), " has no raw data file associated with it!"); return(ret)}
        
        
        # Read raw data (MOLAR)
        headers=read_excel(paste0(inPath_,"/",inputFile), col_names = TRUE)
        colTypes=c("text","date",rep("numeric",dim(headers)[2]-2))
        InputData_=read_excel(paste0(inPath_,"/",inputFile),col_types=colTypes,skip=2,col_names = FALSE);# colnames(RawData)[2:3]=c("Data.File","Date")
        colnames(InputData_)=colnames(headers)
        colnames(InputData_)[1:2]=c("Sample","RunOrder")
        
        # Read compound data
        compNotesDF_=read_excel(  paste0(inPath_,"/",cmpdFile),skip =13) #,col_names = TRUE)
        compNotesDF_=as.data.frame(compNotesDF_)
        
        # Read metaData
        ## Generates sample names using the Study factors from the metadata file
        metaData_=sFactors2Names_Molar(paste0(inPath_,"/",mtdFile))
        
        inMetaFile_<<-paste0(inPath_,"/", prefix_,"_MOLAR_MetadataE.csv")
        write.csv(metaData_,inMetaFile_,row.names = F,quote=F)
        
        metaData_=as.data.frame(metaData_)
        rownames(metaData_)=metaData_[ ,"Replicate"]
        
        return(list(InputData=InputData_, compNotesDF=compNotesDF_, metaData=metaData_, prefix=prefix_, prefix_List=prefix_List_))
}


# MIDA & MOLAR specific functions ----

#'Reads and compares the sample names list between input data and metadata 
#'@param InputData_ input data frame
#'@param metaData_ the metadata data frame
#'@param currModule_ Module currently running (MOLAR or MIDA)
#'
get_sampleNameList_Mida_Molar<-function(InputData_, metaData_, currModule_){
        # testing lines
        # InputData_=cD$InputData; metaData_=cD$metaData; currModule_=currModule
        
        names1=InputData_$Sample
        names2=metaData_$CoreID
        
        ret=""
        if (length(names1)!=length(names2)){
                ret=" Error: the number of sample names differ between the input data and metadata. Please check your files.";return(ret)
        }
        
        if (identical(sort(names1),sort(names2))==FALSE){
                ret=" Error: the list of sample names is not identical in  raw input data and metadata file. Please check your files.";return(ret)
        }
        
        slabeling=NULL
        if (currModule_=="MIDA"){slabeling=metaData_$Labeling}
        
        ret=list(snames=names2,new_snames=metaData_$Replicate, slabeling=slabeling)
        
        return(ret)
}


#'Gets sample names, sample groups, study factors etc
#'@param ss list/vector of sample names
#'
get_splitSampleNames_df_Mida_Molar<-function(ss){
        #ss=sampleNames
        
        ## remove the replicated number
        ss2=sub("_\\d{1,3}$", "",ss)
        group=ss2
        
        sFactors<-lapply(1:length(ss2), function(i){
                s=ss2[i]
                s_=strsplit(s,"_")[[1]]
                return(s_)
        })
        
        df=do.call(rbind,sFactors)
        
        sfCount=apply(df,2,function(xx){length(unique(xx))})
        imxx=which(sfCount==1)
        if (length(imxx)>1){df=df[,-imxx]}
        colnames(df)=paste0("Study_factor_",1:ncol(df))
        
        df=cbind(sample_names=ss, group, df)
        
        df=as.data.frame(df)
        return(df)
}


#### Section 2: ======   Main execution function of the MOLAR program    ==========

## main function

#'Main execution of the MOLAR pipeline
#'@param prefix_  prefix of the current output file  e.g "MC00028"
#'@param inPath_ the path of the input directory
#'@param QC.Blank.thr_ threshold for qc/blank-based filtration
#'@param NA.thr_  threshold for missing values-based filtration in a particular group (compounds with more than NA.thr_% missing values in any group are filtered out ) 
#'@param imput.method_ imputation method ("half of the minimum","minimum", "median" or  "KNN")
#'@param corr.method_ machine learning method for QC-base signal correction ("Loess Regression","random forest classifier)
#'@param QC.span_= Smoothing parameter (default = 0)
#'@param P.degree_ polynomial degree (default =2)

#'@param preproSignalCorr.Flg_ determines whether to skip or perform pre-processing (data filtration and imputation) and QC-correction
#'@param signalCorr.Flg_ determines whether to skip or perform QC-based signal correction 
#'@param qc.blank.Flg_ determines whether to skip or perform qc/blank-based filtration
#'@param NA.thr.Flg_ the threshold of proportion of NA at which a compound must be filtered out
#'@param imputation.Flg_ determines whether to perform imputation or not

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
mainFunMolar<-function(
        
                  prefix_,
                  inPath_,
                  
                  QC.Blank.thr_,
                  NA.thr_,
                  imput.method_,
                  corr.method_,
                  QC.span_=0,
                  P.degree_=2,
                  
                  preproSignalCorr.Flg_,
                  signalCorr.Flg_,
                  qc.blank.Flg_,
                  NA.thr.Flg_,
                  imputation.Flg_,
                  
                  topMetaboFDR_,
                  normMethod_,
                  transMethod_,
                  scalingMethod_,
                  
                  KeggPathwayAnalysisCheck_,
                  KeggPathwayEnrichMethod_,
                  KeggPathwayOrganism_,
                  KeggPathwayFDR_,
                  
                  run.Type_="profiling"
                  
){
     
        # *** Testing lines ----
        # ****** Start executing this block below to systematically run each line/section of the code in order to debug the program.
        testing=FALSE
        if (testing==TRUE){
                
                app.Path<<-"/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/"
                app.Path<<-"/Volumes/projects_secondary/mnp/tool.dev/EMeRGe1.0"
                #app.Path<<-"E:/RHome/Programs/EMeRGe1.0/"
                
                methods_regex <<- "HILIC_[AB]_|LIPID_|T3pos_|PHpos_|PHneg_|PHTCA_"
                methods_regex2<<-gsub("_\\|","\\|",methods_regex); methods_regex2=gsub("_$","",methods_regex2)
                idx_PB_Blanks_pattern<<-"Blank|pb|qc|pos.id|neg.id|sample_"
                        
                setwd(app.Path)
                
                source("midaUtils.R")
                source("molarUtils.R")
                source("idxUtils.R")
                
                
                inputDataFile<<-"/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MOLAR_test/LCQQQ_RawData.xlsx"
                inputDataFile<<-"/Volumes/projects_secondary/mnp/tool.dev/EMeRGe1.0/testData/MOLAR_test/LCQQQ_RawData.xlsx"
                inputDataFile<<-"/Users/Brejnev.Muhire/Downloads/LCQQQ_MOLAR_CmpdNotes.xlsx"
                inputDataFile<<-"/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/MOLAR_test/MC00119_LCQQQ_MOLAR_Metadata.xlsx"
                inputDataFile<<-"/Users/Brejnev.Muhire/Downloads/MC00142_LCQQQ_MOLAR_RawData.xlsx"
                inputDataFile<<-"/Users/Brejnev.Muhire/Downloads/test/GC_MC00137_MOLAR_Metadata.xlsx"
                
                
                inPath_<<-dirname(inputDataFile)
                
                prefix_<<-gsub(pattern="_MOLAR_RawData\\.xlsx$|_MOLAR_CmpdNotes\\.xlsx$|_MOLAR_Metadata\\.xlsx$" ,"", basename(inputDataFile))
                
                QC.Blank.thr_=FALSE
                NA.thr_=5
                imput.method_="minHalf"
                corr.method_="QCRLSC"
                QC.span_=0
                P.degree_=2
                
                preproSignalCorr.Flg_=FALSE
                signalCorr.Flg_=FALSE
                qc.blank.Flg_=5
                NA.thr.Flg_=0.35
                imputation.Flg_=FALSE
                
                topMetaboFDR_=0.1
                normMethod_="None"
                transMethod_="LogNorm"
                scalingMethod_="MeanCenter"
                
                KeggPathwayAnalysisCheck_=TRUE
                KeggPathwayEnrichMethod_="QEA"
                KeggPathwayOrganism_="mmu"
                KeggPathwayFDR_=0.1
                
                run.Type_="profiling"
                
                currModule<<-"MOLAR"
                
                cD<<-read_INPUT_files_Molar(inPath_=inPath_,  prefix_=prefix_)
             
                
        }
        
        
        message( "\n#------------- Running MOLAR - QC-based Signal Correction and Statistical Analysis  ------------#")
        message( "                     *---      ", format(Sys.time(), "%a, %b %Y %d, %X"),"      ---*\n")
        
        options(stringsAsFactors = F, scipen=999)
        
        
        ## Check whether study settings are selected
        
        message("\n#---------------------------------  Current settings  ----------------------------------------#\n")
        message("    * Directory settings:")
        message("      Working directory        : ", inPath_)
        message("      Input files              : ", paste0(unlist(cD$prefix_List), collapse=", "))
       
        message("\n    * Preprocessing and QC-based signal correction : ",preproSignalCorr.Flg_)
        message("       Signal correction ->")
        message("       Sample QC / Blank filteration threshold.      :  > ", QC.Blank.thr_)  # The mean QC sample should be FC.filter times the blank1
        message("       Group missing value, filteration threshold    :  ", NA.thr_)
        message("       Imputation method                             :  ", imput.method_)
        message("       QC-based signal correction method             :  ", corr.method_)
        message("       Smoothing parameter - QC.span                 :  ", QC.span_)
        message("       Polynomial degree                             :  ", P.degree_)
        
        message("\n     * Statistical Analyses settings :")
        message("       Top metabolites FDR                           :  ", topMetaboFDR_)
        message("       Normalization method                          :  ", normMethod_)
        message("       Data transformation method                    :  ", transMethod_)
        message("       Data scaling method                           :  ", scalingMethod_)
        
        message("\n     * Pathway Analyses settings")
        message("       Organism                                      :  ", KeggPathwayOrganism_)
        message("       Pathway enrichment method                     :  ", KeggPathwayEnrichMethod_)
        message("       FDR cutt-off                                  :  ", KeggPathwayFDR_)
        
       
        message("\n#--------------------------------------------------------------------------------------------#\n")
        
        message(paste0("\n*-------------  MOLAR : Analysing ", prefix ,"  data set ------------- ",date(), "\n") )
        
        ### Step 0: Reading input files ... ----
        message("\n** Step 0: Reading input files ...")
        
        # Creating the otuput directory (remove it if it exists)
        subDir=paste0(inPath_,"/",prefix_,"_out/")
        if (!dir.exists(subDir)){dir.create(subDir)} else {unlink(subDir, recursive = TRUE); dir.create(subDir)}
        setwd(subDir)
        
        # loading study settings:
        load(file=paste0(inPath_,"/",prefix_,"_studySettings.RData"))
        
        
        
        ### Step 1: Preprocessing (Reading and validating inputdata) ----
        
        #  Get input table for raw data, metadata and compound note from the cD list of tables made in the server() function
        InputData<<-cD$InputData
        compNotesDF<<-cD$compNotesDF
        metaData<<-cD$metaData
      
        # Reformat compound names from Molar raw input table and compound Note table 
        colnames(InputData)=gsub("\\s|'",".",colnames(InputData))
        colnames(InputData)=gsub("-","_",colnames(InputData))
        colnames(InputData)=RenameStart(colnames(InputData))
        
        compNotesDF$Compound=gsub("\\s|'",".",compNotesDF$Compound)
        compNotesDF$Compound=gsub("-","_",compNotesDF$Compound)
        compNotesDF$Compound=RenameStart(compNotesDF$Compound) # any starting digit get a "X." in front of it (5-Glucose -> X.5-Glucose)
        rownames(compNotesDF)=compNotesDF$Compound
        
        
        # Check whether the compound names match between raw data and compound note
        if (!identical( sort(colnames(InputData)[3:dim(InputData)[2]]),sort(compNotesDF$Compound))){
                
                message("   - The compounds names in the raw data table do not match those in the compound note table.\nPlease check the files and start again.")
                message("   - Unmatching compound names InputData vs CompoundNotes:")
                
                v11=sort(colnames(InputData)[3:dim(InputData)[2]])
                v22=sort(compNotesDF$Compound)
                for (i in 1:length(v11)){
                        if (v11[i]!=v22[i]){
                                message(v11[i], " <> ",v22[i])
                        }
                }
                message(paste0(inFile_," is not processed \n------------------"))
                setwd(getCurrentFileLocation()) 
                return()
        }
        
        
        
        # Filter out the compound marked in the compound note file
        
        message(paste("   * Starting with ",dim(InputData)[2]-2," compounds"))
        
        XsID=which(tolower(compNotesDF$Response)=="x")
        if (length(XsID)>0){                   #if there are compound marked by x remove them
                cmpL=compNotesDF$Compound[-XsID]
                message(paste0("   * ",length(XsID)," compounds marked by 'x'  in CompoundNote are filtered out"))
        }else{
                cmpL=compNotesDF$Compound
        }
        message(paste("   * Remaining with ",length(cmpL)," compounds"))
        
        ix=which(colnames(InputData)%in%cmpL)
        InputData=as.data.frame(InputData[ ,c(1,2,ix)])
        
        
        
        # Check for duplicated names & Get sample names from metaData table
        if (length(which(duplicated(metaData$Replicate)))>0){
                message("   - Error: Duplicated sample names in metaData table detected, 'Replicate colunm'. Please check  your input data!")
                message(metaData$Replicate[which(duplicated(metaData$Replicate))])
                setwd(getCurrentFileLocation()) 
                return()
        } 
        
        if (length(which(duplicated(InputData$Sample)))>0){
                message("   - Error: Duplicated sample name(s) detected in  InputData table, 'Sample' colunm. Please check your input data!")
                message(InputData$Sample[which(duplicated(InputData$Sample))])
                setwd(getCurrentFileLocation()) 
                return()
        }
        
        # Check whether the sample names match between metaData  and InputData
        if (!identical( sort(InputData$Sample),sort(metaData$CoreID))){
                meta=setdiff(metaData$CoreID,InputData$Sample)
                raw=setdiff(InputData$Sample,metaData$CoreID)
                
                if (length(meta) > 0){
                        message (paste("   - Metedata unmatched names : ",setdiff(InputData$Sample,metaData$CoreID)))
                }
                if (length(raw)>0){
                        message ("   - InputData unmatched names: ",setdiff(InputData$Sample,metaData$CoreID))
                }
                message("   - The sample names in the InputData table do not match those in the metaData table.\nPlease check the files and start again.")
                setwd(getCurrentFileLocation()) 
                return()
        }
        
        Sample.Names=unlist(lapply(InputData$Sample,function(rn){ metaData$Replicate[which(metaData$CoreID==rn)]}))
        InputData$Sample=Sample.Names
        rownames(InputData)=Sample.Names
        
        
        ### Step 2-4: QC/Blank filtration, signal correction ----
        
        if (preproSignalCorr.Flg_==TRUE){ # Perform filtration/imputation and QC-based signal corrections
                
                if(qc.blank.Flg_==FALSE){
                        message("** Step 1,  QCs / Blank filtration is be skipped")
                }
                if(qc.blank.Flg_==TRUE){
                        message("   - Performing QCs / Black filtration starts ...")
                        InputData = Filter_qc_blnk(df=InputData,cutoff=QC.Blank.thr_)  #QC.Blank.thr=0.1
                        message(paste("   - Remaining with ",dim(InputData)[2]-2," compounds"))
                }
                
                # Remove blanks from InputData1 & metaData1
                InputData1 = InputData
                metaData1=metaData
                
                #mcDbPos=grep("mc-db",tolower(rownames(InputData))) # remove proce double blanks if there is any
                mcDbPos=which(startsWith(tolower(rownames(InputData)),"mc-db"))
                
                if (length(mcDbPos)>0){
                        InputData1 = InputData[-mcDbPos,]
                        metaData1 = metaData[-mcDbPos,] 
                        InputData1 = as.data.frame(InputData1)
                }
                
                #mcPbPos=grep("mc-pb",tolower(rownames(InputData1)))  # remove processed blanks if there is any
                mcPbPos=which(startsWith(tolower(rownames(InputData1)),"mc-pb"))
                if (length(mcPbPos)>0){
                        metaData1 = metaData1[-mcPbPos,]
                        InputData1=InputData1[-mcPbPos,]
                }
                
                
                colx=unlist(lapply(c("UserSampleName","Replicate"),function(coln){grep(tolower(coln),tolower(colnames(metaData1)))}))
                if (length(colx)!=2){message("   - Input error: Please check the metadata column names 'UserSampleName','Replicate'"); setwd(getCurrentFileLocation()) ;return()}
                
                metaData1=metaData1[,colx]
                metaData1[,"GroupName"]= unlist(lapply(metaData1[,"Replicate"],function(nm){strsplit(nm,"_")[[1]][1]}))
                
                
                ##sorting row
                #o=order(metaData1[,"GroupName"])
                #metaData1=metaData1[o,]
                #InputData1=InputData1[o,]
                OrigSampleOrder<<-metaData1[,"Replicate"]
                
                
                ## For InputData2 and metaData2 the QCs are put at the end and the GroupNumber is added
                qcPos=which(startsWith(tolower(metaData1[,"GroupName"]),"qc"))
                qcPos2=which(startsWith(tolower(metaData1[,"Replicate"]),"qc_")) ## check the format
            
                if (length(qcPos)==0){
                        message("   - Input error: No QC samples where dected. Please check your input file or set the program to skip QC-based signal correction.\nQC-sample names must start with 'QC_' or 'qc_'\t")
                        return()
                } 
                        
                if(length(qcPos2)==0){
                        message("   - Input error: Please check the  QC-sample names, they must start with 'QC_' or 'qc_' eg: QC_1, QC_2 etc\n")
                        return()
                }
                metaData2=rbind(metaData1[-qcPos,],metaData1[qcPos,])
                InputData2=rbind(InputData1[-qcPos,],InputData1[qcPos,])
                   
                ##add group number
                metaData2[,"GroupNum"]= unlist(lapply(1:length(table(metaData2[,"GroupName"])),function(j){
                        gName=unique(metaData2[,"GroupName"])[j]
                        if ( length(grep("qc",tolower(gName)))>0){
                                ret=rep(NA,table(metaData2[,"GroupName"])[gName])
                        }else{
                                ret=rep(j,table(metaData2[,"GroupName"])[gName])
                        }
                        return(ret)
                }))
                
                
                
                # Get the runOder from data formated order, and sort by order
                a <- as.POSIXlt( InputData2$RunOrder, format = "%m/%d/%Y %I:%M %p")
                InputData2=InputData2[order(a),]
                InputData2$RunOrder=1:dim(InputData2)[1]
                
                #Ordering the metaData based on InputData sample order / runOrder
                metaData2=metaData2[match(InputData2$Sample,metaData2$Replicate),]
                
                # Testing lines: show the run order of samples
                #showOrder(t(InputData2))
                
                ## checking the starting  QC samples (3>=)
                sq=grep("qc",tolower(rownames(InputData2)))
                if (length(intersect(c(1,2,3),sq))>2){
                        r=NULL
                        for (k in 1:(length(sq)-1)){
                                if (abs(sq[k]-sq[(k+1)])>1){
                                        r=k
                                        break
                                }
                        } 
                        
                        if (is.null(r)){message("   - Input error: QC sample names"); 
                                setwd(getCurrentFileLocation()) 
                                return()  #???? meaning message 
                        }
                        else{
                                InputData2=InputData2[r:dim(InputData2)[1],]
                                metaData2=metaData2[r:dim(metaData2)[1],]
                        }
                }else{
                        message("   - The input does not contain the first 3 QC samples.\n Please check input and try again!\n")
                        setwd(getCurrentFileLocation()) 
                        return()
                }
                
                #showOrder(t(InputData2))
                ### Step 3: QC correction StatTarget, Loess, Random forest, get corrected values in corTable
                #corTable<-shiftCor_bm (metaDF=metaData, samDF=InputData, Frule = 0.8, MLmethod = "QCRFSC",ntree = 500, QCspan = 0, loess.degree = 2, imputeM = "KNN", plot = TRUE)
                
                corTable<-shiftCor_bm (metaDF=metaData2, 
                                       samDF=InputData2, 
                                       Frule = 1-NA.thr_, 
                                       imputeM = imput.method_ ,
                                       MLmethod = corr.method_, 
                                       ntree = 500, 
                                       QCspan = QC.span_, 
                                       loess.degree = P.degree_, 
                                       plot = TRUE,
                                       
                                       signalCorr.Flg_,
                                       NA.thr.Flg_,
                                       imputation.Flg_)
                
                if (signalCorr.Flg_==TRUE){
                        RSD1path=paste0(subDir,"/RSDresult/RSD_all_before.cor.csv")
                        RSD2path=paste0(subDir,"/RSDresult/RSD_all_after.cor.csv")
                        pathDRSD=paste0(subDir,"/DeltaRSD.pdf")
                        deltaRSD(RSD1path=RSD1path,RSD2path=RSD2path,pathDRSD=pathDRSD)
                }      
                corTable<-corTable[,-c(1,4)]; colnames(corTable)[1:2]=c("sample","group")
                ix=grep("qc",tolower(corTable$group)); corTable<-corTable[-ix,]
                
                ## Put back the InputData in the original sample order
                OrigSampleOrder= OrigSampleOrder[-grep("qc",tolower(OrigSampleOrder))]
                
                ixx=match(OrigSampleOrder,corTable$sample)
                
                ### Step 4: RSD based filter??, and saving the temp csv table for MetaboAnalyst
                #write.csv(corTable[ixx,],paste0(subDir,"/MolarInput_allgroups.csv"),row.names = FALSE)
                write.csv(corTable[ixx,],paste0(subDir, prefix_,"_PoolSize.csv"),row.names = FALSE)
                
        }else { # Skip filtration/imputation and QC-based signal corrections
                
                message("\n** Steps 2-4, QCs / Blank filtration, inputation, & QC based signal correction are skipped")
                
                # Remove blanks from InputData1 & metaData1
                InputData1 = InputData
                metaData1=metaData
                
                #mcDbPos=grep("mc-db",tolower(rownames(InputData))) # remove proce double blanks if there is any
                mcDbPos=which(startsWith(tolower(rownames(InputData)),"mc-db"))
                
                if (length(mcDbPos)>0){
                        InputData1 = InputData[-mcDbPos,]
                        metaData1 = metaData[-mcDbPos,] 
                        InputData1 = as.data.frame(InputData1)
                }
                
                #mcPbPos=grep("mc-pb",tolower(rownames(InputData1)))  # remove processed blanks if there is any
                mcPbPos=which(startsWith(tolower(rownames(InputData1)),"mc-pb"))
                if (length(mcPbPos)>0){
                        metaData1 = metaData1[-mcPbPos,]
                        InputData1=InputData1[-mcPbPos,]
                }
                
                
                # Get group names
                colx=unlist(lapply(c("UserSampleName","Replicate"),function(coln){grep(tolower(coln),tolower(colnames(metaData1)))}))
                if (length(colx)!=2){message("   - Input error: Please check the metadata column names 'UserSampleName','Replicate'"); setwd(getCurrentFileLocation()) ;return()}
                
                metaData1=metaData1[,colx]
                metaData1[,"GroupName"]= unlist(lapply(metaData1[,"Replicate"],function(nm){strsplit(nm,"_")[[1]][1]}))
                
                #qcPos=grep("qc",tolower(metaData1[,"GroupName"])) # remove qc samples if there is any
                qcPos=which(startsWith(tolower(metaData1[,"GroupName"]),"qc"))
                if (length(qcPos) > 0) {
                        metaData1=metaData1[-qcPos,]
                        InputData1=InputData1[-qcPos,]
                }
                
                ##sorting row
                #o=order(metaData1[,"GroupName"])
                #metaData2=metaData1[o,]
                #InputData2=InputData1[o,]
                
                #o=order(metaData1[,"GroupName"])
                
                metaData2=metaData1
                InputData2=InputData1
                
                
                ##******put qc ad the end
                InputData2$RunOrder=metaData2$GroupName
                corTable=InputData2
                colnames(corTable)[1:2]=c("sample","group")
                write.csv(corTable,paste0(subDir,prefix_,"_PoolSize.csv"),row.names = FALSE)
                
                
        }
        
        ### Step 5: Reformating the tables  ----
        message("\n** Step 5: Reformating the table") 
        

        ## Reformatting the MID_FE, PoolSize, and TotalLabel_FE ----
        #  Use the study settings for (sample removal, group exclusion, or study factor scaling)
        #  Adding a third replicates for the sample with Two group
        r=format_Tb_all(subDir_=subDir, prefix1__=prefix_, suffix_="_PoolSize.csv", studySettings_=studySettings, mod_="MOLAR"); if (!is.data.frame(r)){message(r); return()}
        PoolSize_global<<-r
        
        ##MetaboAnalyst
        ### Step 6: Statistical analysis and reporting ... (with MetaboAnalyst) ----
        message("\n** Step 6: Statistical analysis and reporting ... (with MetaboAnalyst)") 
        
        
        outPathDirs <- statistical_Analyses_all ( subDirCurr=subDir, 
                                                  
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
                                                  
                                                  lbl_="6"
                                                )
        
                                                 
        ##MetaboAnalyst
        ### Step 7: Generating the final report ----
        message("\n** Step 7: Generating the final report") 
        
        generateFinalReport(outPathDirz=outPathDirs, subDir_=subDir)
        
        setwd(app.Path)
        cat(app.Path)
        
        
        message(paste0("\n#-----------------------------   MOLAR run ", prefix_," completed successfully   --------------------------#"))
        message(paste0("                        *---      ", format(Sys.time(), "%a, %b %Y %d, %X"),"      ---*\n"))
        
        
}
