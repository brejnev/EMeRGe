#### idxUtil.R : contains function that are specific to idxMOLAR pipeline and functions shared between idxMOLAR, MOLAR and MIDA pipelines ##### ----  

#### Section 1: ====== Library import ====== 

source("./external.functions/general_data_utils.R")
source("./external.functions/general_lib_utils.R")
source("./external.functions/general_load_libs.R")
source("./external.functions/general_misc_utils.R")
source("./external.functions/general_proc_utils.R")
source("./external.functions/general_anot_utils.R")
source("./external.functions/general_norm_utils.R")

source("./external.functions/batch_effect_utils.R")

source("./external.functions/enrich_name_match.R")
source("./external.functions/enrich_mset.R")
source("./external.functions/enrich_stats.R")
source("./external.functions/enrich_graphics.R")

source("./external.functions/stats_chemometrics.R")
source("./external.functions/stats_plot3d.R")
source("./external.functions/stats_univariates.R")


suppressPackageStartupMessages({
  
      library(png)
      library(FactoMineR)
      library(factoextra)
      library(gridExtra)
      library(grid)
      library(readxl)
      library(readr)
      library(pathview)
      library(KEGGREST)
      library(curl)
      library(qs)
      library(rlang)
      library(plotly)
      library(tidyverse)
      library(shiny)
      library(shinyjs)
      library(RColorBrewer)
      library(colorRamps)
      library(reshape)
      library(readxl)
      library(dplyr)
      library(kableExtra)
      library(ggplot2)
      library(gridExtra)
      library(pheatmap)
      library(Matrix)
      library(crayon)
      library(igraph)
      library(magick)
      library(pdftools)
  
})


#*********
## This code we have used some MetaboanalystR functions for statistical analyses and ploting.
## Functions with names ending with "_bm" are those which Brejnev Muhire changed from their original versions in MetaboanalystR R-packages.
## Functions with names ending with "_all" are those share by all metabolomics pipelines (idxMOLAR, MOLAR and MIDA)
 #*********


## Drawing function
iArrows <- igraph:::igraph.Arrows


#### Section 2: ====== Definition Functions ======

# Modified Metaboanalyst common functions for MIDA / MOLAR / IDX-MOLAR  ----

#'Two plot summary plot: Feature View of before and after normalization
#'@description For each plot, the top is a box plot, bottom is a density plot
#'@usage PlotNormSummary(mSetObj, imgName, format, dpi, width)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images,
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}, Jasmine Chong
#'McGill University, Canada
#'@export
#'
PlotNormSummary_bm <- function(mSetObj=NA, imgName, title, format="png", dpi=72, width=NA){
        
        mSetObj <- .get.mSet(mSetObj);
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        
        if(is.na(width)){
                w <- 10.5; h <- 12.5;
        }else if(width==0){
                w = 7.2
                h = 9.5
        }else if(width>0){
                w = width
                h = width*1.25
                # w <- 7.2; h <- 9;
        }
        
        mSetObj$imgSet$norm <- imgName
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        layout(matrix(c(1,2,2,2,3,4,4,4), 4, 2, byrow = FALSE))
        
        # since there may be too many compounds, only plot a subsets (50) in box plot
        # but density plot will use all the data
        
        pre.inx<-GetRandomSubsetIndex(ncol(mSetObj$dataSet$proc), sub.num=50);
        namesVec <- colnames(mSetObj$dataSet$proc[,pre.inx, drop=FALSE]);
        
        # only get common ones
        nm.inx <- namesVec %in% colnames(mSetObj$dataSet$norm)
        namesVec <- namesVec[nm.inx];
        pre.inx <- pre.inx[nm.inx];
        
        norm.inx<-match(namesVec, colnames(mSetObj$dataSet$norm));
        namesVec <- substr(namesVec, 1, 12); # use abbreviated name
        
        rangex.pre <- range(mSetObj$dataSet$proc[, pre.inx, drop=FALSE], na.rm=T);
        rangex.norm <- range(mSetObj$dataSet$norm[, norm.inx, drop=FALSE], na.rm=T);
        
        x.label<-GetAbundanceLabel(mSetObj$dataSet$type);
        y.label<-GetVariableLabel(mSetObj$dataSet$type);
        
        # fig 1
        if(anal.type == "roc" & mSetObj$dataSet$roc_cols == 1){
                op<-par(mar=c(4,7,4,0), xaxt="s");
                plot.new()
        }else{
                op<-par(mar=c(4,7,4,0), xaxt="s");
                plot(density(apply(mSetObj$dataSet$proc, 2, mean, na.rm=TRUE)), col='darkblue', las =2, lwd=2, main="", xlab="", ylab="");
                mtext("Density", 2, 5,cex=.7);
                mtext(paste0(title,"\nBefore transformation"),3, 1, cex=.7)
        }
        
        # fig 2
        op<-par(mar=c(7,7,0,0), xaxt="s");
        boxplot(mSetObj$dataSet$proc[,pre.inx, drop=FALSE], names=namesVec, ylim=rangex.pre, las = 2, col="lightgreen", horizontal=T, show.names=T);
        mtext(x.label, 1, 5, cex=.7);
        
        # fig 3
        if(anal.type == "roc" & mSetObj$dataSet$roc_cols == 1){
                op<-par(mar=c(4,7,4,2), xaxt="s");
                plot.new()
        }else{
                op<-par(mar=c(4,7,4,2), xaxt="s");
                plot(density(apply(mSetObj$dataSet$norm, 2, mean, na.rm=TRUE)), col='darkblue', las=2, lwd =2, main="", xlab="", ylab="");
                mtext(paste0(title,"\nAfter transformation"),3, 1, cex=.7);
        }
        
        # fig 4
        op<-par(mar=c(7,7,0,2), xaxt="s");
        boxplot(mSetObj$dataSet$norm[,norm.inx, drop=FALSE], names=namesVec, ylim=rangex.norm, las = 2, col="lightgreen", horizontal=T, show.names=T);
        mtext(paste("Normalized",x.label),1, 5, cex=.7);
        
        dev.off();
        
        return(.set.mSet(mSetObj));
}

#'Two plot summary plot: Sample View of before and after normalization
#'@description For each plot, the top is a density plot and the bottom is a box plot.
#'@usage PlotSampleNormSummary(mSetObj=NA, imgName, format="png", dpi=72, width=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", of "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images,
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}, Jasmine Chong
#'McGill University, Canada
#'@export
#'
PlotSampleNormSummary_bm <- function(mSetObj=NA, imgName, title, format="png", dpi=72, width=NA){
        
        mSetObj <- .get.mSet(mSetObj);
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 10.5; h <- 12.5;
        }else if(width == 0){
                w <- 7.2;h <- 9.5;
        }else if(width>0){
                w = width
                h = width*1.25
                # w <- 7.2; h <- 9;
        }
        
        mSetObj$imgSet$summary_norm <-imgName;
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, title=title, width=w, height=h, type=format, bg="white");
        layout(matrix(c(1,1,1,2,3,3,3,4), 4, 2, byrow = FALSE))
        
        # since there may be too many samples, only plot a subsets (50) in box plot
        # but density plot will use all the data
        
        pre.inx<-GetRandomSubsetIndex(nrow(mSetObj$dataSet$proc), sub.num=50);
        namesVec <- rownames(mSetObj$dataSet$proc[pre.inx, , drop=FALSE]);
        
        # only get common ones
        nm.inx <- namesVec %in% rownames(mSetObj$dataSet$norm)
        namesVec <- namesVec[nm.inx];
        pre.inx <- pre.inx[nm.inx];
        
        norm.inx<-match(namesVec, rownames(mSetObj$dataSet$norm));
        namesVec <- substr(namesVec, 1, 12); # use abbreviated name
        
        rangex.pre <- range(mSetObj$dataSet$proc[pre.inx, , drop=FALSE], na.rm=T);
        rangex.norm <- range(mSetObj$dataSet$norm[norm.inx, , drop=FALSE], na.rm=T);
        
        x.label<-GetAbundanceLabel(mSetObj$dataSet$type);
        y.label<-"Samples";
        
        # fig 1
        op<-par(mar=c(5.75,8,4,0), xaxt="s");
        boxplot(t(mSetObj$dataSet$proc[pre.inx, , drop=FALSE]), names= namesVec, ylim=rangex.pre, las = 2, col="lightgreen", horizontal=T);
        mtext(paste0(title,"\nBefore Normalization"), 3,1,cex=.7)
        
        # fig 2
        op<-par(mar=c(6.5,7,0,0), xaxt="s");
        plot(density(apply(mSetObj$dataSet$proc, 1, mean, na.rm=TRUE)), col='darkblue', las =2, lwd=2, main="", xlab="", ylab="");
        mtext(x.label, 1, 4, cex=.7);
        mtext("Density", 2, 5,cex=.7);
        
        # fig 3
        op<-par(mar=c(5.75,8,4,2), xaxt="s");
        boxplot(t(mSetObj$dataSet$norm[norm.inx, , drop=FALSE]), names=namesVec, ylim=rangex.norm, las = 2, col="lightgreen", ylab="", horizontal=T);
        mtext(paste0(title,"\nAfter Normalization"), 3, 1,cex=.7);
        
        # fig 4
        op<-par(mar=c(6.5,7,0,2), xaxt="s");
        plot(density(apply(mSetObj$dataSet$norm, 1, mean, na.rm=TRUE)), col='darkblue', las=2, lwd =2, main="", xlab="", ylab="");
        mtext(paste("Normalized",x.label),1, 4,cex=.7)
        
        dev.off();
        return(.set.mSet(mSetObj));
}


#'Plot PCA pair summary, format image in png, tiff, pdf, ps, svg
#'@description Rotate PCA analysis
#'@usage PlotPCAPairSummary(mSetObj=NA, imgName, format="png", dpi=72, width=NA, pc.num)
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@param pc.num Numeric, input a number to indicate the number of principal components to display in the pairwise score plot.
#'@export
#'
PlotPCAPairSummary_bm <- function(mSetObj=NA, imgName, title, format="png", dpi=72, width=NA, pc.num){
        
        mSetObj <- .get.mSet(mSetObj);
        pclabels <- paste("PC", 1:pc.num, "\n", round(100*mSetObj$analSet$pca$variance[1:pc.num],1), "%");
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 10;
        }else if(width == 0){
                w <- 8;
        }else{
                w <- width;
        }
        
        mSetObj$imgSet$pca.pair <- imgName;
        
        h <- w;
        Cairo::Cairo(file = imgName, unit="in", title=title, dpi=dpi, width=w, height=h, type=format, bg="white");
        if(mSetObj$dataSet$cls.type == "disc"){
                pairs(mSetObj$analSet$pca$x[,1:pc.num], col=GetColorSchema(mSetObj$dataSet$cls),main=title, pch=as.numeric(mSetObj$dataSet$cls)+1, labels=pclabels);
        }else{
                pairs(mSetObj$analSet$pca$x[,1:pc.num],main=title , labels=pclabels);
        }
        dev.off();
        return(.set.mSet(mSetObj));
}


#'Plot PCA scree plot
#'@description Rotate PCA analysis
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@param scree.num Numeric, input a number to indicate the number of principal components to display in the scree plot.
#'@usage PlotPCAScree(mSetObj=NA, imgName, format="png", dpi=72, width=NA, scree.num)
#'@export
#'
PlotPCAScree_bm <- function(mSetObj=NA, imgName, title, format="png", dpi=72, width=NA, scree.num){
        
        mSetObj <- .get.mSet(mSetObj);
        
        stds <-mSetObj$analSet$pca$std[1:scree.num];
        pcvars<-mSetObj$analSet$pca$variance[1:scree.num];
        cumvars<-mSetObj$analSet$pca$cum.var[1:scree.num];
        
        ylims <- range(c(pcvars,cumvars));
        extd<-(ylims[2]-ylims[1])/10
        miny<- ifelse(ylims[1]-extd>0, ylims[1]-extd, 0);
        maxy<- ifelse(ylims[2]+extd>1, 1.0, ylims[2]+extd);
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 10;
        }else if(width == 0){
                w <- 8;
        }else{
                w <- width;
        }
        h <- w*2/3;
        
        mSetObj$imgSet$pca.scree <- imgName;
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        par(mar=c(5,5,6,3));
        plot(pcvars, type='l', col='blue', main=paste0('Scree plot - ',title), xlab='PC index', ylab='Variance explained', ylim=c(miny, maxy), axes=F)
        text(pcvars, labels =paste(100*round(pcvars,3),'%'), adj=c(-0.3, -0.5), srt=45, xpd=T)
        points(pcvars, col='red');
        
        lines(cumvars, type='l', col='green')
        text(cumvars, labels =paste(100*round(cumvars,3),'%'), adj=c(-0.3, -0.5), srt=45, xpd=T)
        points(cumvars, col='red');
        
        abline(v=1:scree.num, lty=3);
        axis(2);
        axis(1, 1:length(pcvars), 1:length(pcvars));
        dev.off();
        return(.set.mSet(mSetObj));
}

#'Create 2D PCA score plot
#'@description Rotate PCA analysis
#'@usage PlotPCA2DScore(mSetObj=NA, imgName, format="png", dpi=72, width=NA, pcx, pcy, reg = 0.95, show=1, grey.scale = 0)
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param format Select the image format, "png", or "pdf".
#'@param title The title to add on top of the graph 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@param pcx Specify the principal component on the x-axis
#'@param pcy Specify the principal component on the y-axis
#'@param reg Numeric, input a number between 0 and 1, 0.95 will display the 95 percent confidence regions, and 0 will not.
#'@param show Display sample names, 1 = show names, 0 = do not show names.
#'@param grey.scale Use grey-scale colors, 1 = grey-scale, 0 = not grey-scale.
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotPCA2DScore_bm <- function(mSetObj=NA, imgName, title, format="png", dpi=72, width=NA, pcx, pcy, reg = 0.95, show=1, grey.scale = 0){
        
        mSetObj <- .get.mSet(mSetObj);
        
        xlabel = paste("PC",pcx, "(", round(100*mSetObj$analSet$pca$variance[pcx],1), "%)");
        ylabel = paste("PC",pcy, "(", round(100*mSetObj$analSet$pca$variance[pcy],1), "%)");
        pc1 = mSetObj$analSet$pca$x[, pcx];
        pc2 = mSetObj$analSet$pca$x[, pcy];
        text.lbls<-substr(names(pc1),1,14) # some names may be too long
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 9;
        }else if(width == 0){
                w <- 7.2;
        }else{
                w <- width;
        }
        h <- w;
        
        mSetObj$imgSet$pca.score2d <- imgName;
        
        Cairo::Cairo(file = imgName, unit="in", title=title, dpi=dpi, width=w, height=h, type=format, bg="white");
        op<-par(mar=c(5,5,3,3));
        
        if(mSetObj$dataSet$cls.type == "disc"){
                # obtain ellipse points to the scatter plot for each category
                
                if(mSetObj$dataSet$type.cls.lbl=="integer"){
                        cls <- as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls]);
                }else{
                        cls <- mSetObj$dataSet$cls;
                        #cls <- mSetObj$dataSet$cls_bm;
                }
                
                lvs <- levels(cls);
                pts.array <- array(0, dim=c(100,2,length(lvs)));
                for(i in 1:length(lvs)){
                        inx <-mSetObj$dataSet$cls == lvs[i];
                        groupVar<-var(cbind(pc1[inx],pc2[inx]), na.rm=T);
                        groupMean<-cbind(mean(pc1[inx], na.rm=T),mean(pc2[inx], na.rm=T));
                        pts.array[,,i] <- ellipse::ellipse(groupVar, centre = groupMean, level = reg, npoints=100);
                }
                
                xrg <- range(pc1, pts.array[,1,]);
                yrg <- range(pc2, pts.array[,2,]);
                x.ext<-(xrg[2]-xrg[1])/12;
                y.ext<-(yrg[2]-yrg[1])/12;
                xlims<-c(xrg[1]-x.ext, xrg[2]+x.ext);
                ylims<-c(yrg[1]-y.ext, yrg[2]+y.ext);
                
                cols <- GetColorSchema(cls, grey.scale==1);
                uniq.cols <- unique(cols);
                
                plot(pc1, pc2, xlab=xlabel, xlim=xlims, ylim=ylims, ylab=ylabel, type='n', main=paste0('Scores plot - ',title),
                     col=cols, pch=as.numeric(mSetObj$dataSet$cls)+1); ## added
                grid(col = "lightgray", lty = "dotted", lwd = 1);
                
                # make sure name and number of the same order DO NOT USE levels, which may be different
                legend.nm <- unique(as.character(sort(cls)));
                ## uniq.cols <- unique(cols);
                
                ## BHAN: when same color is choosen; it makes an error
                if ( length(uniq.cols) > 1 ) {
                        names(uniq.cols) <- legend.nm;
                }
                
                # draw ellipse
                for(i in 1:length(lvs)){
                        if (length(uniq.cols) > 1) {
                                polygon(pts.array[,,i], col=adjustcolor(uniq.cols[lvs[i]], alpha=0.2), border=NA);
                        } else {
                                polygon(pts.array[,,i], col=adjustcolor(uniq.cols, alpha=0.2), border=NA);
                        }
                        if(grey.scale) {
                                lines(pts.array[,,i], col=adjustcolor("black", alpha=0.5), lty=2);
                        }
                }
                
                pchs <- GetShapeSchema(mSetObj, show, grey.scale);
                if(grey.scale) {
                        cols <- rep("black", length(cols));
                }
                if(show == 1){
                        text(pc1, pc2, label=text.lbls, pos=4, xpd=T, cex=0.75);
                        points(pc1, pc2, pch=pchs, col=cols);
                }else{
                        if(length(uniq.cols) == 1){
                                points(pc1, pc2, pch=pchs, col=cols, cex=1.0);
                        }else{
                                if(grey.scale == 1 | (exists("shapeVec") && all(shapeVec>=0))){
                                        my.cols <- adjustcolor(cols, alpha.f = 0.4);
                                        my.cols[pchs == 21] <- "black";
                                        points(pc1, pc2, pch=pchs, col=my.cols, bg=adjustcolor(cols, alpha.f = 0.4), cex=1.8);
                                }else{
                                        points(pc1, pc2, pch=21, bg=adjustcolor(cols, alpha.f = 0.4), cex=2);
                                }
                        }
                }
                uniq.pchs <- unique(pchs);
                if(grey.scale) {
                        uniq.cols <- "black";
                }
                
                if(length(lvs) < 6){
                        legend("topright", legend = legend.nm, pch=uniq.pchs, col=uniq.cols);
                }else if (length(lvs) < 10){
                        legend("topright", legend = legend.nm, pch=uniq.pchs, col=uniq.cols, cex=0.75);
                }else{
                        legend("topright", legend = legend.nm, pch=uniq.pchs, col=uniq.cols, cex=0.5);
                }
                
        }else{
                plot(pc1, pc2, xlab=xlabel, ylab=ylabel, type='n', main=paste0('Scores plot - ',title),);
                points(pc1, pc2, pch=15, col="magenta");
                text(pc1, pc2, label=text.lbls, pos=4, col ="blue", xpd=T, cex=0.8);
        }
        par(op);
        dev.off();
        return(.set.mSet(mSetObj));
}

#'Create 3D PCA score plot
#'@description Rotate PCA analysis
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@usage PlotPCA3DScore(mSetObj=NA, imgName, format="json", inx1, inx2, inx3)
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf". 
#'@param inx1 Numeric, indicate the number of the principal component for the x-axis of the loading plot.
#'@param inx2 Numeric, indicate the number of the principal component for the y-axis of the loading plot.
#'@param inx3 Numeric, indicate the number of the principal component for the z-axis of the loading plot.
#'@export
#'
PlotPCA3DScore_bm <- function(mSetObj=NA, imgName,title, format="json", inx1, inx2, inx3){
        mSetObj <- .get.mSet(mSetObj);
        
        pca <-  mSetObj$analSet$pca;
        pca3d <- list();
        pca3d$score$axis <- paste("PC", c(inx1, inx2, inx3), " (", 100*round( mSetObj$analSet$pca$variance[c(inx1, inx2, inx3)], 3), "%)", sep="");
        coords <- data.frame(t(signif(pca$x[,c(inx1, inx2, inx3)], 5)));
        colnames(coords) <- NULL;
        pca3d$score$xyz <- coords;
        pca3d$score$name <- rownames(mSetObj$dataSet$norm);
        
        if(mSetObj$dataSet$type.cls.lbl=="integer"){
                cls <- as.character(sort(as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls])));
        }else{
                cls <- as.character(mSetObj$dataSet$cls);
                #cls <- as.character(mSetObj$dataSet$cls_bm);
                
        }
        
        if(all.numeric(cls)){
                cls <- paste("Group", cls);
        }
        
        pca3d$score$facA <- cls;
        
        # now set color for each group
        cols <- unique(GetColorSchema(mSetObj$dataSet$cls));
        #cols <- unique(GetColorSchema(mSetObj$dataSet$cls_bm));
        rgbcols <- col2rgb(cols);
        cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")})
        pca3d$score$colors <- cols;
        imgName = paste(imgName, ".", format, sep="");
        json.obj <- RJSONIO::toJSON(pca3d, .na='null');
        sink(imgName);
        cat(json.obj);
        sink();
        
        if(!.on.public.web){
                return(.set.mSet(mSetObj));
        }
}

#'@export
PlotPCA3DLoading_bm <- function(mSetObj=NA, imgName,title, format="json", inx1, inx2, inx3){
        mSetObj <- .get.mSet(mSetObj);
        pca = mSetObj$analSet$pca
        coords<-signif(as.matrix(cbind(pca$rotation[,inx1],pca$rotation[,inx2],pca$rotation[,inx3])),5);
        pca3d <- list();
        
        pca3d$loading$axis <- paste("Loading ", c(inx1, inx2, inx3), sep="");
        coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
        
        colnames(coords) <- NULL; 
        pca3d$loading$xyz <- coords;
        pca3d$loading$name <- rownames(pca$rotation);
        pca3d$loading$entrez <-rownames(pca$rotation); 
        
        if(mSetObj$dataSet$type.cls.lbl=="integer"){
                cls <- as.character(sort(as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls])));
        }else{
                cls <- as.character(mSetObj$dataSet$cls);
        }
        
        if(all.numeric(cls)){
                cls <- paste("Group", cls);
        }
        
        pca3d$cls = cls;
        # see if there is secondary
        
        require(RJSONIO);
        imgName = paste(imgName, ".", format, sep="");
        json.mat <- toJSON(pca3d, .na='null');
        sink(imgName);
        cat(json.mat);
        sink();
        current.msg <<- "Annotated data is now ready for PCA 3D visualization!";
        
        if(.on.public.web){
                return(1);
        }else{
                return(.set.mSet(mSetObj));
        }
        
}


#'Plot PCA loadings and also set up the matrix for display
#'@description Rotate PCA analysis
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@usage PlotPCALoading(mSetObj=NA, imgName, format="png", dpi=72, width=NA, inx1, inx2, plotType, lbl.feat=1)
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@param inx1 Numeric, indicate the number of the principal component for the x-axis of the loading plot.
#'@param inx2 Numeric, indicate the number of the principal component for the y-axis of the loading plot.
#'@export
#'
PlotPCALoading_bm <- function(mSetObj=NA, imgName,title, format="png", dpi=72, width=NA, inx1, inx2){
        
        mSetObj <- .get.mSet(mSetObj);
        
        loadings<-as.matrix(cbind(mSetObj$analSet$pca$rotation[,inx1],mSetObj$analSet$pca$rotation[,inx2]));
        # sort based on absolute values of 1, 2 
        ord.inx <- order(-abs(loadings[,1]), -abs(loadings[,2]));
        loadings <- signif(loadings[ord.inx,],5);
        
        ldName1<-paste("Loadings", inx1);
        ldName2<-paste("Loadings", inx2);
        colnames(loadings)<-c(ldName1, ldName2);
        load.x.uniq <- jitter(loadings[,1]);
        names(load.x.uniq) <- rownames(loadings);
        mSetObj$analSet$pca$load.x.uniq <- load.x.uniq;
        mSetObj$analSet$pca$imp.loads<-loadings; # set up the loading matrix
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 9;
        }else if(width == 0){
                w <- 7.2;
        }else{
                w <- width;
        }
        h <- w;
        
        mSetObj$imgSet$pca.loading <- imgName;
        plotType <- mSetObj$analSet$pca$loading.type;
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        
        par(mar=c(6,5,2,6));
        
        plot(loadings[,1],loadings[,2], las=2, xlab=ldName1, ylab=ldName2, main=paste0("Loadings plot - ",title));
        
        mSetObj$pca.axis.lims <- par("usr"); # x1, x2, y1 ,y2
        grid(col = "lightgray", lty = "dotted", lwd = 1);
        points(loadings[,1],loadings[,2], pch=19, col=adjustcolor("magenta", alpha.f = 0.4));
        
        if(plotType=="all"){
                text(loadings[,1],loadings[,2], labels=substr(rownames(loadings), 1, 16), pos=4, col="blue", xpd=T);
        }else if(plotType == "custom"){
                if(length(mSetObj$custom.cmpds) > 0){
                        hit.inx <- colnames(mSetObj$dataSet$norm) %in% mSetObj$custom.cmpds;
                        text(loadings[hit.inx,1],loadings[hit.inx,2], labels=rownames(loadings)[hit.inx], pos=4, col="blue", xpd=T);
                }
        }else{
                # do nothing
        }
        
        dev.off();
        return(.set.mSet(mSetObj));
        
}

#'Create PCA Biplot, set xpd = T to plot outside margin
#'@description Rotate PCA analysis
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@usage PlotPCABiplot(mSetObj=NA, imgName, format="png", dpi=72, width=NA, inx1, inx2)
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@param inx1 Numeric, indicate the number of the principal component for the x-axis of the loading plot.
#'@param inx2 Numeric, indicate the number of the principal component for the y-axis of the loading plot.
#'@export
#'
PlotPCABiplot_bm <- function(mSetObj=NA, imgName,title, format="png", dpi=72, width=NA, inx1, inx2){
        
        mSetObj <- .get.mSet(mSetObj);
        choices = c(inx1, inx2);
        scores <- mSetObj$analSet$pca$x;
        lam <- mSetObj$analSet$pca$sdev[choices]
        n <- NROW(scores)
        lam <- lam * sqrt(n);
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 9;
        }else if(width == 0){
                w <- 7.2;
        }else{
                w <- width;
        }
        h <- w;
        
        mSetObj$imgSet$pca.biplot<-imgName;
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        biplot(t(t(scores[, choices]) / lam), t(t(mSetObj$analSet$pca$rotation[, choices]) * lam), main=title, xpd =T, cex=0.9);
        dev.off();
        return(.set.mSet(mSetObj));
}

### R script for MetaboAnalyst
### Description: 3D scatter plot for PDF report (not used for web display)
### Author: Jeff Xia, jeff.xia@mcgill.ca
### McGill University, Canada
### License: GNU GPL (>= 2)
Plot3D_bm <- function(x, y = NULL, z = NULL, color = par("col"), pch = NULL,
                      main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
                      xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40,
                      axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE,
                      x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL,
                      y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"),
                      lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE,
                      mar = c(5, 3, 4, 3) + 0.1, col.axis = par("col.axis"),
                      col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"),
                      cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"),
                      font.axis = par("font.axis"), font.lab = par("font.lab"),
                      lty.axis = par("lty"), lty.grid = 2, lty.hide = 1,
                      lty.hplot = par("lty"), log = "", ...)
# log not yet implemented
{
        ## Uwe Ligges <ligges@statistik.tu-dortmund.de>,
        ## http://www.statistik.tu-dortmund.de/~ligges
        ##
        ## For MANY ideas and improvements thanks to Martin Maechler!!!
        ## Parts of the help files are stolen from the standard plotting functions in R.
        
        mem.par <- par(mar = mar)
        x.scal <- y.scal <- z.scal <- 1
        xlabel <- if (!missing(x)) deparse(substitute(x))
        ylabel <- if (!missing(y)) deparse(substitute(y))
        zlabel <- if (!missing(z)) deparse(substitute(z))
        
        ## color as part of `x' (data.frame or list):
        if(!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 4))
                color <- x[,4]
        else if(is.list(x) && !is.null(x$color))
                color <- x$color
        
        ## convert 'anything' -> vector
        xyz <- xyz.coords(x=x, y=y, z=z, xlab=xlabel, ylab=ylabel, zlab=zlabel,
                          log=log)
        if(is.null(xlab)) { xlab <- xyz$xlab; if(is.null(xlab)) xlab <- "" }
        if(is.null(ylab)) { ylab <- xyz$ylab; if(is.null(ylab)) ylab <- "" }
        if(is.null(zlab)) { zlab <- xyz$zlab; if(is.null(zlab)) zlab <- "" }
        
        if(length(color) == 1)
                color <- rep(color, length(xyz$x))
        else if(length(color) != length(xyz$x))
                stop("length(color) ", "must be equal length(x) or 1")
        
        angle <- (angle %% 360) / 90
        yz.f <- scale.y * abs(if(angle < 1) angle else if(angle > 3) angle - 4 else 2 - angle)
        yx.f <- scale.y * (if(angle < 2) 1 - angle else angle - 3)
        if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
                temp <- xlab;  xlab <- ylab;   ylab <- temp
                temp <- xlim;  xlim <- ylim;   ylim <- temp
        }
        angle.1 <- (1 < angle && angle < 2) || angle > 3
        angle.2 <- 1 <= angle && angle <= 3
        dat <- cbind(as.data.frame(xyz[c("x","y","z")]), col = color)
        
        n <- nrow(dat);
        y.range <- range(dat$y[is.finite(dat$y)])
        
        ### 3D-highlighting / colors / sort by y
        if(type == "p" || type == "h") {
                y.ord <- rev(order(dat$y))
                dat <- dat[y.ord, ]
                if(length(pch) > 1)
                        if(length(pch) != length(y.ord))
                                stop("length(pch) ", "must be equal length(x) or 1")
                else pch <- pch[y.ord]
                daty <- dat$y
                daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
                if(highlight.3d && !(all(diff(daty) == 0)))
                        dat$col <- rgb(seq(0, 1, length = n) * (y.range[2] - daty) / diff(y.range), g=0, b=0)
        }
        
        ### optim. axis scaling
        p.lab <- par("lab")
        ## Y
        y.range <- range(dat$y[is.finite(dat$y)], ylim)
        y.prty <- pretty(y.range, n = lab[2],
                         min.n = max(1, min(.5 * lab[2], p.lab[2])))
        y.scal <- round(diff(y.prty[1:2]), digits = 12)
        y.add <- min(y.prty)
        dat$y <- (dat$y - y.add) / y.scal
        y.max <- (max(y.prty) - y.add) / y.scal
        
        x.range <- range(dat$x[is.finite(dat$x)], xlim)
        x.prty <- pretty(x.range, n = lab[1],
                         min.n = max(1, min(.5 * lab[1], p.lab[1])))
        x.scal <- round(diff(x.prty[1:2]), digits = 12)
        dat$x <- dat$x / x.scal
        x.range <- range(x.prty) / x.scal
        x.max <- ceiling(x.range[2])
        x.min <-   floor(x.range[1])
        if(!is.null(xlim)) {
                x.max <- max(x.max, ceiling(xlim[2] / x.scal))
                x.min <- min(x.min,   floor(xlim[1] / x.scal))
        }
        x.range <- range(x.min, x.max)
        ## Z
        z.range <- range(dat$z[is.finite(dat$z)], zlim)
        z.prty <- pretty(z.range, n = lab.z,
                         min.n = max(1, min(.5 * lab.z, p.lab[2])))
        z.scal <- round(diff(z.prty[1:2]), digits = 12)
        dat$z <- dat$z / z.scal
        z.range <- range(z.prty) / z.scal
        z.max <- ceiling(z.range[2])
        z.min <-   floor(z.range[1])
        if(!is.null(zlim)) {
                z.max <- max(z.max, ceiling(zlim[2] / z.scal))
                z.min <- min(z.min,   floor(zlim[1] / z.scal))
        }
        z.range <- range(z.min, z.max)
        
        ### init graphics
        plot.new()
        if(angle.2) {x1 <- x.min + yx.f * y.max; x2 <- x.max}
        else        {x1 <- x.min; x2 <- x.max + yx.f * y.max}
        plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
        temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
        if(angle.2) x1 <- x1 - temp - y.margin.add
        else        x2 <- x2 + temp + y.margin.add
        plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
        if(angle > 2) par("usr" = par("usr")[c(2, 1, 3:4)])
        usr <- par("usr") # we have to remind it for use in closures
        title(main, sub, ...)
        
        ### draw axis, tick marks, labels, grid, ...
        xx <- if(angle.2) c(x.min, x.max) else c(x.max, x.min)
        if(grid) {
                ## grids
                ###################
                # XY wall
                i <- x.min:x.max;
                segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min,
                         col = col.grid, lty = lty.grid);
                
                i <- 0:y.max;
                segments(x.min + (i * yx.f), i * yz.f + z.min,
                         x.max + (i * yx.f), i * yz.f + z.min,
                         col = col.grid, lty = lty.grid);
                
                ######################
                # XZ wall
                # verticle lines
                temp <- yx.f * y.max;
                temp1 <- yz.f * y.max;
                i <- (x.min + temp):(x.max + temp);
                segments(i, z.min + temp1, i, z.max + temp1,
                         col = col.grid, lty = lty.grid);
                
                # horizontal lines
                i <- (z.min + temp1):(z.max + temp1);
                segments(x.min + temp, i, x.max + temp, i,
                         col = col.grid, lty = lty.grid)
                
                
                ##################
                # YZ wall
                # horizontal lines
                i <- xx[2]:x.min;
                mm <- z.min:z.max;
                segments(i, mm, i + temp, mm + temp1,
                         col = col.grid, lty = lty.grid);
                # verticle lines
                i <- 0:y.max;
                segments(x.min + (i * yx.f), i * yz.f + z.min,
                         xx[2] + (i * yx.f), i * yz.f + z.max,
                         col = col.grid, lty = lty.grid)
                
                
                # make the axis into solid line
                segments(x.min, z.min, x.min + (yx.f * y.max), yz.f * y.max + z.min,
                         col = col.grid, lty = lty.hide);
                segments(x.max, z.min, x.max + (yx.f * y.max), yz.f * y.max + z.min,
                         col = col.axis, lty = lty.hide);
                segments(x.min + (y.max * yx.f), y.max * yz.f + z.min,
                         x.max + (y.max* yx.f), y.max * yz.f + z.min,
                         col = col.grid, lty = lty.hide);
                segments(x.min + temp, z.min + temp1, x.min + temp, z.max + temp1,
                         col = col.grid, lty = lty.hide);
                segments(x.max + temp, z.min + temp1, x.max + temp, z.max + temp1,
                         col = col.axis, lty = lty.hide);
                segments(x.min + temp, z.max + temp1, x.max + temp, z.max + temp1,
                         col = col.axis, lty = lty.hide);
                segments(xx[2], z.max, xx[2] + temp, z.max + temp1,
                         col = col.axis, lty = lty.hide);
        }
        if(axis) {
                if(tick.marks) { ## tick marks
                        xtl <- (z.max - z.min) * (tcl <- -par("tcl")) / 50
                        ztl <- (x.max - x.min) * tcl / 50
                        mysegs <- function(x0,y0, x1,y1)
                                segments(x0,y0, x1,y1, col=col.axis, lty=lty.axis)
                        ## Y
                        i.y <- 0:y.max
                        mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min,
                               yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
                        ## X
                        i.x <- x.min:x.max
                        mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
                        ## Z
                        i.z <- z.min:z.max
                        mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)
                        
                        if(label.tick.marks) { ## label tick marks
                                las <- par("las")
                                mytext <- function(labels, side, at, ...)
                                        mtext(text = labels, side = side, at = at, line = -.5,
                                              col=col.lab, cex=cex.axis, font=font.lab, ...)
                                ## X
                                if(is.null(x.ticklabs))
                                        x.ticklabs <- format(i.x * x.scal)
                                mytext(x.ticklabs, side = 1, at = i.x)
                                ## Z
                                if(is.null(z.ticklabs))
                                        z.ticklabs <- format(i.z * z.scal)
                                mytext(z.ticklabs, side = if(angle.1) 4 else 2, at = i.z,
                                       adj = if(0 < las && las < 3) 1 else NA)
                                ## Y
                                temp <- if(angle > 2) rev(i.y) else i.y ## turn y-labels around
                                if(is.null(y.ticklabs))
                                        y.ticklabs <- format(y.prty)
                                else if (angle > 2)
                                        y.ticklabs <- rev(y.ticklabs)
                                text(i.y * yx.f + xx[1],
                                     i.y * yz.f + z.min, y.ticklabs,
                                     pos=if(angle.1) 2 else 4, offset=1,
                                     col=col.lab, cex=cex.axis/par("cex"), font=font.lab)
                        }
                }
                
                ## axis and labels
                
                mytext2 <- function(lab, side, line, at)
                        mtext(lab, side = side, line = line, at = at, col = col.lab,
                              cex = cex.lab, font = font.axis, las = 0)
                ## X
                lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, lty = lty.axis)
                mytext2(xlab, 1, line = 1.5, at = mean(x.range))
                ## Y
                lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + z.min),
                      col = col.axis, lty = lty.axis)
                mytext2(ylab, if(angle.1) 2 else 4, line= 0.5, at = z.min + y.max * yz.f)
                
                ## Z
                lines(xx[c(2,2)], c(z.min, z.max), col = col.axis, lty = lty.axis)
                mytext2(zlab, if(angle.1) 4 else 2, line= 1.5, at = mean(z.range))
                
        }
        
        ### plot points
        x <- dat$x + (dat$y * yx.f)
        z <- dat$z + (dat$y * yz.f)
        col <- as.character(dat$col)
        if(type == "h") {
                z2 <- dat$y * yz.f + z.min
                segments(x, z, x, z2, col = col, cex = cex.symbols, lty = lty.hplot, ...)
                points(x, z, type = "p", col = col, pch = pch, cex = cex.symbols, ...)
        }
        else points(x, z, type = type, col = col, pch = pch, cex = cex.symbols, ...)
        
        ### box-lines in front of points (overlay)
        if(axis && box) {
                lines(c(x.min, x.max), c(z.max, z.max),
                      col = col.axis, lty = lty.axis)
                lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max,
                      col = col.axis, lty = lty.axis)
                lines(xx[c(1,1)], c(z.min, z.max), col = col.axis, lty = lty.axis)
        }
        
        
        # par(mem.par) # we MUST NOT set the margins back
        ### Return Function Object
        ob <- ls() ## remove all unused objects from the result's enviroment:
        rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", "y.scal", "z.scal", "yx.f",
                                "yz.f", "y.add", "z.min", "z.max", "x.min", "x.max", "y.max",
                                "x.prty", "y.prty", "z.prty")])
        rm(ob)
        invisible(list(
                xyz.convert = function(x, y=NULL, z=NULL) {
                        xyz <- xyz.coords(x, y, z)
                        if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
                        }
                        y <- (xyz$y - y.add) / y.scal
                        return(list(x = xyz$x / x.scal + yx.f * y,
                                    y = xyz$z / z.scal + yz.f * y))
                },
                points3d = function(x, y = NULL, z = NULL, type = "p", ...) {
                        xyz <- xyz.coords(x, y, z)
                        if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
                        }
                        y2 <- (xyz$y - y.add) / y.scal
                        x <- xyz$x / x.scal + yx.f * y2
                        y <- xyz$z / z.scal + yz.f * y2
                        mem.par <- par(mar = mar, usr = usr)
                        on.exit(par(mem.par))
                        if(type == "h") {
                                y2 <- z.min + yz.f * y2
                                segments(x, y, x, y2, ...)
                                points(x, y, type = "p", ...)
                        }
                        else points(x, y, type = type, ...)
                },
                plane3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                                   lty = "dashed", lty.box = NULL, ...){
                        if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
                        if(is.null(lty.box)) lty.box <- lty
                        if(is.null(x.coef) && length(Intercept) == 3){
                                x.coef <- Intercept[if(angle > 2) 3 else 2]
                                y.coef <- Intercept[if(angle > 2) 2 else 3]
                                Intercept <- Intercept[1]
                        }
                        mem.par <- par(mar = mar, usr = usr)
                        on.exit(par(mem.par))
                        x <- x.min:x.max
                        ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
                        x.coef <- x.coef * x.scal
                        z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
                        z2 <- (Intercept + x * x.coef +
                                       (y.max * y.scal + y.add) * y.coef) / z.scal
                        segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
                        y <- 0:y.max
                        ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
                        y.coef <- (y * y.scal + y.add) * y.coef
                        z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
                        z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
                        segments(x.min + y * yx.f, z1 + y * yz.f,
                                 x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
                },
                
                wall3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                                  lty = "dashed", lty.box = NULL, ...){
                        if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
                        if(is.null(lty.box)) lty.box <- lty
                        if(is.null(x.coef) && length(Intercept) == 3){
                                x.coef <- Intercept[if(angle > 2) 3 else 2]
                                y.coef <- Intercept[if(angle > 2) 2 else 3]
                                Intercept <- Intercept[1]
                        }
                        mem.par <- par(mar = mar, usr = usr)
                        on.exit(par(mem.par))
                        x <- x.min:x.max
                        ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
                        x.coef <- x.coef * x.scal
                        z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
                        z2 <- (Intercept + x * x.coef +
                                       (y.max * y.scal + y.add) * y.coef) / z.scal
                        segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
                        y <- 0:y.max
                        ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
                        y.coef <- (y * y.scal + y.add) * y.coef
                        z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
                        z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
                        segments(x.min + y * yx.f, z1 + y * yz.f,
                                 x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
                },
                box3d = function(...){
                        mem.par <- par(mar = mar, usr = usr)
                        on.exit(par(mem.par))
                        lines(c(x.min, x.max), c(z.max, z.max), ...)
                        lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, ...)
                        lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
                        lines(c(x.max, x.max), c(z.min, z.max), ...)
                        lines(c(x.min, x.min), c(z.min, z.max), ...)
                        lines(c(x.min, x.max), c(z.min, z.min), ...)
                }
        ))
}

#'Create 3D PCA score plot
#'@description This function creates both a static 3D PCA score plot as well as an interactive 3D PCA score plot
#'using the plotly R package. The 3D PCA score plot is stored in the mSetObj (mSetObj$imgSet$pca.3d). To view
#'the plot, if your mSetObj is named mSet, type "mSet$imgSet$pca.3d" inro your R console, and the 3D plot will appear.
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@usage PlotPCA3DScoreImg(mSetObj=NA, imgName, format="png", dpi=72, width=NA, inx1, inx2, inx3, angl)
#'@param mSetObj Input name of the created mSet Object.
#'@param imgName Input a name for the plot.
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@param inx1 Numeric, indicate the number of the principal component for the x-axis of the loading plot.
#'@param inx2 Numeric, indicate the number of the principal component for the y-axis of the loading plot.
#'@param inx3 Numeric, indicate the number of the principal component for the z-axis of the loading plot.
#'@param angl Input the angle 
#'@usage mSet <- PlotPCA3DScore(mSetObj=NA, imgName, format="json", dpi=72, width=NA, inx1, inx2, inx3, angl)
#'@export
#'@importFrom plotly plot_ly add_markers layout
#'
PlotPCA3DScoreImg_bm <- function(mSetObj=NA, imgName,  title, format="png", dpi=72, width=NA, inx1, inx2, inx3, angl){
        
        mSetObj <- .get.mSet(mSetObj);
        
        xlabel = paste("PC",inx1, "(", round(100*mSetObj$analSet$pca$variance[inx1],1), "%)");
        ylabel = paste("PC",inx2, "(", round(100*mSetObj$analSet$pca$variance[inx2],1), "%)");
        zlabel = paste("PC",inx3, "(", round(100*mSetObj$analSet$pca$variance[inx3],1), "%)");
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        
        if(is.na(width)){
                w <- 9;
        }else if(width == 0){
                w <- 7.2;
        }else{
                w <- width;
        }
        h <- w;
        
        mSetObj$imgSet$pca.score3d <- imgName;
        
        Cairo::Cairo(file = imgName, unit="in", title=title, dpi=dpi, width=w, height=h, type=format, bg="white");
        
        pchs <- as.numeric(mSetObj$dataSet$cls)+1;
        uniq.pchs <- unique(pchs);
        
        if(mSetObj$dataSet$cls.type == "disc"){
                
                cols <- GetColorSchema(mSetObj$dataSet$cls);
                legend.nm <- unique(as.character(mSetObj$dataSet$cls));
                
                #cols <- GetColorSchema(mSetObj$dataSet$cls_bm);
                #legend.nm <- unique(as.character(mSetObj$dataSet$cls_bm));
                
                
                uniq.cols <- unique(cols);
                Plot3D_bm(mSetObj$analSet$pca$x[, inx1], mSetObj$analSet$pca$x[, inx2], mSetObj$analSet$pca$x[, inx3], main=title, xlab= xlabel, ylab=ylabel,
                          zlab=zlabel, angle =angl, color=cols, pch=pchs, box=F);
                legend("topleft", legend =legend.nm, pch=uniq.pchs, col=uniq.cols);
                dev.off();
                
                if(!.on.public.web){
                        # 3D View using plotly
                        if(length(uniq.pchs) > 3){
                                col <- RColorBrewer::brewer.pal(length(uniq.pchs), "Set3")
                        }else{
                                col <- c("#1972A4", "#FF7070")
                        }
                        
                        p <- plotly::plot_ly(x = mSetObj$analSet$pca$x[, inx1], y = mSetObj$analSet$pca$x[, inx2], z = mSetObj$analSet$pca$x[, inx3],
                                             color = mSetObj$dataSet$cls, colors = col) 
                        p <- plotly::add_markers(p, sizes = 5)
                        p <- plotly::layout(p, scene = list(xaxis = list(title = xlabel),
                                                            yaxis = list(title = ylabel),
                                                            zaxis = list(title = zlabel)))
                        mSetObj$imgSet$pca.3d <- p;
                        print("The Interactive 3D PCA plot has been created, please find it in mSet$imgSet$pca.3d.")
                }
                
        }else{
                
                Plot3D_bm(mSetObj$analSet$pca$x[, inx1], mSetObj$analSet$pca$x[, inx2], mSetObj$analSet$pca$x[, inx3], xlab= xlabel, ylab=ylabel,
                          zlab=zlabel, angle =angl, pch=pchs, box=F);
                dev.off();
                
                if(!.on.public.web){
                        # 3D View using plotly
                        col <- c("#C61951", "#1972A4")
                        p <- plotly::plot_ly(x = mSetObj$analSet$pca$x[, inx1], y = mSetObj$analSet$pca$x[, inx2], z = mSetObj$analSet$pca$x[, inx3],
                                             color = pchs, colors = col, marker = list(colorbar = list(len = 1, tickmode = array, tickvals = range(unique(pchs)),
                                                                                                       ticktext = levels(mSetObj$dataSet$cls)))); 
                        p <- plotly::add_markers(p, sizes = 1000);
                        p <- plotly::layout(p, scene = list(xaxis = list(title = xlabel),
                                                            yaxis = list(title = ylabel),
                                                            zaxis = list(title = zlabel)));
                        mSetObj$imgSet$pca.3d <- p;
                        print("The Interactive 3D PCA plot has been created, please find it in mSet$imgSet$pca.3d.")
                }
        }
        return(.set.mSet(mSetObj));
}


#'Create Sub Heat Map Plot
#'@description Plot a sub heatmap based on results from t-tests/ANOVA, VIP or randomforest
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@param dataOpt Set data options
#'@param scaleOpt Set the image scale
#'@param smplDist Input the sample distance method
#'@param clstDist Input the clustering distance method
#'@param palette Input color palette choice
#'@param method.nm Input the method for sub-heat map
#'@param top.num Input the top number
#'@param viewOpt Set heatmap options, default is set to "detail"
#'@param rowV Default is set to T
#'@param colV Default is set to T
#'@param border Indicate whether or not to show cell-borders, default is set to T
#'@param grp.ave Logical, default is set to F
#'@param top.FDR  FDR threshold determining which compound are significantly regulated between groups
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotSubHeatMapFDR_bm <- function(mSetObj=NA, imgName, title, format="png", dpi=72, width=NA, dataOpt, scaleOpt, 
                                 smplDist, clstDist, palette, method.nm, viewOpt, rowV=T, colV=T, border=T, grp.ave=F,top.FDR){
        
        
        ## Testing lines
        test=FALSE
        if (test == TRUE){
                mSetObj=mSet; imgName=paste0("6_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_all_samples_");
                title=paste0(nameMultiG_,"\nMetabolites vs all samples (FDR< ",topMetaboFDR__,")");
                
                mSetObj=mSet; imgName=paste0(fn,"_",nams2_[ii],"_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_all_samples_"); 
                title=paste0(nameG_,"\nMetabolites vs all samples (FDR< ",topMetaboFDR__,")");
                
                
                format="pdf"; dpi=72; width=NA; dataOpt="norm"; scaleOpt="row"; smplDist="euclidean"; clstDist="ward.D"; palette="bwm"
                method.nm="tanova"; viewOpt="overview"; rowV=F; colV=T; border=T; grp.ave=F; top.FDR=topMetaboFDR__
                
        }                           
        
        
        mSetObj <- .get.mSet(mSetObj);
        var.nms = colnames(mSetObj$dataSet$norm);
        
        
        
        #if(top.num < length(var.nms)){
        if(top.FDR <1.1){                         #* the line does nothing, is just to replate the previous line that include top.num
                if(method.nm == 'tanova'){
                        
                        if(GetGroupNumber(mSetObj) == 2){
                                if(is.null(mSetObj$analSet$tt)){
                                        Ttests.Anal(mSetObj);
                                        mSetObj <- .get.mSet(mSetObj);
                                }
                                
                                if (mSetObj$analSet$tt$sig.num==0){message("   * No significant feature found for this SubHeatmap");return(.set.mSet(mSetObj))}
                                
                                options(scipen = 999)  
                                var.nms <- rownames(mSetObj$analSet$tt$sig.mat)[which(mSetObj$analSet$tt$sig.mat[,"FDR"]<top.FDR)]
                                
                                if (length(var.nms)==0){message("   * No significant feature found for this SubHeatmap"); return(.set.mSet(mSetObj))}
                                
                                
                        }else{
                                if(is.null(mSetObj$analSet$aov)){
                                        ANOVA.Anal(mSetObj);
                                        mSetObj <- .get.mSet(mSetObj);
                                }
                                
                                if (mSetObj$analSet$aov$sig.num==0){message("   * No significant feature found for this SubHeatmap");return(.set.mSet(mSetObj))}
                                
                                options(scipen = 999)  
                                var.nms <- rownames(mSetObj$analSet$aov$sig.mat)[which(mSetObj$analSet$aov$sig.mat$FDR<top.FDR)]
                                if (length(var.nms)==0){message("   * No significant feature found for this SubHeatmap"); return(.set.mSet(mSetObj))}
                        }
                }else if(method.nm == 'cor'){
                        
                        if(is.null(mSetObj$analSet$cor.res)){
                                Match.Pattern(mSetObj);
                                mSetObj <- .get.mSet(mSetObj);
                        }
                        
                        # re-order for pretty view
                        cor.res <- mSetObj$analSet$cor.res;
                        
                        ord.inx<-order(cor.res[,3]);
                        cor.res <- cor.res[ord.inx, ];
                        
                        ord.inx<-order(cor.res[,1]);
                        cor.res <- cor.res[ord.inx, ];
                        
                        var.nms <- rownames(cor.res)[1:top.num];
                        
                }else if(method.nm == 'vip'){
                        if(is.null(mSetObj$analSet$plsda)){
                                PLSR.Anal(mSetObj);
                                PLSDA.CV(mSetObj);
                                mSetObj <- .get.mSet(mSetObj);
                        }
                        vip.vars <- mSetObj$analSet$plsda$vip.mat[,1];# use the first component
                        var.nms <- names(rev(sort(vip.vars)))[1:top.num];
                }else if(method.nm == 'rf'){
                        if(is.null(analSet$rf)){
                                RF.Anal(mSetObj);
                                mSetObj <- .get.mSet(mSetObj);
                        }
                        var.nms <- GetRFSigRowNames()[1:top.num];
                }
        }
        
        var.inx <- match(var.nms, colnames(mSetObj$dataSet$norm));
        
        PlotHeatMapStar_bm(mSetObj, imgName, title, format, dpi, width, dataOpt, scaleOpt, smplDist, clstDist, palette, viewOpt, rowV, colV, var.inx, border, grp.ave, star=FALSE, top.FDR =top.FDR);
        
}


#'Create Heat Map Plot
#'@description Plot a heatmap based on results from t-tests/ANOVA, VIP or randomforest
#'@param mSetObj Input name of the created mSet Object
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width. 
#'@param dataOpt Set data options
#'@param scaleOpt Set the image scale
#'@param smplDist Input the sample distance method
#'@param clstDist Input the clustering distance method
#'@param palette Input color palette choice
#'@param viewOpt Set heatmap options, default is set to "detail"
#'@param rowV Default is set to T
#'@param colV Default is set to T
#'@param var.inx Default is set to NA
#'@param border Indicate whether or not to show cell-borders, default is set to T
#'@param grp.ave Logical, default is set to F
#'@param top.FDR  FDR threshold determining which compound are significantly regulated between groups
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotHeatMapStar_bm <- function(mSetObj=NA, imgName, title, format="png", dpi=72, width=NA, dataOpt, scaleOpt, smplDist, 
                               clstDist, palette, viewOpt="detail", rowV=T, colV=T, var.inx=NULL, border=T, grp.ave=F, star=FALSE, top.FDR){
        
        ##testing lines 
        test=FALSE
        if (test==TRUE){
                mSetObj=mSet; imgName="testStar.pdf"; star=FALSE; var.inx=NULL; title=title=paste0("All metabolites vs group av.(* FDR<",topMetaboFDR__,")"); format="pdf"; dpi=72; width=NA; dataOpt= "norm"; scaleOpt="row"; smplDist="euclidean"; clstDist="ward.D"
                palette="bwm"; viewOpt="overview";rowV= T; colV=T; var.inx=var.inx;border=T;grp.ave=F; top.FDR=topMetaboFDR__
        }
        
        
        filenm = paste0(imgName, ".json")
        
        mSetObj <- .get.mSet(mSetObj);
        
        # record the paramters
        mSetObj$analSet$htmap <- list(dist.par=smplDist, clust.par=clstDist);
        
        # set up data set
        if(dataOpt=="norm"){
                my.data <- mSetObj$dataSet$norm;
        }else{
                my.data <- qs::qread("prenorm.qs");
        }
        
        # if(is.null(var.inx)){
        #   hc.dat<-as.matrix(my.data);
        # }else{
        #   hc.dat<-as.matrix(my.data[,var.inx]);
        # }
        
        if(length(var.inx)>1) {hc.dat<-as.matrix(my.data[,var.inx])} else { hc.dat<-as.matrix(my.data)}
        
        
        hc.dat.copy=hc.dat; #*bm
        tmpNames=colnames(hc.dat);names(tmpNames)=tmpNames   #*bm
        colnames(hc.dat) <- substr(colnames(hc.dat),1,18) # some names are too long
        tmpNames2=colnames(hc.dat); names(tmpNames2)=tmpNames2 #*bm
        
        
        if(mSetObj$dataSet$type.cls.lbl=="integer"){
                hc.cls <- as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls]);
        }else{
                hc.cls <- mSetObj$dataSet$cls;
                #hc.cls <- mSetObj$dataSet$cls_bm;
        }
        
        if(grp.ave){ # only use group average
                lvs <- levels(hc.cls);
                my.mns <- matrix(ncol=ncol(hc.dat),nrow=length(lvs));
                for(i in 1:length(lvs)){
                        inx <-hc.cls == lvs[i];
                        my.mns[i,]<- apply(hc.dat[inx, ], 2, mean);
                }
                rownames(my.mns) <- lvs;
                colnames(my.mns) <- colnames(hc.dat);
                hc.dat <- my.mns;
                hc.cls <- as.factor(lvs);
        }
        
        # set up colors for heatmap
        if(palette=="gbr"){
                colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
        }else if(palette == "heat"){
                colors <- heat.colors(256);
        }else if(palette == "topo"){
                colors <- topo.colors(256);
        }else if(palette == "gray"){
                colors <- colorRampPalette(c("grey90", "grey10"), space="rgb")(256);
        }else{
                colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256));
        }
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        
        if(is.na(width)){
                minW <- 630;
                myW <- nrow(hc.dat)*18 + 150;
                
                if(myW < minW){
                        myW <- minW;
                }   
                w <- round(myW/72,2);
        }else if(width == 0){
                w <- 7.2;
        }else{
                w <- 7.2;
        }
        
        mSetObj$imgSet$heatmap <- imgName;
        
        myH <- ncol(hc.dat)*18 + 150;
        h <- round(myH/72,2);
        
        if(viewOpt == "overview"){
                if(is.na(width)){
                        if(w > 9){
                                w <- 9;
                        }
                }else if(width == 0){
                        if(w > 7.2){
                                w <- 7.2;
                        }
                        
                }else{
                        w <- 7.2;
                }
                if(h > w){
                        h <- w;
                }
                
                mSetObj$imgSet$heatmap <- imgName;
        }
        
        # make the width smaller fro group average
        if(grp.ave){
                w <- nrow(hc.dat)*25 + 300;
                w <- round(w/72,2);
        }
        
        if(border){
                border.col<-"grey60";
        }else{
                border.col <- NA;
        }
        if(format=="pdf"){
                pdf(file = imgName, width=w, height=h, bg="white", onefile=FALSE);
        }else{
                Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        }
        if(mSetObj$dataSet$cls.type == "disc"){
                annotation <- data.frame(class = hc.cls);
                rownames(annotation) <- rownames(hc.dat); 
                
                
                ### add a start on the names
                if (star==TRUE){
                        
                        if(!is.null(mSetObj$analSet$tt)){
                                
                                newNames= rep("",length(tmpNames))       
                                for (i in 1:length(tmpNames)){       
                                        cl_=tmpNames[i]
                                        if (cl_ %in% rownames(mSetObj$analSet$tt$sig.mat)){
                                                ix=which(rownames(mSetObj$analSet$tt$sig.mat)==cl_)
                                                newNames[i]=tmpNames2[i]
                                                if (length(ix)>0){
                                                        if (mSetObj$analSet$tt$sig.mat[ix,"FDR"] < top.FDR){
                                                                newNames[i]=paste0(tmpNames2[i],"*")
                                                        }
                                                }
                                        }else{
                                                newNames[i]=tmpNames2[i]
                                        }
                                }
                                
                        }
                        
                        if(!is.null(mSetObj$analSet$aov)){
                                
                                newNames= rep("",length(tmpNames))       
                                for (i in 1:length(tmpNames)){  
                                        cl_=tmpNames[i]
                                        if (cl_ %in% rownames(mSetObj$analSet$aov$sig.mat)){
                                                ix=which(rownames(mSetObj$analSet$aov$sig.mat)==cl_)
                                                newNames[i]=tmpNames2[i]
                                                if (length(ix)>0){
                                                        if (mSetObj$analSet$aov$sig.mat$FDR[ix] < top.FDR){
                                                                newNames[i]=paste0(tmpNames2[i],"*")
                                                        }
                                                }
                                        }else{
                                                newNames[i]=tmpNames2[i]
                                        }
                                }
                        }
                        #View(colnames(hc.dat))
                        #View(newNames)
                        
                        colnames(hc.dat)=newNames
                }    
                
                # set up color schema for samples
                cols <- GetColorSchema(mSetObj$dataSet$cls, palette == "gray");
                uniq.cols <- unique(cols);
                
                if(mSetObj$dataSet$type.cls.lbl=="integer"){
                        cls <- as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls]);
                }else{
                        cls <- mSetObj$dataSet$cls;
                        #cls <- mSetObj$dataSet$cls_bm;
                }
                
                names(uniq.cols) <- unique(as.character(sort(cls)));
                
                ann_colors <- list(class= uniq.cols);
                
                ph<-pheatmap::pheatmap(t(hc.dat), 
                                       annotation=annotation, 
                                       main=title,
                                       fontsize=8, fontsize_row=8, 
                                       clustering_distance_rows = smplDist,
                                       clustering_distance_cols = smplDist,
                                       clustering_method = clstDist, 
                                       border_color = border.col,
                                       cluster_rows = colV, 
                                       cluster_cols = rowV,
                                       scale = scaleOpt, 
                                       color = colors,
                                       annotation_colors = ann_colors);
                
                
                
                if (endsWith(imgName,".pdf")){
                        write.csv(t(hc.dat.copy)[ph$tree_row$order,],file=paste0(strsplit(imgName,".pdf")[[1]][1],".csv"))
                }
                if (endsWith(imgName,".png")){
                        write.csv(t(hc.dat.copy)[ph$tree_row$order,],file=paste0(strsplit(imgName,".png")[[1]][1] ,".csv"))
                }
                
                dat = t(hc.dat)
                if(scaleOpt == "row"){
                        res <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
                }else{
                        res <- t(apply(dat, 2, function(x){as.numeric(cut(x, breaks=30))}));
                }
                colnames(dat) = NULL
                netData <- list(data=res, annotation=annotation, smp.nms = colnames(t(hc.dat)), met.nms = rownames(t(hc.dat)), colors = colors);
                sink(filenm);
                cat(RJSONIO::toJSON(netData));
                sink();
        }else{
                heatmap(hc.dat, Rowv = rowTree, Colv=colTree, col = colors, scale="column");
                dat = t(hc.dat)
                if(scaleOpt == "row"){
                        res <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
                }else{
                        res <- t(apply(dat, 2, function(x){as.numeric(cut(x, breaks=30))}));
                }
                colnames(dat) = NULL
                netData <- list(data=res, annotation="NA", smp.nms = colnames(t(hc.dat)), met.nms = rownames(t(hc.dat)), colors = colors);
                sink(filenm);
                cat(RJSONIO::toJSON(netData));
                sink();
        }
        dev.off();
        return(.set.mSet(mSetObj));
}


#'Plot QEA overview
#'@description Plot QEA overview
#'@usage PlotQEA.Overview(mSetObj=NA, imgName, imgOpt, format="png", dpi=72, width=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param title The title to add on top of the graph
#'@param imgOpt "net" 
#'@param format Select the image format, "png", or "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'@param top.FDR  FDR threshold determining which pathways are significantly regulated in two groups 
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotQEA.Overview_bm<-function(mSetObj=NA, imgName, tittle=tittle, imgOpt, format="png", dpi=72, width=NA, top.FDR){
        
        # Testing line
        #mSetObj=mSetMET; imgName="12_metabolite_set_qea_";imgOpt="net"; format="pdf"; dpi=72;  width=NA; top.FDR=topMetaboFDR__
        
        
        mSetObj <- .get.mSet(mSetObj);
        
        #calculate the enrichment fold change
        folds <- mSetObj$analSet$qea.mat[,3]/mSetObj$analSet$qea.mat[,4];
        names(folds) <- GetShortNames(rownames(mSetObj$analSet$qea.mat));
        pvals <- mSetObj$analSet$qea.mat[,6];
        fdrs<- mSetObj$analSet$qea.mat[,7];
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        
        if(is.na(width)){
                w <- 9;
        }else if(width == 0){
                w <- 7;
        }else{
                w <- width;
        }
        h <- w;
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        PlotMSEA.Overview_bm(folds, pvals,fdrs,top.FDR,tittle=tittle);
        dev.off();
        
        if(.on.public.web){
                mSetObj$analSet$enrich.net <- PlotEnrichNet.Overview(folds, pvals);
        }
        
        mSetObj$imgSet$qea <-imgName;
        mSetObj$imgSet$current.img <- imgName;
        return(.set.mSet(mSetObj));
}

#'Plot MSEA overview
#'@description Barplot height is enrichment fold change
#'color is based on p values, used in higher functions
#'@usage PlotMSEA.Overview(folds, pvals)
#'@param folds Input the fold-change values
#'@param pvals Input the p-values
#'@param tittle tittle for the plot
#'@param top.FDR  FDR threshold determining which pathways are significantly regulated in two groups 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotMSEA.Overview_bm<-function(folds, pvals,fdrs,top.FDR,tittle=""){
        
        ## Test lines
        # folds=folds; pvals=pvals; fdrs=fdrs; top.FDR=top.FDR;
        
        # due to space limitation, plot top 50 if more than 50 were given
        title <- paste0("MSE-",tittle,"\n(* FDR <",top.FDR,")");
        if(length(folds) > 50){
                folds <- folds[1:50];
                pvals <- pvals[1:50];
                fdrs  <- fdrs[1:50];
                title <- paste0("MSE-",tittle ," (top 50)\n(* FDR < ",top.FDR,")");
        }
        
        op <- par(mar=c(5,20,4,6), oma=c(0,0,0,4));
        ht.col <- rev(heat.colors(length(folds)));
        
        br=barplot(rev(folds), horiz=T, col=ht.col,
                   xlab="Fold Enrichment", las=1, cex.name=0.75, space=c(0.5, 0.5),
                   main= title);
        
        starsIx=which(format(rev(fdrs),scientific = F)<top.FDR)
        
        for (ix in starsIx){
                #text(rev(folds)[ix]+.02, br[ix,1],labels="*")
                text(max(folds)/2, br[ix,1],labels="*")
        }
        
        minP <- min(pvals);
        maxP <- max(pvals);
        medP <- (minP+maxP)/2;
        
        axs.args <- list(at=c(minP, medP, maxP), labels=format(c(maxP, medP, minP), scientific=T, digit=1), tick = F);
        image.plot(legend.only=TRUE, zlim=c(minP, maxP), col=ht.col,
                   axis.args=axs.args, legend.shrink = 0.4, legend.lab="P value");
        par(op);
}

#'Plot over-representation analysis (ORA)
#'@description Plot over-representation analysis (ORA)
#'@usage PlotORA(mSetObj=NA, imgName, imgOpt, format="png", dpi=72, width=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param tittle tittle for the plot
#'@param imgOpt "net" 
#'@param format Select the image format, "png", or "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width. 
#'@param top.FDR  FDR threshold determining which pathways are significantly regulated in two groups  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotORA_bm<-function(mSetObj=NA, imgName,tittle="" , imgOpt, format="png", dpi=72, width=NA,top.FDR){
        
        mSetObj <- .get.mSet(mSetObj);
        
        #calculate the enrichment fold change
        folds <- mSetObj$analSet$ora.mat[,3]/mSetObj$analSet$ora.mat[,2];
        names(folds) <- GetShortNames(rownames(mSetObj$analSet$ora.mat));
        pvals <- mSetObj$analSet$ora.mat[,4];
        fdrs<- mSetObj$analSet$ora.mat[,6];
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        
        if(is.na(width)){
                w <- 9;
        }else if(width == 0){
                w <- 7;
        }else{
                w <-width;
        }
        h <- w;
        
        #record img
        mSetObj$imgSet$ora <- imgName
        mSetObj$imgSet$current.img <- imgName;
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        
        PlotMSEA.Overview_bm(folds, pvals=pvals,fdrs=fdrs,top.FDR,tittle=tittle);
        dev.off();
        
        if(.on.public.web){
                mSetObj$analSet$enrich.net <- PlotEnrichNet.Overview(folds, pvals);
        }
        
        return(.set.mSet(mSetObj));
}

#'Plot MSEA Dot Plot
#'@description Dot plot of enrichment analysis results.
#'@usage PlotEnrichDotPlot(mSet)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param enrichType Input whether the enrichment analysis was over-respresentation
#'analysis (ora) or quantitative enrichment analysis (qea).
#'@param imgName Input a name for the plot
#'@param tittle tittle for the plot
#'@param format Select the image format, "png", or "pdf". 
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images, 
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.  
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotEnrichDotPlot_bm <- function(mSetObj=NA, enrichType = "ora", imgName, tittle="", format="png", dpi=72, width=NA){
        
        mSetObj <- .get.mSet(mSetObj);
        
        if(.on.public.web){
                load_ggplot()
        }
        
        if(enrichType == "ora"){
                results <- mSetObj$analSet$ora.mat
                my.cols <- GetMyHeatCols(nrow(results));
                if(nrow(results) > 25){
                        results <- results[1:25,]
                        my.cols <- my.cols[1:25];
                }
                
                df <- data.frame(Name = factor(row.names(results), levels = rev(row.names(results))),
                                 rawp = results[,4],
                                 logp = -log10(results[,4]),
                                 folds = results[,3]/results[,2])
                
        }else if(enrichType == "qea"){
                results <- mSetObj$analSet$qea.mat
                my.cols <- GetMyHeatCols(nrow(results));   
                if(nrow(results) > 25){
                        results <- results[1:25,]
                        my.cols <- my.cols[1:25];
                }
                df <- data.frame(Name = factor(row.names(results), levels = rev(row.names(results))),
                                 rawp = results[,5],
                                 logp = -log10(results[,5]),
                                 folds = results[,3]/results[,4]
                )
                
        }
        
        maxp <- max(df$rawp);
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        
        if(is.na(width)){
                w <- 12;
                h <- 9;
        }else if(width == 0){
                h <- w <- 7;
        }else{
                h <- w <- width;
        }
        
        p <- ggplot(df, 
                    aes(x = logp, y = Name)) + 
                geom_point(aes(size = folds, color = rawp)) + scale_size_continuous(range = c(2, 8)) +
                theme_bw(base_size = 14.5) +
                scale_colour_gradient(limits=c(0, maxp), low=my.cols[1], high = my.cols[length(my.cols)]) +
                ylab(NULL) + xlab("-log10 (p-value)") + 
                ggtitle(paste0("Overview of Enriched Metabolite Sets\n",tittle," (Top 25)")) +
                theme(legend.text=element_text(size=14),
                      legend.title=element_text(size=15))
        
        p$labels$colour <- "P-value"
        p$labels$size <- "Enrichment Ratio"
        
        ggsave(p, filename = imgName, dpi=dpi, width=w, height=h)
        
        if(enrichType == "ora"){
                mSetObj$imgSet$ora <- imgName
        }else if(enrichType == "qea"){
                mSetObj$imgSet$qea <-imgName;
        }
        
        mSetObj$imgSet$current.img <- imgName;
        return(.set.mSet(mSetObj));
}


#'Create volcano plot
#'@description For labelling interesting points, it is defined by the following rules:
#'need to be signficant (sig.inx) and or 2. top 5 p, or 2. top 5 left, or 3. top 5 right.
#'@usage PlotVolcano(mSetObj=NA, imgName, plotLbl, format="png", dpi=72, width=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param plotLbl Logical, plot labels, 1 for yes and 0 for no.
#'@param tittle tittle for the plot
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images,
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotVolcano_bm <- function(mSetObj=NA, imgName, plotLbl, title, format="png", dpi=72, width=NA){
        
        mSetObj <- .get.mSet(mSetObj);
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        
        if(is.na(width)){
                w <- 10;
        }else if(width == 0){
                w <- 8;
        }else{
                w <- width;
        }
        h <- w*6/10;
        
        mSetObj$imgSet$volcano <- imgName;
        
        Cairo::Cairo(file = imgName, unit="in", title=title, dpi=dpi, width=w, height=h, type=format, bg="white");
        par(mar=c(5,5,3,4));
        vcn <- mSetObj$analSet$volcano;
        MyGray <- rgb(t(col2rgb("black")), alpha=40, maxColorValue=255);
        MyHighlight <- rgb(t(col2rgb("magenta")), alpha=80, maxColorValue=255);
        
        if(vcn$paired){
                xlim<-c(-nrow(mSetObj$dataSet$norm)/2, nrow(mSetObj$dataSet$norm)/2)*1.2;
                
                # merge fc.all two rows into one, bigger one win
                fc.all <- apply(vcn$fc.all, 2, function(x){ if(x[1] > x[2]){return(x[1])}else{return(-x[2])}})
                
                hit.inx <- vcn$inx.p & (vcn$inx.up | vcn$inx.down);
                plot(fc.all, vcn$p.log, main=title, xlim=xlim, pch=20, cex=ifelse(hit.inx, 1.2, 0.8),
                     col = ifelse(hit.inx, MyHighlight, MyGray),
                     xlab="Count of Significant Pairs", ylab="-log10(p)");
                
                sig.upInx <- vcn$inx.p & vcn$inx.up;
                p.topInx <- GetTopInx(vcn$p.log, 5, T) & vcn$inx.up;
                fc.rtInx <- GetTopInx(vcn$fc.all[1,], 5, T);
                lblInx <- p.topInx & sig.upInx & fc.rtInx;
                if(plotLbl & sum(lblInx, na.rm=T) > 0){
                        text.lbls<-substr(colnames(mSetObj$dataSet$norm)[lblInx],1,14) # some names may be too long
                        text(vcn$fc.all[1,lblInx], vcn$p.log[lblInx],labels=text.lbls, pos=4, col="blue", srt=30, xpd=T, cex=0.8);
                }
                
                sig.dnInx <- vcn$inx.p & vcn$inx.down;
                p.topInx <- GetTopInx(vcn$p.log, 5, T) & vcn$inx.down;
                fc.leftInx <- GetTopInx(vcn$fc.all[2,], 5, T) & vcn$inx.down;
                lblInx <-p.topInx & sig.dnInx & fc.leftInx;
                if(plotLbl & sum(lblInx, na.rm=T) > 0){
                        text.lbls<-substr(colnames(mSetObj$dataSet$norm)[lblInx],1,14) # some names may be too long
                        text(-vcn$fc.all[2,lblInx], vcn$p.log[lblInx],labels=text.lbls, pos=2, col="blue", srt=-30, xpd=T, cex=0.8);
                }
                
        }else{
                imp.inx<-(vcn$inx.up | vcn$inx.down) & vcn$inx.p;
                plot(vcn$fc.log, vcn$p.log, main=title, pch=20, cex=ifelse(imp.inx, 1.2, 0.7),
                     col = ifelse(imp.inx, MyHighlight, MyGray),
                     xlab="log2 (FC)", ylab="-log10(p)");
                
                sig.inx <- imp.inx;
                p.topInx <- GetTopInx(vcn$p.log, 5, T) & (vcn$inx.down);
                fc.leftInx <- GetTopInx(vcn$fc.log, 5, F);
                lblInx <-  sig.inx & (p.topInx | fc.leftInx);
                if(plotLbl &  sum(lblInx, na.rm=T) > 0){
                        text.lbls<-substr(colnames(mSetObj$dataSet$norm)[lblInx],1,14) # some names may be too long
                        text(vcn$fc.log[lblInx], vcn$p.log[lblInx],labels=text.lbls, pos=2, col="blue", srt=-30, xpd=T, cex=0.8);
                }
                
                p.topInx <- GetTopInx(vcn$p.log, 5, T) & (vcn$inx.up);
                fc.rtInx <- GetTopInx(vcn$fc.log, 5, T);
                lblInx <- sig.inx & (p.topInx | fc.rtInx);
                if(plotLbl & sum(lblInx, na.rm=T) > 0){
                        text.lbls<-substr(colnames(mSetObj$dataSet$norm)[lblInx],1,14) # some names may be too long
                        text(vcn$fc.log[lblInx], vcn$p.log[lblInx],labels=text.lbls, pos=4, col="blue", srt=30, xpd=T, cex=0.8);
                }
        }
        
        abline (v = vcn$max.xthresh, lty=3);
        abline (v = vcn$min.xthresh, lty=3);
        abline (h = vcn$thresh.y, lty=3);
        axis(4); # added by Beomsoo
        dev.off();
        return(.set.mSet(mSetObj));
}


#'Plot fold change
#'@description Plot fold change analysis
#'@usage PlotFC(mSetObj=NA, imgName, format="png", dpi=72, width=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param tittle tittle for the plot
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images,
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotFC_bm <- function(mSetObj=NA, imgName,title ,format="png", dpi=72, width=NA){
        
        mSetObj <- .get.mSet(mSetObj);
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 8;
        }else if(width == 0){
                w <- 7;
        }else{
                w <- width;
        }
        h <- w*6/8;
        
        mSetObj$imgSet$fc <- imgName;
        
        Cairo::Cairo(file = imgName, unit="in", title=title, dpi=dpi, width=w, height=h, type=format, bg="white");
        
        par(mar=c(5,5,2,3));
        
        fc = mSetObj$analSet$fc;
        if(fc$paired){
                ylim <- c(-nrow(mSetObj$dataSet$norm)/2, nrow(mSetObj$dataSet$norm)/2);
                xlim <- c(0, ncol(mSetObj$dataSet$norm));
                plot(NULL, main=title,xlim=xlim, ylim=ylim, xlab = GetVariableLabel(mSetObj$dataSet$type),
                     ylab=paste("Count with FC >=", fc$max.thresh, "or <=", fc$min.thresh));
                for(i in 1:ncol(fc$fc.all)){
                        segments(i,0, i, fc$fc.all[1,i], col= ifelse(fc$inx.up[i],"magenta", "darkgrey"),
                                 lwd= ifelse(fc$inx.up[i], 2, 1));
                        segments(i,0, i, -fc$fc.all[2,i], col= ifelse(fc$inx.down[i], "magenta", "darkgrey"),
                                 lwd= ifelse(fc$inx.down[i], 2, 1));
                }
                abline(h=fc$max.thresh, lty=3);
                abline(h=fc$min.thresh, lty=3);
                abline(h=0, lwd=1);
        }else{
                if(fc$raw.thresh > 0){
                        # be symmetrical
                        topVal <- max(abs(fc$fc.log));
                        ylim <- c(-topVal, topVal);
                        plot(fc$fc.log, main=title, ylab="Log2 (FC)", ylim = ylim, xlab = GetVariableLabel(mSetObj$dataSet$type), pch=19, axes=F,
                             col= ifelse(fc$inx.imp, "magenta", "darkgrey"));
                        axis(2);
                        axis(4); # added by Beomsoo
                        abline(h=log(fc$max.thresh,2), lty=3);
                        abline(h=log(fc$min.thresh,2), lty=3);
                        abline(h=0, lwd=1);
                }else{ # plot side by side
                        
                        dat1 <- mSetObj$dataSet$norm[as.numeric(mSetObj$dataSet$cls) == 1, ];
                        dat2 <- mSetObj$dataSet$norm[as.numeric(mSetObj$dataSet$cls) == 2, ];
                        
                        mns1 <- apply(dat1, 2, mean);
                        mn1 <- mean(mns1);
                        sd1 <- sd(mns1);
                        msd1.top <- mn1 + 2*sd1;
                        msd1.low <- mn1 - 2*sd1;
                        
                        mns2 <- apply(dat2, 2, mean);
                        mn2 <- mean(mns2);
                        sd2 <- sd(mns2);
                        msd2.top <- mn2 + 2*sd2;
                        msd2.low <- mn2 - 2*sd2;
                        
                        ylims <- range(c(mns1, mns2, msd1.top, msd2.top, msd1.low, msd2.low));
                        new.mns <- c(mns1, rep(NA, 5), mns2);
                        cols <- c(rep("magenta", length(mns1)), rep(NA, 5), rep("blue", length(mns2)));
                        pchs <- c(rep(15, length(mns1)), rep(NA, 5), rep(19, length(mns2)));
                        plot(new.mns, main=title, ylim=ylims, pch = pchs, col = cols, cex = 1.25, axes=F, ylab="");
                        axis(2);
                        axis(4); # added by Beomsoo
                        abline(h=mn1, col="magenta", lty=3, lwd=2);
                        abline(h=msd1.low, col="magenta", lty=3, lwd=1);
                        abline(h=msd1.top, col="magenta", lty=3, lwd=1);
                        abline(h=mn2, col="blue", lty=3, lwd=2);
                        abline(h=msd2.low, col="blue", lty=3, lwd=1);
                        abline(h=msd2.top, col="blue", lty=3, lwd=1);
                        # abline(h=mean(all.mns), col="darkgrey", lty=3);
                        axis(1, at=1:length(new.mns), labels=c(1:length(mns1),rep(NA, 5),1:length(mns2)));
                }
        }
        dev.off();
        
        return(.set.mSet(mSetObj));
}

#'Plot t-test
#'@description Plot t-test
#'@usage PlotTT(mSetObj=NA, imgName, format="png", dpi=72, width=NA)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param imgName Input a name for the plot
#'@param tittle tittle for the plot
#'@param format Select the image format, "png", or "pdf".
#'@param dpi Input the dpi. If the image format is "pdf", users need not define the dpi. For "png" images,
#'the default dpi is 72. It is suggested that for high-resolution images, select a dpi of 300.
#'@param width Input the width, there are 2 default widths, the first, width = NULL, is 10.5.
#'The second default is width = 0, where the width is 7.2. Otherwise users can input their own width.
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotTT_bm <- function(mSetObj=NA, imgName,title, format="png", dpi=72, width=NA){
        
        mSetObj <- .get.mSet(mSetObj);
        
        imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
        if(is.na(width)){
                w <- 8;
        }else if(width == 0){
                w <- 7;
        }else{
                w <- width;
        }
        h <- w*6/8;
        
        mSetObj$imgSet$tt <- imgName;
        
        Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
        plot(mSetObj$analSet$tt$p.log, ylab="-log10(p)", xlab=GetVariableLabel(mSetObj$dataSet$type), main=paste0( title ," - ",mSetObj$analSet$tt$tt.nm ), pch=19,
             col= ifelse(mSetObj$analSet$tt$inx.imp, "magenta", "darkgrey"));
        abline(h=mSetObj$analSet$tt$thresh, lty=3);
        axis(4);
        dev.off();
        return(.set.mSet(mSetObj));
}

#'Set current user selected metset library for search
#'@description if enrichment analysis, also prepare lib by
#'creating a list of metabolite sets
#'@usage SetCurrentMsetLib(mSetObj=NA, lib.type, excludeNum)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param lib.type Input user selected name of library, "self", "kegg_pathway",
#'"smpdb_pathway", "blood", "urine", "csf", "snp", "predicted", "location", and "drug".
#'@param excludeNum Users input the mimimum number compounds within selected metabolite sets (metabolitesets < excludeNum)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#'@export
#'
SetCurrentMsetLib_bm <- function(mSetObj=NA, libname, excludeNum=0){
        
        mSetObj <- .get.mSet(mSetObj);
        
        mSetObj$analSet$msetlibname <- libname;
        
        if(libname=="self"){
                ms.list <- mSetObj$dataSet$user.mset;
                ms.list <- lapply(ms.list, function(x) unlist(strsplit(x, "; ", fixed=TRUE)))
                current.msetlib <- vector("list", 3)
                names(current.msetlib) <- c("name", "member", "reference")
        }else{
                
                if(!.on.public.web & grepl("kegg", libname)){ # api only for KEGG msets
                        mSetObj$api$libname <- libname
                        mSetObj$api$excludeNum = excludeNum
                        return(.set.mSet(mSetObj));
                }
                
                if(!exists("current.msetlib") || mSetObj$analSet$msetlibname != libname) {
                        destfile <- paste(libname, ".qs", sep = "");
                        if(.on.public.web){
                                my.qs  <- paste("../../libs/msets/", destfile, sep="");
                                current.msetlib <- qs::qread(my.qs);
                        }else{
                                my.qs <- paste("https://www.metaboanalyst.ca/resources/libs/msets/", destfile, sep="");
                                if(!file.exists(destfile)){
                                        download.file(my.qs, destfile);
                                }
                                current.msetlib <- qs::qread(destfile);
                        }
                        mSetObj$analSet$msetlibname <- libname;
                }
                # create a named list, use the ids for list names
                ms.list <- strsplit(current.msetlib[,3],"; ", fixed=TRUE);
                names(ms.list) <- current.msetlib[,2];
        }
        
        if(excludeNum > 0){
                cmpd.count <- lapply(ms.list, length);
                sel.inx <- cmpd.count >= excludeNum;
                ms.list <- ms.list[sel.inx];
                
                if(libname!="self"){
                        current.msetlib <- current.msetlib[sel.inx,];
                }
        }
        
        # total uniq cmpds in the mset lib
        mSetObj$dataSet$uniq.count <- length(unique(unlist(ms.list, use.names = FALSE)));
        
        # update current.mset and push to global env
        current.msetlib$member <- ms.list;
        
        if(libname=="self"){
                current.msetlib$name <- names(ms.list)
                current.msetlib$reference <- rep("User-uploaded", length(ms.list))
        }
        
        current.msetlib <<- current.msetlib;
        return(.set.mSet(mSetObj));
}

#'Perform ANOVA analysis
#'@description ANOVA analysis
#'@usage ANOVA.Anal(mSetObj=NA, nonpar=F, thresh=0.05, post.hoc="fisher")
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param nonpar Logical, use a non-parametric test (T) or not (F)
#'@param thresh Numeric, from 0 to 1, indicate the p-value threshold
#'@param post.hoc Input the name of the post-hoc test, "fisher" or "tukey"
#'@param all_results Logical, if TRUE, it will output the ANOVA results for all compounds 
#'with no post-hoc tests performed.
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
ANOVA.Anal_bm<-function(mSetObj=NA, nonpar=F, thresh=0.05, post.hoc="fisher", all_results=FALSE){
        
        mSetObj <- .get.mSet(mSetObj);
        
        sig.num <- 0;
        if(nonpar){
                aov.nm <- "Kruskal Wallis Test";
                anova.res <- apply(as.matrix(mSetObj$dataSet$norm), 2, kwtest, cls=mSetObj$dataSet$cls);
                
                #extract all p values
                res <- unlist(lapply(anova.res, function(x) {c(x$statistic, x$p.value)}));
                res <- data.frame(matrix(res, nrow=length(anova.res), byrow=T), stringsAsFactors=FALSE);
                
                fstat <- res[,1];
                p.value <- res[,2];
                
                names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
                fdr.p <- p.adjust(p.value, "fdr");
                
                #inx.imp <- p.value <= thresh;
                inx.imp <- fdr.p <= thresh;
                sig.num <- sum(inx.imp);
                
                if(sig.num > 0){
                        sig.f <- fstat[inx.imp];
                        sig.p <- p.value[inx.imp];
                        fdr.p <- fdr.p[inx.imp];
                        
                        sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(fdr.p,5), 'NA');
                        rownames(sig.mat) <- names(sig.p);
                        colnames(sig.mat) <- c("chi.squared", "p.value", "-log10(p)", "FDR", "Post-Hoc");
                        
                        # order the result simultaneously
                        ord.inx <- order(sig.p, decreasing = FALSE);
                        sig.mat <- sig.mat[ord.inx,,drop=F];
                        
                        fileName <- "kw_posthoc.csv";
                        my.mat <- sig.mat[,1:4];
                        colnames(my.mat) <- c("chi_squared", "pval_KW", "-log10(p)", "FDR");
                }
        }else{
                aov.nm <- "One-way ANOVA";
                if(.on.public.web & RequireFastUnivTests(mSetObj)){
                        res <- PerformFastUnivTests(mSetObj$dataSet$norm, mSetObj$dataSet$cls);
                }else{
                        aov.res <- apply(as.matrix(mSetObj$dataSet$norm), 2, aof, cls=mSetObj$dataSet$cls);
                        anova.res <- lapply(aov.res, anova);
                        
                        #extract all p values
                        res <- unlist(lapply(anova.res, function(x) { c(x["F value"][1,], x["Pr(>F)"][1,])}));
                        res <- data.frame(matrix(res, nrow=length(aov.res), byrow=T), stringsAsFactors=FALSE);
                }
                fstat <- res[,1];
                p.value <- res[,2];
                names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
                
                fdr.p <- p.adjust(p.value, "fdr");
                
                if(all_results==TRUE){
                        all.mat <- data.frame(signif(p.value,5), signif(-log10(p.value),5), signif(fdr.p,5));
                        rownames(all.mat) <- names(p.value);
                        colnames(all.mat) <- c("p.value", "-log10(p)", "FDR");
                        o=order(all.mat$p.value)
                        write.csv(all.mat[o,], "anova_all.csv")
                }
                
                # do post-hoc only for signficant entries
                # inx.imp <- p.value <= thresh;
                inx.imp <- fdr.p <= thresh;
                sig.num <- sum(inx.imp);
                if(sig.num > 0){
                        # note aov obj is not avaible using fast version
                        # need to recompute using slower version for the sig ones
                        if(.on.public.web & RequireFastUnivTests(mSetObj)){
                                aov.imp <- apply(as.matrix(mSetObj$dataSet$norm[,inx.imp,drop=FALSE]), 2, aof, cls=mSetObj$dataSet$cls);
                        }else{
                                aov.imp <- aov.res[inx.imp];
                        }
                        sig.f <- fstat[inx.imp];
                        sig.p <- p.value[inx.imp];
                        fdr.p <- fdr.p[inx.imp];
                        cmp.res <- NULL;
                        post.nm <- NULL;
                        if(post.hoc=="tukey"){
                                tukey.res<-lapply(aov.imp, TukeyHSD, conf.level=1-thresh);
                                cmp.res <- unlist(lapply(tukey.res, parseTukey, cut.off=thresh));
                                post.nm = "Tukey's HSD";
                        }else{
                                fisher.res<-lapply(aov.imp, FisherLSD, thresh);
                                cmp.res <- unlist(lapply(fisher.res, parseFisher, cut.off=thresh));
                                post.nm = "Fisher's LSD";
                        }
                        
                        # create the result dataframe,
                        # note, the last column is string, not double
                        
                        sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(fdr.p,5), cmp.res);
                        rownames(sig.mat) <- names(sig.p);
                        colnames(sig.mat) <- c("f.value", "p.value", "-log10(p)", "FDR", post.nm);
                        
                        # order the result simultaneously
                        ord.inx <- order(sig.p, decreasing = FALSE);
                        sig.mat <- sig.mat[ord.inx,,drop=F];
                        fileName <- "anova_posthoc.csv";
                }
        }
        
        AddMsg(paste(c("A total of", sum(inx.imp), "significant features were found."), collapse=" "));
        if(sig.num> 0){
                res <- 1;
                write.csv(sig.mat,file=fileName);
                aov<-list (
                        aov.nm = aov.nm,
                        sig.num = sig.num,
                        sig.nm = fileName,
                        raw.thresh = thresh,
                        thresh = -log10(thresh), # only used for plot threshold line
                        p.value = p.value,
                        p.log = -log10(p.value),
                        inx.imp = inx.imp,
                        post.hoc = post.hoc,
                        sig.mat = sig.mat
                );
        }else{
                res <- 0;
                aov<-list (
                        aov.nm = aov.nm,
                        sig.num = sig.num,
                        raw.thresh = thresh,
                        thresh = -log10(thresh), # only used for plot threshold line
                        p.value = p.value,
                        p.log = -log10(p.value),
                        inx.imp = inx.imp
                );
        }
        mSetObj$analSet$aov <- aov;
        
        if(.on.public.web){
                .set.mSet(mSetObj);
                return(res);
        }else{
                return(.set.mSet(mSetObj));
        }
}



#'Fold change analysis, unpaired
#'@description Perform fold change analysis, method can be mean or median
#'@usage FC.Anal(mSetObj, fc.thresh=2, cmp.type = 0)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param fc.thresh Fold-change threshold, numeric input
#'@param cmp.type Comparison type, 0 for group 1 minus group 2, and 1 for group 
#'1 minus group 2
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
FC.Anal.unpaired_bm<-function(mSetObj=NA, fc.thresh=2, cmp.type = 0){
        
        mSetObj <- .get.mSet(mSetObj);
        
        # make sure threshold is above 1
        fc.thresh = ifelse(fc.thresh>1, fc.thresh, 1/fc.thresh);
        max.thresh = fc.thresh;
        min.thresh = 1/fc.thresh;
        
        res <- GetFC(mSetObj, F, cmp.type);
        fc.all <- res$fc.all;
        fc.log <- res$fc.log;
        
        imp.inx <- fc.all > max.thresh | fc.all < min.thresh;
        sig.mat <- cbind(fc.all[imp.inx, drop=F], fc.log[imp.inx, drop=F]);
        colnames(sig.mat) <- c("Fold Change", "log2(FC)");
        
        # order by absolute log value (since symmetrical in pos and neg)
        inx.ord <- order(abs(sig.mat[,2]), decreasing=T);
        sig.mat <- sig.mat[inx.ord,,drop=F];
        
        fileName <- "fold_change.csv";
        write.csv(sig.mat,file=fileName);
        all.mat=cbind(fc.all, fc.log); colnames(all.mat) <- c("Fold Change", "log2(FC)"); o=order(all.mat[,1])
        fileName <- "fold_change_all.csv";
        write.csv(all.mat[o,],file=fileName);
        
        # create a list object to store fc
        mSetObj$analSet$fc<-list (
                paired = FALSE,
                raw.thresh = fc.thresh,
                max.thresh = max.thresh,
                min.thresh = min.thresh,
                fc.all = fc.all, # note a vector
                fc.log = fc.log,
                inx.imp = imp.inx,
                sig.mat = sig.mat
        );
        return(.set.mSet(mSetObj));
}

#'Perform Volcano Analysis
#'@description Perform volcano analysis
#'@usage Volcano.Anal(mSetObj=NA, paired=FALSE, fcthresh, cmpType, nonpar=F, threshp, equal.var=TRUE, pval.type="raw")
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param paired Logical, T if data is paired, F if data is not.
#'@param fcthresh Numeric, input the fold change threshold
#'@param cmpType Comparison type, 0 indicates group 1 vs group 2, and 1 indicates group 2 vs group 1
#'@param percent.thresh Only for paired data, numeric, indicate the significant count threshold 
#'@param nonpar Logical, indicate if a non-parametric test should be used (T or F)
#'@param threshp Numeric, indicate the p-value threshold
#'@param equal.var Logical, indicates if the group variance is equal (T) or unequal (F)
#'@param pval.type To indicate raw p-values, use "raw". To indicate FDR-adjusted p-values, use "fdr".  
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Volcano.Anal_bm<-function(mSetObj=NA, paired=FALSE, fcthresh, cmpType, percent.thresh, nonpar=F, threshp, equal.var=TRUE, pval.type="raw"){
        
        mSetObj <- .get.mSet(mSetObj);
        
        #### t-tests
        # check to see if already done by microservice
        if(.on.public.web & !nonpar & RequireFastUnivTests(mSetObj)){
                res <- PerformFastUnivTests(mSetObj$dataSet$norm, mSetObj$dataSet$cls, var.equal=equal.var);
        }else{
                res <- GetTtestRes(mSetObj, paired, equal.var, nonpar);
        }
        p.value <- res[,2];
        if(pval.type == "fdr"){
                p.value <- p.adjust(p.value, "fdr");
        }
        
        
        inx.p <- p.value <= threshp;
        p.log <- -log10(p.value);
        
        ### fold change analysis
        # make sure threshold is above 1
        fcthresh = ifelse(fcthresh>1, fcthresh, 1/fcthresh);
        max.xthresh <- log(fcthresh,2);
        min.xthresh <- log(1/fcthresh,2);
        
        res <- GetFC(mSetObj, paired, cmpType);
        
        # create a named matrix of sig vars for display
        fc.log <- res$fc.log;
        fc.all <- res$fc.all;
        
        inx.up <- fc.log > max.xthresh;
        inx.down <- fc.log < min.xthresh;
        
        if(paired){
                count.thresh<-round(nrow(mSetObj$dataSet$norm)/2*percent.thresh);
                mat.up <- res >= max.xthresh;
                mat.down <- res <= min.xthresh;
                
                count.up <- apply(mat.up, 2, sum);
                count.down <- apply(mat.down, 2, sum);
                fc.all <- rbind(count.up, count.down);
                
                inx.up <- count.up>=count.thresh;
                inx.down <- count.down>=count.thresh;
                
                colnames(fc.all) <- colnames(mSetObj$dataSet$norm);
                rownames(fc.all) <- c("Count (up)", "Count (down)");
                
                # replace the count.thresh for plot
                max.xthresh <- count.thresh;
                min.xthresh <- -count.thresh;
        }
        
        # create named sig table for display
        inx.imp <- (inx.up | inx.down) & inx.p;
        if(paired){
                sig.var <- cbind(fc.all[1,][inx.imp,drop=F], fc.all[2,][inx.imp, drop=F], p.value[inx.imp, drop=F], p.log[inx.imp, drop=F]);
                if(pval.type == "fdr"){
                        colnames(sig.var)<-c("Counts (up)","Counts (down)", "p.adjusted", "-log10(p)");
                }else{
                        colnames(sig.var)<-c("Counts (up)","Counts (down)", "raw.pval", "-log10(p)");
                }
                # first order by count difference, then by log(p)
                dif.count<-abs(sig.var[,1]-sig.var[,2]);
                ord.inx<-order(dif.count, sig.var[,4], decreasing=T);
                sig.var<-sig.var[ord.inx,,drop=F];
                sig.var[,c(3,4)]<-signif(sig.var[,c(3,4)],5);
        }else{
                sig.var <- cbind(fc.all[inx.imp,drop=F], fc.log[inx.imp,drop=F], p.value[inx.imp,drop=F], p.log[inx.imp,drop=F]);
                if(pval.type == "fdr"){
                        colnames(sig.var) <- c("FC", "log2(FC)", "p.ajusted", "-log10(p)");
                }else{
                        colnames(sig.var) <- c("FC", "log2(FC)", "raw.pval", "-log10(p)");
                }
                # first order by log(p), then by log(FC)
                ord.inx <- order(sig.var[,4], abs(sig.var[,2]), decreasing=T);
                sig.var <- sig.var[ord.inx,,drop=F];
                sig.var <- signif(sig.var,5);
        }
        
        fileName <- "volcano_sig.csv";
        write.csv(signif(sig.var,5), file=fileName);
        
        volcano <- list (
                raw.threshx = fcthresh,
                raw.threshy = threshp,
                paired = paired,
                max.xthresh = max.xthresh,
                min.xthresh = min.xthresh,
                thresh.y = -log10(threshp),
                fc.all = fc.all,
                fc.log = fc.log,
                fc.log.uniq = jitter(fc.log),
                inx.up = inx.up,
                inx.down = inx.down,
                p.log = p.log,
                inx.p = inx.p,
                sig.mat = sig.var
        );
        
        ## All retulst
        
        all.mat=cbind(volcano$fc.all, volcano$fc.log, 10^(-1*volcano$p.log), volcano$p.log )
        colnames(all.mat)=c("FC","log2(FC)","raw.pval","-log10(p)")
        fileName <- "volcano_all.csv";
        o=order(all.mat[,3])
        write.csv(all.mat[o,], file=fileName)
        
        mSetObj$analSet$volcano <- volcano;
        return(.set.mSet(mSetObj));
}



#'Retrieves T-test p-values
#'@description Utility method to get p values
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param paired Default set to FALSE
#'@param equal.var Default set to TRUE
#'@param nonpar Use non-parametric tests, default is set to FALSE
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
GetTtestRes_bm <- function(mSetObj=NA, paired=FALSE, equal.var=TRUE, nonpar=F){
  
        #mSetObj=mSet; paired=FALSE; equal.var=TRUE; nonpar=F
        
        mSetObj <- .get.mSet(mSetObj);
        inx1 <- which(mSetObj$dataSet$cls==levels(mSetObj$dataSet$cls)[1]);
        inx2 <- which(mSetObj$dataSet$cls==levels(mSetObj$dataSet$cls)[2]);
        if(length(inx1) ==1 || length(inx2) == 1){
          equal.var <- TRUE; # overwrite use option if one does not have enough replicates
        }
        univ.test <- function(x){t.test(x[inx1], x[inx2], paired = paired, var.equal = equal.var)};
        if(nonpar){
          univ.test <- function(x){wilcox.test(x[inx1], x[inx2], paired = paired)};
        }
        
        my.fun <- function(x) {
              if (length(table(x[inx1]))<3 || length(table(x[inx2]))<3){
                    return(c(NA,NA))
              }else{ 
                  tmp <- try(univ.test(x));
                  if(class(tmp) == "try-error") {
                         return(c(NA, NA));
                  }else{
                         return(c(tmp$statistic, tmp$p.value));
                  }  
              }
        }
        res <- apply(as.matrix(mSetObj$dataSet$norm), 2, my.fun);
        return(t(res));
}

#'Perform ANOVA analysis
#'@description ANOVA analysis
#'@usage ANOVA.Anal(mSetObj=NA, nonpar=F, thresh=0.05, post.hoc="fisher")
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param nonpar Logical, use a non-parametric test (T) or not (F)
#'@param thresh Numeric, from 0 to 1, indicate the p-value threshold
#'@param post.hoc Input the name of the post-hoc test, "fisher" or "tukey"
#'@param all_results Logical, if TRUE, it will output the ANOVA results for all compounds 
#'with no post-hoc tests performed.
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
ANOVA.Anal_all_bm<-function(mSetObj=NA, nonpar=F, thresh=0.05, post.hoc="fisher", all_results=FALSE, name_=""){
  
      mSetObj <- .get.mSet(mSetObj);
      
      sig.num <- 0;
      if(nonpar){
            aov.nm <- "Kruskal Wallis Test";
            anova.res <- apply(as.matrix(mSetObj$dataSet$norm), 2, kwtest, cls=mSetObj$dataSet$cls);
            
            #extract all p values
            res <- unlist(lapply(anova.res, function(x) {c(x$statistic, x$p.value)}));
            res <- data.frame(matrix(res, nrow=length(anova.res), byrow=T), stringsAsFactors=FALSE);
            
            fstat <- res[,1];
            p.value <- res[,2];
            
            names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
            fdr.p <- p.adjust(p.value, "fdr");
            
            #inx.imp <- p.value <= thresh;
            inx.imp <- fdr.p <= thresh;
            sig.num <- sum(inx.imp);
            
            if(sig.num > 0){
                    sig.f <- fstat[inx.imp];
                    sig.p <- p.value[inx.imp];
                    fdr.p <- fdr.p[inx.imp];
                    
                    sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(fdr.p,5), 'NA');
                    rownames(sig.mat) <- names(sig.p);
                    colnames(sig.mat) <- c("chi.squared", "p.value", "-log10(p)", "FDR", "Post-Hoc");
                    
                    # order the result simultaneously
                    ord.inx <- order(sig.p, decreasing = FALSE);
                    sig.mat <- sig.mat[ord.inx,,drop=F];
                    
                    fileName <- paste0(name_,"_kw_posthoc.csv");
                    my.mat <- sig.mat[,1:4];
                    colnames(my.mat) <- c("chi_squared", "pval_KW", "-log10(p)", "FDR");
            }
      }else{
            aov.nm <- "One-way ANOVA";
            if(.on.public.web & RequireFastUnivTests(mSetObj)){
                    res <- PerformFastUnivTests(mSetObj$dataSet$norm, mSetObj$dataSet$cls);
            }else{
                    aov.res <- apply(as.matrix(mSetObj$dataSet$norm), 2, aof, cls=mSetObj$dataSet$cls);
                    anova.res <- lapply(aov.res, anova);
                    
                    #extract all p values
                    res <- unlist(lapply(anova.res, function(x) { c(x["F value"][1,], x["Pr(>F)"][1,])}));
                    res <- data.frame(matrix(res, nrow=length(aov.res), byrow=T), stringsAsFactors=FALSE);
            }
            fstat <- res[,1];
            p.value <- res[,2];
            names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
            
            fdr.p <- p.adjust(p.value, "fdr");
            
            if(all_results==TRUE){
                    all.mat <- data.frame(signif(p.value,5), signif(-log10(p.value),5), signif(fdr.p,5));
                    rownames(all.mat) <- names(p.value);
                    colnames(all.mat) <- c("p.value", "-log10(p)", "FDR");
                    o=order(all.mat$p.value)
                    write.csv(all.mat[o,], paste0(name_,"_anova_all.csv"))
            }
            
            # do post-hoc only for signficant entries
            # inx.imp <- p.value <= thresh;
            inx.imp <- fdr.p <= thresh;
            sig.num <- sum(inx.imp);
            if(sig.num > 0){
                  # note aov obj is not avaible using fast version
                  # need to recompute using slower version for the sig ones
                  if(.on.public.web & RequireFastUnivTests(mSetObj)){
                          aov.imp <- apply(as.matrix(mSetObj$dataSet$norm[,inx.imp,drop=FALSE]), 2, aof, cls=mSetObj$dataSet$cls);
                  }else{
                          aov.imp <- aov.res[inx.imp];
                  }
                  sig.f <- fstat[inx.imp];
                  sig.p <- p.value[inx.imp];
                  fdr.p <- fdr.p[inx.imp];
                  cmp.res <- NULL;
                  post.nm <- NULL;
                  if(post.hoc=="tukey"){
                          tukey.res<-lapply(aov.imp, TukeyHSD, conf.level=1-thresh);
                          cmp.res <- unlist(lapply(tukey.res, parseTukey, cut.off=thresh));
                          post.nm = "Tukey's HSD";
                  }else{
                          fisher.res<-lapply(aov.imp, FisherLSD, thresh);
                          cmp.res <- unlist(lapply(fisher.res, parseFisher, cut.off=thresh));
                          post.nm = "Fisher's LSD";
                  }
                  
                  # create the result dataframe,
                  # note, the last column is string, not double
                  
                  sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(fdr.p,5), cmp.res);
                  rownames(sig.mat) <- names(sig.p);
                  colnames(sig.mat) <- c("f.value", "p.value", "-log10(p)", "FDR", post.nm);
                  
                  # order the result simultaneously
                  ord.inx <- order(sig.p, decreasing = FALSE);
                  sig.mat <- sig.mat[ord.inx,,drop=F];
                  fileName <- paste0(name_,"_anova_posthoc.csv");
            }
      }
      
      AddMsg(paste(c("A total of", sum(inx.imp), "significant features were found."), collapse=" "));
      if(sig.num> 0){
              res <- 1;
              write.csv(sig.mat,file=fileName);
              aov<-list (
                      aov.nm = aov.nm,
                      sig.num = sig.num,
                      sig.nm = fileName,
                      raw.thresh = thresh,
                      thresh = -log10(thresh), # only used for plot threshold line
                      p.value = p.value,
                      p.log = -log10(p.value),
                      inx.imp = inx.imp,
                      post.hoc = post.hoc,
                      sig.mat = sig.mat
              );
      }else{
            res <- 0;
            aov<-list (
                  aov.nm = aov.nm,
                  sig.num = sig.num,
                  raw.thresh = thresh,
                  thresh = -log10(thresh), # only used for plot threshold line
                  p.value = p.value,
                  p.log = -log10(p.value),
                  inx.imp = inx.imp
            );
      }
      mSetObj$analSet$aov <- aov;
      
      if(.on.public.web){
            .set.mSet(mSetObj);
            return(res);
      }else{
            return(.set.mSet(mSetObj));
      }
}


#'Perform PCA analysis
#'@description Perform PCA analysis, obtain variance explained, store item to PCA object
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'@param mSetObj Input name of the created mSet Object
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PCA.Anal_bm <- function(mSetObj=NA,name_){
  
      mSetObj <- .get.mSet(mSetObj);
      
      pca <- prcomp(mSetObj$dataSet$norm, center=TRUE, scale=F);
      
      # obtain variance explained
      sum.pca <- summary(pca);
      imp.pca <- sum.pca$importance;
      std.pca <- imp.pca[1,]; # standard devietation
      var.pca <- imp.pca[2,]; # variance explained by each PC
      cum.pca <- imp.pca[3,]; # cummulated variance explained
      
      # store the item to the pca object
      mSetObj$analSet$pca<-append(pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca));
      fast.write.csv(signif(mSetObj$analSet$pca$x,5), file=paste0(name_,"_pca_score.csv"));
      fast.write.csv(signif(mSetObj$analSet$pca$rotation,5), file=paste0(name_,"_pca_loadings.csv"));
      mSetObj$analSet$pca$loading.type <- "all";
      mSetObj$custom.cmpds <- c();
      return(.set.mSet(mSetObj));
}


## Performs volcano analysis
# (NOTE, not being used at the moment)
Volcano.Anal_all_bm<-function(mSetObj=NA, paired=FALSE, fcthresh, cmpType, percent.thresh, nonpar=F, threshp, equal.var=TRUE, pval.type="raw", name_){
  
      mSetObj <- .get.mSet(mSetObj);
      
      #### t-tests
      # check to see if already done by microservice
      if(.on.public.web & !nonpar & RequireFastUnivTests(mSetObj)){
            res <- PerformFastUnivTests(mSetObj$dataSet$norm, mSetObj$dataSet$cls, var.equal=equal.var);
      }else{
            res <- GetTtestRes_bm(mSetObj, paired, equal.var, nonpar);
      }
      p.value <- res[,2];
      if(pval.type == "fdr"){
            p.value <- p.adjust(p.value, "fdr");
      }
      
      
      inx.p <- p.value <= threshp;
      p.log <- -log10(p.value);
      
      ### fold change analysis
      # make sure threshold is above 1
      fcthresh = ifelse(fcthresh>1, fcthresh, 1/fcthresh);
      max.xthresh <- log(fcthresh,2);
      min.xthresh <- log(1/fcthresh,2);
      
      res <- GetFC(mSetObj, paired, cmpType);
      
      # create a named matrix of sig vars for display
      fc.log <- res$fc.log;
      fc.all <- res$fc.all;
      
      inx.up <- fc.log > max.xthresh;
      inx.down <- fc.log < min.xthresh;
      
      if(paired){
              count.thresh<-round(nrow(mSetObj$dataSet$norm)/2*percent.thresh);
              mat.up <- res >= max.xthresh;
              mat.down <- res <= min.xthresh;
              
              count.up <- apply(mat.up, 2, sum);
              count.down <- apply(mat.down, 2, sum);
              fc.all <- rbind(count.up, count.down);
              
              inx.up <- count.up>=count.thresh;
              inx.down <- count.down>=count.thresh;
              
              colnames(fc.all) <- colnames(mSetObj$dataSet$norm);
              rownames(fc.all) <- c("Count (up)", "Count (down)");
              
              # replace the count.thresh for plot
              max.xthresh <- count.thresh;
              min.xthresh <- -count.thresh;
      }
      
      # create named sig table for display
      inx.imp <- (inx.up | inx.down) & inx.p;
      if(paired){
              sig.var <- cbind(fc.all[1,][inx.imp,drop=F], fc.all[2,][inx.imp, drop=F], p.value[inx.imp, drop=F], p.log[inx.imp, drop=F]);
              if(pval.type == "fdr"){
                      colnames(sig.var)<-c("Counts (up)","Counts (down)", "p.adjusted", "-log10(p)");
              }else{
                      colnames(sig.var)<-c("Counts (up)","Counts (down)", "raw.pval", "-log10(p)");
              }
              # first order by count difference, then by log(p)
              dif.count<-abs(sig.var[,1]-sig.var[,2]);
              ord.inx<-order(dif.count, sig.var[,4], decreasing=T);
              sig.var<-sig.var[ord.inx,,drop=F];
              sig.var[,c(3,4)]<-signif(sig.var[,c(3,4)],5);
      }else{
            sig.var <- cbind(fc.all[inx.imp,drop=F], fc.log[inx.imp,drop=F], p.value[inx.imp,drop=F], p.log[inx.imp,drop=F]);
            if(pval.type == "fdr"){
                  colnames(sig.var) <- c("FC", "log2(FC)", "p.ajusted", "-log10(p)");
            }else{
                  colnames(sig.var) <- c("FC", "log2(FC)", "raw.pval", "-log10(p)");
            }
            # first order by log(p), then by log(FC)
            ord.inx <- order(sig.var[,4], abs(sig.var[,2]), decreasing=T);
            sig.var <- sig.var[ord.inx,,drop=F];
            sig.var <- signif(sig.var,5);
      }
      
      fileName <- paste0(name_, "_volcano_sig.csv");
      write.csv(signif(sig.var,5), file=fileName);
      
      volcano <- list (
            raw.threshx = fcthresh,
            raw.threshy = threshp,
            paired = paired,
            max.xthresh = max.xthresh,
            min.xthresh = min.xthresh,
            thresh.y = -log10(threshp),
            fc.all = fc.all,
            fc.log = fc.log,
            fc.log.uniq = jitter(fc.log),
            inx.up = inx.up,
            inx.down = inx.down,
            p.log = p.log,
            inx.p = inx.p,
            sig.mat = sig.var
      );
      
      ## All result
      
      all.mat=cbind(volcano$fc.all, volcano$fc.log, 10^(-1*volcano$p.log), volcano$p.log )
      colnames(all.mat)=c("FC","log2(FC)","raw.pval","-log10(p)")
      fileName <- paste0(name_,"_volcano_all.csv");
      o=order(all.mat[,3])
      write.csv(all.mat[o,], file=fileName)
      
      mSetObj$analSet$volcano <- volcano;
      return(.set.mSet(mSetObj));
}


#'Perform t-test analysis
#'@description This function is used to perform t-test analysis.
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param nonpar Logical, use a non-parametric test, T or F. False is default. 
#'@param threshp Numeric, enter the adjusted p-value (FDR) cutoff
#'@param paired Logical, is data paired (T) or not (F).
#'@param equal.var Logical, evaluates if the group variance is equal (T) or not (F). 
#'@param all_results Logical, if TRUE, returns T-Test analysis results
#'for all compounds. 
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Ttests.Anal_all_bm <- function(mSetObj=NA, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, all_results=FALSE, name_){
        
       # test lines
       # mSetObj=mSet; nonpar=F; threshp=0.05; paired=FALSE;  equal.var=TRUE;  all_results=TRUE; name_=nams2_[[ii]]   
        
        mSetObj <- .get.mSet(mSetObj);
        
        if(.on.public.web & !nonpar & RequireFastUnivTests(mSetObj)){
                res <- PerformFastUnivTests(mSetObj$dataSet$norm, mSetObj$dataSet$cls, var.equal=equal.var);
        }else{
                res <- GetTtestRes_bm(mSetObj, paired, equal.var, nonpar);
        }
        
        t.stat <- res[,1];
        p.value <- res[,2];
        names(t.stat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
        
        p.log <- -log10(p.value);
        fdr.p <- p.adjust(p.value, "fdr");
        
        if(all_results==TRUE){
          
              all.mat <- data.frame(signif(t.stat,5), signif(p.value,5), signif(p.log,5), signif(fdr.p,5));
              
              if(nonpar){
                      tt.nm = "Wilcoxon Rank Test";  
                      file.nm <- paste0(name_,"_wilcox_rank_all.csv")
                      colnames(all.mat) <- c("V", "p.value", "-log10(p)", "FDR");
              }else{
                      tt.nm = "T-Tests";
                      file.nm <- paste0(name_,"_t_test_all.csv");
                      colnames(all.mat) <- c("t.stat", "p.value", "-log10(p)", "FDR");
              }
              fast.write.csv(all.mat, file=file.nm);
        }
        
        inx.imp <- fdr.p <= threshp;
        # if there is no sig cmpds, it will be errors, need to improve
        
        AddMsg(paste("A total of", sum(inx.imp, na.rm=T), "significant features were found."));
        sig.num <- sum(inx.imp,  na.rm=T);    ## removing NA from the results table
        
        if(sig.num > 0){
              sig.t <- t.stat[inx.imp];
              sig.p <- p.value[inx.imp];
              lod<- -log10(sig.p);
              sig.q <-fdr.p[inx.imp];
              
              sig.mat <- cbind(sig.t, sig.p, lod, sig.q);
              colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", "FDR");
              ord.inx <- order(sig.p);
              sig.mat <- sig.mat[ord.inx,,drop=F];
              sig.mat <- signif(sig.mat, 5);
              
              if(nonpar){
                      tt.nm = "Wilcoxon Rank Test";  
                      file.nm <- paste0(name_,"_wilcox_rank.csv")
                      colnames(sig.mat) <- c("V", "p.value", "-log10(p)", "FDR");
              }else{
                      tt.nm = "T-Tests";
                      file.nm <- paste0(name_,"_t_test_sig.csv");
                      colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", "FDR");
              }
              fast.write.csv(sig.mat, file=file.nm);
              
              tt <- list (
                      tt.nm = tt.nm,
                      sig.nm = file.nm,
                      sig.num = sig.num,
                      paired = paired,
                      raw.thresh = threshp,
                      t.score = sort(t.stat),
                      p.value = sort(p.value),
                      p.log = p.log,
                      thresh = -log10(threshp), # only used for plot threshold line
                      inx.imp = inx.imp,
                      sig.mat = sig.mat
              );
        }else{
                tt <- list (
                        sig.num = sig.num,
                        paired = paired,
                        raw.thresh = threshp,
                        t.score = sort(t.stat),
                        p.value = sort(p.value),
                        p.log = p.log,
                        thresh = -log10(threshp), # only used for plot threshold line
                        inx.imp = inx.imp
                );
        }
        
        mSetObj$analSet$tt <- tt;
        
        if(.on.public.web){
          .set.mSet(mSetObj);
          return(as.numeric(sig.num));
        }
        return(.set.mSet(mSetObj));
}


#'Fold change analysis, unpaired
#'@description Perform fold change analysis, method can be mean or median
#'@usage FC.Anal.unpaired(mSetObj, fc.thresh=2, cmp.type = 0)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param fc.thresh Fold-change threshold, numeric input
#'@param cmp.type Comparison type, 0 for group 1 minus group 2, and 1 for group 
#'1 minus group 2
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
FC.Anal.unpaired_all_bm <- function(mSetObj=NA, fc.thresh=2, cmp.type = 0,name_){
  
      mSetObj <- .get.mSet(mSetObj);
      
      # make sure threshold is above 1
      fc.thresh = ifelse(fc.thresh>1, fc.thresh, 1/fc.thresh);
      max.thresh = fc.thresh;
      min.thresh = 1/fc.thresh;
      
      res <- GetFC(mSetObj, F, cmp.type);
      fc.all <- res$fc.all;
      fc.log <- res$fc.log;
      
      imp.inx <- fc.all > max.thresh | fc.all < min.thresh;
      sig.mat <- cbind(fc.all[imp.inx, drop=F], fc.log[imp.inx, drop=F]);
      colnames(sig.mat) <- c("Fold Change", "log2(FC)");
      
      # order by absolute log value (since symmetrical in pos and neg)
      inx.ord <- order(abs(sig.mat[,2]), decreasing=T);
      sig.mat <- sig.mat[inx.ord,,drop=F];
      
      fileName <- paste0(name_,"_fold_change.csv");
      fast.write.csv(sig.mat,file=fileName);
      
      # create a list object to store fc
      mSetObj$analSet$fc<-list (
            paired = FALSE,
            raw.thresh = fc.thresh,
            max.thresh = max.thresh,
            min.thresh = min.thresh,
            fc.all = fc.all, # note a vector
            fc.log = fc.log,
            inx.imp = imp.inx,
            sig.mat = sig.mat
      );
      
      return(.set.mSet(mSetObj));
}

#'Gets the current file location
getCurrentFileLocation <-  function(){
        this_file <- commandArgs() %>% 
                tibble::enframe(name = NULL) %>%
                tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
                dplyr::filter(key == "--file") %>%
                dplyr::pull(value)
        if (length(this_file)==0)
        {
                this_file <- rstudioapi::getSourceEditorContext()$path
        }
        return(dirname(this_file))
}


# idxMOLAR specific functions ----

#'Renames duplicated compound names in the input table ADP, ADPa , APDb etc
#' @param  compData a data frame with compound concentration data 
#' 
renameDuplicateCompound<-function(compData__){
      
      ## Testing lines
         #compData__ =compData.B_    
      
      chkList=list()
      for (v in 1:length(compData__$Name)){
            nm=compData__$Name[v]  
            
            if (!is.na(nm) & nm!="FALSE"){
                    if (nm %in% names(chkList)){
                          chkList[[nm]]<- chkList[[nm]]+1
                          compData__$Name[v]=paste0(compData__$Name[v],letters[chkList[[nm]]])
                    }else{
                          chkList[[nm]]<- 1         
                    }      
            }
      }
      return(compData__)
}

#'Performs sanity check, ie. checks if the input file  matches the  study type : profiling or tracer
#'Checks whether there is enough input data for the study type selected: OneMethod, MultiMethod
#'Checks whether proper/allowed method names are used in the sample name "HILIC_A" "HILIC_B" "LIPID" "T3
#'Prompts for error messages
#'@param compData.L_ list of data frames containing input data for several methods i.e HILIC_A, HILIC_B, ... 
#'@param prefixM_L_ list of input file prefixes  "MC00111_HILIC_A", "MC00111_HILIC_B ..."
#'@param studyType_ study type, tracer or profiling 
#'@param studyMethod_ study method, module or pipeline i.e idxMOLAR MOALAR MIDA
#' 
check_idxMETHOD_files<-function(compData.L_, prefixM_L_, studyType_, studyMethod_){
        
            ##  Testing lines
              # compData.L_=cD$compData_List ; prefixM_L_=cD$prefixM_List ; studyType_=run.Type_; studyMethod_=run.Method_
  
            compData.A_<-compData.L_[[1]]
            prefixM.A_<-prefixM_L_[[1]]
            
            nrows=500; if (dim(compData.A_)[1]<nrows){nrows=dim(compData.A_)[1]}
        
            LblngA=grep("exchange rate \\[%\\]",tolower(compData.A_[2:nrows,]))
            
            
            if (length(LblngA)>0 & studyType_=="profiling"){
              message("           Error: The input file '",prefixM.A_, "' contains labeling data.\n                  Please consider setting the study type 'Tracer analysis'.\n\n")
              return("Error")
            }
            
            if(length(LblngA)==0 & studyType_=="tracer"){
                    message("           Error: The input file '",prefixM.A_, "' contains no labeling data.\n                  Please consider setting the study type to 'Metabolite profiling'.\n\n")
                    return("Error")
            }  
            
     
            baseNames=sapply(compData.A_$FileName_[grep("\\.raw$",compData.A_$FileName_)],function(nm){tail(strsplit(nm,"\\\\")[[1]],n=1)})
            
            if (length(grep(methods_regex, baseNames))==0){
                    message("           Error: No proper method name is detected in the input sample names.\n                  Allowed method names are: 'HILIC_A', 'HILIC_B' 'LIPID' 'T3pos', 'PHpos' , 'PHneg' and 'PHTCA'.\n\n")
                    return("Error")
            }
        
            message("           File '",prefixM.A_,"' is ok")
            
            if (studyMethod_=="MultiMethods"){
                    
                    
                    for (ij in 2:length(compData.L_)){
                                
                                compData.B_=compData.L_[[ij]]
                                prefixM.B_<-prefixM_L_[[ij]]
                                
                                nrows=500; if (dim(compData.B_)[1]<nrows){nrows=dim(compData.B_)[1]}
                                
                                LblngB=grep("exchange rate \\[%\\]",tolower(compData.B_[2:nrows,]))
                                
                         
                                if ((length(LblngA)>0 & length(LblngB)==0 )|| (length(LblngA)==0 & length(LblngB)>0) ){
                                        message("           Error: The input file '",prefixM.B_, "' contains no ",studyType_, " data as in file '",prefixM.A_,"'\n                  Please check whether they  all contain the same data i.e. labeling or profiling.\n\n")
                                        return("Error")
                                } 
                                
                                ## Checking for the method names
                                baseNames=sapply(compData.B_$FileName_[grep("\\.raw$",compData.B_$FileName_)],function(nm){tail(strsplit(nm,"\\\\")[[1]],n=1)})
                                
                                if (length(grep(methods_regex, baseNames))==0){
                                        message("           Error: No proper method name is detected in the input sample names.\n                  Allowed names are: 'HILIC_A', 'HILIC_B' 'LIPID' and 'T3'.\n\n")
                                        return("Error")
                                }
                                message("           File '",prefixM.B_,"' is ok")
                                

                    }
             
            }
            

            message("           The input files passed sanity-check! They are OK for ", studyType_," analysis.")
            return("OK")
  }


#'Verifies whether excel file name is the right format eg: 'MC00099_HILIC_A.xlsx'
#'@param prefix_ is the input file name without file extension "MC00099_HILIC_A" 
#'
check_idxFileName_Format<-function(prefix__){
        methods_regex2=gsub("_\\|","\\|",methods_regex); methods_regex2=gsub("_$","",methods_regex2)
        currPrefix_=gsub("\\.xlsx$","",prefix__)
        
        ret="OK"
        if (length(grep(methods_regex2,currPrefix_))==0){
                ret=paste0("No method found in the file name '",currPrefix_,"'\nPlease use allowed method name i.e. ",gsub("\\|",", ",methods_regex2))
                return(ret)
        }
        
        if (length(grep("^MC\\d{5}",currPrefix_))==0){
                ret=paste0("No proper study number detected at the begining of the file name '",currPrefix_,"'\n Plase use MC{5digits} format:  'MC00012'")
                return(ret)
        }
        return("OK")
}

#'Reads the input file, to extract raw table from compound discover from each detection method used, it also returns the list of input file prefixes
#'Check if multiple input file are available in case study type is selected as "MultiMethods"
#'Re-formats sample names in the tables
#'@param inPath_ path to the input files
#'@param prefix_ input file name prefix, name without .xlsx extension : MC00111_HILIC_A )
#'@param studyType_ study type, tracer or profiling
#'@param studyMethod_ run method, OneMethod or MultiMethods
#'
read_INPUT_files_Idx<-function (inPath_, prefix_, studyType_, studyMethod_){

      ## Testing lines        
      # inPath_=inPath_;  prefix_=prefix_;  studyType_="tracer"; studyMethod_=run.Method_
      
        
      ## Checking if the input file names format is respected:
      ret=check_idxFileName_Format(prefix_); if (ret!="OK"){return(ret)}        
        
        
      prefixM_List<-compData_List<-list()        
      prefixM_List[[1]]<-prefix_
      
      
      inFiles=list.files(path=inPath_, pattern=paste0("^",strsplit(prefix_,"_")[[1]][1],"_.{1,8}xlsx$"))   # Matching the pattern get the MC00084 pattern and files with method name and ending with .xlxs
      
      ## Finding the names of other method's input files found in the same directory
      if (studyMethod_=="MultiMethods"){
            
              inFiles=inFiles[-grep(prefix_,inFiles)]
              if (length(inFiles)==0){return(paste0("No .xlsx file associtated  with ", prefix_," was found,\nPlease change the settings to 'one method'"))}      
              
              
              counter=2
              for (File in inFiles){
                      currPrefix_=gsub("\\.xlsx","",File)
                      
                      # Check if the input file names format is respected:
                      ret=check_idxFileName_Format(currPrefix_); if (ret!="OK"){return(ret)}        
                      
                      prefixM_List[[counter]] <-currPrefix_
                      counter=counter+1
              }
              
              if (length(prefixM_List)<2){return(paste0("No .xlsx file associtated  with ", prefix_," is found,\nPlease change the settings to 'one method'"))}
             
      }
      
      #### File one
      p1=paste0(inPath_,"/",prefixM_List[[1]],".xlsx" )
      if (!file.exists(p1)){ return(paste0("* The input file missing '", basename(p1),"' is missing"))}
      
      cat(paste0("\nReading : ", p1, ", please wait ...\n"))
      compData.1_<-read_excel(path=p1,col_types="text")
      
      if (!identical(tolower(colnames(compData.1_)[1:3]),c("checked","name","formula")) & !identical(tolower(colnames(compData.1_)[1:4]),c("checked","tags","name","formula")) ){ return(paste0("*  Error: wrong info for first three column names.\n         There is              :  ", paste0(colnames(compData.1_)[1:3],collapse=", "),"\n         While expecting       :  Checked, Name, Formula"))}
      
  
      cat("Reformatting sample names in current file, please wait ...\n")
      
      ## Find and rename columns that contain sample file names, and peak area
      compData.1_=as.data.frame(compData.1_)
      
      ## Sample file names
      fn.col1<<- which.max( unlist(apply(compData.1_, 2, function(cl){length(grep("File Name",cl,ignore.case = T))})))   ## column indeces that will be used to retrieve file names
      colnames(compData.1_)[fn.col1]<-"FileName_"
      
      ## Removal of date stamps from sample names
      cat("Removing any date stamp (14 digit format) from sample names.\n")
      compData.1_$FileName_=gsub( pattern="_\\d{14}","", compData.1_$FileName_)   
      
      ## Rename the duplicated compound names( eg. ADP, ADPb, ADPc)
      cat("Renaming duplicated coumpound names ( eg. ADP, ADPb, ADPc), please wait ...\n")   
      compData.1_=renameDuplicateCompound(compData.1_)
      
      compData_List[[1]]<-compData.1_
      
    
      if (studyMethod_=="MultiMethods"){
              
              for (ij in  2:length(prefixM_List)) { # form second file to the last excluding the comb
                      
                        cat("\n",ij)
                        #### File two
                        p2=paste0(inPath_,"/",prefixM_List[[ij]],".xlsx" )
                        if (!file.exists(p2)){ return(paste0("* The input file missing '", basename(p2),"' is missing"))}
                        
                        cat(paste0("\nReading : ", p2, ", please wait ...\n"))
                        compData.2_<-read_excel(path=p2,col_types="text")
                        
                        cat("Reformatting sample names in current file, please wait ...\n")
                        
                        ## Find and rename columns that contain sample file names, and peack area
                        compData.2_=as.data.frame(compData.2_)
                        
                        ## Sample file names
                        fn.col2= which.max( unlist(apply(compData.2_, 2, function(cl){length(grep("File Name",cl,ignore.case = T))})))
                        colnames(compData.2_)[fn.col2]<-"FileName_"
                        
                        ## Removal of date stamps from sample names
                        cat("Removing any date stamp (14 digit format) from sample names.\n")
                        compData.2_$FileName_=gsub( pattern="_\\d{14}","", compData.2_$FileName_)
                        
                        ## Rename the duplicated compund names( eg. ADP, ADPb, ADPc)
                        cat("Renaming duplicated coumpound names ( eg. ADP, ADPb, ADPc), please wait ...\n")   
                        compData.2_=renameDuplicateCompound(compData.2_)
                        
                        compData_List[[ij]]<- compData.2_
                }
                
                
      }
      
      cat(paste0("Reading input files completed\n"))
      
      return (list(compData_List=compData_List, prefixM_List=prefixM_List)) 
     
}


#'Extracts the Metabolomics Core experiment ID MCxxxxx
#'@param sampleNames_ list of sample names
#'
get_sample_MCddddd<-function(sampleNames_){
      string=sampleNames_
      pattern="MC\\d{5}"
      rr=regmatches(string, regexpr(pattern, string)) 
      return(rr)
}


#'Extracts the Method name either HILIC_B, or HILIC_C ... from sample names
#'@param sampleNames_ list of sample names
#'
get_sample_idxMETHOD <-function(sampleNames_){
      string=sampleNames_
      pattern= methods_regex
      rr=regmatches(string, regexpr(pattern, string,ignore.case = T)) 
      return(rr)
}


#'Extracts sample names, sample groups, study factors etc and returns a dataframe showing, sample names, group name, method name, study factors
#'@param ss sample names
#'
get_splitSampleNames_df_idx<-function(ss){
        #ss=sampleNames
        
        # get idx method
        idxMETHOD<-get_sample_idxMETHOD(ss)
      
        #remove the idxMethod
        ss2=sub(pattern=methods_regex, "", ss, ignore.case = T)
       
        #remove the study number
        ss2=sub(pattern="MC\\d{5}_", "", ss2)
        
        sample_names=ss2
        
        ## remove the sample identifier
        ss2=sub(pattern="S\\d{1,3}_", "", ss2)
        ss2=sub(pattern="_S\\d{1,3}", "", ss2)
        ss2=sub(pattern="S\\d{1,3}", "", ss2)
        
        
        ## remove the replicated number
        group=sub("_\\d{1,3}$", "",ss2)
        
        
        sFactors<-lapply(1:length(ss2), function(i){
           
                s=ss2[i]
                
                s_=strsplit(s,"_")[[1]]
                s_=s_[which(s_!="")]
                s_=s_[1:(length(s_)-1)]
                return(s_)
        })
        
        df=do.call(rbind,sFactors)
        
        sfCount=apply(df,2,function(xx){length(unique(xx))}) # removes non study factor part of the sample name
        
        if (length(which(sfCount>1))>1){
                imxx=which(sfCount==1)
                if (length(imxx)>0){df=df[,-imxx]}
        }else{
                df=as.data.frame(rep("None",length(sample_names)))
        }
        
        colnames(df)=paste0("Study_factor_",1:ncol(df))
        
        df=cbind(sample_names, group, idxMETHOD, df)
        
        df=as.data.frame(df)
        return(df)
}


#'Reformats names if that do not have replicate numbers. It will add numbering as follow: _1.raw, _2.raw, _3.raw
#'@param compData_ input dataframe for compound concentration data
#'
reformat_sampleNames_PRISM<-function(compData_){
        
        ## Testing lineß
        # compData_=TB_TPOOL_
 
        compData_$sample=gsub(pattern="_S\\d{1,3}|_S\\d{1,3}_|^S\\d{1,3}_|\\*$","",compData_$sample) 
    
        counterRep<-list()
        
        for (x in 1:length(compData_$sample)){
          
                ## List counter initialization at every compound data section     
          
                currName=compData_$sample[x]
                
                if (length(grep(currName ,pattern="_\\d{1,3}$"))!=0){
                        next
                }else{
                        if (currName %in% names(counterRep)){
                                counterRep[[currName]]=counterRep[[currName]]+1
                        }else{
                                counterRep[[currName]]=1
                        }
                        newName=paste0(currName,"_",counterRep[[currName]])
                        compData_$sample[x]=newName
                }
                
        }
        return(compData_)
}


#'Gets the list of sample names from the raw sample names as it is in the excel file
#'@param compData_ input data frame containing compound concentration data
#'
get_sampleNameList_idx<-function(compData_){
  
      ## Test lines
          # compData_=cD$compData_List[[2]]        
          # compData_=compData.1

      sampleNames_=list()

      for (x in 1:length(compData_$FileName_)){
            currText=compData_$FileName_[x]
            
            if (is.na(currText)){next}
            
            if (tolower(currText)=="file name"){
                  
                  ##pick the sample name,  reformat the path
                  currName=compData_$FileName_[(x+1)]
                  currName<-strsplit(currName,"\\\\")[[1]][length(strsplit(currName,"\\\\")[[1]]) ]
                  currName=strsplit(currName,".raw")[[1]][1]

                  sampleNames_[[currName]]<-currName
            }
      }
     
      return(sampleNames_)
}


#'Gets a list of sub data frames, each containing information for an individual compound 
#'Each list element is named after the compound 
#' For unkown compound detected (Checked==FALSE), we adopt naming of 'unk mass @ retention time' ie. unk_3232.11_32 -> 3232.11 is mass and 32 is RT)
#'@param compDta_ input row excel table from Compound Discover
#'@param cpIx_ indices each corresponding the position starting block containing data of a specific compound 
#' 
get_compoundSubDF_list<-function(compData_, cpIx_){
        
        ## Test line
        # compData_=compData_; cpIx_=cpIx
  
        compSubDF_list_<-lapply(1:(length(cpIx_)),function(i){
                df=NULL
                nm=NULL
                
                ## dataFrame
                if (i!=length(cpIx_)){
                        df=compData_[cpIx_[i]:(cpIx_[(i+1)]-1),] 
                }else{
                        df=compData_[cpIx_[i]:dim(compData_)[1],]        
                }
                return(df)
        }); 
        
        
        ## Generating compound name 
        ## For the check and using the "GuessName_Mass@RT, "Unk_mass@RT"
        cmpNames_<-sapply(cpIx_,function(i){
                nm=NULL
                if (is.na(compData_$Name[i])){
                        #mwIx=grep("Molecular Weight",colnames(compData_),ignore.case = T)
                        mwIx=grep(pattern="Calc\\. MW",colnames(compData_),ignore.case = T)
                        rtIx=grep("RT \\[min\\]",colnames(compData_),ignore.case = T)
                        nm=paste0("Unk_",round(as.numeric(compData_[i,mwIx ]),5),"_",round(as.numeric(compData_[i,rtIx]),5))
                }else if (compData_$Checked[i]==FALSE){
                        #mwIx=grep("Molecular Weight",colnames(compData_),ignore.case = T)
                        mwIx=grep(pattern="Calc\\. MW",colnames(compData_),ignore.case = T)
                        rtIx=grep("RT \\[min\\]",colnames(compData_),ignore.case = T)
                        nm=paste0(compData_$Name[i],"_",round(as.numeric(compData_[i,mwIx ]),5),"_",round(as.numeric(compData_[i,rtIx]),5))
                }else{
                        nm<-compData_$Name[i]
                }
                return(nm)
        })
        
        
        names(compSubDF_list_)=cmpNames_
        
        return(compSubDF_list_)
        
}


#'Extracts the table of MID FE (M+0, M+1, M+2 ...,M+n) from input data frame in Compound Discover format
#'@param compSubDF_list_ list of sub-dataframe containing data of each compound
#'@param sampleNames_ list of sample names for which to extract the lebelling data
#'@param outputDF_ a one column dataframe containing a list of unique sample names. Each compound MID data is appended to this DF and it is returned at the end
#'

get_MID_tbl<-function(compSubDF_list_, sampleNames_, outputDF_){
        
        # test line
        # compSubDF_list_=compSubDF_list; sampleNames_=sampleNames; outputDF_=outputDF
        
        message("   * Generating the table for MID")
        
        total=length(compSubDF_list_)
        counter=1
        for (compName in names(compSubDF_list_)){
                
                #cat(compName,"\n-------------\n")
                
                ## Compound level
                compSubDF<-compSubDF_list_[[compName]]
                
                LL<-list()
                
                pattern_="Exchange Rate"
                
                ixCol=grep(pattern=pattern_,as.character(compSubDF[2,]))
                nmCol=as.character(compSubDF[2,])[ixCol]
                
                ## Sample level
                for (sampleName in sampleNames_){
                        
                        #cat(sampleName,"\n")
                        
                        #checking the sample raw
                        ix=grep(paste0(sampleName,".raw"), compSubDF$FileName_)
                        
                        if (length(ix)>1){message(paste0("         --- * Error ::: sample name  '",sampleName,"' is duplicated for ", compName, ". Quiting now..."));
                                        Sys.sleep(5)
                                        stop()
                                        }
                        
                        # getting the Isopopoloyte M0, +1, +2 ... values
                        subVec=round(as.numeric(compSubDF[(ix-2),ixCol]),5)
                        names(subVec)=nmCol
                        
                        subVec=subVec[!is.na(subVec)]
                        
                        names(subVec)= sapply(names(subVec),function(n){paste0(compName,"_M+",strsplit(n,": ")[[1]][2])})
                        
                        LL[[sampleName]]<-subVec
                }
                
                ## Missing data for some samples / compounds  filled with zeros (# the sample is not available for the particular compound (there are some ampty lists in some cases))
                lens=sapply(LL,length)
                ix=which(lens==0)
                if (length(ix)>0){
                        for (ii in ix){
                                LL[[names(LL)[ii]]]<-rep(0,max(lens))
                        }
                }
                
                DF=do.call(rbind,LL)
                
                outputDF_=cbind(outputDF_, DF)
                
        }
        return(outputDF_)
}


#'Extracts the  total peak area pool size and total label, total enrichment FE : sum  of M+1, M+2, ..., M+n
#'@param compSubDF_list_ list of sub-dataframe containing data of each compound in compound discover format
#'@param sampleNames list of unique sample names for which data will be extracted
#'@param outputDF_ a one column dataframe containing a list of unique sample names. Each compound data is appended to this DF and it is returned at the end
#'@param what determins whether to extract peak area or labelling data i.e AREA, LABELLING
#'
get_AREA_LABEL_tbl<-function(compSubDF_list_, sampleNames_, outputDF_, what="AREA"){
        
        # Test line
        # compSubDF_list_=compSubDF_list; sampleNames_=sampleNames; outputDF_=outputDF; what="AREA"
        
   
        message("   * Generating the table for total ", what)
        
        total=length(compSubDF_list_)
        counter=1
        for (compName in names(compSubDF_list_)){
                
                # Test lines
                # cat(compName,"\n-------------\n")
                # compName=names(compSubDF_list_)[1]
          
                ## Compound level
                compSubDF<-compSubDF_list_[[compName]]
                compSubDF=as.data.frame(compSubDF)
                
                #rea index
                area.ix= which.max( unlist(apply(compSubDF, 2, function(cl){length(grep("AREA$",cl,ignore.case = T))})))
                
                
                VV<-c()
                
                ## Sample level
                for ( i in 1:length(sampleNames_)){
                        
                        sampleName=sampleNames_[[i]]
                        #cat(sampleName,"\n")
                        
                        #checking the sample raw
                        ix=grep(paste0(sampleName,".raw"), compSubDF$FileName_)
                        val=NULL
                        
                        
                        if (length(ix)==0){   #  In case the sample is not available for the particular compound
                                val=0
                        }else{
                              
                                if (what=="AREA"){
                                        
                                        val=as.numeric(compSubDF[(ix-2),area.ix])
                        
                                        if (length(val)>1){
                                              message(paste0("         --- * Warning ::: in compound ", compName," ", sampleName , " is dublicated now using  max (",paste0(val,collapse=", "),", max=",max(val),")"))
                                              val=max(val)
                                        
                                        } ### make sure you peak the max area if more  samples file appears more that  once  in a compound
                                        
                                }else if(what=="LABELLING"){
                                        
                                        iy=grep( "Exchange Rate \\[%\\]: 0",as.character(compSubDF[2,]))  ## 100 - "M+0"   100 - 12C
                                        val=100-as.numeric(compSubDF[(ix-2),iy])
                                }else{
                                        message("  * Please specify whether to find total pool size or total labeling / enrichment")
                                }
                                
                                if (is.na(val)){val=0}
                        }      
                        
                        VV=c(VV,round(as.numeric(val),5))
                }
                
                DF=as.data.frame(VV); rownames(DF)=sampleNames_; colnames(DF)=compName
                
                outputDF_=cbind(outputDF_, DF)
                
        }
        return(outputDF_)
}


#'From Compound Discover table format, it extracts a table the for all compounds, the table is in sample x compounds format 
#'NB: avoid that the last compound subDataFrame include all rows at the bottom
#'@param compData_ row input data frame after proprocessing to sample names
#'@param subDir_ the output path where the table will be saved
#'@param sampleNames_ list of unique sample names from the input excluding QC, Blanks and other unwanted sample names
#'@param prefix__ prefix of current method data.
#'@param outputDF_ the one column dataframe containing sample names, to which compound data will be append vertically/columnwise, these become the output dataframe 
#'@param  run.Type__ the run type or experiment type i.e tracer, profiling
#'
ExtractTBL_all_compound<-function(compData_, subDir_, sampleNames_,prefix__, outputDF_, run.Type__){

        # Testing lines
        # compData_=compData; subDir_=subDir; sampleNames_=sampleNames; prefix__=prefixTmp_; outputDF_=outputDF; run.Type__=run.Type_
  
  
        cpIx=which(compData_$Checked=="TRUE" | compData_$Checked=="FALSE") ## Get compound indices,  for all 
  
        compSubDF_list<-get_compoundSubDF_list(compData_,cpIx) ## slice the whole table into subData frame contains data for each compound samples files
        

        if (run.Type__=="profiling"){
                AREA_table<-get_AREA_LABEL_tbl(compSubDF_list, sampleNames_, outputDF_, what="AREA")
                write.csv( AREA_table, file=paste0(subDir_,prefix__,"_All_PoolSize.csv"))
        }
        
        if (run.Type__=="tracer"){
                ## Get the MID , AREA, LABALLING tables
                MID_table<-get_MID_tbl(compSubDF_list, sampleNames_, outputDF_)
                write.csv( MID_table, file=paste0(subDir_,prefix__,"_All_MID_FE.csv"))
                
                AREA_table<-get_AREA_LABEL_tbl(compSubDF_list, sampleNames_, outputDF_, what="AREA")
                write.csv( AREA_table, file=paste0(subDir_,prefix__,"_All_PoolSize.csv"))
                
                LABEL_table<-get_AREA_LABEL_tbl(compSubDF_list, sampleNames_, outputDF_, what="LABELLING" )
                write.csv( LABEL_table, file=paste0(subDir_,prefix__,"_All_TotalLabel_FE.csv"))
        }
}


#'From Compound Discover table format, it extracts a poosize/labbeling tables the for selected compounds, the table is in sample x compounds format 
#'The function works for idxMOLAR tracer and profiling data formats
#'NB: avoid that the last compound subDataFrame include all rows at the bottom
#'@param compData_ row input data frame after proprocessing to sample names
#'@param subDir_ the output path where the table will be saved
#'@param sampleNames_ list of unique sample names from the input excluding QC, Blanks and other unwanted sample names
#'@param prefix__ prefix of current method data.
#'@param outputDF_ the one column dataframe containing sample names, to which compound data will be append vertically/columnwise, these become the output dataframe 
#'@param  run.Type__ the run type or experiment type i.e tracer, profiling
#'
ExtractTBL_selected_compound<-function(compData_, subDir_, sampleNames_,prefix__, outputDF_, run.Type__){
        
        #testing line
        #compData_=compData; subDir_=subDir; sampleNames_=sampleNames; prefix__=prefixTmp_; outputDF_=outputDF;  run.Type__= run.Type_
        
        cpIx1=which(compData_$Checked==TRUE | compData_$Checked==FALSE)                     ## Get compound indeces,  for all 
        compSubDF_list1<-get_compoundSubDF_list(compData_,cpIx1)
        
        # indices to include
        cpIx2=which(compData_$Checked==TRUE)                     ## Get compound indices,  for selected TRUE
        compSubDF_list2<-get_compoundSubDF_list(compData_,cpIx2)
        
        compSubDF_list=compSubDF_list1[names(compSubDF_list2)]  ##  Get the right indices and corresponding compound subdataframes
        
        
        if (run.Type__=="profiling"){
                AREA_table<-get_AREA_LABEL_tbl(compSubDF_list, sampleNames_, outputDF_,what="AREA")
                write.csv( AREA_table, file=paste0(subDir_,prefix__,"_PoolSize.csv"))
        }
        if (run.Type__=="tracer"){
                # Get the MID , AREA, LABALLING tables
                MID_table<-get_MID_tbl(compSubDF_list, sampleNames_, outputDF_)
                write.csv( MID_table, file=paste0(subDir_,prefix__,"_MID_FE.csv"))
                
                AREA_table<-get_AREA_LABEL_tbl(compSubDF_list, sampleNames_, outputDF_,what="AREA")
                write.csv( AREA_table, file=paste0(subDir_,prefix__,"_PoolSize.csv"))
                
                LABEL_table<-get_AREA_LABEL_tbl(compSubDF_list, sampleNames_, outputDF_,what="LABELLING")
                write.csv(LABEL_table, file=paste0(subDir_, prefix__,"_TotalLabel_FE.csv"))
        }
}     




#'Merges the data table generated after reformatting
#'From the file name removes the method names 'HILIC_A', 'HILIC_B' etc to Comb 
#'returns the final table that will be used for statistical analyses
#'@param subDir_ path to the output directory
#'@param prefix.L_ list of file name prefixes of tables that will be merged
#'@param suffix_  the suffix that is added generated table "_PoolSize.csv", "_MID_FE.csv" ...
#'@param studySettings_ list of dataframes containing sample info, group info, study factor info etc
#'
merge_idxMETHOD_Multi<-function(subDir_,prefix.L_, suffix_, studySettings_){
        
        # Test lines (debug)
        # selected
        # subDir_=subDir; prefix.L_=cD$prefixM_List; suffix_="_PoolSize.csv"; studySettings_=studySettings 
        # All  
        # suffix_="_All_PoolSize.csv"
        
        
        message( "        - Merging   ",paste0(unlist(prefix.L_),collapse=",  "),"  :  ",gsub("^_","",suffix_))
        
        prefixComb__<-strsplit(prefix.L_[[1]],"_")[[1]][1]
        fComb=paste0(subDir_, prefixComb__,"_Comb",suffix_)
        
        tmpL<-lapply(1:length(prefix.L_),function(x_){f=paste0(subDir_,prefix.L_[[x_]],suffix_); read.csv(f,stringsAsFactors = F)})
        
      
        
        ## Check if we have the same sample names between different files 
        diff<-c()
        for (ix in 2:length(tmpL)){
                #if (!identical(tmpL[[1]]$sample,tmpL[[ix]]$sample)==TRUE){
                if (!identical(tmpL[[1]]$Sample_Name,tmpL[[ix]]$Sample_Name)==TRUE){
                        diff<-c(diff, ix)
                }
        }
        
        if (length(diff)>0){
                return(paste0("  * Error messages, sample names different methods do not macth. \n     Please check the input data for: ",paste0(unlist(prefix.L_)[diff],collapse=", "),"\n\n"))
        }  
        
        
        tmpDF=do.call(cbind, tmpL)
        
        ix=which(duplicated(colnames(tmpDF)))   ## remove repeated colnames sample group
        tmpDF=tmpDF[,-ix]
        
           
        ## Writing the final table        
        write.csv(tmpDF, fComb, row.names = F)
        
        ## Removing the temporally files
        lapply(1:length(prefix.L_),function(x_){f=paste0(subDir_,prefix.L_[[x_]],suffix_); unlink(f)})
        
        return(tmpDF)
}



# Common MIDA / MOLAR / IDX-MOLAR functions ----
#'Transposes MATRIX columns for prism formatting        
#'@param sampNames_ list of sample names
#'@param Mat1 input matrix (MID, MID intensity)
#'
transColsPrism<-function(sampNames_, Mat1){
        
        # Test lines
        # sampNames_=sampNames; Mat1= MatMID_I
        # sampNames_=sampNames; Mat1= MatMID
        # sampNames_ <- sampNames; Mat1 <- FE_MID
        
        list2= lapply(sampNames_,function(sampName_){ 
                ## added "^" for MC00189 sample
                # subMat2=Mat1[grep(pattern=paste0(sampName_, "_\\d$"),rownames(Mat1)),] ##brej's code
                subMat2=Mat1[grep(pattern=paste0("^", sampName_, "_\\d$"),rownames(Mat1)),]
                
                my_row=c()
                for (ic in 1:ncol(subMat2)){
                        rr=subMat2[,ic]
                        names(rr)=rep(colnames(subMat2)[ic],nrow(subMat2))
                        my_row=c(my_row,subMat2[,ic])
                }
                return(my_row)
                
        }) ## list2 end
        
        MatOut=do.call(rbind,list2)
        return(MatOut)
}


## Removes the indicated groups from the Group order table, sample normalization table, order and pairwise comparison settings table
#' @param studySettings_ list of dataframes containing sample info, group info, study factor info etc
#' 
updateSettings<-function(studySettings_){
        
        ## testing lines
        
        gTB=studySettings_[["groupsDF"]]
        oTB=studySettings_[["orderDF"]]
        sTB=studySettings_[["samplesDF"]]
        #fTB=studySettings_[["sfactors"]]
        gix=grep("n|no",tolower(oTB[,"Included"]))
        
        if (length(gix)>0){
                gnm=rownames(oTB)[gix]
                ic=grep(gnm,colnames(gTB))
                ir=grep(gnm,rownames(gTB))
                gTB=gTB[-ir,-ic]
                studySettings_[["groupsDF"]]=gTB
                sTB=sTB[-grep(pattern=paste0(gnm,collapse="|"),sTB$Sample_names),]
                studySettings_[["samplesDF"]]<-sTB
                
        }
        
        return(studySettings_) 
}

#'Format digits/numbers with leading zeros 01, 02, ... (result must be 2 digits)
#' @param nmbr_ the number to format
#' 
f00<-function(nmbr_){
        res=formatC(nmbr_, width=2, flag="0")
        return(res)
}


#'Matches the compound names with the compound Kegg ID for pathway analyses
#' @param tmp_DF the data frame contain Total Poolsize or Total Labeling data
#' @param compData_ the current data frame containng the compound note data (compound name and Kegg ids)
#' 
sub_CmpNames_with_KeggIDs<-function(tmp_DF, compData_){
      ## Test lines
      # tmp_DF=tmp_Pool; compData_=compNotesDF
 
      cmpds_=compData_[,"Compound"]
      ids_=compData_[,"Kegg ID"]
      
      for ( k in 3: dim(tmp_DF)[2]){
            cmp_name=colnames(tmp_DF)[k]
            ix=which(cmpds_==cmp_name)
            
            if (length(ix)>0){
                  colnames(tmp_DF)[k]=gsub(" ","",ids_[ix])
            }else{
                  colnames(tmp_DF)[k]=NA
            }
            cat(ids_[ix],"\n")
      }
      
      nas<-length(which(is.na(colnames(tmp_DF))))
      if (nas>0){
                message("            ",nas," compound(s) were not matched to KEGG ids in the compound notes")
      }else if(nas==0){
                message("           All compound(s) were matched to KEGG ids in the compound notes")
      }
      
      return(tmp_DF)
}


#'Reset global variable used to display tables
remove_tb_globals<-function(){
      MID_FE_global<<-NULL
      TotalLabel_FE_global<<-NULL
      PoolSize_global<<-NULL
      MID_PoolSize_global<<-NULL 
      TotalLabel_PoolSize_global<<-NULL
}

#'Computes Z-scores
#'@param v a vector or matrix containing values to compute z-score for
#' 
zscore<-function(v){(v-mean(v,na.rm=T)) / sqrt(var(v,na.rm=T))}      # z-score: (x-μ)/σ

#'Converts a dataframe into a matrix of numeric values
#'@param input dataframe
#'
data.matrix_bm<-function(df_){
        
        # dataDF[r.ix,c.ix]
        mat <- apply(as.matrix(df_[,1:dim(df_)[2]]),2, as.numeric)
        rownames(mat)=rownames(df_)
        return(mat)
}

#Removes the indicated groups or Samples from the save tables
#'@param tmpDF_ input dataframe
#'@param studySettings__ list of dataframes containing sample info, group info, study factor info etc
#'
remove_Samples_Groups<-function(tmpDF_,studySettings__){
          ## testing line
          # tmpDF_=tmpDF; studySettings__=studySettings_
          
        tbS=studySettings__$samplesDF; 
        six=grep("n|no",tolower(tbS[,"Included"]),ignore.case = T)
        if (length(six)>0){
                snms=tbS[six,"Sample_names"]
                sixRmv=sapply(snms,function(snm){which(tmpDF_[,"sample"]==snm)})
                tmpDF_=tmpDF_[-sixRmv,]
        }
        
        # groups removal
        tbG=studySettings__$orderDF; 
        gix=grep("n|no",tolower(tbG[,"Included"]),ignore.case = T)
        if (length(gix)>0){
              gnms=rownames(tbG)[gix]
              gixRmv=sapply(gnms,function(gnm){which(tmpDF_[,"group"]==gnm)})
              tmpDF_=tmpDF_[-gixRmv,]
        }
          
        return(tmpDF_)
}   


## Scales values using normalization factors (e.g. cell number)
#' @param tmpDF_ the data frame to scale
#' @param studySettings__ list of dataframes containing sample info, group info, study factor info etc
#' 
normalizeWithFactors<-function(tmpDF_,studySettings__){
        
  ## Test lines
  # tmpDF_=tmp1; studySettings__=studySettings_    
        
  sDF<-studySettings__$samplesDF
  
  if (length(unique(sDF$Normalization_factors))!=1){
         ix=match(tmpDF_$sample,sDF$Sample_names)
            
         nFactors=sDF$Normalization_factors[ix]
         for (cI in 1:dim(tmpDF_)[1]){
                tmpDF_[cI,3:dim(tmpDF_)[2]]<-tmpDF_[cI,3:dim(tmpDF_)[2]]/nFactors[cI]       
         }
            message("          Scaling sample(s) with normalization factor(s) completed")
  }
  return(tmpDF_)
}


#'Disables or enabling widgets for idxMOLAR, MOLAR, MIDA
#'@param currModule_ designates the current program running (idxMOLAR, MOLAR, MIDA)
#'
Dis.EnableWidgets_all<-function(currModule_){
        
        if (currModule=="MIDA" | currModule_=="MOLAR"){
               shinyjs::disable("idxMolarRun.Type")
               shinyjs::disable("idxMolarRun.Method")
        }else{
               shinyjs::enable("idxMolarRun.Type")
               shinyjs::enable("idxMolarRun.Method")
        }
  
  
        if (currModule=="idxMIDA" | currModule_=="MOLAR"){
                shinyjs::disable("midaAnalysis")
                shinyjs::disable("midaBaseLDet")
                shinyjs::disable("midaCMPDList")
        
          }else{
                shinyjs::enable("midaAnalysis")
                shinyjs::enable("midaBaseLDet")
                shinyjs::enable("midaCMPDList")
        }
  
  
        if (currModule=="MIDA" | currModule_=="idxMOLAR"){
                
                shinyjs::disable("molarPreprocSignalCorrection")
                shinyjs::disable("molarQC.Blank.thresh") 
                shinyjs::disable("molarNA.thresh") 
                shinyjs::disable("molarImputationMethod") 
                shinyjs::disable("molarCorrectionMethod")  
                shinyjs::disable("molarQC.span")  
                shinyjs::disable("molarP.degree") 
                
                shinyjs::disable("molarQC.Blank.threshCheck")
                shinyjs::disable("molarImputationMethodCheck")
                shinyjs::disable("molarNA.threshCheck")
                shinyjs::disable("molarSignalCorrCheck")
                
          }else{
            
                shinyjs::enable("molarPreprocSignalCorrection")
                shinyjs::enable("molarQC.Blank.thresh")
                shinyjs::enable("molarNA.thresh") 
                shinyjs::enable("molarImputationMethod") 
                shinyjs::enable("molarCorrectionMethod")  
                shinyjs::enable("molarQC.span")  
                shinyjs::enable("molarP.degree") 
                
                shinyjs::enable("molarQC.Blank.threshCheck")
                shinyjs::enable("molarImputationMethodCheck")
                shinyjs::enable("molarNA.threshCheck")
                shinyjs::enable("molarSignalCorrCheck")
            
        }
}


#' Generates the list of dataframes containing sample info, group info, study factor info etc  (also refered to studySettings dataframe)
#' @param sgroups_sfactors_ data frame containing sampleName, group name, method name, study factor 1, 2, ..., n
#' 
get_Group_Samples_Factors_DF_all<-function(sgroups_sfactors_){
          ## Test line
          #sgroups_sfactors_=sgroups_sfactors
          
          # df_S sample names and normalization factors
          MAT=matrix(NA, nrow=length(sgroups_sfactors_$sample_names),ncol=3, dimnames=list(NULL,c("Sample_names","Included","Normalization_factors")))
          MAT[,"Sample_names"]<-as.character(sgroups_sfactors_$sample_names)
          MAT=MAT[order(MAT[,"Sample_names"]),]
          
          df_S=data.frame(MAT)
          rownames(df_S)=NULL
          df_S$`Normalization_factors`=rep(1,nrow(df_S))   
          df_S$Included=rep("Y",nrow(df_S))
          
          gns=unique(sgroups_sfactors_$group); 
          
          MAT=matrix(NA, nrow=length(gns),ncol=length(gns), dimnames=list(gns,gns))
          
          # add FALSE to the upper triangular side
          MAT[upper.tri(MAT)]=FALSE
          MAT[lower.tri(MAT)]=FALSE
          
          df_G=data.frame(MAT)
          colnames(df_G)=c(colnames(MAT))
          rownames(df_G)=rownames(MAT)
          
          ## Getting the ordering
          MAT=matrix(cbind(rep(TRUE,length(gns)),1:length(gns)), nrow=length(gns),ncol=2,dimnames=list(gns,c("Included","Order")))
          
          orderDF=data.frame(MAT); 
          orderDF$Included=rep("Y",dim(MAT)[1])
          
          
          ## Getting the factor dataframe
          sfactors_<-colnames(sgroups_sfactors_)[grep("Study_factor",colnames(sgroups_sfactors_))]
          
          sfactors=sapply(sfactors_, function(factor){ paste0(unique(sgroups_sfactors_[,factor]),collapse="|")})
          
          
          MAT=matrix(NA, nrow=2,ncol=length(sfactors), dimnames=list(c("Separate","Compare"),sfactors))
          MAT[c(1,2),]=FALSE
          df_F=data.frame(MAT)
          colnames(df_F)=sfactors
          rownames(df_F)=c("Separate","Compare")
          
          ## Return the two data frames in a list
          #ordering sample names
          df_S=df_S[order(df_S$Sample_names),]
          
          DF_list<-list(samplesDF=df_S, groupsDF=df_G, factorsDF=df_F,orderDF=orderDF, sfactors=sfactors)
          return(DF_list)
  
}

#'Generates PRISM tables for Fractional Enrichment, Intensity, MID intentity for each compound 
#'@param p_MID_ path to the MID FE table
#'@param p_TPOOL_ path to the Total PoolSize table
#'@param subDir_ current outpout directory where the PRISM folder will be placed 
#'@param studySettings_ list of dataframes containing sample info, group info, study factor info etc

make_PRISM_Tables<-function(p_MID_, p_TPOOL_, subDir_, studySettings_ ){
        # Test line
        # p_MID_=p_MID; p_TPOOL_=p_TPOOL; subDir_=subDir;  studySettings_=studySettings
        # p_MID_=p_MID; p_TPOOL_=p_TPOOL; subDir_=subDir1; studySettings_=studySettings ## brej's code
        
        oRowNames=c()
        
        oD=studySettings_$orderDF
        
        gix=grep("n|no",tolower(oD[,"Included"]))
        if (length(gix)>0){oD=oD[-gix,]}
        oRowNames=rownames(oD)[order(oD$Order)]
        
        
        dirout.FE_MID=paste0(subDir_,"PRISM/FE_MID/")
        if (!file.exists(dirout.FE_MID)) {dir.create(dirout.FE_MID,r=T)}
        
        dirout.Intensity_MID=paste0(subDir_,"PRISM/Intensity_MID/")
        if (!file.exists(dirout.Intensity_MID)) {dir.create(dirout.Intensity_MID,r=T)}
        
        dirout.FE_12C_13C=paste0(subDir_,"PRISM/FE_12C_13C/")
        if (!file.exists(dirout.FE_12C_13C)) {dir.create(dirout.FE_12C_13C,r=T)}
        
        dirout.Intensity_12C_13C=paste0(subDir_,"PRISM/Intensity_12C_13C/")
        if (!file.exists(dirout.Intensity_12C_13C)) {dir.create(dirout.Intensity_12C_13C,r=T)}
        
        
        TB_MID_=read.csv(p_MID_,stringsAsFactors = F)
        TB_TPOOL_=read.csv(p_TPOOL_,stringsAsFactors = F)
        
        TB_MID_=reformat_sampleNames_PRISM(TB_MID_)
        TB_TPOOL_=reformat_sampleNames_PRISM(TB_TPOOL_)
        
        
        compNamesTPOOL_ALL_<-colnames(TB_TPOOL_)[3:dim(TB_TPOOL_)[2]]
        compNamesTPOOL<-unique(compNamesTPOOL_ALL_); # not needed #str(compNamesTPOOL) 
        
        compNamesMID_ALL_<-sub("_M\\.\\d{1,3}$","",colnames(TB_MID_)[3:dim(TB_MID_)[2]])
        compNamesMID<-unique(compNamesMID_ALL_); #str(compNamesMID)
        
        
        compNames=c()
        
        ## For debug purpose
        # View(compNamesTPOOL)
        # View(compNamesMID)
        
        if (!identical(compNamesTPOOL,compNamesMID)){
                message("    Warning: cannot match the compound names from the MID table with those from the total pool size")
        }else{
                compNames=compNamesTPOOL
        }
        
        total=length(compNames)
        ctr=1
        for (compName in compNames){
                
                ## Line for debug purpose
                # cat(compName,"\n")
                
                # compName <- compNames[1]
                # MID compound columns
                cTbl1=TB_MID_[,c(1,2+which(compNamesMID_ALL_==compName))]
                MatMID=data.matrix(cTbl1[,2:ncol(cTbl1)]); rownames(MatMID)=cTbl1$sample
                
                
                ### MID heatmap for every compound ** to be moved the point of generating stats & heatmaps
                ht_MID=MatMID
                colnames(ht_MID)=str_extract(colnames(ht_MID),pattern="M\\.\\d{1,2}$")
                ht_MID=ht_MID[,-grep("M\\.0",colnames(ht_MID))]
                
                if (sum(ht_MID)>0){
                        pdf(paste0( dirout.FE_MID,compName,".pdf"),width=dim(ht_MID)[2]*.8, height=dim(ht_MID)[1]*.3)
                        pheatmap::pheatmap(ht_MID, cluster_rows =F, cluster_cols=F, main=compName)
                        dev.off()
                }
                
                ## FE_MID and Int_MID
                FE_MID=MatMID/100
                
                #TOTAL PoolSize compund columns
                cTbl2=TB_TPOOL_[,c(1,2+which(compNamesTPOOL_ALL_==compName))]
                MatTPOOL=data.matrix(cTbl2[,2:ncol(cTbl2)]); 
                
                MatTPOOL=apply(MatTPOOL, 2,as.numeric)
                rownames(MatTPOOL)=cTbl2$sample
                
                
                if (!identical(rownames(FE_MID),rownames(MatTPOOL))){
                        message("    Warning: Missmatch between sample names for MID values and total pool size values")
                }else{
                        Int_MID=apply(FE_MID,2,function(col){col*MatTPOOL})
                        rownames(Int_MID)=rownames(MatTPOOL)
                }
                
                
                ######## 1-2. MID Intensities & Fractional enrichment 
                ## getting the row names (sample unique names)
                sampNamesALL<-sub("_\\d{1,3}$|_\\d{1,3}\\*$","",rownames(Int_MID)) # \\*$ catches any other letter at the end of replicate name
                sampNames=unique(sampNamesALL)
                names(sampNames)=sampNames
                
                #check whether the same number of replicate for each sample
                if (length(table(sampNamesALL))!=length(which(duplicated(table(sampNamesALL))))+1){
                        message("       *  Error not the same number of replicates are found: Please ckec input \n",paste0("         ",rownames(MatTPOOL),collaspe="\n"))
                        message("         Creation of PRISM tables is skipped ...")
                        return()
                }
                
                Intensity_MID=transColsPrism(sampNames_=sampNames, Mat1= Int_MID)
                Intensity_MID=cbind(rownames(Intensity_MID), Intensity_MID) ## adding sample names
                
                ## error
                FractionalEnrichment_MID=transColsPrism(sampNames_=sampNames, Mat1= FE_MID)
                FractionalEnrichment_MID=FractionalEnrichment_MID ## why?
                FractionalEnrichment_MID=cbind(rownames(FractionalEnrichment_MID), FractionalEnrichment_MID) ## adding sample names
                
                repNum=table(sampNamesALL)[1]
                
                colN=c()
                for (i in seq(0,ncol(FE_MID)-1))
                  {colN=c(colN,rep(i,repNum))}
                
                colN=c("Sample",paste0("M+",colN,collaspe="")) 
                
                ## An error 2021.12.23
                ## ncol(FractionalEnrichment_MID) != colN
                colnames(FractionalEnrichment_MID)=colN ## adding M+ column lebels
                
                write.table(FractionalEnrichment_MID[oRowNames,],file=paste0( dirout.FE_MID, compName,".csv"),col.names = TRUE,row.names=FALSE,quote=FALSE,sep=",")
                
                
                colnames(Intensity_MID)=colN ## adding M+ column lebels
                
                write.table(Intensity_MID[oRowNames,],file=paste0( dirout.Intensity_MID, compName,".csv"),col.names = TRUE,row.names=FALSE,quote=FALSE,sep=",")
                
                ######## done MID Intensity & Fractional enrichment ***
                
                
                ######## 3. Intensity using C12 C13 Total
                Int_MID<-t(apply(Int_MID,1,function(r){c(r[1],sum(r[2:length(r)]))}))
                colnames(Int_MID)=c("C12",rep("C13",ncol(Int_MID)-1))                    #adding all c13
                
                ## Transposing the replicates
                Intensity_12C_13C=transColsPrism(sampNames_=sampNames, Mat1=Int_MID)
                
                ## Transform and adding appropriate heading
                Intensity_12C_13C=cbind(rownames(Intensity_12C_13C), Intensity_12C_13C)
                colnames(Intensity_12C_13C)=c("Sample",rep("C12",table(sampNamesALL)[1]),rep("C13",table(sampNamesALL)[1]))
                
                write.table(Intensity_12C_13C[oRowNames,],file=paste0(  dirout.Intensity_12C_13C, compName,".csv"),col.names = TRUE,row.names=FALSE,quote=FALSE,sep=",")
                
                ######## done with Intensity***
                
                
                
                ######## 4. Fractional Enrichment (C12, C12, C12, C13, C13, C13) total
                FE_12C_13C<-t(apply(FE_MID,1,function(r){c(r[1],sum(r[2:length(r)]))}))
                colnames(FE_12C_13C)=c("C12",rep("C13",ncol(FE_12C_13C)-1))                    #adding all c13
                
                ## Transposing the replicates
                FractionalEnrichment_12C_13C=transColsPrism(sampNames_=sampNames, Mat1=FE_12C_13C)
                
                ## Transform and adding appropriate heading
                FractionalEnrichment_12C_13C=cbind(rownames(FractionalEnrichment_12C_13C), FractionalEnrichment_12C_13C)
                colnames(FractionalEnrichment_12C_13C)=c("Sample",rep("C12",table(sampNamesALL)[1]),rep("C13",table(sampNamesALL)[1]))
                
                
                write.table(FractionalEnrichment_12C_13C[oRowNames,],file=paste0( dirout.FE_12C_13C,compName,".csv"),col.names = TRUE,row.names=FALSE,quote=FALSE,sep=",")
                
                ######## done with Fractional Enrichment ***
                
                if (ctr%%20==0) message("                " , round((ctr*100/total),1),"% completed")
                ctr=ctr+1
        }
        message("                " , round(100,1),"% completed")
        
}

#'Generates the MID_PoolSize table formula :  MID_FE * MID_PoolSize
#'@param subDir_  the output path where the table will be saved
#'@param currPrefix_ prefix of current method data
#'@param suffix_MID_ suffix of the MID FE table 
#'@param suffix_PS_ suffix of the PoolSize table
#'@param suffix_Out_ Suffix of the MID_PoolSize  output file
#'@param mod_ program or module i.e idxMOLAR, MOLAR, MIDA
#'
make_MID_PoolSize_tb<-function(subDir_, currPrefix_, suffix_MID_, suffix_PS_, suffix_Out_, mod_){
        ## Testing lines
        
        # idxMOLAR
        # subDir_=subDir;  currPrefix_=prefix.1Multi; suffix_MID_="_MID_FE.csv";  suffix_PS_="_PoolSize.csv"; suffix_Out_="_MID_PoolSize"
        # MIDA
        # subDir_=subDir1; currPrefix_=prefix_; suffix_MID_="_MID_FE.csv"; suffix_PS_="_PoolSize.csv"; suffix_Out_="_MID_PoolSize"
  
  
        outname=paste0(currPrefix_, suffix_Out_)
        
        message( "        - Generating table :    ",outname)
        
        MID_FE = read.csv(file=paste0(subDir_, currPrefix_, suffix_MID_), stringsAsFactors = F)
        PoolSize = read.csv(file=paste0(subDir_, currPrefix_, suffix_PS_), stringsAsFactors = F)
        
        colnzPS=colnames(PoolSize)[3:length(colnames(PoolSize))]
        
        colnzFE=gsub(pattern="_M\\.\\d{1,2}$" ,"",colnames(MID_FE))
        
        for (coln in colnzPS ){
          
              clPS=which(colnames(PoolSize)==coln)
              
              clFE=which(colnzFE==coln )
              
              PoolSizeCol = PoolSize [,clPS]
              MID_FE[,clFE] = round((MID_FE[,clFE]/100)*PoolSizeCol,2)
          
        }
        
        MID_PoolSize <-MID_FE
        
        outfile=paste0(subDir_, outname)
        
        write.csv( MID_PoolSize, outfile, row.names = F)
        
        return(MID_PoolSize)
}


## Generates the TotalLabeling PoolSize table formula :  PoolSize  * TotalLabel_FE
#'Generates the MID_PoolSize table formula :  MID_FE * MID_PoolSize
#'@param subDir_  the output path where the table will be saved
#'@param currPrefix_ prefix of current method data
#'@param suffix_PS_ suffix of the PoolSize table
#'@param suffix_TlFE_ suffix of the Total FE table 
#'@param suffix_Out_ Suffix of the TotalLabel_FE  output file
#'@param mod_ program or module i.e idxMOLAR, MOLAR, MIDA
#'
make_TotalLabel_PoolSize_tb<-function(subDir_, currPrefix_, suffix_PS_, suffix_TlFE_, suffix_Out_ ){
        ## Testing lines
        # subDir_=subDir;  currPrefix_=prefix.1Multi; suffix_PS_="_PoolSize.csv" ; suffix_TlFE_="_TotalLabel_FE.csv"; 
        
        outname=paste0(currPrefix_, suffix_Out_)
        message( "        - Generating table :    ",outname)
        
        f11=paste0(subDir_, currPrefix_, suffix_PS_)
        f22=paste0(subDir_, currPrefix_, suffix_TlFE_)
        
        PoolSize = read.csv(file=f11, stringsAsFactors = F)
        TotalLabelFE= read.csv(file=f22, stringsAsFactors = F)
        
        colnzPS=colnames(PoolSize)[3:length(colnames(PoolSize))]
        colnzTLFE=colnames(TotalLabelFE)[3:length(colnames(TotalLabelFE))]
        
        if (!identical(colnzPS,colnzTLFE)){return( paste0("   Error colnames do not match between PoolSize and TotalLabelFE tables :\n",f11,"\n",f22))}
        
        PoolSize[,colnzPS]=round(PoolSize[,colnzPS]*TotalLabelFE[,colnzTLFE],2)
        
        TotalLabel_PoolSize<-PoolSize
        
        outfile=paste0(subDir_, outname)
        write.csv(TotalLabel_PoolSize, outfile, row.names = F)
        
        return(TotalLabel_PoolSize)
}

#'Reformats names and remove QC and Blank samples and reformat the names to factor1-factor2-factor3_1 by removing MC00000 (experiment ID)
#'@param inDF input dataframe
#'@param mod__ program or module i.e idxMOLAR, MOLAR, MIDA
#'
reformat_sNames_rmv_QC_Blnk<-function (inDF, mod__){
        
        ## Testing lines
        #inDF=tmp1;  mod_=mod
        
        colnames(inDF)[1]<-"sample"
        
        ix_=grep(pattern=idx_PB_Blanks_pattern,inDF[,"sample"], ignore.case = T)
        if (length(ix_)>0){
                outDF=inDF[-ix_,]
        }else{ outDF=inDF}
        
        ## remove the method eg: HELIC_A_
        if (mod__ =="idxMOLAR"){
                outDF[,1]=gsub(pattern=methods_regex,"",outDF[,1])
                outDF[,1]=gsub(pattern=gsub("_","",methods_regex),"",outDF[,1])
                
                outDF[,1]=gsub(pattern="MC\\d{5}_","",outDF[,1])
                outDF[,1]=gsub(pattern="MC\\d{5}","",outDF[,1])
        }        
        
        ## Add group names
        groups=sub("_\\d{1,3}$", "",outDF[,1])
        
        for (pat in c("_S\\d{1,3}_", "S\\d{1,3}_", "_S\\d{1,3}")){
                groups=gsub(pat,"", groups)
        }
        
        colnames(outDF)[1]<-"sample"
        outDF[,2]<-groups; colnames(outDF)[2]<-"group"
        
        
        return(outDF)
}

#'The function reformat generated MID_FE, PoolSise and TotalLalbel tables (all modules: MIDA, MOLAR, idxMOLAR)
#'Removes the method name eg.'HILIC_A' sample from names and group name, 
#'Adds a method prefix to compound names (Method abbreviation) for (idxMOLAR)
#'Removes the indicated groups
#'Creates new replicates for group with less than 3 reps
#'Scales final values with normalization factors
#'Writes the final table - which could be merged, and will be used for statistical analyses
#'@param subDir_  path to the output directory
#'@param prefix1__ file name prefix of tables
#'@param suffix_ the suffix that is added generated table
#'@param studySettings_ list of dataframes containing sample info, group info, study factor info etc
#'@param mod_ module/pipeline name can be idxMOLAR, MOLAR or MIDA
#'
format_Tb_all<-function(subDir_, prefix1__, suffix_, studySettings_, mod_){
        
        ## Test lines (debug)
                #idxMOLAR
                # subDir_=subDir; prefix1__=prefix.1; suffix_="_MID_FE.csv"; suffix_="_TotalLabel_FE.csv" ; studySettings_=studySettings ; mod_="idxMOLAR"
                # MOLAR
                # subDir_=subDir; prefix1__=prefix_; suffix_="_PoolSize.csv"; studySettings_=studySettings; mod_="MOLAR"
                # MIDA
                # subDir_=subDir; prefix1__=prefix_; suffix_="_PoolSize.csv"; studySettings_=studySettings; mod_="MIDA"
                # subDir_=subDir; prefix1__=prefix_; suffix_="_MID_FE.csv"; studySettings_=studySettings; mod_="MIDA"                
                
        message("        - Reformating ", prefix1__, suffix_)
        f1=paste0(subDir_,prefix1__,suffix_)
        
        tmp1=read.csv(f1,stringsAsFactors = F)
        
        
        ## Reformats names and remove QC and Blank samples and reformat the names to factor1-factor2-factor3_1 by removing MC00000
        tmp1=reformat_sNames_rmv_QC_Blnk(inDF=tmp1, mod__=mod_)   ## This has been done while exrtacting the tables
        
        ## Replace any .., ..., .... introduced by spaces with "."
        
        colnames(tmp1)=gsub("\\.{2,4}",".",colnames(tmp1))
        
        ## For debug purpose
        #View(tmp1$Sample_Name)
        
        if (mod_=="idxMOLAR"){
                sfx.1=""; if(length(grep("HILIC_A",prefix1__))>0){sfx.1="a"}; if(length(grep("HILIC_B",prefix1__))>0){sfx.1="b"};  if(length(grep("LIPID",prefix1__))>0){sfx.1="l"};  if(length(grep("T3",prefix1__))>0){sfx.1="t"}
                if(length(grep("T3pos",prefix1__))>0){sfx.1="tp"};  if(length(grep("PHneg",prefix1__))>0){sfx.1="pn"};  if(length(grep("PHpos",prefix1__))>0){sfx.1="pp"};  if(length(grep("PHTCA",prefix1__))>0){sfx.1="pt"}
                
                for (v in 3:length(colnames(tmp1))){   colnames(tmp1)[v]=paste0(sfx.1,"_",colnames(tmp1)[v]) }
        }        
        
        
        ## Scaling final values with normalization factors, 
        if (suffix_=="_All_PoolSize.csv" | suffix_=="_PoolSize.csv"){
                tmp1<-normalizeWithFactors(tmpDF_=tmp1, studySettings__=studySettings_)
        }
        
        
        tmp1$group=gsub("S\\d{1,3}_","",tmp1$group)
        
        
        ## Remove the indicated groups
        tmp1=remove_Samples_Groups(tmpDF_=tmp1, studySettings__=studySettings_)
        
        studySettings_=updateSettings(studySettings_)
        
        ## Create new replicates for group with less than 3 reps
        groupNamesALL<-tmp1[,"group"]
        
        if (max(table(groupNamesALL))!=min(table(groupNamesALL))){
                message("          Warning: some sample have fewer replicates \n")
                message("          Using average for available replicates to create a new replciate (eg repName_3*)\n")
                
                gprs=names(which(table(groupNamesALL)==2))
                
                for (gpr in gprs){
                        m=tmp1[ grep(gpr,tmp1[,"group"]),3:dim(tmp1)[2]]
                        m=data.matrix_bm(m)
                        sim=colMeans(m)
                        snames=tmp1[ grep(gpr,tmp1[,"group"]),"sample"]
                        
                        repN=max(as.numeric(regmatches(snames, regexpr(pattern="\\d{1,3}$", snames)))) 
                        simNam=paste0(gpr,"_",repN+1,"*")
                        tmp1=rbind(tmp1,c(simNam, gpr,sim))
                }
                
        }
        
      
        ## Writing the final table
        ret_DF<-tmp1[order(tmp1[,"group"]),]
        
        #ret_DF[,3:dim(ret_DF)[2]]<-round(ret_DF[,3:dim(ret_DF)[2]],2) ## round tp 2 decimals
        
        ret_DF[,3:dim(ret_DF)[2]]=apply(ret_DF[,3:dim(ret_DF)[2]], 2, function(v){round(as.numeric(v),2)})
        
        write.csv(ret_DF, f1, row.names = F)
   
        return(ret_DF)
}

#'Formats Data Table for visualization on the application interface
#'@param DF_ input table / dataframe
#' 
format_DT<-function(DF_){
        DF_=t(DF_)
        DF_=cbind(rownames(DF_),DF_)
        rownames(DF_)<-NULL
        colnames(DF_)<-DF_[1,]
        DF_=DF_[-c(1,2),] 
        colnames(DF_)[1]="Compound"

        DF_=as.data.frame(DF_)
        
       for (p in 2:dim(DF_)[2]){
                  DF_[ , p] <-as.numeric(DF_[, p])
       }
       
        return(DF_)
}
  

#'Quantitative enrichment analysis with globaltest
#'@description Various enrichment analysis algorithms
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
CalculateGlobalTestScore_bm<-function(mSetObj=NA){
        #browser()
        mSetObj <- .get.mSet(mSetObj);
        
        if(.on.public.web){
                .prepare.globaltest.score(mSetObj);
                .perform.computing();   
                .save.globaltest.score(mSetObj);  
                return(.set.mSet(mSetObj));
        }else{
                mSetObj <- .prepare.globaltest.score(mSetObj);
                .perform.computing();   
                mSetObj <- .save.globaltest.score(mSetObj);  
        } 
        
        return(.set.mSet(mSetObj));
}


## Generates the final reports
#' @param outPathDirz list of paths of all pdf figures
#' @param subDir_ the output path where report will be written
#' 
generateFinalReport <-function(outPathDirz, subDir_){
        
        ## outPathDirz=outPathDirs; subDir_=subDir
        
        for (pName in names(outPathDirz)){
                message("        -  Generating pdf report for ",pName, " ~")
                pathPDF=paste0(outPathDirz[[pName]],"/")
                pdfs=list.files(pattern = "\\.pdf$", path = pathPDF)
                
                rmv=grep(pattern="^N_",pdfs)
                if (length(rmv)>0){pdfs=pdfs[- rmv]}
                
                
                ## Skip including damaged PDFs from the report
                size=sapply(pdfs,function(f){file.info(paste0(pathPDF,f))$size})
                
                ix=which(size==0)
                if (length(ix>0)){
                        message("        *  Damaged PDF figure(s) (plese check them):\n     ",paste0(pdfs[ix],collapse="\n     "))
                        pdfs=pdfs[-ix]
                }
                
                ## Ordering of pdfs
                o=order(as.numeric(sapply(pdfs,function(p){strsplit(p,"_")[[1]][1]})))
                pdfs<-pdfs[o]
                
                print(pdfs)
                
                ## Save a combined pdf
                out <- tryCatch(
                        pdf_combine(input=paste0(pathPDF,pdfs),output=paste0(subDir_,pName,".pdf")),
                        
                        error=function(cond) {
                                message(cond)
                                return(NA)
                        },
                        warning=function(cond) {
                                message(cond)
                                # Choose a return value in case of warning
                                return(NULL)
                        },
                        finally={
                                next
                        }
                )
                
        }  
}


#'Reformats the table for MetaboAnalyst use, returns the path of save tmp.csv .file (does group ordering/numbering for heatmaps and PCA plot)
#'@param input_ full path of the input table that needs to be reformated
#'@param tmp_name a temporally name used  
#'@param  studySettings_ list of dataframes containing sample info, group info, study factor info etc
#'
reSaveTbl<-function(input_,tmp_name, studySettings_ ){
        
        # Test lines
                # input_=paste0(subDir,prefix.1Multi,"_MID13C_only.csv"); tmp_name="tmpTPool.csv"; studySettings_=studySettings
                # input_=paste0(subDirCurr,prefix.1Multi,"_PoolSize.csv"); tmp_name="tmpTPool.csv";  studySettings_=studySettingz  
                # input_=paste0(subDirCurr,prefixCurr, "_PoolSize.csv"); tmp_name="tmpTPool.csv"; inPath__=inPath_;  studySettings_=studySettingz
                
        
        inDF=read.csv(input_,stringsAsFactors = F)
        
        ## Re-saving table with order number   
        oRowNames=c()
        oD=studySettings_[["orderDF"]]
        oix=grep("n|no",tolower(oD$Included))
        if (length(oix)>0){oD=oD[-oix,]} 
        
        mmax=max(as.numeric(oD[,2]))  # for fromating with leading zeros
        
        groupsN=unlist(lapply(inDF$group,function(gn){
                ii=which(rownames(oD)==gn)
                
                nmbr=oD[ii,2]
                if (mmax>9){nmbr=f00(nmbr)}             # for fromating with leading zeros
                return(paste0(nmbr,"_",gn))
        }))
        
        inDF[,"group"]=groupsN
        inDF=inDF[order(inDF$group),]
        
        
        outf=paste0(getwd(),"/",tmp_name)
        
        ## rename the psdseudo
        ijj=grep("\\*$",inDF$sample)
        inDF$sample[ijj]<-paste0("ps-",gsub("\\*$","",inDF$sample[ijj]))
        
        write.csv(inDF, file=outf,row.names = FALSE)
        return(outf)
}

## Prepared the input files for the multiple groups statistical analyses
#' @param run.TypeCurr_ variable indicating whether it is a tracer or profiling study
#' @param subDirCurr_ the path to the current output directory
#' @param prefixCurr_ the current prefix, it is the name of the output directory
#' @param studySettingz_ list of dataframes containing sample info, group info, study factor info etc
#' 
prep.MetaboAllG_Stats_all<-function(run.TypeCurr_, subDirCurr_, prefixCurr_, studySettingz_){
        
        #Debug/test lines
                # run.TypeCurr_=run.TypeCurr; subDirCurr_=subDirCurr, 
                # prefixCurr_=prefixCurr, studySettingz_=studySettingz
        
                        
        ## Preparing files
        files<-nams<-nams2<-c()
        
        if (run.TypeCurr_=="profiling"){
                filePool=reSaveTbl(paste0(subDirCurr_,prefixCurr_,"_PoolSize.csv"),"tmpTPool.csv", studySettings_=studySettingz_)   
                files=c(filePool)
                nams=c("Total pool size")
                nams2=c("T.pool")  ### label for figure file name
        }
        
        if (run.TypeCurr_=="tracer"){
                filePool=reSaveTbl(paste0(subDirCurr_,prefixCurr_,"_PoolSize.csv"),"tmpTPool.csv",  studySettings_=studySettingz_)   
                fileLabel=reSaveTbl(paste0(subDirCurr_,prefixCurr_,"_TotalLabel_FE.csv"),"tmpTLabel.csv",  studySettings_=studySettingz_)   
                files=c(filePool,fileLabel)
                nams=c("Total pool size", "Total labelling")
                nams2=c("T.pool","T.label")  ### label for figure file name
        }
        
        return(list("files"=files,"nams"=nams,"nams2"=nams2))
}


#'Prepared the input files for the two groups statistical analyses
#'@param run.TypeCurr_ variable indicating whether it is a tracer or profiling study
#'@param subDirCurr_ the path to the current output directory
#'@param prefixCurr_ the current prefix, it is the name of the output directory
#'@param studySettingz_ list of dataframes containing sample info, group info, study factor info etc
#'@param files_PoolLabel full paths to the input Total PoolSize and Total Labelling tables
#'@param g1 the name of the first group
#'@param g2 the name of the second group
#'@param dirout.2G path of the directory where group comparison results will be saved
#' 
prep.Metabo2G_Stats_all<-function(run.TypeCurr_, subDirCurr_, prefixCurr_, studySettingz_,files_PoolLabel, g1 , g2, dirout.2G){
        
        ## Lines for testing
        # run.TypeCurr_=run.TypeCurr; subDirCurr_=subDirCurr; prefixCurr_=prefixCurr; studySettingz_=studySettingz; files_PoolLabel=outL.AllGroups$files; g1=g1 ; g2=g2
        
        files<-nams<-nams2<-c()
        
        if (run.TypeCurr_=="profiling"){ 
                filePool=files_PoolLabel[1]
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                ix1=c(grep(paste0("_",g1,"$"),tmp_Pool$group),grep(paste0("_",g2,"$"),tmp_Pool$group))
                
                tmp_Pool=tmp_Pool[ix1,]
                
                filePool2G=paste0(dirout.2G,"/tmpTPool.csv")
                write.csv(tmp_Pool,file=filePool2G, row.names = FALSE)
                
                files=c(filePool2G)
                nams=c("Total pool size")
                nams2=c("T.pool")  ### label for figure file name
        }
        
        if (run.TypeCurr_=="tracer"){ 
                
                filePool=files_PoolLabel[1]
                fileLabel=files_PoolLabel[2]
                
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                tmp_Label=read.csv(fileLabel,stringsAsFactors = F)
                
                ix1=c(grep(paste0("_",g1,"$"),tmp_Pool$group),grep(paste0("_",g2,"$"),tmp_Pool$group))
                
                tmp_Pool=tmp_Pool[ix1,]
                ix2=c(grep(paste0("_",g1,"$"),tmp_Label$group),grep(paste0("_",g2,"$"),tmp_Label$group))
                
                
                tmp_Label=tmp_Label[ix2,]
                
                filePool2G=paste0(dirout.2G,"/tmpTPool.csv")
                fileLabel2G=paste0(dirout.2G,"/tmpTLabel.csv") 
                write.csv(tmp_Pool,file=filePool2G, row.names = FALSE)
                write.csv(tmp_Label,file=fileLabel2G, row.names = FALSE)
                
                files=c(filePool2G,fileLabel2G)
                nams=c("Total pool size","Total labelling")
                nams2=c("T.pool","T.label")  ### label for figure file name
        }
        
        return(list("files"=files,"nams"=nams,"nams2"=nams2))
}

## Prepared the input files for two  group pathway analyses
#'@param run.TypeCurr_ variable indicating whether it is a tracer or profiling study
#'@param subDirCurr_ the path to the current output directory
#'@param prefixCurr_ the current prefix, it is the name of the output directory
#'@param studySettingz_ list of dataframes containing sample info, group info, study factor info etc
#'@param compNotesDF_ the CompoundNotes table (MIDA- and MOLAR-specific), used to get the Kegg IDs
#'@param files_PoolLabel full paths to the input Total PoolSize and Total Labelling tables
#'@param g1 the name of the first group
#'@param g2 the name of the second group
#'@param dirout.2G path to the output directory where group comparison results will be saved
#' 
prep.Metabo2G_Pathway_all<-function(run.TypeCurr_, subDirCurr_, prefixCurr_, compNotesDF_, studySettingz_, files_PoolLabel, g1 , g2, dirout.2G){
        
        ## Lines for testing
        # run.TypeCurr_=run.TypeCurr; subDirCurr_=subDirCurr; prefixCurr_=prefixCurr; studySettingz_=studySettingz; files_PoolLabel=outL.AllGroups$files; g1=g1 ; g2=g2
        
        files<-nams<-nams2<-c()
        
        if (run.TypeCurr_=="profiling"){ 
                
                filePool=files_PoolLabel[1]
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                
                ix1=c(grep(paste0("_",g1,"$"),tmp_Pool$group),grep(paste0("_",g2,"$"),tmp_Pool$group))
                
                tmp_Pool=tmp_Pool[ix1,]
                
                filePool2G_kegg=paste0(dirout.2G,"/tmpTPool_kegg.csv")
                
                tmp_Pool_kegg<-sub_CmpNames_with_KeggIDs(tmp_Pool, compNotesDF_)
                write.csv(tmp_Pool_kegg,file=filePool2G_kegg, row.names = FALSE)
                
                files_kegg=c(filePool2G_kegg)
                nams=c("Total pool size")
                nams2=c("T.pool")  ### label for figure file name
        }
        
        if (run.TypeCurr_=="tracer"){ 
                
                filePool=files_PoolLabel[1]
                fileLabel=files_PoolLabel[2]
                
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                tmp_Label=read.csv(fileLabel,stringsAsFactors = F)
                
                ix1=c(grep(paste0("_",g1,"$"),tmp_Pool$group),grep(pattern=paste0("_",g2,"$"),tmp_Pool$group))
                
                tmp_Pool=tmp_Pool[ix1,]
                ix2=c(grep(paste0("_",g1,"$"),tmp_Label$group),grep(paste0("_",g2,"$"),tmp_Label$group))
                tmp_Label=tmp_Label[ix2,]
                
                
                filePool2G_kegg=paste0(dirout.2G,"/tmpTPool_kegg.csv")
                fileLabel2G_kegg=paste0(dirout.2G,"/tmpTLabel_kegg.csv") 
                
                
                tmp_Pool_kegg<-sub_CmpNames_with_KeggIDs(tmp_Pool, compNotesDF_)
                tmp_Label_kegg<-sub_CmpNames_with_KeggIDs(tmp_Label, compNotesDF_)
                write.csv(tmp_Pool_kegg,file=filePool2G_kegg, row.names = FALSE)
                write.csv(tmp_Label_kegg,file=fileLabel2G_kegg, row.names = FALSE)
                
                files_kegg=c(filePool2G_kegg, fileLabel2G_kegg)
                nams=c("Total pool size","Total labelling")
                nams2=c("T.pool","T.label")  ### label for figure file name
        }
        return(list("files_kegg"=files_kegg, "nams"=nams,"nams2"=nams2))
}


#'Prepared the input files for statistical analyses for a specific study factor
#'@param run.TypeCurr_ variable indicating whether it is a tracer or profiling study
#'@param subDirCurr_ the path to the current output directory
#'@param prefixCurr_ the current prefix, it is the name of the output directory
#'@param studySettingz_ list of dataframes containing sample info, group info, study factor info etc
#'@param files_PoolLabel full paths to the input Total PoolSize and Total Labelling tables
#'@param sFactor  study factor to study
#'@param dirout.sF.sep path to the output directory wheter results will be save
#'
prep.Metabo1F_Stats_all<-function(run.TypeCurr_, subDirCurr_, prefixCurr_, studySettingz_,files_PoolLabel, f1, dirout.sF.sep){
        
        
        ## Lines for testing
        # run.TypeCurr_=run.TypeCurr; subDirCurr_=subDirCurr; prefixCurr_=prefixCurr; studySettingz_=studySettingz; files_PoolLabel=outL.AllGroups$files; sFactor=sF; dirout.sF.sep=dirout.1F
        
        
        f1ptrn=paste0("^",f1,"_|_",f1,"_")
        
        files<-nams<-nams2<-c()
        
        if (run.TypeCurr_=="profiling"){ 
                filePool=files_PoolLabel[1]
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                
                ix1=grep(f1ptrn, tmp_Pool$sample)
                tmp_Pool=tmp_Pool[ix1,]
                
                groupN=length(table(tmp_Pool$group))
                
                filePool2G=paste0(dirout.sF.sep,"/tmpTPool.csv")
                write.csv(tmp_Pool,file=filePool2G, row.names = FALSE)
                
                files=c(filePool2G)
                nams=c("Total pool size")
                nams2=c("T.pool")  ### label for figure file name
        }
        
        if (run.TypeCurr_=="tracer"){ 
                
                filePool=files_PoolLabel[1]
                fileLabel=files_PoolLabel[2]
                
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                tmp_Label=read.csv(fileLabel,stringsAsFactors = F)
                
                ix1=grep(f1ptrn, tmp_Pool$sample)
                tmp_Pool=tmp_Pool[ix1,]
                
                ix2=grep(f1ptrn, tmp_Label$sample)
                tmp_Label=tmp_Label[ix2,]
                
                groupN=c(length(table(tmp_Pool$group)), length(table(tmp_Label$group)))
                
                filePool1F=paste0(dirout.sF.sep,"/tmpTPool.csv")
                fileLabel1F=paste0(dirout.sF.sep,"/tmpTLabel.csv") 
                write.csv(tmp_Pool,file=filePool1F, row.names = FALSE)
                write.csv(tmp_Label,file=fileLabel1F, row.names = FALSE)
                
                files=c(filePool1F,fileLabel1F)
                nams=c("Total pool size","Total labelling")
                nams2=c("T.pool","T.label")  ### label for figure file name
        }
        
        return(list("files"=files,"nams"=nams,"nams2"=nams2,"groupN"=groupN))
        
}


#' Prepares the input files for two study factors statistical analyses (two factors comparison)
#'@param run.TypeCurr_ variable indicating whether it is a tracer or profiling study
#'@param subDirCurr_ the path to the current output directory
#'@param prefixCurr_ the current prefix, it is the name of the output directory
#'@param studySettingz_ list of dataframes containing sample info, group info, study factor info etc
#'@param files_PoolLabel full paths to the input Total PoolSize and Total Labelling tables
#'@param f1 the first factor 
#'@param f2 the second facotr
#'@param dirout.2F path to the output directory where results will be saved
#'
prep.Metabo2F_Stats_all<-function(run.TypeCurr_, subDirCurr_, prefixCurr_, studySettingz_, files_PoolLabel, f1 , f2, dirout.2F){
        
        ## Lines for testing
        # run.TypeCurr_=run.TypeCurr; subDirCurr_=subDirCurr; prefixCurr_=prefixCurr; studySettingz_=studySettingz; files_PoolLabel=outL.AllGroups$files; f1=sF1 ; f2=sF2 ;  dirout.2F=dirout.2F
        
        
        f1ptrn=paste0("^",f1,"_|_",f1,"_")
        f2ptrn=paste0("^",f2,"_|_",f2,"_")
        
        files<-nams<-nams2<-c()
        
        if (run.TypeCurr_=="profiling"){ 
                
                filePool=files_PoolLabel[1]
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)

                ix11=grep(f1ptrn,tmp_Pool$sample)
                ix12=grep(pattern=f2ptrn,tmp_Pool$sample)
                
                tmp_Pool[ix11,"group"]=f1
                tmp_Pool[ix12,"group"]=f2
                tmp_Pool=tmp_Pool[c(ix11,ix12),]
                
                groupN=length(table(tmp_Pool$group))
                
                filePool2F=paste0(dirout.2F,"/tmpTPool.csv")
                write.csv(tmp_Pool,file=filePool2F, row.names = FALSE)
                
                files=c(filePool2F)
                nams=c("Total pool size")
                nams2=c("T.pool")  ### label for figure file name
        }
        
        if (run.TypeCurr_=="tracer"){ 
                
                filePool=files_PoolLabel[1]
                fileLabel=files_PoolLabel[2]
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                tmp_Label=read.csv(fileLabel,stringsAsFactors = F)
                
                ix11=grep(f1ptrn,tmp_Pool$sample)
                ix12=grep(f2ptrn,tmp_Pool$sample)
                tmp_Pool[ix11,"group"]=f1
                tmp_Pool[ix12,"group"]=f2
                tmp_Pool=tmp_Pool[c(ix11,ix12),]
                
                ix21=grep(f1ptrn,tmp_Label$sample)
                ix22=grep(f2ptrn,tmp_Label$sample)
                
                tmp_Label[ix21,"group"]=f1
                tmp_Label[ix22,"group"]=f2
                tmp_Label=tmp_Label[c(ix21,ix22),]
                
                groupN=c(length(table(tmp_Pool$group)), length(table(tmp_Label$group)))
                
                filePool2F=paste0(dirout.2F,"/tmpTPool.csv")
                fileLabel2F=paste0(dirout.2F,"/tmpTLabel.csv") 
                write.csv(tmp_Pool,file=filePool2F, row.names = FALSE)
                write.csv(tmp_Label,file=fileLabel2F, row.names = FALSE)
                
                files=c(filePool2F,fileLabel2F)
                nams=c("Total pool size","Total labelling")
                nams2=c("T.pool","T.label")  ### label for figure file name
        }
        
        return(list("files"=files,"nams"=nams,"nams2"=nams2, "groupN"=groupN ))
}


#'Prepares the input files for two study factors pathway analyses 
#'@param run.TypeCurr_ variable indicating whether it is a tracer or profiling study
#'@param subDirCurr_ the path to the current output directory
#'@param prefixCurr_ the current prefix, it is the name of the output directory
#'@param compNotesDF_ the CompoundNotes table frame (MIDA- and MOLAR-specific)
#'@param studySettingz_ list of dataframes containing sample info, group info, study factor info etc
#'@param files_PoolLabel full paths to the input Total PoolSize and Total Labeling tables
#'@param f1 the first factor 
#'@param f2 the second facotr
#'@param dirout.2F path to the output directory where results will be saved
#'
prep.Metabo2F_Pathway_all<-function(run.TypeCurr_, subDirCurr_, prefixCurr_, compNotesDF_, studySettingz_, files_PoolLabel, f1 , f2, dirout.2F){
        
        ## Lines for testing
        # run.TypeCurr_=run.TypeCurr; subDirCurr_=subDirCurr; prefixCurr_=prefixCurr; studySettingz_=studySettingz; files_PoolLabel=outL.AllGroups$files; f1=sF1 ; f2=sF2
        
        
        f1ptrn=paste0("^",f1,"_|_",f1,"_")
        f2ptrn=paste0("^",f2,"_|_",f2,"_")
        
        files<-nams<-nams2<-c()
        
        if (run.TypeCurr_=="profiling"){
                
                filePool=files_PoolLabel[1]
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)

                ix11=grep(f1ptrn,tmp_Pool$sample)
                ix12=grep(f2ptrn,tmp_Pool$sample)
                tmp_Pool[ix11,"group"]=f1
                tmp_Pool[ix12,"group"]=f2
                
                groupN=length(table(tmp_Pool$group))
                
                tmp_Pool=tmp_Pool[c(ix11,ix12),]
                
                filePool2F_kegg=paste0(dirout.2F,"/tmpTPool_kegg.csv")
                
                tmp_Pool_kegg<-sub_CmpNames_with_KeggIDs(tmp_Pool, compNotesDF_)
                write.csv(tmp_Pool_kegg,file=filePool2F_kegg, row.names = FALSE)
                
                files_kegg=c(filePool2F_kegg)
                nams=c("Total pool size")
                nams2=c("T.pool")  ### label for figure file name
        }
        
        if (run.TypeCurr_=="tracer"){ 
                
                filePool=files_PoolLabel[1]
                fileLabel=files_PoolLabel[2]
                tmp_Pool=read.csv(filePool,stringsAsFactors = F)
                tmp_Label=read.csv(fileLabel,stringsAsFactors = F)
                
                ix11=grep(f1ptrn,tmp_Pool$sample)
                ix12=grep(f2ptrn,tmp_Pool$sample)
                tmp_Pool[ix11,"group"]=f1
                tmp_Pool[ix12,"group"]=f2
                tmp_Pool=tmp_Pool[c(ix11,ix12),]

                ix21=grep(f1ptrn,tmp_Label$sample)
                ix22=grep(f2ptrn,tmp_Label$sample)
                tmp_Label[ix21,"group"]=f1
                tmp_Label[ix22,"group"]=f2
                tmp_Label=tmp_Label[c(ix21,ix22),]

                groupN=c(length(table(tmp_Pool$group)), length(table(tmp_Label$group)))
                
                filePool2F_kegg=paste0(dirout.2F,"/tmpTPool_kegg.csv")
                fileLabel2F_kegg=paste0(dirout.2F,"/tmpTLabel_kegg.csv") 
                
                tmp_Pool_kegg<-sub_CmpNames_with_KeggIDs(tmp_Pool, compNotesDF_)
                tmp_Label_kegg<-sub_CmpNames_with_KeggIDs(tmp_Label, compNotesDF_)
                write.csv(tmp_Pool_kegg,file=filePool2F_kegg, row.names = FALSE)
                write.csv(tmp_Label_kegg,file=fileLabel2F_kegg, row.names = FALSE)
                
                files_kegg=c(filePool2F_kegg, fileLabel2F_kegg)
                nams=c("Total pool size","Total labelling")
                nams2=c("T.pool","T.label")  ### label for figure file name
        }
        return(list("files_kegg"=files_kegg, "nams"=nams,"nams2"=nams2,"groupN"=groupN))
}


# Runs statistical analyses, generates tables and several figures including PCA plots and Heatmaps for a multi-group analysis
#'@param files_ full paths to the input Total PoolSize and Total Labeling tables
#'@param nams_ labels corresponding to the input files e.g. "Total pool size" "Total labeling"
#'@param nams2_  short name corresponding to the input files e.g. "T.pool"  "T.label"
#'@param prefix__  the name of given to current analysis, e.g. MC00134_Comb,
#'@param nameData_  the type of the statistical analysis performed e.g "AllGroups"
#'@param topMetaboFDR__ FDR threshold to determining metabolites with significant difference between groups 
#' 
run.MetaboAllG_Stats_all<-function(files_, 
                             nams_, 
                             nams2_, 
                             prefix__, 
                             nameData_, 
                             topMetaboFDR__){
        
        #Testing lines
        # IdxMOLAR
        # files_=outL.AllGroups$files; nams_=outL.AllGroups$nams; nams2_=outL.AllGroups$nams2; prefix__=prefixCurr; nameData_="AllGroups";  topMetaboFDR__=topMetaboFDR_
        # idxMOLAR
        # prefix__=prefix.1Multi
        
        fn=1
        
        for (ii in 1:length(files_)){
                
                fileMetabo=files_[ii]
                
                nameMultiG_=paste0(nams_[ii],"\n",nameData_)
                
                
                options( scipen = 0,digits = 6 )
                
                
                mSet<-InitDataObjects("pktable", "stat", FALSE)
                mSet<-Read.TextData(mSet, fileMetabo, "rowu", "disc")
                
                mSet<-SanityCheckData(mSet)
                mSet<-ReplaceMin(mSet)
                mSet<-PreparePrenormData(mSet)
                
                mSet<-Normalization(mSet, rowNorm="None", transNorm="LogNorm", scaleNorm="AutoNorm", ref=NULL, ratio=FALSE, ratioNum=20)
                mSet<-PlotNormSummary_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_norm_summary_"), title=nameMultiG_,"pdf",100, width=8.5)           ;fn=fn+1 
                mSet<-PlotNormSummary_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_norm_summary_"), title=nameMultiG_,"png",100, width=8.5)           
                mSet<-PlotSampleNormSummary_bm(mSet, paste0(f00(fn), "_", nams2_[ii], "_snorm_summary_"), title=nameMultiG_ ,"pdf", 100, width=8.5)                 ;fn=fn+1 
                mSet<-PlotSampleNormSummary_bm(mSet, paste0(f00(fn), "_", nams2_[ii], "_snorm_summary_"), title=nameMultiG_ ,"png", 100, width=8.5)
                
                
                if (length(mSet$dataSet$prenorm.feat.nms)<5){
                        message("           PCA plots skipped for : ", prefix__," ~ ",nams_[ii], ", due to less than 5 compounds")
                }else{

                        message("           Generating PCA plots for : ", prefix__," ~ ",nams_[ii])
                        
                        mSet<-PCA.Anal_bm(mSet,nams2_[ii])
                        mSet<-PlotPCAPairSummary_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_pair_"), title=nameMultiG_ , "pdf", dpi=72, width=NA, 5)           ;fn=fn+1 
                        mSet<-PlotPCAPairSummary_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_pair_"), title=nameMultiG_ , "png", 72, width=NA, 5)         
                        mSet<-PlotPCAScree_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_scree_"), title=nameMultiG_ , "pdf", 72, width=NA, 5)           ;fn=fn+1 
                        mSet<-PlotPCAScree_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_scree_"), title=nameMultiG_ , "png", 72, width=NA, 5)           
                        mSet<-PlotPCA2DScore_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_score2d_"), title=nameMultiG_ , "pdf", 72, width=NA,  1,2,0.95,0,0)           ;fn=fn+1 
                        mSet<-PlotPCA2DScore_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_score2d_"), title=nameMultiG_ , "png", 72, width=NA,  1,2,0.95,0,0)         
                        mSet<-PlotPCALoading_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_loading_"), title=nameMultiG_ , "pdf", 72, width=NA, 1,2)           ;fn=fn+1 
                        mSet<-PlotPCALoading_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_loading_"), title=nameMultiG_ , "png", 72, width=NA, 1,2)           
                        mSet<-PlotPCABiplot_bm(mSet, paste0("N_",nams2_[ii],"_pca_biplot_"), title=nameMultiG_ , "pdf", 72, width=NA, 1,2)           ;fn=fn+1 
                        mSet<-PlotPCA3DScoreImg_bm(mSet, paste0(f00(fn), "_",nams2_[ii],"_pca_score3d_"), title=nameMultiG_ , "pdf", 72, width=NA, 1,2,3, 40)           ;fn=fn+1 
                }       
                
                mSet<-ANOVA.Anal_all_bm(mSet,nonpar=FALSE,thresh=0.05,all_results=TRUE, name_=nams2_[ii])
                
                if (mSet$analSet$aov$sig.num==0){
                        message("           ANOVA shows no significantly different compounds for ", nameMultiG_)
                        message("           Now generating heatmaps for all compound.\n")
                        
                        ## HEATMAP of all metabolites
                        # All samples
                        mSet<-PlotHeatMapStar_bm(mSet, imgName=paste0(f00(fn), "_",nams2_[ii],"_htmp_all_metabolites_vs_all_samples_"),
                                                 title=paste0(nameMultiG_,"\n All metabolites vs all samples (* FDR<",topMetaboFDR__,")"),
                                                 format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, F, star=TRUE, topMetaboFDR__)           ;fn=fn+1 
                        # Group Averages
                        mSet<-PlotHeatMapStar_bm(mSet,paste0(f00(fn), "_",nams2_[ii],"_htmp_all_metabolite_vs_group_av_"), 
                                                 title=paste0(nameMultiG_,"\n All metabolites vs group av. (* FDR<",topMetaboFDR__,")"),
                                                 format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, T, star=TRUE, topMetaboFDR__)           ;fn=fn+1 
                        
                }else{
                
                        # HEATMAP of top sgnificant metabolistes (using P-values)
                        message("           Generating heatmaps for  : ",nameData_," ~ ", nams_[ii])
                        
                        
                        #All samples
                        mSet<-PlotSubHeatMapFDR_bm(mSet, imgName=paste0(f00(fn), "_",nams2_[ii],"_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_all_samples_"), 
                                                   title=paste0(nameMultiG_,"\nMetabolites vs all samples (FDR< ",topMetaboFDR__,")"),
                                                   format="pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", "overview", F, T, T, F, topMetaboFDR__)           ;fn=fn+1 
                        
                        # Group Average
                        mSet<-PlotSubHeatMapFDR_bm(mSet, imgName=paste0(f00(fn), "_",nams2_[ii],"_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_group_av_"), 
                                                   title=paste0(nameMultiG_,"\nMetabolites vs group av. (FDR< ",topMetaboFDR__,")"),
                                                   format="pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", "overview", F, T, T, T, topMetaboFDR__)           ;fn=fn+1 
                        
                        ## HEATMAP of all metabolites
                        # All samples
                        mSet<-PlotHeatMapStar_bm(mSet, imgName=paste0(f00(fn), "_",nams2_[ii],"_htmp_all_metabolites_vs_all_samples_"),
                                                 title=paste0(nameMultiG_,"\n All metabolites vs all samples (* FDR<",topMetaboFDR__,")"),
                                                 format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, F, star=TRUE, topMetaboFDR__)           ;fn=fn+1 
                        # Group Averages
                        mSet<-PlotHeatMapStar_bm(mSet,paste0(f00(fn), "_",nams2_[ii],"_htmp_all_metabolite_vs_group_av_"), 
                                                 title=paste0(nameMultiG_,"\n All metabolites vs group av. (* FDR<",topMetaboFDR__,")"),
                                                 format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, T, star=TRUE, topMetaboFDR__)           ;fn=fn+1 
                }      
                
        }
        
        unlink( list.files(pattern = "\\.json$"))
}


#'Runs statistical analyses, generates tables and several figures including PCA plots and Heatmaps for a two-group comparison
#'@param files_ full paths to the input Total PoolSize and Total Labeling tables
#'@param nams_ labels corresponding to the input files e.g. "Total pool size" "Total labeling"
#'@param nams2_  short name corresponding to the input files e.g. "T.pool"  "T.label"
#'@param nameData_ name of the output directory for this pairwise comparison "Stats_2Groups_exp1_ARID1A_EV_veh_v_exp1_ARID1A_OE_insulin"
#'@param prefix__  the name of given to current analysis, e.g. MC00134_Comb,
#'@param topMetaboFDR__ FDR threshold to determining metabolites with significant difference between groups 
#'@param  normMethod__ normalization method ("none", "quantile normalization", "normalization by reference feature" etc.),default set to "none"
#'@param transMethod__ data transformation method ("log normalization", "cubic-root normalization")
#'@param  scalingMethod__ data scaling method (mean-centering, auto-scaling, pareto-scaling or range-scaling)
#'
run.Metabo2G_Stats_all<-function(   files_,
                                   nams_, 
                                   nams2_, 
                                   nameData_, 
                                   prefix__, 
                                   topMetaboFDR__,
                                   normMethod__,
                                   transMethod__,
                                   scalingMethod__
                                  
){
     
        ## Testing lines
        test=FALSE
        if (test==TRUE){
                
                #IDX
                files_=files; nams_=nams; nams2_=nams2; prefix__=prefix.1Multi; nameData_=name2G; topMetaboFDR__=topMetaboFDR_
                #MOLAR MDIA
                files_=files; nams_=nams; nams2_=nams2; prefix__=prefixCurr; nameData_=name2G; topMetaboFDR__=topMetaboFDR_
                files_=files; nams_=nams; nams2_=nams2; prefix__=prefixCurr; nameData_="AllGroups"; topMetaboFDR__=topMetaboFDR_ 
                
                normMethod__=normMethod_; transMethod__=transMethod_;  scalingMethod__=scalingMethod_
        }
        
     
        for (ii in 1:length(files_)){
                
                fn.2=1
                
                fileMetabo=files_[ii]
                nameG_=paste0(nams_[ii],"\n",nameData_)
                
                options( scipen = 0,digits = 6 )
                
                mSet<-InitDataObjects("pktable", "stat", FALSE)
                mSet<-Read.TextData(mSet, fileMetabo, "rowu", "disc")
                
                mSet<-SanityCheckData(mSet)
                mSet<-ReplaceMin(mSet)
                
                mSet<-PreparePrenormData(mSet)
                mSet<-Normalization(mSet, rowNorm="None", transNorm="LogNorm", scaleNorm="AutoNorm", ref=NULL, ratio=FALSE, ratioNum=20)
                mSet<-PlotNormSummary_bm(mSet, paste0(f00(ii), "_", f00(fn.2),"_", nams2_[ii], "_norm_summary_"), title=nameG_,"pdf",100, width=8.5)     ;fn.2=fn.2+1 
                mSet<-PlotNormSummary_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_norm_summary_"), title=nameG_,"png",100, width=8.5)            
                mSet<-PlotSampleNormSummary_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_snorm_summary_"), title=nameG_,"pdf", 100, width=8.5)      ;fn.2=fn.2+1 
                mSet<-PlotSampleNormSummary_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_snorm_summary_"), title=nameG_,"png", 100, width=8.5)    
          
                if (length(mSet$dataSet$prenorm.feat.nms)<5){
                        message("           PCA plots skipped for : ",nams_[ii], ", due to less than 5 compounds")
                }else{
                
                        message("           Generating PCA plots for : ",nams_[ii])
                           
                        mSet<-PCA.Anal_bm(mSet,nams2_[[ii]])
                        mSet<-PlotPCAPairSummary_bm(mSet, paste0(f00(ii), "_", f00(fn.2) ,"_", nams2_[ii],"_pca_pair_"), title=nameG_ , "pdf", dpi=72, width=NA, 5)      ;fn.2=fn.2+1 
                        mSet<-PlotPCAPairSummary_bm(mSet, paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_pca_pair_"), title=nameG_ , "png", 72, width=NA, 5)  
                        mSet<-PlotPCAScree_bm(mSet, paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_pca_scree_"), title=nameG_ , "pdf", 72, width=NA, 5)       ;fn.2=fn.2+1 
                        mSet<-PlotPCAScree_bm(mSet, paste0(f00(ii), "_", f00(fn.2) ,"_", nams2_[ii],"_pca_scree_"), title=nameG_ , "png", 72, width=NA, 5)
                        mSet<-PlotPCA2DScore_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_pca_score2d_"), title=nameG_ , "pdf", 72, width=NA,  1,2,0.95,0,0)      ;fn.2=fn.2+1 
                        mSet<-PlotPCA2DScore_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_pca_score2d_"), title=nameG_ , "png", 72, width=NA,  1,2,0.95,0,0)      
                        mSet<-PlotPCALoading_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_pca_loading_"), title=nameG_ , "pdf", 72, width=NA, 1,2)      ;fn.2=fn.2+1 
                        mSet<-PlotPCALoading_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_pca_loading_"), title=nameG_ , "png", 72, width=NA, 1,2)      
                        mSet<-PlotPCABiplot_bm(mSet, paste0("N_",nams2_[ii],"_pca_biplot_"), title=nameG_ , "pdf", 72, width=NA, 1,2)  
                        mSet<-PlotPCA3DScoreImg_bm(mSet, paste0("N_",nams2_[ii],"_pca_score3d_"), title=nameG_ , "pdf", 72, width=NA, 1,2,3, 40)     
                }
                
                mSet<-FC.Anal.unpaired_all_bm(mSet, 2.0, 0, nams2_[ii] )                                                                         ## Fold Change Analysist
                mSet<-PlotFC_bm(mSet, paste0("N_",nams2_[ii],"_fc_"), title=nameG_, "pdf", 72, width=NA)    
                mSet<-Volcano.Anal_all_bm(mSet, FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw", nams2_[ii])                                         ## Volcano Analysist
                
                mSet<-PlotVolcano_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_volcano_"),1, title=nameG_,  "pdf", 72, width=NA)      ;fn.2=fn.2+1 
                mSet<-PlotVolcano_bm(mSet, paste0(f00(ii), "_", f00(fn.2) , "_", nams2_[ii],"_volcano_"),1, title=nameG_,  "png", 72, width=NA)    
                mSet<-Ttests.Anal_all_bm(mSet, F, 0.05, FALSE, TRUE, TRUE,nams2_[[ii]])                           ## FDR by Benjamini-Hochberg                                      ## T-test
                mSet<-PlotTT_bm(mSet, paste0("N_",nams2_[ii],"_t-test_"), title=nameG_, "pdf", 72, width=NA)      ;fn.2=fn.2+1 
                
                
                if (mSet$analSet$tt$sig.num==0){
                        
                        message(paste0("           T-test shows no significantly different compounds for ",nams_[ii]))
                        
                        message("           Now generating heatmaps for all compound.")
                        mSet<-PlotHeatMapStar_bm(mSet, imgName=paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_htmp_all_metabolites_vs_all_samples_"),
                                                 title=paste0(nameG_,"\n All metabolites vs all samples (* FDR<",topMetaboFDR__,")"),
                                                 format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, F, star=TRUE, topMetaboFDR__)      ;fn.2=fn.2+1 
                        # Group Averages
                        #mSet<-PlotHeatMapStar_bm(mSet,paste0(f00(ii), "_", f00(fn.2) ,"_", nams2_[ii],"_htmp_all_metabolite_vs_group_av_"), 
                        #                                  title=paste0(nameG_,"\n All metabolites vs group av. (* FDR<",topMetaboFDR__,")"),
                        #                                  format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, T, star=TRUE, topMetaboFDR__)      ;fn.2=fn.2+1 
                        
                
                }else{
                
                        
                        # HEATMAP of top significant metabolites (using P-values)
                        message("           Generating heatmaps for  : ", nameData_," ~ ",nams_[ii])
                        
                        #All samples
                        mSet<-PlotSubHeatMapFDR_bm(mSet, imgName=paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_all_samples_"), 
                                                   title=paste0(nameG_,"\nMetabolites vs all samples (FDR< ",topMetaboFDR__,")"),
                                                   format="pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", "overview", F, T, T, F, topMetaboFDR__)      ;fn.2=fn.2+1 
                        mSet<-PlotSubHeatMapFDR_bm(mSet, imgName=paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_all_samples_"), 
                                                   title=paste0(nameG_,"\nMetabolites vs all samples (FDR< ",topMetaboFDR__,")"),
                                                   format="png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", "overview", F, T, T, F, topMetaboFDR__)     
                        
                        # Group Average
                        # mSet<-PlotSubHeatMapFDR_bm(mSet, imgName=paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_group_av_"), 
                        #                            title=paste0(nameG_,"\nMetabolites vs group av. (FDR< ",topMetaboFDR__,")"),
                        #                            format="pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", "overview", F, T, T, T, topMetaboFDR__)      ;fn.2=fn.2+1 
                        # mSet<-PlotSubHeatMapFDR_bm(mSet, imgName=paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_htmp_FDR_",topMetaboFDR__,"_metabolite_vs_group_av_"), 
                        #                            title=paste0(nameG_,"\nMetabolites vs group av. (FDR< ",topMetaboFDR__,")"),
                        #                            format="png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", "overview", F, T, T, T, topMetaboFDR__) 
                        
                        
                        ## HEATMAP of all metabolites
                        # All samples
                        mSet<-PlotHeatMapStar_bm(mSet, imgName=paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_htmp_all_metabolites_vs_all_samples_"),
                                                 title=paste0(nameG_,"\n All metabolites vs all samples (* FDR<",topMetaboFDR__,")"),
                                                 format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, F, star=TRUE, topMetaboFDR__)      ;fn.2=fn.2+1 
                        mSet<-PlotHeatMapStar_bm(mSet, imgName=paste0(f00(ii), "_", f00(fn.2) ,"_",nams2_[ii],"_htmp_all_metabolites_vs_all_samples_"),
                                                 title=paste0(nameG_,"\n All metabolites vs all samples (* FDR<",topMetaboFDR__,")"),
                                                 format="png", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, F, star=TRUE, topMetaboFDR__)      ;fn.2=fn.2+1 
                        
                        # Group Averages
                        # mSet<-PlotHeatMapStar_bm(mSet,paste0(f00(ii), "_", f00(fn.2) ,"_", nams2_[ii],"_htmp_all_metabolite_vs_group_av_"), 
                        #                          title=paste0(nameG_,"\n All metabolites vs group av. (* FDR<",topMetaboFDR__,")"),
                        #                          format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, T, star=TRUE, topMetaboFDR__)      ;fn.2=fn.2+1 
                        # mSet<-PlotHeatMapStar_bm(mSet,paste0(f00(ii), "_", f00(fn.2) ,"_", nams2_[ii],"_htmp_all_metabolite_vs_group_av_"), 
                        #                          title=paste0(nameG_,"\n All metabolites vs group av. (* FDR<",topMetaboFDR__,")"),
                        #                          format="pdf", dpi=72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, T, star=TRUE, topMetaboFDR__)    
                }        
                unlink( list.files(pattern = "\\.json$"))
                
                message(paste0("           Statistical analyses for two sample groups [ ",gsub(pattern="\n"," : ",nameG_)," ] is completed"))
        }
        
        return(fn.2)
             
}


#'Runs pathway analyses for a two two-group comparison
#'@param files_kegg_ contains paths to file containing PoolSize and Total Labeling data
#'@param nams_ labels corresponding to the input files e.g. "Total pool size" "Total labeling"
#'@param nams2_ short name corresponding to the input files e.g. "T.pool"  "T.label"
#'@param nameData_ !!! name of the output directory for this pairwise comparison " ... Stats_2Groups_exp1_ARID1A_EV_veh_v_exp1_ARID1A_OE_insulin"
#'@param prefix__ the name of given to current analysis, e.g. MC00134_Comb,
#'@param topMetaboFDR__ FDR threshold to determining metabolites with significant difference between groups 
#'@param normMethod__ normalization method ("none", "quantile normalization", "normalization by reference feature" etc.),default set to "none"
#'@param transMethod__ data transformation method ("log normalization", "cubic-root normalization")
#'@param scalingMethod__ data scaling method (mean-centering, auto-scaling, pareto-scaling or range-scaling)
#'@param KeggPathwayAnalysisCheck_ determines whether or not do perform pathway analysis
#'@param KeggPathwayEnrichMethod__  Kegg pathway enrichment method (QEA or ORA)
#'@param KeggPathwayOrganism__ the organism being studied
#'@param KeggPathwayFDR__ FDR threshold for significantly up-regulated pathways between two groups
#'@param fn_2 the figure number that follows from statistical analyses, it helps numbering the figures such one can control their order in the report 
#'@param lbl_ the number corresponding the pipeline execution step

run.Metabo2G_Pathway_all<-function( files_kegg_,
                                   nams_, 
                                   nams2_, 
                                   nameData_, 
                                   prefix__, 
                                 
                                   topMetaboFDR__,
                                   normMethod__,
                                   transMethod__,
                                   scalingMethod__,
                                 
                                   KeggPathwayAnalysisCheck_,
                                   KeggPathwayEnrichMethod__,
                                   KeggPathwayOrganism__,
                                   KeggPathwayFDR__,
                                   fn_2=fn_2,
                                   
                                   lbl_="8"
){
        
        ## Testing lines
        test=FALSE
        if (test==TRUE){
                #IDX
                files_kegg_=files_kegg; nams_=nams; nams2_=nams2; prefix__=prefix.1Multi; nameData_=name2G; topMetaboFDR__=topMetaboFDR_
                #MOLAR MIDA
                
                files_kegg_=outL.2G.Pathway$files_kegg; nams_=outL.2G.Pathway$nams; nams2_=outL.2G.Pathway$nams2
                prefix__=prefixCurr; nameData_=name2G; topMetaboFDR__=topMetaboFDR_
                
                normMethod__=normMethod_; transMethod__=transMethod_;  scalingMethod__=scalingMethod_; KeggPathwayEnrichMethod__=KeggPathwayEnrichMethod_; KeggPathwayAnalysisCheck_=TRUE; KeggPathwayOrganism__=KeggPathwayOrganism_; KeggPathwayFDR__=KeggPathwayFDR_
                
                cat(paste(c(files_kegg_,
                            nams_,
                            topMetaboFDR__,
                            normMethod__,
                            transMethod__,
                            scalingMethod__,
                            KeggPathwayEnrichMethod__,
                            KeggPathwayAnalysisCheck_,
                            KeggPathwayOrganism__,
                            KeggPathwayFDR__),collabse="\n"))
        }
        
        
        if(is.null(fn_2)){fn_2=1}
       
        for (ii in 1:length(files_kegg_)){
                fn.2=fn_2
                
                fileMetabo_kegg=files_kegg_[ii]
                
                nameG_=paste0(nams_[ii],"\n",nameData_)
                
                options( scipen = 0,digits = 6 )
                
                message( "           Kegg pathway analysis based on FC enrichment [ ",gsub(pattern="\n"," : ",nameG_)," ] starts ...")
                
                namesTopPathways<-c()
                
                if  (KeggPathwayEnrichMethod__=="ORA"){ # Over Representation Analysis (CalculateHyperScore)
                        
                        mSetORA<-InitDataObjects("conc", "msetora", FALSE)
                        mSetORA<-Read.TextData(mSetORA, fileMetabo_kegg, "rowu", "disc")
                        
                        mSetORA<-SanityCheckData(mSetORA)
                        mSetORA<-ReplaceMin(mSetORA)
                        mSetORA<-PreparePrenormData(mSetORA)
                        #mSetORA<-Normalization(mSetORA, "NULL", "LogNorm", "MeanCenter", ratio=FALSE, ratioNum=20)
                        mSetORA<-Normalization(mSetORA, rowNorm="None", transNorm="LogNorm", scaleNorm="AutoNorm", ref=NULL, ratio=FALSE, ratioNum=20)
                        
                        mSetORA<-FC.Anal.unpaired_all_bm(mSetORA, 2.0, 0, nams2_[ii] ) 
                        mSetORA<-Ttests.Anal_all_bm(mSetORA, F, 0.05, FALSE, TRUE, TRUE,nams2_[[ii]])  
                        
                        
                        ### quit the test if there are not significant features
                        if (mSetORA$analSet$tt$sig.num==0){message("           Warning -  no significantly different features are found. ") ;next}
                        
                        
                        Log2FC_up=mSetORA$dataSet$cmpd[which(mSetORA$analSet$fc$fc.log>0 & mSetORA$analSet$tt$inx.imp==TRUE)]
                        Log2FC_dwn=mSetORA$dataSet$cmpd[which(mSetORA$analSet$fc$fc.log<0 & mSetORA$analSet$tt$inx.imp==TRUE)]
                        
                        if (length(Log2FC_up)==0){message("           Error -  ORA failed, no significanly different commpounds were found ");return() }
                        
                        mSetORA<-Setup.MapData(mSetORA, Log2FC_up);
                        mSetORA<-CrossReferencing(mSetORA, "kegg");
                        
                        # check if mapping worked and return
                        if (is.numeric(mSetORA)){message("           Error -  ORA failed due to used unconvetional compound names"); return()}
                        
                        mSetORA<-CreateMappingResultTable(mSetORA)
                     
                        mSetORA<-SetMetabolomeFilter(mSetORA, F);
                        mSetORA<-SetCurrentMsetLib_bm(mSetORA, "kegg_pathway", 2);
                        
                        
                        message("           Metabolite sets enrichment - Over Represationation Analysis (MSE-QEA) starts ...")
                        message("           MSE-QEA using the significantly upreagulated compounds")
                        

                        result = tryCatch({
                                mSetORA<-CalculateHyperScore(mSetORA)
                        }, warning = function(w) {
                                message(w)
                                message("           Error -  MSE-ORA failed via MetaboAnalyst web-server, please check the code ---")
                                return()
                        }, error = function(e) {
                                message(e)
                                message("           Error -  MSE-ORA failed via MetaboAnalyst web-server, please check the code ---")
                                return()
                        }, finally = {
                                 message("           Computation of Metabolite set (ORA) compled ")
                        })
                        
                        
                        mSetORA<-PlotORA_bm(mSetORA, paste0(f00(ii),"_",f00(fn.2),"_metabolite_set_ora_"), tittle=paste0("ORA (",nams_[ii],")"), "net", "pdf", 72, width=NA,KeggPathwayFDR__); fn.2=fn.2+1
                        mSetORA<-PlotEnrichDotPlot_bm(mSetORA, "ora", paste0(f00(ii),"_",f00(fn.2),"_Enrich_Anal_ORA_"),tittle=paste0("ORA (",nams_[ii],")"), "pdf", 72, width=NA); fn.2=fn.2+1
                        
                        
                        options(scipen = 999)        
                        namesTopPathways <-rownames(mSetORA$analSet$ora.mat)[ which(mSetORA$analSet$ora.mat[,"FDR"]< KeggPathwayFDR__)]
                        
                        if (length(namesTopPathways)>0){
                                message("   - Getting the significantly enriched pathways ")
                                message("     Found ", length(namesTopPathways)," significant pathway(s) for ",gsub(pattern="\n"," : ",nameG_)," at FDR<",KeggPathwayFDR__)
                                message(paste0(capture.output(as.data.frame(mSetORA$analSet$ora.mat[ which( mSetORA$analSet$ora.mat[,"FDR"] < KeggPathwayFDR__),])), collapse = "\n"))
                        }
                        
                        #The other method table
                        mSetORA<-CalculateGlobalTestScore(mSetORA)
                         
                }
                        
                if(KeggPathwayEnrichMethod__=="QEA"){   ## Quantitative Enrichment Analysis (CalculateGlobalTestScore)
                        
                        mSetQEA<-InitDataObjects("conc", "msetqea", FALSE)
                        mSetQEA<-Read.TextData(mSetQEA, fileMetabo_kegg, "rowu", "disc");
                        
                        mSetQEA<-SanityCheckData(mSetQEA)
                        mSetQEA<-ContainMissing(mSetQEA)
                        mSetQEA<-ReplaceMin(mSetQEA)
                        mSetQEA<-PreparePrenormData(mSetQEA)
                        #mSetORA<-Normalization(mSetORA, "NULL", "LogNorm", "MeanCenter", ratio=FALSE, ratioNum=20)
                        mSetQEA<-Normalization(mSetQEA, rowNorm="None", transNorm="LogNorm", scaleNorm="AutoNorm", ref=NULL, ratio=FALSE, ratioNum=20)
                       
                        mSetQEA<-FC.Anal.unpaired_all_bm(mSetQEA, 2.0, 0, nams2_[ii] ) 
                        mSetQEA<-Ttests.Anal_all_bm(mSetQEA, F, 0.05, FALSE, TRUE, TRUE,nams2_[[ii]])  
                        
                        
                        ### quit the test if there are not significant features
                        if (mSetQEA$analSet$tt$sig.num==0){message("           Warning -  no significantly different features are found. ") ;next}
                        
                        mSetQEA<-CrossReferencing(mSetQEA, "kegg")
                        
                        # check if mapping worked and return
                        if (is.numeric(mSetQEA)){message("           Error -  QEA failed due to used unconvetional compound names"); return()}
                        
                        mSetQEA<-CreateMappingResultTable(mSetQEA)
                        mSetQEA<-PreparePrenormData(mSetQEA)
                        mSetQEA<-Normalization(mSetQEA, "NULL", "LogNorm", "MeanCenter", ratio=FALSE, ratioNum=20)
                        mSetQEA<-SetMetabolomeFilter(mSetQEA, F);
                        mSetQEA<-SetCurrentMsetLib(mSetQEA, "kegg_pathway", 2);
                        
                        message("           Metabolite sets enrichment - Quantitative Enrichment Analysis (MSE-QEA) starts ...")
                        
                        result = tryCatch({
                                mSetQEA<-CalculateGlobalTestScore(mSetQEA)
                        }, warning = function(w) {
                                message(w)
                                message("           Error -  MSE-QEA failed via MetaboAnalyst web-server, please check the code ---")
                                return()
                        }, error = function(e) {
                                message(e)
                                message("           Error -  MSE-QEA failed via MetaboAnalyst web-server, please check the code ---")
                                return()
                        }, finally = {
                                message("          MSE-QEA computation done")
                        })
                        
                        title=paste0(nameG_,"\nMetabolites vs all samples (FDR< ",topMetaboFDR__,")")
                        
                        mSetQEA<-PlotQEA.Overview_bm(mSetObj=mSetQEA, paste0(f00(ii),"_",f00(fn.2),"_MSet_Enrich_Anal_QEA_bar_"),tittle=paste0("QEA (",nams_[ii],")"), "net", "pdf", 72, width=NA,KeggPathwayFDR__);fn.2=fn.2+1
                        mSetQEA<-PlotEnrichDotPlot_bm(mSetQEA, "qea", paste0(f00(ii),"_",f00(fn.2),"_MSet_Enrich_Anal_QEA_dot_"),tittle=paste0("QEA (",nams_[ii],")"), "pdf", 72, width=NA);fn.2=fn.2+1
                        
                        
                        options(scipen = 999)        
                        namesTopPathways <-rownames(mSetQEA$analSet$qea.mat)[ which(mSetQEA$analSet$qea.mat[,"FDR"]< KeggPathwayFDR__)]
                        
                        if (length(namesTopPathways)>0){
                                message("   - Getting the significantly enriched pathways ")
                                message( "     Found ", length(namesTopPathways)," significant pathways for ",gsub(pattern="\n"," : ",nameG_)," at FDR<",KeggPathwayFDR__)
                                message(paste0(capture.output(as.data.frame(mSetQEA$analSet$qea.mat[ which( mSetQEA$analSet$qea.mat[,"FDR"] < KeggPathwayFDR__),])), collapse = "\n"))
                        }
                        
                        #The method's other table
                        mSetQEA<-CalculateHyperScore(mSetQEA)
                        
                }
                
                
                if (KeggPathwayEnrichMethod__==" Other pathway method not yest implemented !!!! CalculateQeaScore "){
                        #### add code here
                }
                
                
                if (length(namesTopPathways)>0){
                        
                        pathList <- keggList("pathway", KeggPathwayOrganism__); 
                        
                        pathList=gsub("\\)","",pathList); pathList=gsub("\\(","- ",pathList);
                        namesTopPathways=gsub("\\)","",namesTopPathways);
                        namesTopPathways=gsub("\\(","- ",namesTopPathways) 
                        
                        message("   - Singificant pathways IDs:")
                        TopPathwaysIDs=rep("",length(namesTopPathways))
                        for (ip in 1:length( namesTopPathways)){
                                nm=namesTopPathways[ip]
                                
                                if (nm=="Glycolysis or Gluconeogenesis"){nm="Glycolysis / Gluconeogenesis"}
                                ix=grep(nm,pathList)
                                nam=names(pathList)[ix]
                                message(nam)
                                TopPathwaysIDs[ip]=strsplit(nam,"path:")[[1]][2]
                        }
                        
                        #print(as.data.frame(pathList))
                        
                        ## Generate pathways heatmaps based on compound FC enrichement for Top Kegg Pathways based on FDR
                        if (KeggPathwayEnrichMethod__=="ORA"){
                                FC=mSetORA$analSet$fc$fc.log
                                names(FC)=mSetORA$dataSet$cmpd
                        }
                        if (KeggPathwayEnrichMethod__=="QEA"){
                                FC=mSetQEA$analSet$fc$fc.log
                                names(FC)=mSetQEA$dataSet$cmpd
                        }
                        
                        
                        message("\n   - Compound names mapped to Kegg ID")
                        
                        message(paste0("NAs : ",length(grep("^NA",names(FC))), " | Non-NAs : ",length(names(FC))-length(grep("^NA",names(FC)))))
                        
                        message(paste0(names(FC),collabse=" "))
                        
                       
                        topPath=5; if (topPath>length(TopPathwaysIDs)){topPath=length(TopPathwaysIDs)}
                        tryCatch(
                                for ( iz in 1:topPath){
                                        
                                        ID=TopPathwaysIDs[iz]
                                        message("\n     ~",iz,"~ Pathway : ", names(TopPathwaysIDs)[iz]," , ",ID )  
                                        pv.out<-pathview(cpd.data=FC,  pathway.id = ID, species = KeggPathwayOrganism__, out.suffix ="kegg.pathway", low="blue", mid="grey", high="red")
                                        
                                        file.rename(paste0(ID,".kegg.pathway.png"), paste0(f00(ii),"_",f00(fn.2),"_",ID,".kegg.pathway.png") )
                                        
                                        unlink(paste0(ID,".png")); unlink(paste0(ID,".xml"))
                                        
                                        #Make pdf
                                        pdf(paste0(ii,"_",fn.2,"_",ID,".kegg.pathway.pdf"),width=9, height=6)
                                        pngRaster <- readPNG(paste0(f00(ii),"_",f00(fn.2),"_",ID,".kegg.pathway.png"))
                                        grid.raster(pngRaster) 
                                        dev.off()
                                        fn.2=fn.2+1 
                                }
                        )
                        
                        message("           KeggPathway heatmap produced for ",gsub(pattern="\n"," : ",nameG_),"   \n")
                }else{
                        message("           No significantly enriched KeggPathways found for ",gsub(pattern="\n"," : ",nameG_),"  \n")
                }
        }
        unlink( list.files(pattern = "\\.json$"))
}


#'Performs statistical analyses for both tracer and profiling, for all groups, group pairwise comparison etc. 
#'It also performs pathway analysis in case a two-group comparison has been selected
#'It works for the idxMOLAR, MOLAR, MIDA pipelines
#'@param subDirCurr the main output directory path where all results placed e.g. "/Users/Analysis/MC00134/MC00134_Comb_out/"
#'@param prefixCurr prefix of the current output file e.g. MC00134_Comb
#'@param run.TypeCurr indicates whether it is a tracer or profiling experiment
#'@param compNotesDF  compound notes data frame (MIDA- and MOLAR-specific)
#'@param studySettingz list of dataframes containing sample info, group info, study factor info etc
#'@param topMetaboFDR_ FDR threshold to determining metabolites with significant difference between groups 
#'@param normMethod__ normalization method ("none", "quantile normalization", "normalization by reference feature" etc.),default set to "none"
#'@param transMethod__ data transformation method ("log normalization", "cubic-root normalization")
#'@param scalingMethod_ data scaling method (mean-scaling, auto-scaling, pareto-scaling, range-scaling)
#'@param KeggPathwayAnalysisCheck_ determines whether or not do perform pathway analysis
#'@param KeggPathwayEnrichMethod_  Kegg pathway enrichment method (MSEA,QEA or ORA)
#'@param KeggPathwayOrganism_ the organism being studied
#'@param KeggPathwayFDR_ FDR threshold for significantly up-regulated pathways between two groups
#'@param lbl_ the number corresponding the pipeline execution step
#' 
statistical_Analyses_all<-function(subDirCurr, 
                                   
                                   prefixCurr,
                                   run.TypeCurr,
                                   
                                   compNotesDF=NULL,
                                   studySettingz,  
                                   
                                   topMetaboFDR_,
                                   normMethod_,
                                   transMethod_,
                                   scalingMethod_,
                                   
                                   KeggPathwayAnalysisCheck_,
                                   KeggPathwayEnrichMethod_,
                                   KeggPathwayOrganism_,
                                   KeggPathwayFDR_,
                                   
                                   lbl_
){     
        
        
        ## Testing lines
        test=FALSE
        if (test==TRUE){
                
                #MIDA & MOLAR
                subDirCurr=subDir; prefixCurr=prefix_; run.TypeCurr=run.Type_
                compNotesDF=compNotesDF; # Skip this if you are working wiht idxMOLAR pipeline
                studySettingz=studySettings  
                topMetaboFDR_=topMetaboFDR_; normMethod_=normMethod_; transMethod_=transMethod_
                scalingMethod_=scalingMethod_; KeggPathwayAnalysisCheck_=KeggPathwayAnalysisCheck_
                KeggPathwayEnrichMethod_=KeggPathwayEnrichMethod_; KeggPathwayOrganism_=KeggPathwayOrganism_;KeggPathwayFDR_=KeggPathwayFDR_
                lbl_="8"
                #idxMOLAR
                prefixCurr=prefix.1Multi
                compNotesDF=NULL
        }
        
        ## Updating the pairwise comparison setting to remove indicated groups
        studySettingz <-updateSettings(studySettings_=studySettingz)
        
        nGroups=length(grep("y|yes",tolower(studySettingz$orderDF$Included)))
        
        outPathDirs<-list()
        counterD=1
        
        if (nGroups>2){
                
                ## Prepare the output directory for  multiple groups statistical analyses
                dirout.multiG=paste0(subDirCurr,"/Stats_AllGroups")
                if (!file.exists(dirout.multiG)) {dir.create(dirout.multiG)}
                setwd(dirout.multiG); outPathDirs[["Stats_AllGroups"]]<-dirout.multiG
                
                ## Prepare the input files for  multiple groups statistical analyses
                outL.AllGroups<-prep.MetaboAllG_Stats_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr, studySettingz_=studySettingz)
                
                
                ## Performing statistics for  multiple groups
                message("\n    * ",lbl_,".1 Performing statistical analyses for all groups -: ",prefixCurr)
                
                run.MetaboAllG_Stats_all( files_=outL.AllGroups$files, 
                                              nams_=outL.AllGroups$nams,
                                              nams2_=outL.AllGroups$nams2, 
                                              prefix__=prefixCurr, 
                                              nameData_="AllGroups" , 
                                              topMetaboFDR__=topMetaboFDR_
                )
                message("           Statistical analyses for all groups completed")
                
                
                ## Performing pairwise group comparison
                M=studySettingz[["groupsDF"]]
                
                Gs=list()
                counterG=1
                for (y in  1:dim(M)[1]){
                        for (z in 1:dim(M)[2]){
                                if (!is.na(M[y,z])){
                                        if (M[y,z]==TRUE){
                                                Gs[[counterG]]= c(rownames(M)[y], colnames(M)[z])
                                                counterG=counterG+1
                                        }
                                }       
                        } 
                }
                
                message("\n    * ",lbl_,".2 Pairwise group comparison ")
                
                if (length(Gs)>0) {
                        
                        for (i in 1:length(Gs)){
                                
                                ## Prepare the output directory  for two groups statistical analyses
                                
                                g1=Gs[[i]][1];  g2=Gs[[i]][2]
                                
                                oD=studySettingz[["orderDF"]]
                                o1=oD[g1,"Order"]; o2=oD[g2,"Order"]
                                
                                if(o1>o2){ g1=Gs[[i]][2];  g2=Gs[[i]][1]}   ## Performing the ordering for titles and directory
                                
                                name2G=paste0("Stats_2Groups_",g1,"_v_", g2)
                                message("\n       ",counterD," - Pairwise comparison      :   ",name2G)
                                
                                dirout.2G=paste0(subDirCurr,"/",name2G)
                                if (!file.exists(dirout.2G)) {dir.create(dirout.2G)}
                                setwd(dirout.2G); outPathDirs[[name2G]]<-dirout.2G
                                
                                ## Prepare the input files for  two groups statistical analyses
                                outL.2G <-prep.Metabo2G_Stats_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr, studySettingz_=studySettingz,files_PoolLabel=outL.AllGroups$files, g1=g1 , g2=g2, dirout.2G=dirout.2G)
                                
                                ## Run two group statistical analysis
                                fn_2=run.Metabo2G_Stats_all( files_=outL.2G$files,
                                                            nams_=outL.2G$nams,
                                                            nams2_=outL.2G$nams2,
                                                            nameData_=name2G, 
                                                            prefix__=prefixCurr, 
                                                            
                                                            topMetaboFDR__=topMetaboFDR_,  
                                                            normMethod__=normMethod_,
                                                            transMethod__=transMethod_,
                                                            scalingMethod__=scalingMethod_
                                )
                                
                                
                                if (KeggPathwayAnalysisCheck_==TRUE & !is.null(compNotesDF)){
                                        
                                        ## Preparing input files for Pathway Analysis
                                        outL.2G.Pathway <-prep.Metabo2G_Pathway_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr, compNotesDF_=compNotesDF, studySettingz_=studySettingz,files_PoolLabel=outL.AllGroups$files, g1=g1 , g2=g2,dirout.2G=dirout.2G)
                                        
                                        ## Run the KEGG PATHWAY analysis
                                        run.Metabo2G_Pathway_all( files_kegg_=outL.2G.Pathway$files_kegg,
                                                                 nams_=outL.2G.Pathway$nams,
                                                                 nams2_=outL.2G.Pathway$nams2,
                                                                 nameData_=name2G, 
                                                                 prefix__=prefixCurr,
                                                                 
                                                                 topMetaboFDR__=topMetaboFDR_,
                                                                 normMethod__=normMethod_,
                                                                 transMethod__=transMethod_,
                                                                 scalingMethod__=scalingMethod_,
                                                                 
                                                                 KeggPathwayAnalysisCheck_=KeggPathwayAnalysisCheck_,
                                                                 KeggPathwayEnrichMethod_=KeggPathwayEnrichMethod_,
                                                                 KeggPathwayOrganism_=KeggPathwayOrganism_,
                                                                 KeggPathwayFDR_=KeggPathwayFDR_,
                                                                 fn_2=fn_2
                                        )
                                        message("           Statistical/Pathway analyses completed ",name2G)
                                }
                                
                                
                                counterD=counterD+1
                        }
                        
                }else{
                        message("          No groups were selected for pairwise comparison\n")
                }   
                
                
        }else if (nGroups==2){
                
                ## Preparing files
                oD=studySettingz[["orderDF"]]
                
                g1=rownames(oD)[1];  g2=rownames(oD)[2]
                o1=oD[g1,"Order"]; o2=oD[g2,"Order"]
                
                if(o1>o2){ g1=rownames(oD)[2];  g2=rownames(oD)[1]}   ## performing the ordering for titles and directory
                
                name2G=paste0("Stats_2Groups_",g1,"_v_", g2)
                message("\n       ",counterD," - Pairwise comparison      :   ",name2G)
                
                dirout.2G=paste0(subDirCurr,name2G)
                if (!file.exists(dirout.2G)) {dir.create(dirout.2G)}
                setwd(dirout.2G); outPathDirs[[name2G]]<-dirout.2G
                
                ## Prepare the input files ()
                
                outL.AllGroups<-prep.MetaboAllG_Stats_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr, studySettingz_=studySettingz)
                
                outL.2G <-prep.Metabo2G_Stats_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr,  studySettingz_=studySettingz,files_PoolLabel=outL.AllGroups$files, g1=g1 , g2=g2, dirout.2G=dirout.2G)
                
                
                message("\n    * ",lbl_,".1 Performing statistical analyses for the two groups present -: ",prefixCurr)
                
                fn_2=run.Metabo2G_Stats_all( files_=outL.2G$files,
                                            nams_=outL.2G$nams, 
                                            nams2_=outL.2G$nams2,
                                            nameData_=name2G, 
                                            prefix__=prefixCurr, 
                                            
                                            topMetaboFDR__=topMetaboFDR_,
                                            normMethod__=normMethod_,
                                            transMethod__=transMethod_,
                                            scalingMethod__=scalingMethod_
                )  
                
                message("           Statistical analyses for the two groups present completed")
                
                if (KeggPathwayAnalysisCheck_==TRUE & !is.null(compNotesDF)){
                        
                        
                        ## Preparing input files for Pathway Analysis
                        outL.2G.Pathway <-prep.Metabo2G_Pathway_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr, compNotesDF_=compNotesDF, studySettingz_=studySettingz,files_PoolLabel=outL.AllGroups$files, g1=g1 , g2=g2,dirout.2G=dirout.2G)
                        
                        ## Run the KEGG PATHWAY analysis
                        run.Metabo2G_Pathway_all( files_kegg_=outL.2G.Pathway$files_kegg,
                                                 nams_=outL.2G.Pathway$nams,
                                                 nams2_=outL.2G.Pathway$nams2,
                                                 nameData_=name2G, 
                                                 prefix__=prefixCurr,
                                                 
                                                 topMetaboFDR__=topMetaboFDR_,
                                                 normMethod__=normMethod_,
                                                 transMethod__=transMethod_,
                                                 scalingMethod__=scalingMethod_,
                                                 
                                                 KeggPathwayAnalysisCheck_=KeggPathwayAnalysisCheck_,
                                                 KeggPathwayEnrichMethod_=KeggPathwayEnrichMethod_,
                                                 KeggPathwayOrganism_=KeggPathwayOrganism_,
                                                 KeggPathwayFDR_=KeggPathwayFDR_,
                                                 fn_2=fn_2
                        )
                        
                }
                
                message("           Statistical/Pathway analyses completed ",name2G)
                
                counterD=counterD+1
                
        }
        
        
        
        ## Study factor analysis (Separation)
        M=studySettingz[["factorsDF"]]
        
        message("\n    * ",lbl_,".3 Separation samples by specific study factors")
        sFsIx=which(M["Separate",]==TRUE) # Study factor indices
        
        if ( dim(M)[2]==1){
                message("           No study factors available for this analysis \n")
                
        } else if (length(sFsIx)>0 & nGroups>2) {
                
                studyFactors=unlist(lapply(colnames(M)[sFsIx],function(nm){strsplit(nm,"\\|")[[1]]}))
                message("           Statistical analyses for specific study factors:\n           ",paste0(studyFactors,collapse=", "))
                
                for (sF in studyFactors){
                        
                        name1F=paste0("Stats_sFactor_", sF)
                        message("\n       ",counterD," - Study factor separation : ",name1F)
                        
                        #Check if with this study factor will contain atleast two group for this analysis to proceed
                        Groups_=gsub("_\\d{1,3}$","",studySettingz$samplesDF$Sample_names)
                        Groups_sF=Groups_[grep(sF,Groups_)]
                        if (length(unique(Groups_sF))==1){ message("          - Skipping the study factor ",name1F, " because only one group is found obtained"); next}
                        
                        
                        dirout.1F=paste0(subDirCurr,"/",name1F)
                        
                        if (!file.exists(dirout.1F)) {dir.create(dirout.1F)}
                        setwd(dirout.1F); outPathDirs[[name1F]]<-dirout.1F
                        
                     
                        ## Preparing files for factor analysis
                        files<-nams<-nams2<-c()
                        
                        
                        outL.sF.sep <-prep.Metabo1F_Stats_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr,  studySettingz_=studySettingz,files_PoolLabel=outL.AllGroups$files,f1=sF, dirout.sF.sep=dirout.1F)
                        
                        
                        ## Prepare the input files for  multiple groups statistical analyses
                        if (outL.sF.sep$groupN[[1]] > 2){
                                ## Performing statistics for  multiple groups
                                message("\n    * ",lbl_,".1 Performing statistical analyses study factors -: ",prefixCurr)
                                
                                run.MetaboAllG_Stats_all( files_=outL.sF.sep$files,
                                                              nams_=outL.sF.sep$nams,
                                                              nams2_=outL.sF.sep$nams2, 
                                                              prefix__=prefixCurr, 
                                                              nameData_=paste0("Study Factor : ",sF), 
                                                              topMetaboFDR__=topMetaboFDR_
                                )
                                
                                counterD=counterD+1
                                message("           Statistical analyses completed for study factor : ",sF)  
                                
                        }else if (outL.sF.sep$groupN[[1]]==2){
                                ## Performing statistics for two groups
                                fn_2=run.Metabo2G_Stats_all( files_=outL.sF.sep$files,
                                                            nams_=outL.sF.sep$nams, 
                                                            nams2_=outL.sF.sep$nams2,
                                                            nameData_=paste0("Study Factor : ", sF), 
                                                            prefix__=prefixCurr, 
                                                            
                                                            topMetaboFDR__=topMetaboFDR_,
                                                            normMethod__=normMethod_,
                                                            transMethod__=transMethod_,
                                                            scalingMethod__=scalingMethod_
                                )  
                                
                                counterD=counterD+1
                                message("           Statistical analyses completed for study factor : ",sF)
                                
                        }else{
                                
                                message("           No enough groups for study factor : ", sF)
                        } 
                        
                }     
                
                message("           Statistical analyses for specific study factors completed.")       
                
                
                        
        }else{
                message("           No study factors were selected for separation \n")
        } 
        
        
        ## Study factor analysis (Comparison)
        message("\n    * ",lbl_,".4 Comparison of specific study factors ")
        
        M=studySettingz[["factorsDF"]]
        sFsIx2=which(M["Compare",]==TRUE) # Study factor indeces
        
        
        if (dim(M)[2]==1){
                message("           No study factors available for this analysis \n")
                
        } else if (length(sFsIx2)>0) {
                
                
                studyFactors=lapply(colnames(M)[sFsIx2],function(nm){strsplit(nm,"\\|")[[1]]})
                message("           Statistical analyses for comparison of two study factors : \n           ")
                
                
                for (si in 1:length(studyFactors) ){
                        
                        sFactors=studyFactors[[si]]
                        
                        for (s1 in 1:(length(sFactors)-1)){
                                for (s2 in (s1+1):length(sFactors)){
                                
                                        sF1=sFactors[s1]; sF2=sFactors[s2]                        
                                        
                                        name2F=paste0("Stats_sFactor_", sF1, "_vs_",sF2 )
                                        message("\n       ",counterD," - Study factor comparison : ",name2F)
                                        
                                        # Check whether the obtained groups do not overlap
                                        Groups_=gsub("_\\d{1,3}$","",studySettingz$samplesDF$Sample_names)
                                        Group_1=grep(sF1,Groups_)
                                        Group_2=grep(sF2,Groups_)
                                        if (length(length(intersect(Group_1,Group_2)))>0){ message("          - Skipping the study factor ",name2F, " because sample names overlap between two obtained groups"); next}
                                        

                                        dirout.2F=paste0(subDirCurr,"/",name2F)
                                        
                                        if (!file.exists(dirout.2F)) {dir.create(dirout.2F)}
                                        setwd(dirout.2F); outPathDirs[[name2F]]<-dirout.2F
                                        
                                        
                                        ## Preparing files for factor analysis
                                        files<-nams<-nams2<-c()
                                        
                                        outL.AllGroups<-prep.MetaboAllG_Stats_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr, studySettingz_=studySettingz)
                                        outL.sF.comp <-prep.Metabo2F_Stats_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr,  studySettingz_=studySettingz,files_PoolLabel=outL.AllGroups$files, f1=sF1 , f2=sF2 , dirout.2F=dirout.2F)
                                        
                                        
                                        if (outL.sF.comp$groupN==2){
                                        
                                                ## Performing statistics for two groups
                                                fn_2=run.Metabo2G_Stats_all( files_=outL.sF.comp$files,
                                                                            nams_=outL.sF.comp$nams, 
                                                                            nams2_=outL.sF.comp$nams2,
                                                                            nameData_=paste0("Study Factor : ", sF1," vs ", sF2), 
                                                                            prefix__=prefixCurr, 
                                                                            
                                                                            topMetaboFDR__=topMetaboFDR_,
                                                                            normMethod__=normMethod_,
                                                                            transMethod__=transMethod_,
                                                                            scalingMethod__=scalingMethod_
                                                )  
                                                
        
                                                if (KeggPathwayAnalysisCheck_==TRUE & !is.null(compNotesDF)){
                                                        
                                                        
                                                        ## Preparing input files for Pathway Analysis
                                                        outL.2F.Pathway <-prep.Metabo2F_Pathway_all(run.TypeCurr_=run.TypeCurr, subDirCurr_=subDirCurr, prefixCurr_=prefixCurr,  compNotesDF_=compNotesDF, studySettingz_=studySettingz,files_PoolLabel=outL.AllGroups$files, f1=sF1 , f2=sF2,dirout.2F=dirout.2F)
                                                        
                                                        ## Run the KEGG PATHWAY analysis
                                                        run.Metabo2G_Pathway_all( files_kegg_=outL.2F.Pathway$files_kegg,
                                                                                 nams_=outL.2F.Pathway$nams,
                                                                                 nams2_=outL.2F.Pathway$nams2,
                                                                                 nameData_=paste0("Study Factor : ", sF1," vs ",  sF2),  
                                                                                 prefix__=prefixCurr,
                                                                                 
                                                                                 topMetaboFDR__=topMetaboFDR_,
                                                                                 normMethod__=normMethod_,
                                                                                 transMethod__=transMethod_,
                                                                                 scalingMethod__=scalingMethod_,
                                                                                 
                                                                                 KeggPathwayAnalysisCheck_=KeggPathwayAnalysisCheck_,
                                                                                 KeggPathwayEnrichMethod_=KeggPathwayEnrichMethod_,
                                                                                 KeggPathwayOrganism_=KeggPathwayOrganism_,
                                                                                 KeggPathwayFDR_=KeggPathwayFDR_,
                                                                                 fn_2=fn_2
                                                        )
                                                        
                                                }
                                                
                                                
                                                message("           Statistical analyses completed ",name2F)
                                                
                                                counterD=counterD+1
                                        }else{
                                                message("           Wrong number of groups for : ", name2F)
                                        }
                                }
                        } 
                }
                message("           Statistical analyses for study factors comparison completed.")   
        }else{
                message("           No study factors were selected for comparison\n")
        }  
        
        return(outPathDirs)
        
} 


#### Section 3: ======   Main execution function of the idxMOLAR program    ==========

#'Main execution of the idxMOLAR pipeline
#'@param inFile_ full path of the input file
#'@param prefix_  prefix of the current output file  e.g "MC00139_HILIC_A"
#'@param inPath_ the path of the input directory
#'@param topMetaboFDR__ FDR threshold to determining metabolites with significant difference between groups 
#'@param normMethod__ normalization method ("none", "quantile normalization", "normalization by reference feature" etc.),default set to "none"
#'@param transMethod__ data transformation method ("log normalization", "cubic-root normalization")
#'@param scalingMethod__ data scaling method (mean-centering, auto-scaling, pareto-scaling or range-scaling)
#'@param KeggPathwayAnalysisCheck_ determines whether or not do perform pathway analysis
#'@param KeggPathwayEnrichMethod__  Kegg pathway enrichment method (QEA or ORA)
#'@param KeggPathwayOrganism__ the organism being studied
#'@param KeggPathwayFDR__ FDR threshold for significantly up-regulated pathways between two groups
#'@param run.Type_ indicates whether it is a tracer or profiling experiment
#'@param run.Method_ run method i.e "OneMethod" or "MultiMethods"
#'
mainFunIdxMOLAR<-function(     
  
                    inFile_,
                    prefix_,
                    inPath_,
                    
                    topMetaboFDR_,
                    normMethod_,
                    transMethod_,
                    scalingMethod_,
                    
                    KeggPathwayAnalysisCheck_,
                    KeggPathwayEnrichMethod_,
                    KeggPathwayOrganism_,
                    KeggPathwayFDR_,
                    
                    run.Type_,
                    run.Method_
                    
  ){ 
        
       # *** Testing lines for debugging purposes ----
       # ****** Start executing the block below to systematically run each line/section of the code in order to debug the program.
       testing="NO"
       if (testing=="YES") {
               app.Path<<-"/varidata/research/projects/mnp/tool.dev/EMeRGe1.0"
               app.Path<<-"/Volumes/projects_secondary/mnp/tool.dev/EMeRGe1.0"
               setwd(app.Path)
               
               source("midaUtils.R")
               source("molarUtils.R")
               source("idxUtils.R")
               
               .on.public.web<<-FALSE
               api.base <<- "api.metaboanalyst.ca"
               api.base <<- "http://api.xialab.ca"
            
               methods_regex <<- "HILIC_[AB]_|LIPID_|T3pos_|PHpos_|PHneg_|PHTCA_"
               methods_regex2<<-gsub("_\\|","\\|",methods_regex); methods_regex2=gsub("_$","",methods_regex2)
               
               idx_PB_Blanks_pattern<<-"Blank|pb|qc|pos.id|neg.id|sample_"
               
               
               inFile_= "/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00076/MC00076_HILIC_A_plasma.xlsx"
               inFile_="/home/brejnev.muhire/testData/MC00060/MC00060_HILIC_A.xlsx"
               inFile_="/home/brejnev.muhire/testData/MC00063/MC00063_HILIC_B.xlsx"
               inFile_="/home/brejnev.muhire/testData/MC00074/MC00074_HILIC_A.xlsx"
               inFile_="/home/brejnev.muhire/testData/MC00080/MC00080_HILIC_B.xlsx"
               inFile_="/home/brejnev.muhire/testData/MC00080/MC00080_LIPID.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00084/MC00084_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00074_redo/MC00074_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00083/MC00083_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00053/MC00053_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00053/MC00053_T3.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00099/MC00099_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00061/MC00061_HILIC_A.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00082/MC00082_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00094/MC00094_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00098/MC00098_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00098/MC00098_HILIC_A.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00099/MC00099_HILIC_A.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00099/MC00099_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00102/MC00102_HILIC_B.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/MC00106/MC00106_T3pos.xlsx"
               inFile_="/varidata/research/projects/mnp/tool.dev/EMeRGe1.0/testData/idxMOLAR_test/MC00107/MC00107_HILIC_A.xlsx"
               inFile_="/Users/Brejnev.Muhire/My_programs/idxMOLAR/testData/MC00102/MC00102_HILIC_A.xlsx"
               inFile_="/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/MC00102/MC00102_HILIC_A.xlsx"
               inFile_="/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/MC00102/MC00102_HILIC_B.xlsx"
               inFile_="/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/idxMOLAR_test/MC00060/MC00060_HILIC_B.xlsx"
               inFile_="/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/idxMOLAR_test/MC00106/MC00106_PHneg.xlsx"
               inFile_="/Users/Brejnev.Muhire/My_programs/EMeRGe1.0/testData/idxMOLAR_test/MC00111/MC00111_HILIC_A.xlsx"
               inFile_="/Users/Brejnev.Muhire/Downloads/test/MC00118/MC00118_HILIC_A.xlsx"
               inFile_="/Users/Brejnev.Muhire/Downloads/MC00116_HILIC_B.xlsx"
               inFile_="/Users/Brejnev.Muhire/Downloads/MC00128_HILIC_A.xlsx"
               inFile_="/Users/Brejnev.Muhire/Downloads/MC00122_HILIC_A.xlsx"
               inFile_="/Users/Brejnev.Muhire/Downloads/MC00137_HILIC_A.xlsx";inputDataFile=inFile_
               inFile_="/Users/Brejnev.Muhire/Downloads/MC00134_HILIC_A.xlsx";inputDataFile=inFile_
               
               inFile_="/Volumes/projects_secondary/mnp/tool.dev/EMeRGe1.0/testData/idxMOLAR_test/MC00111/MC00111_HILIC_A.xlsx"
               inFile_="/Volumes/projects_secondary/mnp/tool.dev/EMeRGe1.0/testData/idxMOLAR_test/MC00076/MC00076_T3pos.xlsx"
               
               prefix_=strsplit(basename(inFile_),".xlsx")[[1]]; 
               inPath_=dirname(inFile_)
               topMetaboFDR_=0.001; normMethod_="None"; transMethod_="LogNorm"; scalingMethod_="MeanCenter"
               preproSignalCorr.Flg_=FALSE; signalCorr.Flg_=FALSE; qc.blank.Flg_=TRUE; NA.thr.Flg_=TRUE; imputation.Flg_=TRUE
               KeggPathwayEnrichMethod_="QEA"; KeggPathwayAnalysisCheck_=TRUE; KeggPathwayEnrichMethod_="ORA";  
               KeggPathwayOrganism_="mmu"; KeggPathwayFDR_=0.1
               
               
               run.Type_="profiling"
               run.Type_="tracer"
               
               run.Method_="MultiMethods" #
               run.Method_="OneMethod"  #
               
               
               
               ##reading file please wait       
               cD<<-read_INPUT_files_Idx(inPath_=inPath_,  prefix_=prefix_, studyType_=run.Type_, studyMethod_=run.Method_ );str(cD)
               
               
        }
        # ***** End for testing lines
       
       
        # ***** Main execution begins
        message( "\n#------------- Running idxMOLAR - IDX-Orbitrap Mass Isotopologue Distribution Analysis  ------------#")
        message( "                     *---      ", format(Sys.time(), "%a, %b %Y %d, %X"),"      ---*\n")
      
        options(stringsAsFactors = F, scipen=999)
      
      
        ## Check whether study settings are selected
        message("\n#---------------------------------  Current settings  ----------------------------------------#\n")
        message("    * Directory settings")
        
        message("      Working directory                    :     ", inPath_)
        message("      RawData file path                    :     ", inFile_)
        
        message("\n    * Study type                           :     ",run.Type_)
        message("\n    * Study method                         :     ",run.Method_)
      
        message("\n    * Statistical Analyses settings")
        message("      Top metabolites FDR                  :     ", topMetaboFDR_)
        message("      Normalization method                 :     ", normMethod_)
        message("      Data transformation method           :     ", transMethod_)
        message("      Data scaling method                  :     ", scalingMethod_)
      
        message("\n    * Pathway Analyses settings")
        message("      Organism                             :     ", KeggPathwayOrganism_)
        message("      Pathway enrichment method            :     ", KeggPathwayEnrichMethod_)
        message("      FDR cut-off                          :     ", KeggPathwayFDR_)
      
        message("\n#--------------------------------------------------------------------------------------------#\n")
      
        message(paste0("\n*-------------  idxMOLAR : Analysing ", basename(inFile_) ,"  ------------- ",date(), "\n") )
        
        
        message("\n** Step 0: Checking the input file ...")
        
        ## Checking The quality the input data (this is aready done)
        ret=check_idxMETHOD_files(cD$compData_List, cD$prefixM_List, studyType_=run.Type_, studyMethod_=run.Method_)
        if (ret=="Error"){return()}
        
        
        ## Checking that the studySettings are valid
        load(file=paste0(inPath_,"/",prefix_,"_studySettings.RData"))
        
      
        all.Compounds_=TRUE
        
        if (run.Method_=="OneMethod"){
                if (dim(cD$compData_List[[1]])[1]>800000){
                      message("        More than 800,000 were detected in the input file, only tables for selected compounds will be make (not all compound).")
                      all.Compounds_=FALSE
                }
        }else{
              if (max(sapply(cD$compData_List,function(df){dim(df)[1]}))>800000){
                     message("        More than 800,000 were detected in the input file, only tables for selected compounds will be make (not all compound).")
                     all.Compounds_=FALSE
              }
        }
        
        
        ## Creating the output directory ----
        if (run.Method_=="OneMethod"){subDir=paste0(inPath_,"/",cD$prefixM_List[[1]],"_out/")}else{subDir=paste0(inPath_,"/",strsplit(cD$prefixM_List[[1]],"_")[[1]][[1]],"_Comb_out/")}
        if (!dir.exists(subDir)){dir.create(subDir)} else {unlink(subDir, recursive = TRUE); dir.create(subDir)}
        setwd(subDir)
        
  
        message("   * Current output directory : ",subDir)
        
        ## Extracting MID_FE, PoolSize, and TotalLabel_FE data from the xlsx files ----
        for (ij in 1:length(cD$compData_List)){
                
                compData=cD$compData_List[[ij]]
                prefixTmp_=cD$prefixM_List[[ij]]
                
                sampleNames<-get_sampleNameList_idx(compData)                           ## Picking up all sample names, place them in to the list
                
                ## Remove blank and QC and other unwanted samples names
                message(  "   * Removing QC sample names, and other unwanted sample names ", idx_PB_Blanks_pattern)
                ixx=grep(pattern=idx_PB_Blanks_pattern,sampleNames,ignore.case=T)  
                if (length(ixx>0)){sampleNames=sampleNames[-ixx]} 
                
                
                sampleNames<-sampleNames[order(names(sampleNames))]                 ## Sorting by sample name
                
                outputDF<-do.call(rbind,sampleNames); colnames(outputDF)="Sample_Name"  ## Generate the output dataFrame
                outputDF=as.data.frame(outputDF)                                        ## Output dataFrame formatting 
                
                
                # I. Generating tables for all compounds 
                if (all.Compounds_==TRUE){
                      message("\n** Step 1-2: Extracting the tables for all compounds - ",prefixTmp_)
                      ExtractTBL_all_compound(compData_=compData, subDir_=subDir, sampleNames_=sampleNames, prefix__=prefixTmp_, outputDF_=outputDF, run.Type__= run.Type_)
                }
                
                
                # II. Tables for selected compounds  [[NB: avoid that the last compound subDataFrame include all rows at the bottom ]]
                message("\n** Step 1-2: Extracting the tables for selected compounds - ",prefixTmp_)
                ExtractTBL_selected_compound(compData_=compData, subDir_=subDir, sampleNames_=sampleNames, prefix__=prefixTmp_, outputDF_=outputDF,run.Type__= run.Type_)
        }
        
        
        message("\n** Step 3: Reformatitng each individual table ") 
        
        ## Reformating the MID_FE, PoolSize, and TotalLabel_FE
        #  Reformats names, remove the method and project number from the sample names
        #  Use the study settings for (sample removal, group exclusion, or study factor scaling)
        
        
        for (index1 in 1:length(cD$prefixM_List)){
                
                prefix.1=cD$prefixM_List[[index1]]
                
                if (run.Type_=="profiling"){
                  
                        r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1, suffix_="_PoolSize.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}
                        PoolSize_global<<-r
                        
                        if (all.Compounds_==TRUE){
                                r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1,  suffix_="_All_PoolSize.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}
                        }
                }
                
                if (run.Type_=="tracer"){
                        # Table1 MID FE;  Table2 PoolSize;  Table3 Total Label FE
                        r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1, suffix_="_MID_FE.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}
                        MID_FE_global<<-r
                        
                        r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1, suffix_="_PoolSize.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}
                        PoolSize_global<<-r
                        
                        r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1, suffix_="_TotalLabel_FE.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}
                        TotalLabel_FE_global<<-r
                        
                        
                        if (all.Compounds_==TRUE){
                                # Table1 MID FE;  Table2 PoolSize;  Table3 Total Label FE
                                r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1, suffix_="_All_MID_FE.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}
                                r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1, suffix_="_All_PoolSize.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}
                                r=format_Tb_all(subDir_=subDir, prefix1__=prefix.1, suffix_="_All_TotalLabel_FE.csv", studySettings_=studySettings, mod_="idxMOLAR"); if (!is.data.frame(r)){message(r); return()}    
                        }
                }
        }
 

        message("\n** Step 4: Merging Method 1 and Method 2 data ...") 
        
        ## Merging Method in case of multiple methods are used ----
        ##  Change the HILIC_A HILIC_B to Comb
        ##  Get final the tables into global variables for visualization on the interface
        
        
        if (run.Method_=="MultiMethods"){
  
                  if (run.Type_=="profiling"){
                           
                           r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_PoolSize.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                           PoolSize_global<<-r
                           
                           if (all.Compounds_==TRUE){
                                    r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_All_PoolSize.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                           }
                  }
                  
                  if (run.Type_=="tracer"){
                    
                            r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_MID_FE.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                            MID_FE_global<<-r
                            
                            r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_PoolSize.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                            PoolSize_global<<-r
                            
                            r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_TotalLabel_FE.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                            TotalLabel_FE_global<<-r
                          
                            if (all.Compounds_==TRUE){
                                      r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_All_MID_FE.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                                      r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_All_PoolSize.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                                      r=merge_idxMETHOD_Multi(subDir_=subDir, prefix.L_=cD$prefixM_List, suffix_="_All_TotalLabel_FE.csv", studySettings_=studySettings); if (!is.data.frame(r)){message(r); return()}
                            }
                  }
        }else if (run.Method_=="OneMethod"){
                message("\n           Merging skipped since only one method was used.")
        }


        ## Chose a universal variable (prefix.1Multi, for 1 or multiple methods) used to calling Statistical analysis
        if (run.Method_=="OneMethod"){
                  prefix.1Multi=cD$prefixM_List[[1]]
        }else if (run.Method_=="MultiMethods"){
                  prefix.1Multi= paste0(strsplit(cD$prefixM_List[[1]],"_")[[1]][1],"_Comb") 
        }
        
        ## Making MID_Poolsize and TotalLabelPoolSize tables ----  
        
        message("\n           Generating MID PoolSize and Total Labeling PoolSize tables.") 
        
        if (run.Type_=="tracer"){
                # Table4 MID Poolsize; Table6 Total labeling Poolsize 
                
                r=make_MID_PoolSize_tb(subDir_=subDir, currPrefix_=prefix.1Multi,  suffix_MID_="_MID_FE.csv" , suffix_PS_="_PoolSize.csv", suffix_Out_="_MID_PoolSize.csv", mod_="idxMOLAR" )
                MID_PoolSize_global<<-r
                
                r=make_TotalLabel_PoolSize_tb(subDir_=subDir, currPrefix_=prefix.1Multi, suffix_PS_="_PoolSize.csv" , suffix_TlFE_="_TotalLabel_FE.csv",suffix_Out_="_TotalLabel_PoolSize.csv");  if (!is.data.frame(r)){message(r); return()}
                TotalLabel_PoolSize_global<<-r
                
                if (all.Compounds_==TRUE){
                        # Table4 MID Poolsize; Tables Total labeling Poolsize
                        r=make_MID_PoolSize_tb(subDir_=subDir, currPrefix_=prefix.1Multi,  suffix_MID_="_All_MID_FE.csv" , suffix_PS_="_All_PoolSize.csv",suffix_Out_="_All_MID_PoolSize.csv", mod_="idxMOLAR" )
                        r=make_TotalLabel_PoolSize_tb(subDir_=subDir, currPrefix_=prefix.1Multi, suffix_PS_="_All_PoolSize.csv" , suffix_TlFE_="_All_TotalLabel_FE.csv",suffix_Out_="_All_TotalLabel_PoolSize.csv");  if (!is.data.frame(r)){message(r); return()}       
                }
        }
 
        
          
        ## Generating tables for Prism (for only  selected compound) ----
        message("\n** Step 5: Generating tables (.csv) for Prism - only for selected compounds  ...")
        
        if (run.Type_=="profiling"){
                  message("          No PRISM tables generated since this is a profiling study. Function to be impletented for metabolite profiling analysis")  
                
        }else if (run.Type_=="tracer"){
                
                   #Fractional Enrichment  
                   message("           Producing PRISM format fractional enrichment table for each compound, please wait ...")
                   p_MID=paste0(subDir, prefix.1Multi, "_MID_FE.csv")  
                   p_TPOOL=paste0(subDir, prefix.1Multi, "_PoolSize.csv")  
                   make_PRISM_Tables(p_MID_=p_MID,  p_TPOOL_=p_TPOOL, subDir_=subDir , studySettings_=studySettings)
        }
        
        
        ## Statistical analyses (with MetaboAnalyst) ----
        message("\n** Step 6: Statistical analysis (with MetaboAnalyst) ...")
        
        outPathDirs <-statistical_Analyses_all (  subDirCurr=subDir, 
                                                  
                                                  prefixCurr=prefix.1Multi,
                                                  run.TypeCurr=run.Type_,
                                                  
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

        
        if (is.null(outPathDirs)){message("      * Exiting ");return()}
        
       
        ## Generating the final reports ----

        message("\n** Step 7: Generating pdf report")  
        
        generateFinalReport(outPathDirz=outPathDirs, subDir_=subDir)
        
        setwd(app.Path)
        cat(app.Path)
        
        message("\n#-----------------------------   idxMOLAR run completed successfully   --------------------------#")
        message( "                          *---      ", format(Sys.time(), "%a, %b %Y %d, %X"),"      ---*\n")
        
}

