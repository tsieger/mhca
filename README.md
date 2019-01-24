# mhca

Mahalanobis distance-based hierarchical cluster analysis, in which
elliptical clusters get found naturally. 

## Install:

    devtools::install_github("tsieger/mhca")

## Example (Comparison with Classical HCA):

    library(mhca)
    opar<-par(mfrow=c(2,2))
    
    k<-3
    n<-nrow(xy)
    
    # classical HCA
    h<-hclust(dist(xy))
    
    # Mahalanobis HCA
    mh<-mhclust(xy,thresh=.3)
    
    ch<-cutree(h,k=k)
    cmh<-cutree(mh,k=k)
    
    # feature space plots with 3 top clusters
    plot(xy[,1],xy[,2],asp=1,col=ch,main='HCA',frame=FALSE)
    plot(xy[,1],xy[,2],asp=1,col=cmh,main='Mahalanobis HCA',frame=FALSE)
    
    # HCA dendrogram
    plot(h,hang=0,labels=FALSE,main='Dendrogram of HCA')
    y<-min(h$height)-diff(range(h$height))/20
    text(1:n,y,(1:n)[h$order],col=ch[h$order],srt=90)
    
    # MHCA dendrogram
    plot(mh,labels=FALSE,main='Dendrogram of MHCA')
    y<-min(mh$height)-diff(range(mh$height))/10
    text(1:n,y,(1:n)[mh$order],col=cmh[mh$order],srt=90)
    
    par(opar)

![Example](/inst/mhca.png?raw=true "Comparison with classical HCA.")

Find out more at https://github.com/tsieger/mhca.

## Reference

A paper ["Detection and monitoring of normal and leukemic cell 
populations with hierarchical clustering of flow cytometry 
data"](http://dx.doi.org/10.1002/cyto.a.21148) in 
[Cytometry Part A](https://onlinelibrary.wiley.com/journal/15524930).
