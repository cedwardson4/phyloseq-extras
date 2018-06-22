#R code for calculating median absolute devance, variance, and standard deviation around pairwise community distance/dissimilaritys.
#written by Ashley Shade, shade.ashley@gmail.com for calculations in:
#"A meta-analysis of changes in bacterial and archaeal communities with time"
#by A Shade, JG Caporaso, J Handelsman, R Knight, N Fierer
#The ISME Journal 7 (8), 1493-1506

#This code comes with no warranty

#This code has a Creative Commons Universial license.  Use freely, adapt, share, but attribute and note changes.

#Input files
#dist_fp is a tab-delimited pairwise distance or similarity matrix for community structure.  
#Default is column names in row 1 and column 1.
#map_fp is a tab-delimited mapping file with samples in the same order as the distance file. 
#One column name should be labeled "SiteID" to identify the group membership (categories) of the sample names 
#in the distance matrix. Variation, sd, and MAD will be calculated among samples within these designated groups. 

#Function
varDistance.f=function(dist_fp, map_fp){  
  map=read.table(map_fp, header=TRUE, sep="\t", check.names=FALSE)
  dist=read.table(dist_fp, header=TRUE, row.names=1, sep="\t",check.names=FALSE)
  
  u=unique(map[,"SiteID"])
  
  v.out=NULL
  for(i in 1:length(u)){
    dist.tmp=dist[map[,"SiteID"]==u[i],map[,"SiteID"]==u[i]]
    
    v=var(as.dist(dist.tmp))
    stdev=sd(as.dist(dist.tmp))
    m=mad(as.dist(dist.tmp))
    
    v2=c(v,stdev,m)
    
    v.out=rbind(v.out,v2)
    
  }
  row.names(v.out)=u
  colnames(v.out)=c("Var", "SD", "MAD")	
  write.table(v.out, "VarianceTable.txt", sep="\t", quote=FALSE)
  
}

#To run:
varDistance.f(dist_fp, map_fp)