#This script uses an extension of vegan library's bioenv() function and finds the best set of environmental variables with maximum (rank) correlation with community dissimilarities and then plots them as vectors along with the best subset of taxas on the NMDS plot.
# ============================================================
# Tutorial on plotting significant taxa and environmental variables on an NMDS plot using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

library(vegan)
library(ggplot2)
library(grid)

# This R script is an extension of vegan library's bioenv()
# function and uses the bio.env() and bio.step() of
#	http://menugget.blogspot.co.uk/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html
#	The original author suggested these functions to overcome
#	the inflexibility of the bioenv() function which uses
#	a similarity matrix based on normalized "euclidean" distance.
# The new functions are given below and implement the following algorithms: 
# Clarke, K. R & Ainsworth, M. 1993. A method of linking multivariate community structure to environmental variables. Marine Ecology Progress Series, 92, 205-219.
# Clarke, K. R., Gorley, R. N., 2001. PRIMER v5: User Manual/Tutorial. PRIMER-E, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 2001. Changes in Marine Communities: An Approach to Statistical Analysis and Interpretation, 2nd edition. PRIMER-E Ltd, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 1998. Quantifying structural redundancy in ecological communities. Oecologia, 113:278-289. 

bv.step <- function(fix.mat, var.mat, 
                    fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    max.rho=0.95,
                    min.delta.rho=0.001,
                    random.selection=TRUE,
                    prop.selected.var=0.2,
                    num.restarts=10,
                    var.always.include=NULL,
                    var.exclude=NULL,
                    output.best=10
){
  
  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
  require(vegan)
  
  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
  
  fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
  
  #an initial removal phase
  var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
  full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
  var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
  RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
  for(i in 1:dim(var.comb)[2]){
    var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
    temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
    RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
    RES$rho[i] <- temp$estimate
  }
  delta.rho <- RES$rho - full.cor
  exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
  
  if(random.selection){
    num.restarts=num.restarts
    prop.selected.var=prop.selected.var
    prob<-rep(1,ncol(var.mat))
    if(prop.selected.var< 1){
      prob[exclude]<-0
    }
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  } else {
    num.restarts=1
    prop.selected.var=1  
    prob<-rep(1,ncol(var.mat))
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  }
  
  RES_TOT <- c()
  for(i in 1:num.restarts){
    step=1
    RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
    attr(RES$step.dir, "levels") <- c("F","B")
    best.comb <- which.max(RES$rho)
    best.rho <- RES$rho[best.comb]
    delta.rho <- Inf
    selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
    while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
      #forward step
      step.dir="F"
      step=step+1
      var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
      if(RES$n.var[best.comb] == 0){
        var.comb.incl<-1:length(var.comb)
      } else {
        var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
        temp <- NA*1:length(var.comb)
        for(j in 1:length(temp)){
          temp[j] <- all(var.keep %in% var.comb[[j]]) 
        }
        var.comb.incl <- which(temp==1)
      }
      
      RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
      for(f in 1:length(var.comb.incl)){
        var.incl <- var.comb[[var.comb.incl[f]]]
        var.incl <- var.incl[order(var.incl)]
        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES.f$var.incl[f] <- paste(var.incl, collapse=",")
        RES.f$rho[f] <- temp$estimate
      }
      
      last.F <- max(which(RES$step.dir=="F"))
      RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
      best.comb <- which.max(RES$rho)
      delta.rho <- RES$rho[best.comb] - best.rho 
      best.rho <- RES$rho[best.comb]
      
      if(best.comb == step){
        while(best.comb == step & RES$n.var[best.comb] > 1){
          #backward step
          step.dir="B"
          step <- step+1
          var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
          var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
          RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
          for(b in 1:length(var.comb)){
            var.incl <- var.comb[[b]]
            var.incl <- var.incl[order(var.incl)]
            var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES.b$var.incl[b] <- paste(var.incl, collapse=",")
            RES.b$rho[b] <- temp$estimate
          }
          RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
          best.comb <- which.max(RES$rho)
          best.rho<- RES$rho[best.comb]
        }
      } else {
        break()
      }
      
    }
    
    RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
    print(paste(round((i/num.restarts)*100,3), "% finished"))
  }
  
  RES_TOT <- unique(RES_TOT[,3:5])
  
  
  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
  } else {
    order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
  }
  rownames(order.by.best)<-NULL
  
  order.by.i.comb <- c()
  for(i in 1:length(selected.var)){
    f1 <- which(RES_TOT$n.var==i)
    f2 <- which.max(RES_TOT$rho[f1])
    order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
  }
  rownames(order.by.i.comb)<-NULL
  
  if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
  out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
  )
  out
  
}

bio.env <- function(fix.mat, var.mat, 
                    fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    output.best=10,
                    var.max=ncol(var.mat)
){
  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(var.max > dim(var.mat)[2]){stop("var.max cannot be larger than the number of variables (columns) in var.mat")}
  
  require(vegan)
  
  combn.sum <- sum(factorial(ncol(var.mat))/(factorial(1:var.max)*factorial(ncol(var.mat)-1:var.max)))
  
  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
  fix.dist <- vegdist(fix.mat, method=fix.dist.method)
  RES_TOT <- c()
  best.i.comb <- c()
  iter <- 0
  for(i in 1:var.max){
    var.comb <- combn(1:ncol(var.mat), i, simplify=FALSE)
    RES <- data.frame(var.incl=rep(NA, length(var.comb)), n.var=i, rho=0)
    for(f in 1:length(var.comb)){
      iter <- iter+1
      var.dist <- vegdist(as.matrix(var.mat[,var.comb[[f]]]), method=var.dist.method)
      temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
      RES$var.incl[f] <- paste(var.comb[[f]], collapse=",")
      RES$rho[f] <- temp$estimate
      if(iter %% 100 == 0){print(paste(round(iter/combn.sum*100, 3), "% finished"))}
    }
    
    order.rho <- order(RES$rho, decreasing=TRUE)
    best.i.comb <- c(best.i.comb, RES$var.incl[order.rho[1]])
    if(length(order.rho) > output.best){
      RES_TOT <- rbind(RES_TOT, RES[order.rho[1:output.best],])
    } else {
      RES_TOT <- rbind(RES_TOT, RES)
    }
  }
  rownames(RES_TOT)<-NULL
  
  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- order(RES_TOT$rho, decreasing=TRUE)[1:output.best]
  } else {
    order.by.best <- order(RES_TOT$rho, decreasing=TRUE)
  }
  OBB <- RES_TOT[order.by.best,]
  rownames(OBB) <- NULL
  
  order.by.i.comb <- match(best.i.comb, RES_TOT$var.incl)
  OBC <- RES_TOT[order.by.i.comb,]
  rownames(OBC) <- NULL
  
  out <- list(
    order.by.best=OBB,
    order.by.i.comb=OBC,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(OBB$var.incl[1], ",")))], collapse=",") ,
    best.model.rho=OBB$rho[1]
  )
  out
}

abund_table<-read.csv("SPE_pitlatrine.csv",row.names=1,check.names=FALSE)
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
meta_table<-read.csv("ENV_pitlatrine.csv",row.names=1,check.names=FALSE)
#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
meta_table<-meta_table[rownames(abund_table),]
#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
# > head(grouping_info)
# X1 X2 X3
# T_2_1   T  2  1
# T_2_10  T  2 10
# T_2_12  T  2 12
# T_2_2   T  2  2
# T_2_3   T  2  3
# T_2_6   T  2  6


#Parameters
cmethod<-"pearson" #Correlation method to use: pearson, pearman, kendall
fmethod<-"bray" #Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
vmethod<-"bray" #Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
nmethod<-"bray" #NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao


res <- bio.env(wisconsin(abund_table), meta_table,fix.dist.method=fmethod, var.dist.method=vmethod, correlation.method=cmethod,
               scale.fix=FALSE, scale.var=TRUE) 

#Get the 10 best subset of environmental variables
envNames<-colnames(meta_table)
bestEnvFit<-""
for(i in (1:length(res$order.by.best$var.incl)))
{
  bestEnvFit[i]<-paste(paste(envNames[as.numeric(unlist(strsplit(res$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res$order.by.best$rho[i],sep="")
}
bestEnvFit<-data.frame(bestEnvFit)
colnames(bestEnvFit)<-"Best combination of environmental variables with similarity score"

#> bestEnvFit
#Best combination of environmental variables with similarity score
#1                                               Carbo = 0.112933470129642
#2                                  Temp + VS + Carbo = 0.0850146558049758
#3                                       Temp + Carbo = 0.0836568099887996
#4                                               Temp = 0.0772024618792259
#5                                          perCODsbyt = 0.076360226020659
#6              TS + CODs + perCODsbyt + Prot + Carbo = 0.0752573390422113
#7              Temp + CODt + CODs + perCODsbyt + NH4 = 0.0749607646516179
#8  Temp + TS + VFA + perCODsbyt + NH4 + Prot + Carbo = 0.0718412745558854
#9                             pH + Temp + TS + Carbo = 0.0693960850973541
#10                                Temp + NH4 + Carbo = 0.0687961955200212

res.bv.step.biobio <- bv.step(wisconsin(abund_table), wisconsin(abund_table), 
                              fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                              scale.fix=FALSE, scale.var=FALSE, 
                              max.rho=0.95, min.delta.rho=0.001,
                              random.selection=TRUE,
                              prop.selected.var=0.3,
                              num.restarts=10,
                              output.best=10,
                              var.always.include=NULL) 

#Get the 10 best subset of taxa
taxaNames<-colnames(abund_table)
bestTaxaFit<-""
for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl)))
{
  bestTaxaFit[i]<-paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res.bv.step.biobio$order.by.best$rho[i],sep="")
}
bestTaxaFit<-data.frame(bestTaxaFit)
colnames(bestTaxaFit)<-"Best combination of taxa with similarity score"

#> bestTaxaFit
#Best combination of taxa with similarity score
#1                             Bacilli + Bacteroidia + Chrysiogenetes + Clostridia + Dehalococcoidetes + Fibrobacteria + Flavobacteria + Fusobacteria + Gammaproteobacteria + Methanomicrobia + Mollicutes + Opitutae + Synergistia + Thermomicrobia + Unknown = 0.900998476657462
#2                                              Bacilli + Bacteroidia + Clostridia + Dehalococcoidetes + Fibrobacteria + Flavobacteria + Fusobacteria + Gammaproteobacteria + Methanomicrobia + Mollicutes + Opitutae + Synergistia + Thermomicrobia + Unknown = 0.899912718316228
#3                                                        Bacilli + Bacteroidia + Clostridia + Dehalococcoidetes + Fibrobacteria + Flavobacteria + Fusobacteria + Gammaproteobacteria + Methanomicrobia + Mollicutes + Opitutae + Synergistia + Thermomicrobia = 0.896821772937576
#4  Clostridia + Dehalococcoidetes + Deltaproteobacteria + Epsilonproteobacteria + Erysipelotrichi + Flavobacteria + Fusobacteria + Methanobacteria + Mollicutes + Opitutae + Planctomycetacia + Sphingobacteria + Subdivision3 + Synergistia + Thermomicrobia = 0.892670058822226
#5                                                                     Bacilli + Bacteroidia + Clostridia + Dehalococcoidetes + Fibrobacteria + Flavobacteria + Fusobacteria + Gammaproteobacteria + Methanomicrobia + Opitutae + Synergistia + Thermomicrobia = 0.892533063335985
#6                   Clostridia + Dehalococcoidetes + Deltaproteobacteria + Epsilonproteobacteria + Erysipelotrichi + Flavobacteria + Fusobacteria + Methanobacteria + Mollicutes + Opitutae + Planctomycetacia + Sphingobacteria + Subdivision3 + Synergistia = 0.891217789278463
#7                                     Clostridia + Dehalococcoidetes + Deltaproteobacteria + Epsilonproteobacteria + Erysipelotrichi + Flavobacteria + Fusobacteria + Mollicutes + Opitutae + Planctomycetacia + Sphingobacteria + Subdivision3 + Synergistia = 0.888669881927483
#8                                                                                       Bacilli + Bacteroidia + Clostridia + Dehalococcoidetes + Fibrobacteria + Flavobacteria + Fusobacteria + Gammaproteobacteria + Opitutae + Synergistia + Thermomicrobia = 0.887052815492516
#9                                                  Clostridia + Dehalococcoidetes + Deltaproteobacteria + Epsilonproteobacteria + Erysipelotrichi + Flavobacteria + Fusobacteria + Opitutae + Planctomycetacia + Sphingobacteria + Subdivision3 + Synergistia = 0.885880090785632
#10                             Acidobacteria_Gp4 + Actinobacteria + Bacilli + Betaproteobacteria + Clostridia + Dehalococcoidetes + Erysipelotrichi + Fibrobacteria + Lentisphaeria + Methanomicrobia + Opitutae + Sphingobacteria + Thermomicrobia + Unknown = 0.882956559638505

#Generate NMDS plot
MDS_res=metaMDS(abund_table, distance = nmethod, k = 2, trymax = 50)

bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[1], ",")))
bio.fit <- envfit(MDS_res, abund_table[,bio.keep,drop=F], perm = 999)

#> bio.fit
#
#***VECTORS
#
#NMDS1    NMDS2     r2 Pr(>r)    
#Bacilli              0.94632  0.32322 0.0423  0.167    
#Bacteroidia         -0.24478 -0.96958 0.0383  0.218    
#Chrysiogenetes      -0.22674  0.97395 0.0078  0.660    
#Clostridia          -0.98067 -0.19567 0.1092  0.017 *  
#Dehalococcoidetes   -0.79160  0.61103 0.1421  0.009 ** 
#Fibrobacteria       -0.90885 -0.41712 0.1320  0.010 ** 
#Flavobacteria        0.56629  0.82421 0.1746  0.003 ** 
#Fusobacteria         0.26601 -0.96397 0.0284  0.333    
#Gammaproteobacteria  0.67008  0.74229 0.2024  0.001 ***
#Methanomicrobia     -0.99912 -0.04188 0.0602  0.120    
#Mollicutes           0.99885 -0.04798 0.0052  0.769    
#Opitutae             0.55774  0.83002 0.2033  0.002 ** 
#Synergistia         -0.99732  0.07322 0.2079  0.002 ** 
#Thermomicrobia       0.38169  0.92429 0.2514  0.001 ***
#Unknown             -0.99232 -0.12373 0.1570  0.009 ** 
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Permutation: free
#Number of permutations: 999

#use the best set of environmental variables in env.keep
eval(parse(text=paste("env.keep <- c(",res$order.by.best$var.incl[1],")",sep="")))
env.fit <- envfit(MDS_res, meta_table[,env.keep,drop=F], perm = 999) 

#Get site information
df<-scores(MDS_res,display=c("sites"))

#Add grouping information
df<-data.frame(df,Type=grouping_info[rownames(df),1])

#Get the vectors for bioenv.fit
df_biofit<-scores(bio.fit,display=c("vectors"))
df_biofit<-df_biofit*vegan:::ordiArrowMul(df_biofit)
df_biofit<-as.data.frame(df_biofit)

#Get the vectors for env.fit
df_envfit<-scores(env.fit,display=c("vectors"))
df_envfit<-df_envfit*vegan:::ordiArrowMul(df_envfit)
df_envfit<-as.data.frame(df_envfit)

#Draw samples
p<-ggplot()
p<-p+geom_point(data=df,aes(NMDS1,NMDS2,colour=Type))
#Draw taxas
p<-p+geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_biofit*1.1),aes(NMDS1, NMDS2, label = rownames(df_biofit)),color="#808080",alpha=0.5)
#Draw environmental variables
p<-p+geom_segment(data=df_envfit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                  arrow = arrow(length = unit(0.2, "cm")),color="#4C005C",alpha=0.5)

p<-p+geom_text(data=as.data.frame(df_envfit*1.1),aes(NMDS1, NMDS2, label = rownames(df_envfit)),color="#4C005C",alpha=0.5)
p<-p+theme_bw()
pdf("NMDS_bioenv.pdf")
print(p)
dev.off()