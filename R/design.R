newGT<-read.csv('workspace/ikmc_targeting_efficiencies/tables/newGT.csv',sep='\t')
newGT[1:10,]
geenClub=read.csv('workspace/ikmc_targeting_efficiencies/tables/geneClubs.csv',sep='\t')
geneClub=geenClub
rm(geenClub)
rap=function(x){return(which(geneClubGT$EP_PLATE_NAME==x,arr.ind=TRUE))}
#uniqueGTEPplateIndex=lapply(unique(geneClubGT$EP_PLATE_NAME),rap)
rap=function(x){return(which(geneClub$EP_PLATE_NAME==x,arr.ind=TRUE))}
geneClubEPIndex=lapply(unique(geneClub$EP_PLATE_NAME),rap)
geneClub=geneClub[order(geneClub$EP_PLATE_NAME),]
geneClub[1:10,]
rap=function(x){return(which(geneClub$EP_PLATE_NAME==x,arr.ind=TRUE))}
geneClubEPIndex=lapply(unique(geneClub$EP_PLATE_NAME),rap)
geneClubEPIndex
effCalc=function(x){return((x$EPD_DISTRIBUTE+x$TARGETED_TRAP)/x$EPD_WELL_NAME)}
eff=effCalc(geneClub)
eff
eff[which(eff=="NaN")]=0
eff
data.frame(geneClub,efficiencies=eff)
geneClub=data.frame(geneClub,efficiencies=eff)
epEff=function(x){return(sum(geneClubGT$efficiencies[x])/length(x))}
epEffs=unlist(lapply(geneClubEPIndex,epEff))
#epEffs=unlist(lapply(geneClubEPIndex,epEff))
epEff=function(x){return(sum(geneClub$efficiencies[x])/length(x))}
epEffs=unlist(lapply(geneClubEPIndex,epEff))
epEffs
goodEP=unlist(epEffs[which(unlist(epEffs)>0.1)])
rubbishPlates=which(unlist(epEffs)<0.1) 
rubbishPlates
length(rubbishPlates)
rubbishGenes=function(x){return(geneClubEPIndex[x])}
rubGeneseff=geneClub$efficiencies[unlist(rubbishGenes(rubbishPlates))]
rubGeneseff
rubbishPlates
length(rubGeneseff)
unlist(rubbishGenes(rubbishPlates))
rubGenes=geneClub$MARKER_SYMBOL[unlist(rubbishGenes(rubbishPlates))]
identifyGenes=function(x){return(which(newData$MARKER_SYMBOL==x))}
#rubIndex=unlist(lapply(rubGenes,identifyGenes))
identifyGenes=function(x){return(which(geneClub$MARKER_SYMBOL==x))}
rubIndex=unlist(lapply(rubGenes,identifyGenes))
rubIndex
rubGenes
goodGenes=newData[-rubIndex,]
goodGenesEff=effCalc(goodGenes)
#goodGenes=data.frame(goodGenes,efficiencies=goodGenesEff)
goodGenes=geneClub[-rubIndex,]
goodGenesEff=effCalc(goodGenes)
goodGenes=data.frame(goodGenes,efficiencies=goodGenesEff)
goodGenesEff[1:10]
length(goodGenesEff)
unique(goodGenes$MARKER_SYMBOL)
length(unique(goodGenes$MARKER_SYMBOL))
length(geneClub$MARKER_SYMBOL)
length(unique(geneClub$MARKER_SYMBOL))
13196-9521
hist(goodGenesEff,col="lightblue",border="gray")
length(goodGenes$MARKER_SYMBOL)
goodGenes[1:10,]
goodGenes[order(goodGenes$MARKER_SYMBOL)[1:10],]
goodGenes[1@10,]
goodGenes[1:10,]
goodGenes$efficiencies=='NaN'
which(goodGenes$efficiencies=='NaN')
goodGenes[1:10,-8]
goodGenes=goodGenes[,-8]
plot(goodGenes$EP_PLATE_NAME)
plot(goodGenes$EP_PLATE_NAME,xlim=c(1,10))
plot(goodGenes$EP_PLATE_NAME,xlim=c(1,1))
hist(goodGenes$EP_PLATE_NAME,xlim=c(1,1))
plot(goodGenes$EP_PLATE_NAME,xlim=c(1,1))
plot(goodGenes$EP_PLATE_NAME)
plot(goodGenes$efficiencies)
hist(goodGenes$efficiencies)
hist(goodGenes$efficiencies,col="lightblue")
hist(goodGenes$efficiencies,col="lightblue",border="Gray")
hist(goodGenes$efficiencies,col="lightblue",border="Gray",ylim=c(0,3500))
hist(goodGenes$efficiencies,col="lightblue",border="Gray",ylim=c(0,3500),xlab="Targeting Efficiencies")
hist(goodGenes$efficiencies,col="lightblue",border="Gray",ylim=c(0,3500),xlab="Targeting Efficiencies")
history()
hist(goodGenes$efficiencies,col="lightblue",border="Gray",ylim=c(0,3500),xlab="Targeting Efficiencies")
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/edgeGeneFre')
hist(goodGenes$efficiencies,col="lightblue",border="Gray",ylim=c(0,3500),xlab="Targeting Efficiencies")
dev.off()
plot(goodGenes$EP_PLATE_NAME)
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/edgeEPPLATES')
plot(goodGenes$EP_PLATE_NAME)
dev.off()
length(which(goodGenes$MARKER_SYMBOL %in% newGT$genes))
length(which(unique(goodGenes$MARKER_SYMBOL) %in% newGT$genes))
length(which(unique(goodGenes$MARKER_SYMBOL) %in% unique(newGT$genes)))
length(which(goodGenes$MARKER_SYMBOL %in% newGT$genes))
which(goodGenes$MARKER_SYMBOL %in% newGT$genes)[1:2]
goodGenes$MARKER_SYMBOL[4]
newGT$gene[1]
newGT$gene2[1]
newGT$genes[1]
newGT$genes[4]
newGT$genes[1:10]
length(goodGenes$MARKER_SYMBOL)
goodGenes$MARKER_SYMBOL[1:10]
goodInter=goodGenes[which(goodGenes$MARKER_SYMBOL %in% newGT$genes,)
]
goodInter=goodGenes[which(goodGenes$MARKER_SYMBOL %in% newGT$genes),]
goodInter[1:10,]
length(which(newGT$genes %in% goodGenes$MARKER_SYMBOL))
newInter=newGT[(which(newGT$genes %in% goodGenes$MARKER_SYMBOL),]
newInter=newGT[(which(newGT$genes %in% goodGenes$MARKER_SYMBOL),]
newGT[1:5,]
newInter=newGT[which(newGT$genes %in% goodGenes$MARKER_SYMBOL),]
newInter[1:5,]
which(goodGenes$MARKER_SYMBOL==newInter$genes)
which(goodGenes$MARKER_SYMBOL==as.character(newInter$genes))
which(goodGenes$MARKER_SYMBOL==newInter$genes[1])
newInter$genes[1]
which(goodGenes$MARKER_SYMBOL==" 0610007C21Rik")
which(goodGenes$MARKER_SYMBOL=="0610007C21Rik")
which(goodInter$MARKER_SYMBOL=="0610007C21Rik")
g=goodInter[order(goodInter$MARKER_SYMBOL),]
g[1:10,]
rapII=function(x){return(which(g$MARKER_SYMBOL==x,arr.ind=TRUE))}
gIndex=lapply(unique(g$MARKER_SYMBOL),rapII)
gIndex=lapply(unique(g$MARKER_SYMBOL[1:10]),rapII)
gIndex=lapply(unique(g$MARKER_SYMBOL[1:100]),rapII)
#gIndex=lapply(unique(g$MARKER_SYMBOL),rapII)
unique(g$MARKER_SYMBOL)
gIndex=lapply(unique(g$MARKER_SYMBOL),rapII)
gIndex[1]
gIndex[2]
so=function(x){return(g[x])}
so(gIndex[1])
so(gIndex[[1]])
gIndex[2]
g[gIndex[[2]],]
g[gIndex[[2]],]
g$EPD_DISTRIBUTE[gIndex[[2]]]
so=function(x){return(g$EPD_DISTRIBUTE[x])}
lapply(gIndex[1],so)
lapply(gIndex[2],so)
so=function(x){return(sum(g$EPD_DISTRIBUTE[x]))}
uniqEPD=unlist(lapply(gIndex,so))
uniqEPD
mergedData=data.frame(unique(g$MARKER_SYMBOL),uniqEPD)
g$EP_PLATE_NAME[1:5]
g$EP_PLATE_NAME[1:50]
so2=function(x){return(sum(g$TARGETED_TRAP[x]))}
uniqTT=unlist(lapply(gIndex,so2))
mergedData=data.frame(mergedData,TARGETED_TRAP=uniqTT)
mergedData[1:5,]
plateNames=function(x){return(sum(g$EP_PLATE_NAME[x]))}
epPlates=unlist(lapply(gIndex,plateNames))
plateNames=function(x){return(paste(g$EP_PLATE_NAME[x]))}
epPlates=unlist(lapply(gIndex,plateNames))
epPlates[1:10]
epPlates[1:100]
epPlates[1:1000]
length(plateNames)
plateNames
length(epPlates)
length(mergedData$uniqEPD)
length(newInter)
length(newInter$genes)
mergedData$uniqEPD[1:10]
mergedData[1:10,]
newInter[1:10,]
interFinal=data.frame(mergedData,newInter$geneTrap)
#epPlates=unlist(lapply(gIndex,plateNames))
wellNames=function(x){return(sum(g$EPD_WELL_NAME)[x]/g$EPD_WELL_NAME[x])}
wells=unlist(lapply(gIndex,wellNames))
wells
wellNames=function(x){return(sum(g$EPD_WELL_NAME[x])/g$EPD_WELL_NAME[x])}
wells=unlist(lapply(gIndex,wellNames))
wells
wellNames=function(x){return(sum(g$EPD_WELL_NAME[x]))}
wells=unlist(lapply(gIndex,wellNames))
wells
wellNames=function(x){return(g$EPD_WELL_NAME[x])}
wells=unlist(lapply(gIndex,wellNames))
wells
wellNames=function(x){return(sum(g$EPD_WELL_NAME[x]))}
wells=unlist(lapply(gIndex,wellNames))
wells
wells[1:10]
length(gIndex)
wellNamesII==function(x){return(sum(g$EPD_WELL_NAME[x]))}
wellNamesII=function(x){return(sum(g$EPD_WELL_NAME[x]))}
wellsII=unlist(lapply(gIndex,wellNamesII))
wellsII[1:5]
wellsI[1:5]
wells[1:5]
length(wells)
length(wellsII)
wellNamesII=function(x){return(g$EPD_WELL_NAME[x])}
wellsII=unlist(lapply(gIndex,wellNamesII))
length(wellsII)
wellsII[1:5]
wells[1:5]
lapply(gIndex[1:5],print)
lapply(gIndex[1:5],print)
lapply(gIndex[1:5],head)
lapply(gIndex[1:5],floor)
lapply(unlist(gIndex[1:5]),floor)
lapply(gIndex[1:5],floor)
#lapply(gIndex[1:5],floor)
floor(gIndex[[2]])
head(gIndex[[2]])
ceiling(gIndex[[2]])
tail(gIndex[[2]])
initial.digit<-function(x) {substr(as.character(x),1,1)}
lapply(l,initial.digit)
lapply(gIndex[[2]],initial.digit)
f <- function(s) strsplit(s, "_")[[1]][1]
sapply(gIndex[2],f)
#f <- function(s) [[1]][1]
initial.digit<-function(x) {floor(x/10^floor(log(x,10)))}
lapply(gIndex[2],initial.digit)
wellNames=function(x){if (length(x)>1) {return(sum(g$EPD_WELL_NAME[x])/g$EPD_WELL_NAME)} else {return(g$EPD_WELL_NAME)}}
lapply(gIndex[2],wellNames)
#lapply(gIndex[2],wellNames)
wellNames=function(x){if (length(x)>1) {return(sum(g$EPD_WELL_NAME[x])/g$EPD_WELL_NAME[x])} else {return(g$EPD_WELL_NAME[x])}}
lapply(gIndex[2],wellNames)
lapply(gIndex[1],wellNames)
wellNames=function(x){if (length(x)>1) {return(sum(g$EPD_WELL_NAME[x])/length(g$EPD_WELL_NAME[x]))} else {return(g$EPD_WELL_NAME[x])}}
lapply(gIndex[1],wellNames)
lapply(gIndex[2],wellNames)
gIndex[]
gIndex[1:100]
lapply(gIndex[66],wellNames)
g[gIndex[[66]],]
round(lapply(gIndex[66],wellNames))
round(unlist(lapply(gIndex[66],wellNames)))
wells=lapply(gIndex,wellNames)
wells
wellNames=function(x){if (length(x)>1) {return(round(sum(g$EPD_WELL_NAME[x])/length(g$EPD_WELL_NAME[x])))} else {return(g$EPD_WELL_NAME[x])}}
wells=lapply(gIndex,wellNames)
wells
wells=unlist(lapply(gIndex,wellNames))
wells
ls()
length(interFinal$TARGETED_TRAP)
interFinal=data.frame(interFinal,EPD_WELL_NAME=wells)
interFinal[1:10,]
colnames(interFinal)
colnames(interFinal)[1]="MARKER_SYMBOL"
colnames(interFinal)[2]="EPD_DISTRIBUTE"
colnames(interFinal)[4]="Gene Traps"
#colnames(interFinal)[4]="Gene Traps"
effCalc(interFinal)
which(effCalc(interFinal)=='NaN')
interEff=effCalc(interFinal)
which(interEff)=='NaN'
which(interEff=='NaN')
interEff[which(interEff=='NaN')]
interEff[which(interEff=='NaN')]=0
interFinal=data.frame(interFinal,efficiencies=interEff)
interFinal[1:10,]
#cor(interFinal$)
colnames(interFinal)[5]
colnames(interFinal)[4]
colnames(interFinal)[4]="GeneTraps"
cor(interFinal$GeneTraps,interFinal$efficiencies)
plot(interFinal$efficiencies,interFinal$GeneTraps)
length(interFinal$GeneTraps)
maz(interFinal$GeneTraps)
max(interFinal$GeneTraps)
unique(interFinal$GeneTraps)
savehistory('workspace/ikmc_targeting_efficiencies/codes/R/interGrub.R')
length(which(goodGenes$MARKER_SYMBOL %in% newGT$genes))
goodGenes[which(goodGenes$MARKER_SYMBOL %in% newGT$genes),]
goodGenes$MARKER_SYMBOL[which(goodGenes$MARKER_SYMBOL %in% newGT$genes)]
unique(goodGenes$MARKER_SYMBOL[which(goodGenes$MARKER_SYMBOL %in% newGT$genes)])
epClub=read.csv('workspace/ikmc_targeting_efficiencies/tables/epClubs.csv',sep='\t')
epClub[1:10]
epClub[1:10,]
epClub=epClub[order(epClub$MARKER_SYMBOL),]
epClub[1:10,]
length(epClub$MARKER_SYMBOL)
length(which(epClub$MARKER_SYMBOL %in% newGT$genes))
length(which(epClub$MARKER_SYMBOL %in% goodGenes$MARKER_SYMBOL))
length(which(epClub$MARKER_SYMBOL %in% unique(goodGenes$MARKER_SYMBOL)))
eps=epClub[which(epClub$MARKER_SYMBOL %in% goodGenes$MARKER_SYMBOL),]
eps[1:10,]
length(goodGenes$MARKER_SYMBOL)
length(which(goodGenes$MARKER_SYMBOL %in% epClub$MARKER_SYMBOL))
length(unique(goodGenes$MARKER_SYMBOL))
goodGenes[1:10,]
epG=goodGenes[order(goodGenes$MARKER_SYMBOL),]
epG[1:10,]
eps[1:10,]
length(goodInter$MARKER_SYMBOL)
length(which(eps$MARKER_SYMBOL %in% newGT$genes))
length(!which(eps$MARKER_SYMBOL %in% newGT$genes))
length(which(!eps$MARKER_SYMBOL %in% newGT$genes))
length(newInter$genes)
mergedData[1:5,]
length(mergedData$uniqEPD)
length(eps$TARGETED_TRAP)
length(which(eps %in% mergedData$unique.g.MARKER_SYMBOL.))
length(which(eps$MARKER_SYMBOL %in% mergedData$unique.g.MARKER_SYMBOL.))
length(interFinal$EPD_DISTRIBUTE)
length(!which(goodGenes %in% newGT$genes))
length(!which(goodGenes$MARKER_SYMBOL %in% newGT$genes))
length(which(goodGenes$MARKER_SYMBOL %in% newGT$genes))
length(which(!goodGenes$MARKER_SYMBOL %in% newGT$genes))
#eps=epClub[which(epClub$MARKER_SYMBOL %in% goodGenes$MARKER_SYMBOL),]
length(goodGenes$MARKER_SYMBOL)
length(unique(goodGenes$MARKER_SYMBOL))
length(eps$MARKER_SYMBOL)
length(eps$MARKER_SYMBOL %in% newGT$genes)
length(newGT$genes %in% eps$MARKER_SYMBOL)
eps[1:2,]
length(which(eps$MARKER_SYMBOL %in% newGT$genes))
epsGT=eps[which(eps$MARKER_SYMBOL %in% newGT$genes),]
epsGT[1:10,]
length(epsGT$MARKER_SYMBOL)
newGTeps=newGT[which(newGT$genes %in% eps$MARKER_SYMBOL),]
length(newGTeps$genes)
goodGeneswithGT=data.frame(epsGT,newGTeps$geneTrap)
goodGeneswithGT[1,3:]
goodGeneswithGT[1,3,]
goodGeneswithGT[1:3,]
goodGeneswithNOGT=eps[which(!eps$MARKER_SYMBOL %in% newGT$genes),]
length(goodGeneswithNOGT)
length(goodGeneswithNOGT$MARKER_SYMBOL)
goodGeneswithNOGT[1:5,]
goodGeneswithZeroGT=data.frame(goodGeneswithNOGT,geneTRaps=rep(0,length(goodGeneswithNOGT$MARKER_SYMBOL)))
goodGeneswithZeroGT[1:5,]
goodGeneswithZeroGT=data.frame(goodGeneswithNOGT,geneTraps=rep(0,length(goodGeneswithNOGT$MARKER_SYMBOL)))
goodGenes_mixedGT=data.frame(goodGeneswithGT,goodGeneswithNOGT)
goodGenes_mixedGT=rbind(goodGeneswithGT,goodGeneswithNOGT)
goodGenes_mixedGT=rbind(goodGeneswithGT,goodGeneswithZeroGT)
goodGeneswithZeroGT[1:4,]
goodGeneswithGT[1:4,]
colnames(goodGeneswithGT)[6]
colnames(goodGeneswithGT)[7]
colnames(goodGeneswithGT)[7]="geneTraps"
goodGenes_mixedGT=rbind(goodGeneswithGT,goodGeneswithZeroGT)
goodGenes_mixedGT[1:10,]
length(goodGenes_mixedGT$geneTraps)
length(unique(goodGenes_mixedGT$MARKER_SYMBOL))
#plot(goodGenes_mixedGT)
effCalc(goodGenes_mixedGT)
ggMT=effCalc(goodGenes_mixedGT)
which(ggMT=='NaN')
ggMT[which(ggMT=='NaN')]
ggMT[which(ggMT=='NaN')]=0
length(ggMT)
ggMT
sum(ggMT)
sum(ggMT)/length(ggMT)
goodGenes_mixedGTe=data.frame(goodGenes_mixedGT,efficiencies=ggMT)
goodGenes_mixedGTe[1:5,]
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies)
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies,xlim=c(0,100))
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies,xlim=c(0,20))
cor(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies)
cor(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$geneTraps,col=)
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/GTvsEFF_goodGenes_MIXGT')
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies,xlim=c(0,100),col="green",xlab="Gene Trap Counts",ylab="Targeting Efficiencies",main="Correlation is 0.12, Each circle is a gene")
dev.off()
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/GTvsEFF_goodGenes_MIXGT')
dev.off()
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/GTvsEFF_goodGenes_MIXGT')
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies,xlim=c(0,100),col="darkgreen",xlab="Gene Trap Counts",ylab="Targeting Efficiencies",main="Correlation is 0.12, Each circle is a gene")
dev.off()
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/GTvsEFF_goodGenes_MIXGT_100')
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies,xlim=c(0,100),col="darkgreen",xlab="Gene Trap Counts",ylab="Targeting Efficiencies",main="Correlation is 0.12, Each circle is a gene")
dev.off()
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/GTvsEFF_goodGenes_MIXGT')
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies,col="darkgreen",xlab="Gene Trap Counts",ylab="Targeting Efficiencies",main="Correlation is 0.12, Each circle is a gene")
dev.off()
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/GTvsEFF_goodGenes_MIXGT_20')
plot(goodGenes_mixedGTe$geneTraps,goodGenes_mixedGTe$efficiencies,col="darkgreen",xlab="Gene Trap Counts",ylab="Targeting Efficiencies",main="Correlation is 0.12, Each circle is a gene",xlim=c(0,20))
dev.off()
hist(goodGenes_mixedGT$MARKER_SYMBOL)
hist(goodGenes_mixedGTe$efficiencies)
goodGenes_mixedGT=goodGenes_mixedGTe
hist(goodGenes_mixedGT$efficiencies,col="lightblue")
hist(goodGenes_mixedGT$efficiencies,col="lightblue",border="pink")
hist(goodGenes_mixedGT$efficiencies,col="lightblue",border="pink",xlab="Targeting Efficiencies")
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/effvsFre_goodGenes_MIXGT')
hist(goodGenes_mixedGT$efficiencies,col="lightblue",border="pink",xlab="Targeting Efficiencies")
dev.off()
homo=read.csv('workspace/ikmc_targeting_efficiencies/tables/homology_arm_repeats.csv')
length(homo$Design.ID)
homo[1:5,]
homo[1:5,1:5]
homo[1:5,1:6]
homo[1:5,1:7]
homo[1:5,1:6]
homoBuffer=homo[,1:6]
savehistory('workspace/ikmc_targeting_efficiencies/codes/R/interGrub.R')
goodGenes_mixedGT[1:4,]
rm(goodGenes_mixedGTe)
write.table(goodGenes_mixedGT,file='workspace/ikmc_targeting_efficiencies/tables/goodGenes_mixedGT.csv',sep='\t',quote=FALSE,col.names=colnames(goodGenes_mixedGT))
length(which(homoBuffer %in% goodGenes_mixedGT$MARKER_SYMBOL))
length(which(homoBuffer$Marker.Symbol %in% goodGenes_mixedGT$MARKER_SYMBOL))
homoBuf=homoBuffer[which(homoBuffer$Marker.Symbol %in% goodGenes_mixedGT$MARKER_SYMBOL)),]
homoBuf=homoBuffer[which(homoBuffer$Marker.Symbol %in% goodGenes_mixedGT$MARKER_SYMBOL),]
homoBuf[1:10,]
homoBuf=homoBuf[order(homoBuf$Marker.Symbol),]
homoBuf[1:10,]
homoBuf[1:2,]
homoRepeater=function(x){return(x$G5.Homology.Arm.Repeat.Length/x$G5.Homology.Arm.Length)}
homoRepeater(homoBuf)
sum(homoRepeater(homoBuf))
#data.framesum(homoRepeater(homoBuf))
homoBufII=data.frame(homoBuf,repeatEff=sum(homoRepeater(homoBuf)))
homoBufII[1:10,]
homoBufII[1:5,]
homoBufII[1:5,]
goodGenes_mixedGT[1:5,]
homoBuffer=homo[,1:6]
homoBuffer[1:10,]
homoBuf=homoBuffer[order(homoBuffer$Marker.Symbol),]
length(which(homoBuf$Marker.Symbol %in% goodGenes_mixedGT$MARKER_SYMBOL))
homoBufII=homoBuf[which(homoBuf$Marker.Symbol %in% goodGenes_mixedGT$MARKER_SYMBOL),]
homoBufII[1:5,]
goodGenes_mixedGT[1:5,]
#order(homoBufII)
homoBufII[1:5,]
homoBufII=data.frame(homoBufII,repeatEff=sum(homoRepeater(homoBufII)))
homoBufII$Marker.Symbol[1:10]
goodGenes_mixedGT$MARKER_SYMBOL[1:10]
a=goodGenes_mixedGT[which(goodGenes_mixedGT$MARKER_SYMBOL %in% homoBufII$Marker.Symbol),]
b=as.data.frame(homoBufII,a)
#b=as.data.frame(homoBufII,a)
a[1:10,]
b=as.data.frame(homoBufII,a$efficiencies)
length(a$efficiencies)
a=length(which(goodGenes_mixedGT$MARKER_SYMBOL %in% homoBufII$Marker.Symbol))
length(which(goodGenes_mixedGT$MARKER_SYMBOL %in% homoBufII$Marker.Symbol))
length(which(unique(goodGenes_mixedGT$MARKER_SYMBOL) %in% homoBufII$Marker.Symbol))
length(which(homoBufII$Marker.Symbol %in% goodGenes_mixedGT))
length(which(homoBufII$Marker.Symbol %in% goodGenes_mixedGT$MARKER_SYMBOL))
length(which(homoBufII$Marker.Symbol %in% unique(goodGenes_mixedGT$MARKER_SYMBOL)))
length(homo$Marker.Symbol)
homo$Marker.Symbol[1:10]
homo=homo[order(homo$Marker.Symbol),]
homo$Marker.Symbol[1:10]
homo[1:10,]
#homo=homo[order(homo$Marker.Symbol),]
ls()
design<-read.csv('workspace/ikmc_targeting_efficiencies/tables/designIDgenes.csv')
design[1:10,]
length(which(design$MARKER_SYMBOL %in% goodGenes_mixedGT$MARKER_SYMBOL))
designInter=design[which(design$MARKER_SYMBOL %in% goodGenes_mixedGT$MARKER_SYMBOL),]
length(design$MARKER_SYMBOL)
length(designInter$MARKER_SYMBOL)
unique(designInter$MARKER_SYMBOL)
length(unique(designInter$MARKER_SYMBOL))
design[1:10,]
length(which(goodGenes_mixedGT$MARKER_SYMBOL %in% design$MARKER_SYMBOL))
goodGenes[which(goodGenes_mixedGT$MARKER_SYMBOL %in% design$MARKER_SYMBOL),]
which(goodGenes_mixedGT$MARKER_SYMBOL %in% design$MARKER_SYMBOL)[1:10]
which(goodGenes_mixedGT$MARKER_SYMBOL %in% design$MARKER_SYMBOL)[1:100]
length(unique(designInter$DESIGN_ID))
length(unique(designInter$MARKER_SYMBOL))
goodGenes_mixedGT[1:5,]
rap=function(x){return(which(designInter$DESIGN_ID==x,arr.ind=TRUE))}
lapply(unique(designInter$DESIGN_ID[1:100]),rap)
lapply(unique(designInter$DESIGN_ID[1:1000]),rap)
designIndex=lapply(unique(designInter$DESIGN_ID),rap)
designIndex[1:10]
#designInter$DESIGN_ID
ids=function(x){designInter$DESIGN_ID[x]}
lapply(designInter$DESIGN_ID[1:100],ids)
lapply(designIndex[1:100],ids)
designUniq=lapply(designIndex,ids)
designUniq
lapply(designUniq,unique)
unlist(lapply(designUniq,unique))
length(unlist(lapply(designUniq,unique)))
length(unique(designInter$DESIGN_ID))
unique(designInter$DESIGN_ID)[1]
unique(designInter$DESIGN_ID)[1:2]
length(goodGenes_mixedGT$efficiencies)
length(unlist(lapply(designUniq,unique)))
uniqD=lapply(designUniq,unique)
lapply(uniqD,length)
unlist(lapply(uniqD,length))
unique(designInter$DESIGN_ID)[1:2]
savehistory('workspace/ikmc_targeting_efficiencies/codes/R/design.R')
