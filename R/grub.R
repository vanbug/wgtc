# Grub R script for htgt data
# Author : Sukhdeep Singh
# Organization : Sanger Center
################################################################################################

# read in htgtAll.csv - non null marker symbols

# variable declarations
epdLoop<-list()
htgt<-read.csv('htgt.csv')
sym=htgt$MARKER_SYMBOL

epdLength<-function(x,y){return(which(x==y))}
counter<-function(x,y){return(which(x==y))}

for (i in 1:length(genes)){print (which(sym==genes[i]))}

# or
symIndex=which(duplicated(sym,incomparables=FALSE,MARGIN=1)==FALSE)
# index is same
sym[tail(symIndex)]==tail(genes)
epdWells=c(diff(symIndex),1)

# correcting length discrepancy
symIndex[length(symIndex)+1]=length(tab$MARKER_SYMBOL)
grub<-function(x,y,z){return(which(x[which(y==z)]=='yes'))}
# constructing table
tab=cbind(as.character(genes),epdWells)

for (i in 1:length(genes)){
which(sym==genes[i])
}

# fetching numbers for the ditributed and targeted traps -- use apply

tt<-list();epd<-list()
for (i in 1:length(genes)){
epd[[i]]=which(htgt$EPD_DISTRIBUTE[which(sym==genes[i])]=='yes')
tt[[i]]=which(htgt$TARGETED_TRAP[which(sym==genes[i])]=='yes')
print (paste(length(genes)-i,"left"))
}
epdDis=unlist(lapply(epd,length))
tTraps=unlist(lapply(tt,length))


grub<-function(x,y,z){return(which(x[which(y==z)]=='yes'))}

grub2<-function(x,y){return(which(x==y))}

rake<-function(x,y,z){which(x[which(y==z)]=='yes')}

whc<-function(x){return(which(x==0))}

# naming columns of inputed htgt parsed file
colnames(htgtParsed)=c("genes","epWellNames","epdWellNames","epdDistributes","targetedTraps","geneTraps")

################
epClubsSorted=epClubs[order(epClubs$MARKER_SYMBOL),]


fractionAge<-function(x,y,z){return((x+y)/z)}

a<-function(x,y){return(which(x==y))}

# nested apply
yed=function(x){return(which((mapply(match,x,geneClubs$EP_PLATE_NAME))=='1',arr.ind=TRUE))}
red=function(x){return(yed(x))}


# 2nd.Sep.2011
# Targeting Efficiencies vs EP well names
geneFractionage=fractionAge(geneClubs$EPD_DISTRIBUTE,geneClubs$TARGETED_TRAP,geneClubs$EPD_WELL_NAME)
# they has NaN as geneClubs$EPD_DISTRIBUTE is empty for some
geneFractionage[which(geneFractionage=='NaN')]=0

# returns unique value index
# sort geneClubs always as it is not sorted in the csv file
geneClubsEP=geneClubs[order(geneClubs$EP_PLATE_NAME),]
rap=function(x){return(which(x==geneClubsEP$EP_PLATE_NAME,arr.ind=TRUE))}

epPlateNameIndex=lapply(unique(geneClubsEP$EP_PLATE_NAME),rap)
epPlateIndex=unlist(lapply(epPlateNameIndex,length))
# sorts the data.frame with counts of epPlateNames per epPlates
epPlateName=data.frame(epPlateNames=as.character(unique(geneClubsEP$EP_PLATE_NAME)),epCounts=epPlateIndex)

# % efficiencies per epPlateName
oriGeneEPNameIndex=unique(geneClubsEP$EP_PLATE_NAME)

unwindEPDwell=function(x){return(sum(geneClubsEP$EPD_WELL_NAME[unlist(x)]))}
unwindEPDDis=function(x){return(sum(geneClubsEP$EPD_DISTRIBUTE[unlist(x)]))}
unwindEPDTT=function(x){return(sum(geneClubsEP$TARGETED_TRAP[unlist(x)]))}

# storing EPD_wells and EPD_dis and EPD_TTs from geneClubs for each EP plate names - we have to do other way around
geneClubEPDWELL=unlist(lapply(epPlateNameIndex,unwindEPDwell))
geneClubEPDDis=unlist(lapply(epPlateNameIndex,unwindEPDDis))
geneClubTT=unlist(lapply(epPlateNameIndex,unwindEPDTT))
# calculating fractionage
geneClubFractionage=fractionAge(geneClubEPDDis,geneClubTT,geneClubEPDWELL)
# replacing NAN by zero
geneClubFractionage[which(geneClubFractionage=='NaN')]=0

# combining targeting efficiencies per EP Plate name with EP Plate name counts


plot(epF$MGI_GT_COUNT,epF$fractionAge,xlab="Gene trap count",ylab="Fractions of electroporations",main="Gene Traps counts vs Fractions of Electroporations((dis+tarTraps)/total electroporations)",sub="Each point is a gene",col="darkgreen")

plot(epF$MGI_GT_COUNT,epF$fractionAge,xlab="Gene trap count",ylab="Fractions of electroporations",xlim=c(0,100))


# fetches genes per EP plate
naive<-function(x){return(geneClubsEP[x,])}
genesperEP=lapply(epPlateNameIndex,naive)

#
gEPeff=function(x){return((x$EPD_DISTRIBUTE+x$TARGETED_TRAP)/x$EPD_WELL_NAME)}
gEPeffIndex=lapply(genesperEP,gEPeff)

# finding which is NaN as EPD_DISTRIBUTE is 0 and replacing the NaN
replaceNaN=function(x){return(as.numeric(gsub('NaN',0,x)))}
geneEff=lapply(gEPeffIndex,replaceNaN)
















