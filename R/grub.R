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







































