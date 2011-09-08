gtLength<-function(x){length}

# raw gt file 
MRK_GT=read.csv('workspace/ikmc_targeting_efficiencies/tables/MRK_GeneTrap.rpt',sep='\t',header=FALSE)

# counting GeneTraps in raw file
gtCount<-read.table('workspace/ikmc_targeting_efficiencies/tables/MRK_GeneTrapGT.csv',sep='\t',header=FALSE,na.strings='')

# recording the geneTrap counts from the new data
gtCounts<-function(x){return(length(which(!is.na(gtCount[x,]))))}
gt=lapply(seq(1,length(gtCount[,1])),gtCounts)

igo=function(x){return(which(!is.na(x)==TRUE))}


# getting intersected gene trap count
epGenes=epClub$MARKER_SYMBOL
ep=epClub[which(epGenes %in% newGT$genes),]

# these are the genes which are not present in new trap data but present in the ep data (so have eff counts)
nonGT=epClub[which(!epGenes %in% newGT$genes),]

# genes which match epGenes
GT=newGT[which(newGT$genes %in% epGenes),]

# new data with matched genes from both data's and having gene trap and efficiencies values as well
newData=data.frame(ep,geneTraps=GT$geneTrap)

# filling up geneTrap values as 0
epFull=data.frame(nonGT,geneTraps=rep(0,length(nonGT$MGI_GT_COUNT)))
# binding ep+epFull, includes 0 for which geneTrap counts are not available
epFull=rbind(newData,epFull)

# finding efficiencies and putting 0 for NaN's
effCalc=function(x){return((x$EPD_DISTRIBUTE+x$TARGETED_TRAP)/x$EPD_WELL_NAME)}
eff=effCalc(newData)
eff[which(eff=="NaN")]=0

# geneClub merged with EP_PLATE_NAME data with new gene trap data for fetching plate counts
geneClub=read.csv('geneClubs.csv',sep='\t')
geneClubGT=geneClub[which(geneClub$EP_PLATE_NAME %in% epFull$EP_PLATE_NAME),]
geneClubGT=geneClubGT[order(geneClubGT$EP_PLATE_NAME),]
geneClubGTeff=effCalc(geneClubGT)
geneClubGTeff[which(geneClubGTeff=='NaN')]=0
# putting the efficiencies back into the table with no NAN in it
geneClubGT=geneClubGT[,-7]
geneClubGT=data.frame(geneClubGT,efficiencies=geneClubGTeff)

# gives the index of electroporation plates
rap=function(x){return(which(geneClubGT$EP_PLATE_NAME==x,arr.ind=TRUE))}
uniqueGTEPplateIndex=lapply(unique(geneClubGT$EP_PLATE_NAME),rap)

# averaging the efficiencies of genes/per plate to filter out rubbish plates
epEff=function(x){return(sum(geneClubGT$efficiencies[x])/length(x))}
epEffs=unlist(lapply(uniqueGTEPplateIndex,epEff))

# filtering the plates with eff >0.1
goodEP=unlist(epEffs[which(unlist(epEffs)>0.1)])
rubbishPlates=which(unlist(epEffs)<0.1) #133 plates, pull all genes from this plate and discard them
goodPlates=which(unlist(epEffs)>0.1)
# filtering rubbishGenes
rubbishGenes=function(x){return(uniqueGTEPplateIndex[x])}
# calculating rubbish genes eff
rubGeneseff=geneClubGT$efficiencies[unlist(rubbishGenes(rubbishPlates))]
# calculating rubbish gene names
rubGenes=geneClubGT$MARKER_SYMBOL[unlist(rubbishGenes(rubbishPlates))]
# identify rubbish genes in the new data and pull them out for plotting
identifyGenes=function(x){return(which(newData$MARKER_SYMBOL==x))}
rubIndex=unlist(lapply(rubGenes,identifyGenes))
# filtered out good genes and made their data frame from newData
goodGenes=newData[-rubIndex,]
goodGenesEff=effCalc(goodGenes)
goodGenes=data.frame(goodGenes,efficiencies=goodGenesEff)


# plotting histogram of efficiencies
hist(eff,col="lightblue",border="gray",ylim=c(0,2000))

# plotting goodGenes (GT vs EFF)
plot(goodGenes$efficiencies,goodGenes$geneTraps,ylim=c(0,100),xlab="Targeting Efficiencies",ylab="Gene Trap Counts",main="GeneTraps vs Efficiencies",sub=paste(length(goodGenes$efficiencies),"genes"))

