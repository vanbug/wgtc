# inputting data with design ID's and putting NA's fro blank spaces
designII<-read.table('workspace/ikmc_targeting_efficiencies/tables/designIDgenes.csv',blank.lines.skip=NA,sep=',',header=TRUE,na.strings='')

# read file without calling on level using stringsAsFactors
design<-read.table('workspace/ikmc_targeting_efficiencies/tables/htgt/htgtData_designID_organised_NA.csv',sep='\t',header=TRUE,stringsAsFactors=FALSE)

# replace yes in EPD_DISTRIBUTE and TARGETED_TRAP with 1
#epdDis=gsub('yes',as.numeric(1),designII$EPD_DISTRIBUTE)
#tt=gsub('yes',as.numeric(1),designII$TARGETED_TRAP)

# find index and order of geneClub by EP Plates
geneClub=design
geneClub=geneClub[order(geneClub$EP_PLATE_NAME),]
rap=function(x){return(which(geneClub$EP_PLATE_NAME==x,arr.ind=TRUE))}
geneClubEPIndex=lapply(unique(geneClub$EP_PLATE_NAME),rap)

# index of genes per plate
gener=function(x){return(which(geneClub$MARKER_SYMBOL[geneClubEPIndex[[i]]]==x))}

# variable decalarations
epdChar<-list();ttChar<-list();epdWellsChar<-list();epdWells<-list();epd<-list();tt<-list();generIndex<-list()

# gene and plate Indexing function for a superloop to fetch, gene data per plate
# gives design index per gene
geneIndexer=function(x,y){return(geneClubEPIndex[[y]][x])}
# gives gene index per plate
plateIndexer_EPD=function(x){return(geneClub$EPD_DISTRIBUTE[x])}
plateIndexer_TT=function(x){return(geneClub$TARGETED_TRAP[x])}
plateIndexer_wells=function(x){return(geneClub$EPD_WELL_NAME[x])}

yesSir=function(x){return(which(x=='yes'))}
totalWells=function(x){return(which(!is.na(x)))}
for (i in 1:length(geneClubEPIndex)){
generIndex[[i]]=lapply(unique(geneClub$MARKER_SYMBOL[geneClubEPIndex[[i]]]),gener)
epdChar[[i]]=lapply(lapply(generIndex[[i]],y=i,geneIndexer),plateIndexer_EPD)
epd[[i]]=unlist(lapply(lapply(epdChar[[i]],yesSir),length))
ttChar[[i]]=lapply(lapply(generIndex[[i]],y=i,geneIndexer),plateIndexer_TT)
tt[[i]]=unlist(lapply(lapply(ttChar[[i]],yesSir),length))
epdWellsChar[[i]]=lapply(lapply(generIndex[[i]],y=i,geneIndexer),plateIndexer_wells)
epdWells[[i]]=unlist(lapply(lapply(epdWellsChar[[i]],totalWells),length))
print (length(geneClubEPIndex)-i)
}

# we average the tt, epd and wells to find efficiency per plate per gene
effs=function(x,y,z){if (z==0){return (rep(0,length(z)))} else {return(((x+y)/z))}}
eff=mapply(effs,tt,epd,epdWells)

# na remover from efficiency list
for (i in 1:length(eff)){eff[[i]][which(unlist(lapply(eff[[i]],function(x){return(x)}))=='NaN')]=0}

# two directions
# direction I - get all genes per plate after removing rubbish plates and ranking efficiencies of same genes same ID's different plates
# averaging eff
avg=lapply(eff,ave)
effPlate=unlist(lapply(avg,function(x)return(tail(x,1))))

# get plates with rubbish eff (<0.1)
plates=data.frame(epPlates=unique(geneClub$EP_PLATE_NAME),targetingEfficiencies=effPlate)
rubbishPlates=plates$epPlates[which(plates$targetingEfficiencies<0.1)]

# finding unique genes per plate to make the geneClub file - rerun unique indexing
uniqGenes=function(x){return(unique(geneClub$MARKER_SYMBOL[x]))}
uniqGenesperPlate=lapply(geneClubEPIndex,uniqGenes)
uniqID=function(x){return(unique(geneClub$DESIGN_ID[x]))}
uniqIDperPlate=lapply(geneClubEPIndex,uniqID)

# gives same genes on same plate with different efficiencies, remove them manually and then reindex the genes per plate using superloop and uni indexing
sameGeneMultiID=which(unlist(lapply(uniqIDperPlate,length))-unlist(lapply(uniqGenesperPlate,length))!=0)

# removing same genes multi ID's (325,582,628)

# we know Hnf4g has multi id's find the eff and remove the less efficient - for viewing
#geneClub[geneClubEPIndex[[582]][which(geneClub$MARKER_SYMBOL[geneClubEPIndex[[582]]]=='Hnf4g')],]

#removing that from whole club - for removing # follow same steps for new data
geneClub=geneClub[-geneClubEPIndex[[582]][218:226],]
geneClub=geneClub[-geneClubEPIndex[[325]][341],]
geneClub=geneClub[-geneClubEPIndex[[628]][151:160],]
geneClub=geneClub[-geneClubEPIndex[[628]][171:177],]
geneClub=geneClub[-geneClubEPIndex[[628]][183:197],]
# if you run these removals , rerun the unique indxing and eff calclulator (for eff, rerun the superloop)

# QA - quality analysis for previous loop
#length(which(geneClub$EPD_DISTRIBUTE[geneClubEPIndex[[2]][generIndex[[2]][[1]]]]=='yes'))

# storing new values in data frame - we dont have EPplate names here so have to fetch it in next step
geneClubII=data.frame(genes=unlist(uniqGenesperPlate),designID=unlist(uniqIDperPlate),epPlates=unique(geneClub$EP_PLATE_NAME),eff=unlist(eff))

# calculating and adding plate index to the new data.frame
uniqPlate=function(x){return(unique(geneClub$EP_PLATE_NAME[x]))}
uniqPlates=lapply(geneClubEPIndex,uniqPlate)
plateIndexII<-list()
for(i in 1:length(uniqGenesperPlate)){plateIndexII[[i]]=rep(uniqPlates[[i]],length=length(uniqGenesperPlate[[i]]))}

# making data.frame
geneClubII=data.frame(genes=unlist(uniqGenesperPlate),designID=unlist(uniqIDperPlate),plates=unlist(plateIndexII),eff=unlist(eff))
#geneClubII<-data.frame(geneClubII,plates=unlist(plateIndexII))

# removing same gene, same design II on different plates by rankings based on efficiencies#
joy=geneClubII[order(geneClubII$designID),]

# direction - I, getting ranking by max efficiency os same genes same ID's on different plates
# finding gene index again and max efficiency
geneIndexF=function(x){which(joy$genes==x,arr.ind=TRUE)}
geneIndex=lapply(unique(joy$genes),geneIndexF)
maxEff=function(x){which.max(joy$eff[x])}
maxEffIndex=lapply(geneIndex,maxEff)


###################
# Direction -II - get all genes after removing rubbish genes and find same genes on different plates prefoming bad , eff<0.1
worstavgEff=lapply(lapply(worstGeneIndex,worstavgEffF),function(x)return(tail(x,1)))
avgEffGeneallPlates=lapply(lapply(worstGeneIndex,worstavgEffF),tailer)

# making a data.frame of gene by gene efficiencies averaged on diff EP plates
mediumgeneEff=data.frame(unique(mediumGenes$genes),unlist(avgEffGeneallPlates)) 

# returing good genes based on maximum efficiency index and storing joy
yp=function(x){return(geneIndex[[x]][maxEffIndex[[x]]])}
geneIDefficientIndex=unlist(lapply(seq(1,length(geneIndex),by=1),yp))
joy=joy[geneIDefficientIndex,]

# we have rubbishPlates from top and now we will remove them from joy data
# we will observe some values 0 in rubPIndex as we have removed some plates in previous step
rubP=function(x){which(joy$epPlates==x)}
rubbishGenesIndex=lapply(rubbishPlates,rubP) # gives us rubbishGenesIndex

# removing rubbish Plates from all gene data (20K genes)
rubG=function(x){which(geneClubII$epPlates==x)}
rubbishGenesIndexG=lapply(rubbishPlates,rubP) # gives us rubbishGenesIndex
mediumGenes=geneClubII[-unlist(rubbishGenesIndexG),]

# best genes subset of mediumGenes 
bestGenes=mediumGenes[which(mediumGenes$eff>0.9),]

# gene index for medium genes
worstG=function(x){which(mediumGenes$genes==x)}
wGeneIndex=lapply(unique(mediumGenes$genes),worstG) # normal gene Index

# tailer function - return tail of list useful for avg etc.
tailer=function(x)return(tail(x,1))

# we also average the same genes on different EP plates having eff <0.1
worstavgEffF=function(x){mediumGenes$eff[x]}
samegenesdiffPlateEff=lapply(wGeneIndex,worstavgEffF) # gene efficiencies on different plates (subsetted)

# filters out same gene multi plate class
so=function(x){return(which(length(x)>1,arr.ind=TRUE))}
k=lapply(samegenesdiffPlateEff,so)

# filters true values
k2=which(k==1)
sameGeneDiffPlates=mediumGenes[k2,]
sameGeneDiffPlates=sameGeneDiffPlates[order(sameGeneDiffPlates$genes),]
sameIndexF=function(x){return(which(sameGeneDiffPlates$genes==x))}
sameIndex=lapply(unique(sameGeneDiffPlates$genes),sameIndexF)
# index working
#sameGeneDiffPlates[sameIndex[[10]],]

#multiEff=samegenesdiffPlateEff[k2]

# this function returns a value if all the gene efficiencies on all plates <0.1
so2=function(x){b=which(x<0.1);if (length(b)==length(x)){return(b)}}
worstGenes=lapply(sameGeneDiffPlates$eff,so2)
# wg is the genes on diff plate having low efficiency everywhere
wg=lapply(worstGenes,function(x){return(which(!is.null(x)==TRUE))})
wg2=which(wg==1)
totalrubbishGenes=sameGeneDiffPlates[wg2,]

# QA - rubbish efficiency same gene different plates
#rubbishEffsamegenediffPlateEff=samegenesdiffPlateEff[k2[wg2]]

# correlation between best gene and worst gene data and repeat counts
# for best genes
bestInter=bestGenes[which(bestGenes$designID %in% homo$Design.ID),] # bestInter has two plates best
homoBestInter=homo[which(homo$Design.ID %in% bestGenes$designID),]
homoBestInter=homoBestInter[order(homoBestInter$Design.ID),]
bestInter=bestInter[order(bestInter$designID),]


# gene and design is same
#bestInter[237:238,]
# best Inter index
bestIindex=function(x){which(bestInter$designID==x)}
bestInterID=lapply(unique(bestInter$designID),bestIindex)

repeat500Best<-list();repeat1000Best<-list()
for (i in 1:length(homoBestInter$Design.ID)){
repeat500Best[[i]]=rep(homoBestInter$Repeats.500[i]/500,length=length(bestInterID[[i]]))
repeat1000Best[[i]]=rep(homoBestInter$Repeats.1000[i]/1000,length=length(bestInterID[[i]]))
}
besthomoInter=data.frame(bestInter,repeat500=unlist(repeat500Best),repeat1000=unlist(repeat1000Best))
# removing duplicated ID, might not do it, shouldn't do it
#bestInter=bestInter[-238,]

# plotting 5" repeats - 500 KB
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/repeats/bestGenesvsRepeats')
plot(besthomoInter$repeat500,besthomoInter$eff,main=paste("Corelation =",round(cor(besthomoInter$repeat500,besthomoInter$eff),digits=3)),xlab="5\" Homology Arm Repeats_500",ylab="Targeting Efficiency of best perfoming genes eff>0.9")
dev.off()

# plotting 5" repeats - 1000 KB
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/repeats/bestGenesvsRepeats_1000')
plot(besthomoInter$repeat1000,besthomoInter$eff,main=paste("Corelation =",round(cor(besthomoInter$repeat1000,besthomoInter$eff),digits=3)),xlab="5\" Homology Arm Repeats_1000",ylab="Targeting Efficiency of best perfoming genes eff>0.9")
dev.off()

# repeat correlation for worst perfoming genes
homoWorstInter=homo[which(homo$Design.ID %in% totalrubbishGenes$designID),] #1930
worstInter=totalrubbishGenes[which(totalrubbishGenes$designID %in% homo$Design.ID),] #2250
worstInter=worstInter[order(worstInter$designID),]
homoWorstInter=homoWorstInter[order(homoWorstInter$Design.ID),]

# getting good genes
goodGenes=joy[-unlist(rubbishGenesIndex),]

# working with repeats

# reading 5P arm repeats for 500 and 1000 bases
#homo<-read.csv('workspace/ikmc_targeting_efficiencies/tables/repeats/fivep_arm_repeats.csv')

# reading all 3P and 5P repeats
#homoII<-read.csv('workspace/ikmc_targeting_efficiencies/tables/repeats/homology_arm_repeats.csv')
#homoII=homoII[,-c(7:8)] # removes repeat class column for better viewing

goodInter=goodGenes[which(goodGenes$designID %in% homo$Design.ID),]
homoInter=homo[which(homo$Design.ID %in% goodGenes$designID),]
goodInter=goodGenes[order(goodGenes$designID),]
homoInter=homoInter[order(homoInter$Design.ID),]

# getting five prime 500 and 1000 repeat numbers
fiveP5=(homoInter$Repeats.500/500)
fiveP10=(homoInter$Repeats.1000/1000)

# plotting repeats
plot(fiveP5,goodInter$eff,xlab="Five Prime Repeats",ylab="Targeting Efficiencies")
plot(fiveP10,goodInter$eff,xlab="Five Prime Repeats-1000",ylab="Targeting Efficiencies")

# standard error
sd(homoBestInter$Repeats.500/500)/sqrt(length(homoBestInter$Repeats.500/500))

# worst Inter index
worstIindex=function(x){which(worstInter$designID==x)}
worstInterDesignID=lapply(unique(worstInter$designID),worstIindex)

repWorst500<-list();repWorst1000<-list()
for (i in 1:length(worstInterDesignID)){
repWorst500[[i]]=rep(homoWorstInter$Repeats.500[i]/500,length=length(worstInterDesignID[[i]]))
repWorst1000[[i]]=rep(homoWorstInter$Repeats.1000[i]/1000,length=length(worstInterDesignID[[i]]))
}
# new data frame for repeat efficiency intersection
worstInter=data.frame(worstInter,repeat500=unlist(repWorst500),repeat1000=unlist(repWorst1000))
 # 500 repeat vs eff graph
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/repeats/worstGenesvsRepeats500')
plot(worstInter$repeat500,worstInter$eff,xlab="Repeats_500",ylab="Targeting efficiency for worst genes",main=paste("Correlation=",round(cor(worstInter$repeat500,worstInter$eff),digits=3)))
dev.off()

# 1000 repeat vs eff graph
cairo_pdf('workspace/ikmc_targeting_efficiencies/plots/newData/repeats/worstGenesvsRepeats1000')
plot(worstInter$repeat1000,worstInter$eff,xlab="Repeats_1000",ylab="Targeting efficiency for worst genes",main=paste("Correlation=",round(cor(worstInter$repeat1000,worstInter$eff),digits=3)))
dev.off()
########################
# making sheets of best genes, eff, trap numbers and repeat eff
# BEST
gt<-read.csv('workspace/ikmc_targeting_efficiencies/tables/geneTraps/MRK_geneTrap_numbers.csv',sep='\t')
bestGeneTrap=gt[which(gt$genes %in% besthomoInter$genes),]
besthomoInterGT=besthomoInter[which(besthomoInter$genes %in% gt$genes),]
# ordering is very important
bestGeneTrap=bestGeneTrap[order(bestGeneTrap$genes),]
besthomoInterGT=besthomoInterGT[order(besthomoInterGT$genes),]
# finding index
besthomer=function(x){which(besthomoInterGT$genes==x)}
besthomoInterID=lapply(unique(besthomoInterGT$genes),besthomer)

gts<-list()
for (i in 1:length(besthomoInterID)){
gts[[i]]=rep(bestGeneTrap$geneTrap[i],length=length(besthomoInterID[[i]]))
}

bestGeneTrapRepeats=data.frame(besthomoInterGT,geneTraps=unlist(gts))
#writing table
write.table(bestGeneTrapRepeats,'workspace/ikmc_targeting_efficiencies/tables/bestComplete.csv',sep='\t',quote=FALSE,row.names=FALSE)

# WORST
worstInterGT=worstInter[which(worstInter$genes %in% gt$genes),] ##1164
worstGT=gt[which(gt$genes %in% worstInter$genes),]		##1014
worstInterGT=worstInterGT[order(worstInterGT$genes),]
worstGT=worstGT[order(worstGT$genes),]
# finding index
indexIf=function(x){which(worstInterGT$genes==x)}
indexI=lapply(unique(worstInterGT$genes),indexIf)

gtsw<-list()
for(i in 1:length(indexI)){
gtsw[[i]]=rep(worstGT$geneTrap[i],length(indexI[[i]]))
}
worstGeneTrapRepeats=data.frame(worstInterGT,geneTraps=unlist(gtsw))
#writing table
write.table(worstGeneTrapRepeats,'workspace/ikmc_targeting_efficiencies/tables/worstComplete.csv',sep='\t',quote=FALSE,row.names=FALSE)




