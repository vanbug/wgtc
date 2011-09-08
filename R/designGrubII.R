# inputting data with design ID's and putting NA's fro blank spaces
designII<-read.table('workspace/ikmc_targeting_efficiencies/tables/designIDgenes.csv',blank.lines.skip=NA,sep=',',header=TRUE,na.strings='')

# read file without calling on level using stringsAsFactors
design<-read.table('workspace/ikmc_targeting_efficiencies/tables/designSortedEP.csv',sep='\t',header=TRUE,stringsAsFactors=FALSE)

# replace yes in EPD_DISTRIBUTE and TARGETED_TRAP with 1
#epdDis=gsub('yes',as.numeric(1),designII$EPD_DISTRIBUTE)
#tt=gsub('yes',as.numeric(1),designII$TARGETED_TRAP)

# find index and order of EP Plates
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

# we average the tt, epd and wells to find efficiency per plate
effs=function(x,y,z){if (z==0){return (rep(0,length(z)))} else {return(((x+y)/z))}}
eff=mapply(effs,tt,epd,epdWells)

# na remover from efficiency list
for (i in 1:length(eff)){eff[[i]][which(unlist(lapply(eff[[i]],yar))=='NaN')]=0}

# averaging eff
avg=lapply(eff,ave)
effPlate=unlist(lapply(avg,function(x)return(tail(x,1))))

# get plates with rubbish eff (<0.1)
plates=data.frame(epPlates=unique(geneClub$EP_PLATE_NAME),targetingEfficiencies=effPlate)
rubbishPlates=plates$epPlates[which(plates$targetingEfficiencies<0.1)]

# finding unique genes per plate to make the geneClub file
uniqGenes=function(x){return(unique(geneClub$MARKER_SYMBOL[x]))}
uniqGenesperPlate=lapply(geneClubEPIndex,uniqGenes)
uniqID=function(x){return(unique(geneClub$DESIGN_ID[x]))}
uniqIDperPlate=lapply(geneClubEPIndex,uniqID)
sameGeneMultiID=which(unlist(lapply(uniqIDperPlate,length))-unlist(lapply(uniqGenesperPlate,length))!=0)

# finding same genes multi ID's (325,582,628)

# we know Hnf4g has multi id's find the eff and remove the less efficient
#geneClub[geneClubEPIndex[[582]][which(geneClub$MARKER_SYMBOL[geneClubEPIndex[[582]]]=='Hnf4g')],]

#removing that from whole club
#geneClub=geneClub[-geneClubEPIndex[[582]][217:225],]

# QA - quality analysis for previous loop
length(which(geneClub$EPD_DISTRIBUTE[geneClubEPIndex[[2]][generIndex[[2]][[1]]]]=='yes'))

# storing new values in data frame 
geneClubII=data.frame(genes=unlist(uniqGenesperPlate),designID=unlist(uniqIDperPlate),epPlates=unique(geneClub$EP_PLATE_NAME),eff=unlist(eff))

# calculating and adding plate index to the new data.frame
uniqPlate=function(x){return(unique(geneClub$EP_PLATE_NAME[x]))}
uniqPlates=lapply(geneClubEPIndex,uniqPlate)
plateIndexII<-list()
for(i in 1:length(uniqGenesperPlate)){plateIndexII[[i]]=rep(uniqPlates[[i]],length=length(uniqGenesperPlate[[i]]))}
geneClubII<-data.frame(geneClubII,plates=unlist(plateIndexII))

# removing same gene, same design II on different plates by rankings based on efficiencies#
joy=geneClubII[order(geneClubII$designID),]

# finding gene index again and max efficiency
geneIndexF=function(x){which(joy$genes==x,arr.ind=TRUE)}
geneIndex=lapply(unique(joy$genes),geneIndexF)
maxEff=function(x){which.max(joy$eff[x])}
maxEffIndex=lapply(geneIndex,maxEff)

# removing rubbish Plates from all gene data (20K genes)
rubG=function(x){which(geneClubII$epPlates==x)}
rubbishGenesIndexG=lapply(rubbishPlates,rubP) # gives us rubbishGenesIndex
mediumGenes=geneClubII[-unlist(rubbishGenesIndexG),]

# best genes subset of mediumGenes 
bestGenes=mediumGenes[which(mediumGenes$eff>0.9),]

# worst genes - same gene on different epPlates performing bad
worstG=function(x){which(mediumGenes$genes==x)}
worstGeneIndex=lapply(unique(mediumGenes$genes),worstG)

# tailer function - return tail of list useful for avg etc.
tailer=function(x)return(tail(x,1))
# we also average the same genes on different EP plates having eff <0.1
worstavgEffF=function(x){ave(mediumGenes$eff[x])}
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

# getting good genes
goodGenes=joy[-unlist(rubbishGenesIndex),]

# working with repeats

# reading 5P arm repeats for 500 and 1000 bases
homo<-read.csv('workspace/ikmc_targeting_efficiencies/tables/repeats/fivep_arm_repeats.csv')

# reading all 3P and 5P repeats
homoII<-read.csv('workspace/ikmc_targeting_efficiencies/tables/repeats/homology_arm_repeats.csv')
homoII=homoII[,-c(7:8)] # removes repeat class column for better viewing

goodInter=goodGenes[which(goodGenes$designID %in% homo$Design.ID),]
homoInter=homo[which(homo$Design.ID %in% goodGenes$designID),]

# getting five prime 500 and 1000 repeat numbers
fiveP5=(homoInter$Repeats.500/500)
fiveP10=(homoInter$Repeats.1000/1000)

# plotting repeats
plot(fiveP5,goodInter$eff,xlab="Five Prime Repeats",ylab="Targeting Efficiencies")
plot(fiveP10,goodInter$eff,xlab="Five Prime Repeats-1000",ylab="Targeting Efficiencies")

# another gene Indexing
worstG=function(x){which(goodInter$genes==x)}
worstGeneIndex=lapply(unique(goodGenes$genes),worstG)

