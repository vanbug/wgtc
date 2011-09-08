geneClub=read.csv('workspace/ikmc_targeting_efficiencies/tables/geneClubs.csv',sep='\t')

# find index and order of EP Plates
geneClub=geneClub[order(geneClub$EP_PLATE_NAME),]
rap=function(x){return(which(geneClub$EP_PLATE_NAME==x,arr.ind=TRUE))}
geneClubEPIndex=lapply(unique(geneClub$EP_PLATE_NAME),rap)

# find index and order of unique genes in an EP plate
gener=function(x){return(which(geneClub$MARKER_SYMBOL==x,arr.ind=TRUE))}
generII=function(x){return(which(geneClub$MARKER_SYMBOL[x]==unique(geneClub$MARKER_SYMBOL[x])))}
###########
# To find efficiencies, we need to have electroporation numbers
# calculating epd_wells, epd_distributes and targeted traps from the raw data
epd=function(x){return(which(geneClub$EPD_DISTRIBUTE[x]=='yes'))}
tt=function(x){return(which(geneClub$TARGETED_TRAP[x]=='yes'))}
epdWell=function(x){return(which(!geneClub$EPD_WELL_NAME[x]==''))}

# EPD_DISTRIBUTE
epdIndex=lapply(geneClubEPIndex,epd)
epdDis=unlist(lapply(epdIndex,length))

# TARGETED_TRAP
ttIndex=lapply(geneClubEPIndex,tt)
tTrap=unlist(lapply(ttIndex,length))

# EPD_WELL_NAME (Electroporation Wells)
epdwell=lapply(geneClubEPIndex,epdWell)
epdWells=unlist(lapply(epdwell,length))
###########
plates=data.frame(uniqD)

# calculate efficiencies and make data.frame, easy and add them per plate later
effCalc=function(x){return((x$EPD_DISTRIBUTE+x$TARGETED_TRAP)/x$EPD_WELL_NAME)}
eff=effCalc(geneClub)
eff[which(eff=="NaN")]=0
geneClub=data.frame(geneClub,efficiencies=eff)

# averaging all efficiencies and finding efficiencies per plate
epEff=function(x){return(sum(geneClub$efficiencies[x])/length(x))}
epEffs=unlist(lapply(geneClubEPIndex,epEff))

# ep plate efficiency filter - 917 good plates ones
#goodEPeff=unlist(epEffs[which(unlist(epEffs)>0.1)])

# filtering rubbishPlates and rubbishGenes (4421)
rubbishPlates=which(unlist(epEffs)<0.1) 
rubbishGenes=function(x){return(geneClubEPIndex[x])}
rubGenes=geneClub$MARKER_SYMBOL[unlist(rubbishGenes(rubbishPlates))]
#rubGeneseff=geneClub$efficiencies[unlist(rubbishGenes(rubbishPlates))]

# finding rubbish gene Index in original geneClub data
identifyGenes=function(x){return(which(geneClub$MARKER_SYMBOL==x))}
rubIndex=unlist(lapply(rubGenes,identifyGenes))

# removing rubbish genes and getting good genes
goodGenes=geneClub[-rubIndex,]
goodGenesEff=effCalc(goodGenes)
goodGenes=data.frame(goodGenes,efficiencies=goodGenesEff)


###############################################################################


