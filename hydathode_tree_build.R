library(ape)
library(geiger)
library(phytools)
library(dplyr)

setwd()#set wd here)

##Part One: Trimming Testo & Sundue (2016) tree to genus-level sampling

setwd("~/Desktop/Projects/Hydathodes/")

tree<-read.nexus("4000taxon.tre")

tipsToDrop <- read.table(file = "tipsToDrop.txt", header = T, row.names = 1)

name.check(tree,tipsToDrop)

sampledTree<-drop.tip(tree,rownames(tipsToDrop))

write.tree(sampledTree, "generaTree.tre") 

#then manually update names -- change binomials to genus names and update taxonomy


##Part Two: Trimming genus-level tree to remove all genera with unknown hydathode occurrence

tree<-read.tree("generaTree.tre")

scoredTraits<-read.csv("hydathodeScoring.csv")

unknownGenera <- scoredTraits %>% filter(Scoring == "u")

generaToDrop <- unknownGenera$Genus %>% as.vector()

reducedTree <- drop.tip(tree, generaToDrop)

write.tree(reducedTree, "treeReduce.tre")

##Part Three: Generating plots of tip state distribution and ACR

tree<-read.tree("generaTree.tre")

X<-read.csv("hydathodeScoring.csv", row.names = 1)

state<-X[tree$tip.label,]

cols<-setNames(c("#377eb8", "#e41a1c", "#4daf4a", "#ffff33", "#000000"), levels(state))
colsALT<-setNames(c("#377eb8", "#e41a1c", "#ff7f00", "#ffff33", "#000000"), levels(state))
legendNames<-c("Absent", "All Present", "Enlarged Vein Endings", "Present in Some Species", "Unknown")

pdf(file = "Figures/HydathodeTipMap.pdf", width = 8, height = 8)

plotTree(tree,type="fan",fsize=0.3,ftype="i",lwd=0.5, offset = 3)

tiplabels(pie=to.matrix(state,sort(unique(state))),piecol=cols,cex=0.15)

add.simmap.legend(leg= legendNames, colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(0.9*nodeHeights(tree)),fsize=0.5, shape = "circle")
dev.off()


treeReduce<-read.tree("treeReduce.tre")

Y<-read.csv("hydathodeScoringReduced.csv", row.names = 1)

stateReduce<-Y[treeReduce$tip.label,]

cols2<-setNames(c("#377eb8", "#e41a1c", "#4daf4a", "#ffff33"), levels(stateReduce))
cols2ALT<-setNames(c("#377eb8", "#e41a1c", "#ff7f00", "#ffff33"), levels(stateReduce))

legendNamesReduced<-c("Absent", "All Present", "Enlarged Vein Endings", "Present in Some Species")

fitER_Reduce<-ace(stateReduce,treeReduce, type="discrete", method = "ML")

pdf(file = "Figures/Hydathode_ACR", width = 8, height = 8)

plotTree(treeReduce,type="fan",fsize=0.3,ftype="i",lwd=0.5, offset = 3)

tiplabels(pie=to.matrix(stateReduce,sort(unique(stateReduce))),piecol=cols2,cex=0.2)

nodelabels(node=1:treeReduce$Nnode+Ntip(treeReduce),
           pie=fitER_Reduce$lik.anc,piecol=cols2,cex=0.3)

add.simmap.legend(leg=legendNamesReduced, colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(0.9*nodeHeights(treeReduce)),fsize=0.5, shape = "circle")
dev.off()




