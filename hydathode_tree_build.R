library(ape)
library(geiger)
library(phytools)
library(dplyr)

##Part One: Trimming Testo & Sundue (2016) tree to genus-level sampling

setwd("~/Desktop/Projects/Hydathodes/")

tree<-read.nexus("4000taxon.tre")

tipsToDrop <- read.table(file = "tipsToDrop.txt", header = T, row.names = 1)

name.check(tree,tipsToDrop)

sampledTree<-drop.tip(tree,rownames(tipsToDrop))

write.tree(sampledTree, "generaTree_2.tre") 

#then manually update names -- change binomials to genus names and update taxonomy


##Part Two: Trimming genus-level tree to remove all genera with unknown hydathode occurrence

tree<-read.tree("generaTree.tre")

scoredTraits<-read.csv("hydathodeScoring.csv")

unknownGenera <- scoredTraits %>% filter(Scoring == "unknown")

generaToDrop <- unknownGenera$Genus %>% as.vector()

reducedTree <- drop.tip(tree, generaToDrop)

write.tree(reducedTree, "treeReduced.tre")



##Part Three: Generating plots of tip state distribution and ACR

tree<-read.tree("generaTree.tre")

X<-read.csv("hydathodeScoring.csv", row.names = 1)

state<-X[tree$tip.label,]

cols<-setNames(c("#377eb8", "#e41a1c", "#ff7f00", "#ffff33", "#000000"), levels(state))
legendNames<-c("Absent", "All Present", "Enlarged Vein Endings", "Present in Some Species", "Unknown")

pdf(file = "Figures/HydathodeTipMap.pdf", width = 8, height = 8)

plotTree(tree,type="fan",fsize=0.3,ftype="i",lwd=0.5, offset = 3)

tiplabels(pie=to.matrix(state,sort(unique(state))),piecol=cols,cex=0.1)

add.simmap.legend(leg= legendNames, colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(0.9*nodeHeights(tree)),fsize=0.3, shape = "circle")
dev.off()



### ACR for trimmed tree (four states)
treeReduce<-read.tree("treeReduce.tre")

Y<-read.csv("hydathodeScoringReduced.csv", row.names = 1)

stateReduce<-Y[treeReduce$tip.label,]

cols2<-setNames(c("#377eb8", "#e41a1c", "#ff7f00", "#ffff33"), levels(stateReduce))

legendNamesReduced<-c("Absent", "All Present", "Enlarged Vein Endings", "Present in Some Species")

fitER_Reduce<-ace(stateReduce,treeReduce, type="discrete", method = "ML")

pdf(file = "Figures/Hydathode_ACR", width = 8, height = 8)

plotTree(treeReduce,type="fan",fsize=0.3,ftype="i",lwd=0.5, offset = 3)

tiplabels(pie=to.matrix(stateReduce,sort(unique(stateReduce))),piecol=cols2,cex=0.15)

nodelabels(node=1:treeReduce$Nnode+Ntip(treeReduce),
           pie=fitER_Reduce$lik.anc,piecol=cols2,cex=0.3)

add.simmap.legend(leg=legendNamesReduced, colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(0.9*nodeHeights(treeReduce)),fsize=0.5, shape = "circle")
dev.off()


### ACR for trimmed tree (absent/present)
treeReduce<-read.tree("treeReduce.tre")

Y<-read.csv("hydathodeScoringReduced.csv", row.names = 1)

stateReduce<-Y[treeReduce$tip.label,]

stateReduceVector <- as.vector(stateReduce) #transform states to vector

stateReduceVector <- gsub("sp", "ap", stateReduceVector) ##recode presence 

stateReduceVector <- gsub("ev", "ap", stateReduceVector) ##recode presence

names(stateReduceVector) <- treeReduce$tip.label 

levels<-levels(as.factor(stateReduceVector))

cols2<-setNames(c("#377eb8", "#e41a1c"), levels)

legendNamesReduced<-c("Absent", "Present")

fitER_Reduce_Binary <- ace(stateReduce,treeReduce, type="discrete", method = "ML")

pdf(file = "Figures/Hydathode_ACR_Binary", width = 8, height = 8)

plotTree(treeReduce,type="fan",fsize=0.3,ftype="i",lwd=0.5, offset = 3)

tiplabels(pie=to.matrix(stateReduce,sort(unique(stateReduce))),piecol=cols2,cex=0.15)

nodelabels(node=1:treeReduce$Nnode+Ntip(treeReduce),
           pie=fitER_Reduce_Binary$lik.anc,piecol=cols2,cex=0.3)

add.simmap.legend(leg=legendNamesReduced, colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(0.9*nodeHeights(treeReduce)),fsize=0.5, shape = "circle")
dev.off()





### SIMMAP for trimmed tree


treeReduce<-read.tree("treeReduce.tre")

Y<-read.csv("hydathodeScoringReduced.csv", row.names = 1)

stateReduce<-Y[treeReduce$tip.label,]

stateReduceVector <- as.vector(stateReduce) #transform states to vector

stateReduceVector <- gsub("sp", "ap", stateReduceVector) ##recode presence 

stateReduceVector <- gsub("ev", "ap", stateReduceVector) ##recode presence

names(stateReduceVector) <- treeReduce$tip.label 

levels<-levels(as.factor(stateReduceVector))

cols2<-setNames(c("#377eb8", "#e41a1c"), levels)

legendNamesReduced<-c("Absent", "Present")

fitER <- fitMk(treeReduce, stateReduceVector, model = "ER") ##fit model for simmap rates


Q <- matrix(NA, length(levels),length(levels)) #make simmap transition matrix
Q[1,]<-c(-0.005,0.005)
Q[2,]<-c(0.005,-0.005)
colnames(Q) <- rownames(Q) <- fitER$states

simmapReduce <- make.simmap(treeReduce, stateReduceVector, Q = Q, nsim = 500)

simmapReduceSummary <- summary(simmapReduce)

simmapReduceCount <- countSimmap(simmapReduce)


pdf(file = "Figures/Hydathode_SimMap.pdf", width = 8, height = 8)

plot(simmapReduceSummary,colors = cols2,type = "fan",
     fsize = 0.3,cex = c(0.175, 0.175), lwd= 0.25, offset = 10)

nodelabels(node=1:treeReduce$Nnode+Ntip(treeReduce),
           pie=simmapReduceSummary$ace,piecol=cols2,cex=0.3)

add.simmap.legend(leg=legendNamesReduced, colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(0.9*nodeHeights(treeReduce)),fsize=0.5, shape = "circle")

dev.off()

meanAP <- mean(simmapReduceCount$Tr[,3])

meanPA <- mean(simmapReduceCount$Tr[,4])

sdAP <- sd(simmapReduceCount$Tr[,3])

sdPA <- sd(simmapReduceCount$Tr[,4])
