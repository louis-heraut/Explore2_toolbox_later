rm(list=ls(all=TRUE))

library(png)

data_dir = paste0("C:\\Users\\esauquet\\Desktop\\Explore2\\TRACC\\CLIMAT")
setwd(data_dir)

climatZHTot <- read.table("climatZH.csv",sep=";",header=T)
uniqueZH <- unique(climatZHTot[,1])

for (k in 2:2) { #1:length(uniqueZH)){
  
i3 <- which(climatZHTot$NivRechauf=="GWL-30"&climatZHTot$Variable=="RR"&climatZHTot$ZH==uniqueZH[k])
climatZH <- climatZHTot[i3,]

i3DJF <- which(climatZH$Saison=="seas-DJF")
climatZHDJF <- climatZH[i3DJF,]
GCMRCM <- paste(climatZHDJF$GCM,climatZHDJF$RCM,sep="_")

uniqueGCMRCM <- unique(GCMRCM)

deltaRRDJF <- c()
for (j in 1:length(uniqueGCMRCM)){
  iGR <- which(GCMRCM==uniqueGCMRCM[j])
  deltaRRDJF <- rbind(deltaRRDJF,c(uniqueGCMRCM[j],
                      weighted.mean(climatZHDJF$delta[iGR],climatZHDJF$surface[iGR])))
}

i3JJA <- which(climatZH$Saison=="seas-JJA")
climatZHJJA <- climatZH[i3JJA,]
GCMRCM <- paste(climatZHJJA$GCM,climatZHJJA$RCM,sep="_")

deltaRRJJA <- c()
for (j in 1:length(uniqueGCMRCM)){
  iGR <- which(GCMRCM==uniqueGCMRCM[j])
  deltaRRJJA <- rbind(deltaRRJJA,c(uniqueGCMRCM[j],
                                   weighted.mean(climatZHJJA$delta[iGR],
                                                 climatZHJJA$surface[iGR])))
}

#plot(deltaRRDJF[,2],deltaRRJJA[,2],col=1:17,pch=1:17)

####

i3 <- which(climatZHTot$NivRechauf=="GWL-30"&climatZHTot$Variable=="TMm"&climatZHTot$ZH==uniqueZH[k])
climatZH <- climatZHTot[i3,]

i3DJF <- which(climatZH$Saison=="seas-DJF")
climatZHDJF <- climatZH[i3DJF,]
GCMRCM <- paste(climatZHDJF$GCM,climatZHDJF$RCM,sep="_")

uniqueGCMRCM <- unique(GCMRCM)

deltaTMmDJF <- c()
for (j in 1:length(uniqueGCMRCM)){
  iGR <- which(GCMRCM==uniqueGCMRCM[j])
  deltaTMmDJF <- rbind(deltaTMmDJF,c(uniqueGCMRCM[j],
                                   weighted.mean(climatZHDJF$delta[iGR],climatZHDJF$surface[iGR])))
}

i3JJA <- which(climatZH$Saison=="seas-JJA")
climatZHJJA <- climatZH[i3JJA,]
GCMRCM <- paste(climatZHJJA$GCM,climatZHJJA$RCM,sep="_")

deltaTMmJJA <- c()
for (j in 1:length(uniqueGCMRCM)){
  iGR <- which(GCMRCM==uniqueGCMRCM[j])
  deltaTMmJJA <- rbind(deltaTMmJJA,c(uniqueGCMRCM[j],
                                   weighted.mean(climatZHJJA$delta[iGR],
                                                 climatZHJJA$surface[iGR])))
}

#plot(deltaTMmDJF[,2],deltaTMmJJA[,2],col=1:17,pch=1:17)



png(paste(uniqueZH[k],"_climatZones_Hydro.png",sep=""),width = 9, 
    height = 2.5, units = 'in', res = 600)
layout(matrix(c(1,2,3),1,3))
par(mar=c(4, 4.2, 0.2, 1),xpd=FALSE)

#9 MOHC-HadGEM2-ES ICTP-RegCM4-6
#13 MPI-M-MPI-ESM-LR
#15 NCC-NorESM1-M DMI-HIRHAM5

colPt <- rep("grey",14)
colPt[9] <- "#1b9e77"
colPt[12] <- "#d95f02"
colPt[15] <- "#7570b3"

par(fig=c(0,2.5,0,10)/10)
plot(deltaTMmDJF[,2],100*as.numeric(deltaRRDJF[,2]),bg=colPt,pch=21,col="black",
     xlab="Chang. température hiver (°C)",ylab="Chang. précipitations hiver (%)",cex=2.25,cex.axis=1.2,cex.lab=1.2)
par(fig=c(2.5,5,0,10)/10)
par(new=T)
plot(deltaTMmJJA[,2],100*as.numeric(deltaRRJJA[,2]),bg=colPt,pch=21,col="black",
     xlab="Chang. température été (°C)",ylab="Chang. précipitations été (%)",cex=2.25,cex.axis=1.2,cex.lab=1.2)
par(fig=c(5,10,0,10)/10)
par(new=T)
par(mar=c(0.5, 0.5, 0.5, 0.5),xpd=FALSE)
plot(c(0.5,6.5),c(1,10),axes=F,col=NA)

textTableau <- matrix("xxx",7,5)

textTableau[7,1] <- ""
textTableau[6,1] <- "Médiane (ensemble)"
textTableau[5,1] <- "Minimum (ensemble)"
textTableau[4,1] <- "Maximum (ensemble)"
textTableau[3,1] <- uniqueGCMRCM[9]
textTableau[2,1] <- uniqueGCMRCM[12]
textTableau[1,1] <- uniqueGCMRCM[15]

textTableau[7,2] <- "TDJF (°C)"
textTableau[7,3] <- "PDJF (%)"
textTableau[7,4] <- "TJJA (°C)"
textTableau[7,5] <- "PJJA (%)"

textTableau[6,2] <- round(median(as.numeric(deltaTMmDJF[,2])),1)
textTableau[6,3] <- round(median(100*as.numeric(deltaRRDJF[,2])),1)
textTableau[6,4] <- round(median(as.numeric(deltaTMmJJA[,2])),1)
textTableau[6,5] <- round(median(100*as.numeric(deltaRRJJA[,2])),1)

textTableau[5,2] <- round(min(as.numeric(deltaTMmDJF[,2])),1)
textTableau[5,3] <- round(min(100*as.numeric(deltaRRDJF[,2])),1)
textTableau[5,4] <- round(min(as.numeric(deltaTMmJJA[,2])),1)
textTableau[5,5] <- round(min(100*as.numeric(deltaRRJJA[,2])),1)

textTableau[4,2] <- round(max(as.numeric(deltaTMmDJF[,2])),1)
textTableau[4,3] <- round(max(100*as.numeric(deltaRRDJF[,2])),1)
textTableau[4,4] <- round(max(as.numeric(deltaTMmJJA[,2])),1)
textTableau[4,5] <- round(max(100*as.numeric(deltaRRJJA[,2])),1)

textTableau[3,2] <- round((as.numeric(deltaTMmDJF[9,2])),1)
textTableau[3,3] <- round((100*as.numeric(deltaRRDJF[9,2])),1)
textTableau[3,4] <- round((as.numeric(deltaTMmJJA[9,2])),1)
textTableau[3,5] <- round((100*as.numeric(deltaRRJJA[9,2])),1)

textTableau[2,2] <- round((as.numeric(deltaTMmDJF[12,2])),1)
textTableau[2,3] <- round((100*as.numeric(deltaRRDJF[12,2])),1)
textTableau[2,4] <- round((as.numeric(deltaTMmJJA[12,2])),1)
textTableau[2,5] <- round((100*as.numeric(deltaRRJJA[12,2])),1)

textTableau[1,2] <- round((as.numeric(deltaTMmDJF[15,2])),1)
textTableau[1,3] <- round((100*as.numeric(deltaRRDJF[15,2])),1)
textTableau[1,4] <- round((as.numeric(deltaTMmJJA[15,2])),1)
textTableau[1,5] <- round((100*as.numeric(deltaRRJJA[15,2])),1)

colText <- c(rep("black",4),"#1b9e77",
             "#d95f02","#7570b3")
for (ii in 1:1){
  for (jj in 1:7){
    text((ii-0.75),(2+jj),textTableau[jj,ii],col=colText[8-jj],pos=4,cex=1.1)
  }
}


for (ii in 2:5){
  abline(h=ii+2.5,lwd=2)
  for (jj in 1:7){
    text((ii+1),(2+jj),textTableau[jj,ii],col="black",cex=1.1)
  }
}
abline(h=6+2.5,lwd=2)
abline(h=1+2.5,lwd=2)
abline(h=0+2.5,lwd=2)

text(3.5,1.75,"Changements projetés (référence : 1976-2005)",col="black",cex=1.6)

dev.off()

png(paste(uniqueZH[k],"_hydrologieZones_Hydro.png",sep=""),width = 9, 
    height = 2, units = 'in', res = 600)
layout(matrix(c(1),1,1))
par(mar=c(4, 4.2, 0.2, 1),xpd=FALSE)

#9 MOHC-HadGEM2-ES ICTP-RegCM4-6
#13 MPI-M-MPI-ESM-LR
#15 NCC-NorESM1-M DMI-HIRHAM5

colPt <- rep("grey",14)
colPt[9] <- "#1b9e77"
colPt[12] <- "#d95f02"
colPt[15] <- "#7570b3"

par(mar=c(0.5, 0.5, 0.5, 0.5),xpd=FALSE)
plot(c(0.5,10.25),c(1.7,9),axes=F,col=NA)

textTableau <- matrix("xxx",7,9)

textTableau[7,1] <- ""
textTableau[6,1] <- "Médiane (ensemble)"
textTableau[5,1] <- "Minimum (ensemble)"
textTableau[4,1] <- "Maximum (ensemble)"
textTableau[3,1] <- uniqueGCMRCM[9]
textTableau[2,1] <- uniqueGCMRCM[12]
textTableau[1,1] <- uniqueGCMRCM[15]

textTableau[7,2] <- "VCN10(5) (%)"
textTableau[7,3] <- "QA (%)"
textTableau[7,4] <- "Recharge (%)"
textTableau[7,5] <- "QDJF (%)"
textTableau[7,6] <- "QMAM (%)"
textTableau[7,7] <- "QJJA (%)"
textTableau[7,8] <- "QSON (%)"
textTableau[7,9] <- "QJXA10 (%)"

colText <- c(rep("black",4),"#1b9e77",
             "#d95f02","#7570b3")
for (ii in 1:1){
  for (jj in 1:7){
    text((ii-1),(2+jj),textTableau[jj,ii],col=colText[8-jj],pos=4,cex=0.85)
  }
}


for (ii in 2:9){
  for (jj in 1:7){
    text((ii+1),(2+jj),textTableau[jj,ii],col="black",cex=0.85)
  }
}
abline(h=6+2.5,lwd=2)
abline(h=5+2.5,lwd=2)
abline(h=4+2.5,lwd=2)
abline(h=3+2.5,lwd=2)
abline(h=2+2.5,lwd=2)
abline(h=1+2.5,lwd=2)
abline(h=0+2.5,lwd=2)

text(5.5,1.75,"Changements médians projetés (référence : 1976-2005, statistiques spatiales)",col="black",cex=1.1)

dev.off()


}



# write.table(climatZH,"climatZH.csv",sep=";",row.name=F,
#             col.names=c("ZH","Variable","Saison","NivRechauf","GCM","RCM","delta"),quote=F)
