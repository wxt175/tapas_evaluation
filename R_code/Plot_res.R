### Plot multiple results 
### 10/21/2024
### Wen

rm(list=ls())
setwd('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/')
gc()

source('R_code/Utility_func.R')
library(ggplot2)
library(ggsci)
library(ggplot2)
library(gridExtra)


### Read in the result data
downs_Dis<-readRDS('result/pbmc_res_data/GPT5/Downsampling_25res_gpt5.RDS')
Ran_Dis<-readRDS('result/pbmc_res_data/GPT5/Mix_package_25res_gpt5.RDS')
WithinTissue_Dis<-readRDS('result/pbmc_res_data/GPT5/Mix_Within_NK_25res_gpt5.RDS')
WithoutTissue_Dis<-readRDS('result/pbmc_res_data/GPT5/Mix_outside_fat_25res_gpt5_notissue.RDS')

###Get the p-value for plot
shape<-1.9278#Lung:3.10023#Kidney: 3.008188 #Heart:2.582543
rate<-0.0960#Lung:0.161203#Kidney: 0.1486833 #Heart:0.1354775

downs_Dis_p<-get_p_vector(downs_Dis,shape=shape,rate=rate)
Ran_Dis_p<-get_p_vector(Ran_Dis,shape=shape,rate=rate)
WithinTissue_Dis_p<-get_p_vector(WithinTissue4_Dis,shape=shape,rate=rate)
WithoutTissue_Dis_p<-get_p_vector(WithoutTissue3_Dis,shape=shape,rate=rate)



###Set the color 
observable_colors <- scales::hue_pal()(10)

# Assign colors to variables
c1 <- observable_colors[4]
c2 <- observable_colors[2]
c3 <- observable_colors[8]
c4 <- observable_colors[1]
c5 <- observable_colors[7]


### Plot
prop<-c(1:5)
###
pdf('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/manu_fig/fig4_Lung_res1.pdf',width=3.6,height=3.2)
par(mar=c(2,2,0.5,0.5))
plot(x=prop, y=downs_Dis_p[1:5], type = "b",xlim=c(1,5),ylim=c(0,1) ,frame = FALSE, pch = 19, 
     col = c5, xlab = "", ylab = " ",lty=1,xaxt="n",yaxt="n",cex=2,lwd=2,main='')
axis(1, at=c(1:5),labels=c('100%','80%','60%','40%','20%'), las=0,pos = 0,lwd=2,cex.axis=1.3)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c('0','0.2','0.4','0.6','0.8','1.0'),lwd=2,cex.axis=1.3)
lines(x=prop, y=downs_Dis_p[6:10], pch = 5, col = c1, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=downs_Dis_p[11:15], pch = 11, col =c4, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=downs_Dis_p[16:20], pch = 2, col = c3, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=downs_Dis_p[21:25], pch = 8, col = c2, type = "b", lty = 1,cex=1.5,lwd=2)

dev.off()

pdf('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/manu_fig/fig4_Lung_res2.pdf',width=3.6,height=3.2)
par(mar=c(2,2,0.5,0.5))
plot(x=prop, y=Ran_Dis_p[1:5], type = "b",xlim=c(1,5),ylim=c(0,1) ,frame = FALSE, pch = 19, 
     col = c5, xlab = "", ylab = "",lty=1,xaxt="n",yaxt="n",cex=1.5,lwd=2,main='')
axis(1, at=c(1:5),labels=c('100%','80%','60%','40%','20%'), las=0,pos = 0,lwd=2,cex.axis=1.3)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c('0','0.2','0.4','0.6','0.8','1.0'),lwd=2,cex.axis=1.3)
lines(x=prop, y=Ran_Dis_p[6:10], pch = 5, col = c1, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=Ran_Dis_p[11:15], pch = 11, col =c4, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=Ran_Dis_p[16:20], pch = 2, col = c3, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=Ran_Dis_p[21:25], pch = 8, col = c2, type = "b", lty = 1,cex=1.5,lwd=2)

dev.off()


pdf('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/manu_fig/fig4_Lung_res3.pdf',width=3.6,height=3.2)
par(mar=c(2,2,0.5,0.5))
plot(x=prop, y=WithinTissue_Dis_p[1:5], type = "b",xlim=c(1,5),ylim=c(0,1) ,frame = FALSE, pch = 19, 
     col = c5, xlab = "", ylab = "Score",lty=1,xaxt="n",yaxt="n",cex=1.5,lwd=2,main='')
axis(1, at=c(1:5),labels=c('100%','80%','60%','40%','20%'), las=0,pos = 0,lwd=2,cex.axis=1.3)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c('0','0.2','0.4','0.6','0.8','1.0'),lwd=2,cex.axis=1.3)
lines(x=prop, y=WithinTissue_Dis_p[6:10], pch = 5, col = c1, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=WithinTissue_Dis_p[11:15], pch = 11, col =c4, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=WithinTissue_Dis_p[16:20], pch = 2, col = c3, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=WithinTissue_Dis_p[21:25], pch = 8, col = c2, type = "b", lty = 1,cex=1.5,lwd=2)

dev.off()


pdf('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/manu_fig/fig4_Heart_res4.pdf',width=3.6,height=3.2)
par(mar=c(2,2,0.5,0.5))
plot(x=prop, y=WithoutTissue_Dis_p[1:5], type = "b",xlim=c(1,5),ylim=c(0,1) ,frame = FALSE, pch = 19, 
     col = c5, xlab = "", ylab = "Score",lty=1,xaxt="n",yaxt="n",cex=1.5,lwd=2,main='')
axis(1, at=c(1:5),labels=c('100%','80%','60%','40%','20%'), las=0,pos = 0,lwd=2,cex.axis=1.3)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c('0','0.2','0.4','0.6','0.8','1.0'),lwd=2,cex.axis=1.3)
lines(x=prop, y=WithoutTissue_Dis_p[6:10], pch = 5, col = c1, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=WithoutTissue_Dis_p[11:15], pch = 11, col =c4, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=WithoutTissue_Dis_p[16:20], pch = 2, col = c3, type = "b", lty = 1,cex=1.5,lwd=2)
lines(x=prop, y=WithoutTissue_Dis_p[21:25], pch = 8, col = c2, type = "b", lty = 1,cex=1.5,lwd=2)
dev.off()


#legend
pdf('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/manu_fig/fig4_Heart_legend.pdf',width=3,height=3.2)
par(mar=c(0,0,0,0))
plot(x=prop[1], y=0.2, type = "b",xlim=c(1,5),ylim=c(0,1) ,frame = FALSE, pch = 19, 
     col = "white", xlab = "", ylab = "Score",lty=1,xaxt="n",yaxt="n",cex=1.5,lwd=2,main='')
legend(0.7,1,c('Fat Cells','B Cells','Muscle Cells','Endothelial Cells','Fibroblasts'),col=c(c5,c1,c4,c3,c2), pch=c(19,5,11,2,8),pt.cex=2,pt.lwd = 1.5,
       y.intersp=1.2,bty = 'n',x.intersp = 0.6,text.width = 0.3,horiz = F,cex=1.7)

dev.off()
