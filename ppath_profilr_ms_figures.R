rm(list=setdiff(ls(), "filepath"))
#####################################################################################################
# Underlying data (Dubern 1994 Tropical Science) is reproduced in the                               #
# This is non-EpiPv pacakge code used to produce figures in Donnelly Tankam Gilligan 2024 (figs 3-4)# 
# obtained with the present files (AnalyseDubern_Figs.R, model_1.stan, model_2.stan) are presented  #
#####################################################################################################

library(rstudioapi)
library(ggplot2)
library(ggpmisc)
library(ggtext)
library(ggpubr)
library(bayesplot)
library(grid)
library(gridExtra)

setwd("C:/Users/ruair/Desktop/EpiPv")
#setwd("C:/Users/donnelr3/OneDrive - University Of Cambridge/Archive_ZG_WaveCassava/rFiles/rStanFiles/RPanalysis")

### IMPORT CHAIN INFERENCE DATA FOR FIGURE 3 MAIN TEXT
# here's one i made earlier, P EPIDEMICS chains
virus_type='SPT'
SPT_chaindata_table_Pl=as.matrix(read.table((paste('path_',virus_type,'_chains_byF_frPlant.dat'))))
SPT_chaindata_table_Ins=as.matrix(read.table((paste('path_',virus_type,'_chains_byF_frInsect.dat'))))
SPTsz_chns=dim(SPT_chaindata_table_Pl)
SPTburdens=rownames(SPT_chaindata_table_Pl)

virus_type='PT'
PT_chaindata_table_Pl=as.matrix(read.table((paste('path_',virus_type,'_chains_byF_frPlant.dat'))))
PT_chaindata_table_Ins=as.matrix(read.table((paste('path_',virus_type,'_chains_byF_frInsect.dat'))))
PTsz_chns=dim(PT_chaindata_table_Pl)
PTburdens=rownames(PT_chaindata_table_Pl)




# FIGURE 3 # FIGURE 3# FIGURE 3# FIGURE 3# FIGURE 3# FIGURE 3# FIGURE 3# FIGURE 3
################################################################################
###############################################################################
################################################################################
# subplot by subplot
fig_insBurden_indx=1;
ins1PR_frPl=data.frame(PT_chaindata_table_Pl[fig_insBurden_indx,1:PTsz_chns[2]])
ins1PR_frIns=data.frame(PT_chaindata_table_Ins[fig_insBurden_indx,1:(PTsz_chns[2])])
names(ins1PR_frPl)=c('length')
names(ins1PR_frIns)=c('length')
ins1PR_frPl$inoc='plant inoculum'
ins1PR_frIns$inoc='insect inoculum'
inocTypes_1insect_PT=rbind(ins1PR_frPl,ins1PR_frIns)

cmb1=ggplot(inocTypes_1insect_PT, aes(length, fill = inoc)) + geom_density(alpha = 0.2)+
  labs(
    title = "             CMB",
    subtitle = "F=1",
    caption = " ",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none") + annotate(geom="text", x=0.7, y=40, label="infected insect intro",
                                              color="red")+ annotate(geom="text", x=0.45, y=27, label="infected plant intro",
                                                                     color="Dark Cyan")+ annotate(geom="text", x=0.075, y=40-3, label="A",
                                                                                                   color="black", size = unit(5, "pt"))+
  theme(plot.margin = unit(c(0,2.5,0,0), "lines"))
##############################################################################

fig_insBurden_indx=2;
ins2PR_frPl=data.frame(PT_chaindata_table_Pl[fig_insBurden_indx,1:PTsz_chns[2]])
ins2PR_frIns=data.frame(PT_chaindata_table_Ins[fig_insBurden_indx,1:(PTsz_chns[2])])
names(ins2PR_frPl)=c('length')
names(ins2PR_frIns)=c('length')
ins2PR_frPl$inoc='plant inoculum'
ins2PR_frIns$inoc='insect inoculum'
inocTypes_2insect_PT=rbind(ins2PR_frPl,ins2PR_frIns)

cmb2=ggplot(inocTypes_2insect_PT, aes(length, fill = inoc)) + geom_density(alpha = 0.2)+
  labs(
    title = "                ",
    subtitle = "F=2",
    caption = "Source: Dubern '94 AP dataset",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none")+theme(plot.margin = unit(c(0,2.5,0,0), "lines"))+ annotate(geom="text", x=0.075, y=147, label="C",
                                                                                               color="black", size = unit(5, "pt")) 
################################################################################


##############################################################################
fig_insBurden_indx=1;
ins1SPR_frPl=data.frame(SPT_chaindata_table_Pl[fig_insBurden_indx,1:SPTsz_chns[2]])
ins1SPR_frIns=data.frame(SPT_chaindata_table_Ins[fig_insBurden_indx,1:(SPTsz_chns[2])])
names(ins1SPR_frPl)=c('length')
names(ins1SPR_frIns)=c('length')
ins1SPR_frPl$inoc='plant inoculum'
ins1SPR_frIns$inoc='insect inoculum'
inocTypes_1insect=rbind(ins1SPR_frPl,ins1SPR_frIns)

cbsi1=ggplot(inocTypes_1insect, aes(length, fill = inoc)) + geom_density(alpha = 0.2)+
  labs(
    title = "            CBSI",
    subtitle = "F=4",
    caption = " ",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none")+theme(plot.margin = unit(c(0,2.5,0,0), "lines"))+ annotate(geom="text", x=0.35, y=20, label="infected insect intro",
                                                                                                color="red")+ annotate(geom="text", x=0.4, y=5, label="infected plant intro",
                                                                                                                       color="Dark Cyan")+annotate(geom="text", x=0.075, y=35-3, label="B",
                                                                                              color="black", size = unit(5, "pt"))
##############################################################################


fig_insBurden_indx=3;
ins2SPR_frPl=data.frame(SPT_chaindata_table_Pl[fig_insBurden_indx,1:SPTsz_chns[2]])
ins2SPR_frIns=data.frame(SPT_chaindata_table_Ins[fig_insBurden_indx,1:(SPTsz_chns[2])])
names(ins2SPR_frPl)=c('length')
names(ins2SPR_frIns)=c('length')
ins2SPR_frPl$inoc='plant inoculum'
ins2SPR_frIns$inoc='insect inoculum'
inocTypes_2insect=rbind(ins2SPR_frPl,ins2SPR_frIns)

cbsi2=ggplot(inocTypes_2insect, aes(length, fill = inoc)) + geom_density(alpha = 0.2)+
  labs(
    title = "                ",
     subtitle = "F=8",
    caption = "Source: Maruthi et al '17 AP dataset",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none")+theme(plot.margin = unit(c(0,2.5,0,0), "lines"))+ annotate(geom="text", x=0.075, y=21-3, label="D",
                                                                                              color="black", size = unit(5, "pt")) 
################################################################################
################################################################################
# Assemble subplots and save
pdf("Figure_3.pdf", width = 6.5, height = 4.5) # Open a new pdf file
grid.arrange(cmb1,cbsi1,cmb2,cbsi2, ncol=2, nrow =2)
dev.off() # Close the file
################################################################################
################################################################################
################################################################################
################################################################################




# FIGURE 4 # FIGURE 4# FIGURE 4# FIGURE 4# FIGURE 4# FIGURE 4# FIGURE 4# FIGURE 4
################################################################################
################################################################################
################################################################################
# subplot by subplot
fig_insBurden_indx=1;
ins1PR_frPl=data.frame(PT_chaindata_table_Pl[1,1:PTsz_chns[2]])
ins1PR_frPl$density=1
ins2PR_frPl=data.frame(PT_chaindata_table_Pl[2,1:PTsz_chns[2]])
ins2PR_frPl$density=2
ins3PR_frPl=data.frame(PT_chaindata_table_Pl[3,1:PTsz_chns[2]])
ins3PR_frPl$density=3

names(ins1PR_frPl)=c('length','density')
names(ins2PR_frPl)=c('length','density')
names(ins3PR_frPl)=c('length','density')

inocTypes_1insect_PT=rbind(ins1PR_frPl,ins2PR_frPl,ins3PR_frPl)
inocTypes_1insect_PT$density=as.factor(inocTypes_1insect_PT$density)

cmb1=ggplot(inocTypes_1insect_PT, aes(length, fill = density)) + geom_density(alpha = 0.2)+
  labs(
    title = "             CMB",
    subtitle = "Plant inoculum",
    caption = " ",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none")+
  theme(plot.margin = unit(c(0,2.5,0,0), "lines"))

g <- ggplot_build(cmb1)
gu=unique(g$data[[1]]["fill"])

cmbA=cmb1+annotate(geom="text", x=0.625-0.075, y=13, label="F=1",
              color=gu$fill[1])+annotate(geom="text", x=0.825-0.075, y=30, label="F=2",
                                         color=gu$fill[2])+annotate(geom="text", x=0.825-0.075, y=42, label="F=3",
                                                                    color=gu$fill[3])+ annotate(geom="text", x=0.075, y=50-2, label="A",
                                                                                                color="black", size = unit(5, "pt")) 
##########################################################################
fig_insBurden=1;
ins1PR_frIns=data.frame(PT_chaindata_table_Ins[1,1:PTsz_chns[2]])
ins1PR_frIns$density=1
ins2PR_frIns=data.frame(PT_chaindata_table_Ins[2,1:PTsz_chns[2]])
ins2PR_frIns$density=2
ins3PR_frIns=data.frame(PT_chaindata_table_Ins[3,1:PTsz_chns[2]])
ins3PR_frIns$density=3
names(ins1PR_frIns)=c('length','density')
names(ins2PR_frIns)=c('length','density')
names(ins3PR_frIns)=c('length','density')

inocTypes_1insect_PT=rbind(ins1PR_frIns,ins2PR_frIns,ins3PR_frIns)
inocTypes_1insect_PT$density=as.factor(inocTypes_1insect_PT$density)

cmb2=ggplot(inocTypes_1insect_PT, aes(length, fill = density)) + geom_density(alpha = 0.2)+
  labs(
    title = "                ",
    subtitle = "Insect inoculum",
    caption = " ",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none")+
  theme(plot.margin = unit(c(0,2.5,0,0), "lines"))

g <- ggplot_build(cmb2)
gu=unique(g$data[[1]]["fill"])

cmbB=cmb2+ annotate(geom="text", x=0.075, y=147, label="C",
                                                 color="black", size = unit(5, "pt"))
                    #+annotate(geom="text", x=0.675, y=2.25, label="F=1",
           #        color=gu$fill[1])+annotate(geom="text", x=0.825, y=4, label="F=3",
              #                                color=gu$fill[2])+annotate(geom="text", x=0.875, y=7, label="F=5",
                  #                                                       color=gu$fill[3])+ annotate(geom="text", x=0.075, y=8-1, label="C",
                                                                       #                              color="black", size = unit(5, "pt")) 
################### CBSI Plant inoc ############################################
fig_insBurden_indx=1;
ins1SPR_frPl=data.frame(SPT_chaindata_table_Pl[1,1:SPTsz_chns[2]])
ins1SPR_frPl$density=1
ins2SPR_frPl=data.frame(SPT_chaindata_table_Pl[2,1:SPTsz_chns[2]])
ins2SPR_frPl$density=2
ins3SPR_frPl=data.frame(SPT_chaindata_table_Pl[3,1:SPTsz_chns[2]])
ins3SPR_frPl$density=3

names(ins1SPR_frPl)=c('length','density')
names(ins2SPR_frPl)=c('length','density')
names(ins3SPR_frPl)=c('length','density')

inocTypes_1insect_SPT=rbind(ins1SPR_frPl,ins2SPR_frPl,ins3SPR_frPl)
inocTypes_1insect_SPT$density=as.factor(inocTypes_1insect_SPT$density)

cbsi1=ggplot(inocTypes_1insect_SPT, aes(length, fill = density)) + geom_density(alpha = 0.2)+
  labs(
    title = "             CBSI",
    subtitle = "Plant inoculum",
    caption = " ",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  ylim(0, 15)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none")+
  theme(plot.margin = unit(c(0,2.5,0,0), "lines"))

g <- ggplot_build(cbsi1)
gu=unique(g$data[[1]]["fill"])

cbsiA=cbsi1+annotate(geom="text", x=0.05, y=4.25, label="F=4",
                   color=gu$fill[1])+annotate(geom="text", x=0.5-0.15, y=4.25, label="F=6",
                                              color=gu$fill[2])+annotate(geom="text", x=0.75, y=4, label="F=8",
                                                                         color=gu$fill[3])+annotate(geom="text", x=0.075, y=14.5, label="B",
                                                                                                     color="black", size = unit(5, "pt")) 
################### CBSI Insect inoc ###########################################
fig_insBurden_indx=1;
ins1SPR_frIns=data.frame(SPT_chaindata_table_Ins[1,1:SPTsz_chns[2]])
ins1SPR_frIns$density=1
ins2SPR_frIns=data.frame(SPT_chaindata_table_Ins[2,1:SPTsz_chns[2]])
ins2SPR_frIns$density=2
ins3SPR_frIns=data.frame(SPT_chaindata_table_Ins[3,1:SPTsz_chns[2]])
ins3SPR_frIns$density=3

names(ins1SPR_frIns)=c('length','density')
names(ins2SPR_frIns)=c('length','density')
names(ins3SPR_frIns)=c('length','density')

inocTypes_1insect_SPT=rbind(ins1SPR_frIns,ins2SPR_frIns,ins3SPR_frIns)
inocTypes_1insect_SPT$density=as.factor(inocTypes_1insect_SPT$density)

cbsi2=ggplot(inocTypes_1insect_SPT, aes(length, fill = density)) + geom_density(alpha = 0.2)+
  labs(
    title = "                 ",
    subtitle = "Insect inoculum",
    caption = " ",
    x = "Epidemic probability",
    y = "Density"
  )+
  xlim(0, 1)+
  theme_classic() +
  theme(
    plot.title = element_text(color = "#0099F8", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "#00ac42", size = 14, face = "bold"),
    plot.caption = element_text(face = "italic")
  )+ theme(legend.position="none")+
  theme(plot.margin = unit(c(0,2.5,0,0), "lines"))

g <- ggplot_build(cbsi2)
gu=unique(g$data[[1]]["fill"])

cbsiB=cbsi2+annotate(geom="text", x=0.075, y=31, label="D",
                        color="black", size = unit(5, "pt"))
#+annotate(geom="text", x=0.075, y=30, label="F=4",
                   #  color=gu$fill[1])+annotate(geom="text", x=0.1, y=22, label="F=6",
                                    #            color=gu$fill[2])+annotate(geom="text", x=0.125, y=15, label="F=8",
                                                               #            color=gu$fill[3])+ annotate(geom="text", x=0.075, y=34.5, label="D",
                                                                                                    #   color="black", size = unit(5, "pt")) 


cbsi1+annotate(geom="text", x=0.05, y=3, label="F=1",
               color=gu$fill[1])+annotate(geom="text", x=0.5, y=4, label="F=2",
                                          color=gu$fill[2])+annotate(geom="text", x=0.75, y=5, label="F=3",
                                                                     color=gu$fill[3])+annotate(geom="text", x=0.075, y=14.75, label="B",
                                                                                                color="black", size = unit(5, "pt"))
################################################################################
################################################################################
# Assemble subplots and save
pdf("Figure_4.pdf", width = 6.5, height = 4.5) # Open a new pdf file
grid.arrange(cmbA,cbsiA,cmbB,cbsiB, ncol=2, nrow =2)
dev.off() # Close the file
################################################################################
################################################################################
################################################################################
