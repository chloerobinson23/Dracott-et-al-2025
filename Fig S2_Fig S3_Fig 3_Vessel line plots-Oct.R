# Foraging and non-foraging plots against vessel info (type and speed)
# 18 Oct 2024
# Author: Dr. Chloe V. Robinson

library(ggplot2)# for creating graphs
library(gridExtra)#for exporting figure panels
library(ggpubr)#for annotating plots
library(Hmisc)#for error bars

#read in AIS data

hp <- read.csv(("Feeding_DPM_H.csv"),stringsAsFactors = FALSE)
attach(hp)
str(hp)

# Basic scatter plot for sum of total number of tugs by month

plot2 <- ggplot()+
  geom_line(data=hp,aes(x=Month, y=Tugs_2020, color="black"))+
  geom_line(data=hp,aes(x=Month, y=Tugs_2021, color="blue"))+
  geom_line(data=hp,aes(x=Month, y=Tugs_2022, color="red"))+
  scale_color_identity(guide = "legend", name="Year",labels = c("2020", "2021","2022"))+
  ylab("Total number of tugs")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=12, angle = 45, hjust = 1),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(face = "bold", size=12),
    axis.title.y = element_text(face="bold", size=12),
    legend.position = "right",
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size=12))+
    scale_x_continuous(name="Month",limits=c(1, 12), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12),
                       labels=c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot2


#save plot
ggsave("FigureS2_Tugs_monthly.tiff", plot2, dpi=400,width = 7, height = 5)



# Basic scatter plot for sum of total number of vessels going fast by month

plot3 <- ggplot()+
  geom_line(data=hp,aes(x=Month, y=Speed_Fast_2020, color="black"))+
  geom_line(data=hp,aes(x=Month, y=Speed_Fast_2021, color="blue"))+
  geom_line(data=hp,aes(x=Month, y=Speed_Fast_2022, color="red"))+
  scale_color_identity(guide = "legend", name="Year",labels = c("2020", "2021","2022"))+
  ylab("Total number of fast vessels")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=12, angle = 45, hjust = 1),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(face = "bold", size=12),
    axis.title.y = element_text(face="bold", size=12),
    legend.position = "right",
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size=12))+
  scale_x_continuous(name="Month",limits=c(1, 12), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12),
                     labels=c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

plot3


#save plot
ggsave("Figure3_fast_vessels_monthly.tiff", plot3, dpi=400, width = 7, height = 5)



#read in visual data

hp2 <- read.csv(("F_NF_CPOD_vessels.csv"),stringsAsFactors = FALSE)
attach(hp2)
str(hp2)

hp2b <- read.csv(("F_NF_CPOD_vessels_v2.csv"),stringsAsFactors = FALSE)
attach(hp2b)
str(hp2b)

# Basic scatter plot sum of total number of vessels by feeding DPM

hp2b

mydata<-data.frame(
  month = c(1,5,6,7,8,9,10,11,12),
  DPM_F = c(8,24,3,8,6,8,69,0,3),
  DPM_NF = c(11,49,14,16,17,74,142,0,8),
  vessels = c(55,32,249,745,852,84,450,51,81),
  DPM_F_SD = c(4.6,13.9,1.7,3.1,3.6,4.6,34.0,0,1.7),
  DPM_NF_SD = c(6.4,28.3,5.7,5.0,9.8,42.7,76.8,0,4.6),
  vessels_SD= c(31.8,18.5,102.8,266.6,149.9,48.5,104.2,29.4,46.8)
)

mydata


dpmColour<- "blue"
dpmColour2<- "lightblue"
vesselColour<- "black"

plot5 <- ggplot(mydata, aes(x=month))+
  geom_line(data=mydata,aes(x=month, y=DPM_F, color=dpmColour))+
  geom_line(data=mydata,aes(x=month, y=DPM_NF, color=dpmColour2))+
  geom_line(data=mydata,aes(x=month, y=vessels), color=vesselColour)+
  geom_errorbar(aes(ymin=DPM_F-DPM_F_SD, ymax=DPM_F+DPM_F_SD), width=.2)+
  geom_errorbar(aes(ymin=DPM_NF-DPM_NF_SD, ymax=DPM_NF+DPM_NF_SD), width=.2)+
  geom_errorbar(aes(ymin=vessels-vessels_SD, ymax=vessels+vessels_SD), width=.2)+
  scale_color_identity(guide = "legend", name="Data Type",labels = c("Foraging DPM/hr", "Non-Foraging DPM/hr","Number of Vessels"))+
  labs (x="Month", y = "Click Trains/Min")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=14, angle = 45, hjust = 1),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(face = "bold", size=14),
    axis.title.y = element_text(face="bold", size=14),
    axis.title.y.right = element_text(face="bold", size=14, color="black", vjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 14, colour = "black"),
    legend.title = element_text(size=14, face = "bold"))+
  scale_x_continuous(name="Month",limits=c(5, 12), breaks = c(5,6,7,8,9,10,11,12),
                     labels=c("May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))+
  scale_y_continuous(sec.axis = sec_axis(~.*5, name="Total number of vessels"))


plot5



#save plot
ggsave("FigureS3_vessel total_sum_F_CL_v2.tiff", plot5, dpi=400, width = 10, height = 5)



################################################
#### END #######################################
###############################################


