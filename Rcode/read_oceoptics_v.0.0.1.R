########################################################
#  Analisis de Datos Luminaria sotano                  #
#                                                      #
# Gaston Quero - Sebastian Fernandez - Mauro Moré      #
#           03-09-2020                                  #
###################################################

getwd()
# seteo el directorio  #### Aca se cambia segun donde tengas los datos 
# Esta linea se corre una sola vez en cada caso, y se 
# desabilita la linea del directorio de otro usuario
# esto lo usamos asi hasta que trabajemos en Git

#setwd("C:/Users/Usuario/OneDrive/Documentos/Grupo_fotobiologia/R_functions/Funcion_process_spectrum")


# seteo el directorio en MIR ####
#setwd("C:/Gaston/Proyecto_Photosintesis_git")
#setwd("/Users/sebfer/GoogleDrive/relevamientos_espectros/r_espectro")



# Paquetes 
library(lme4)
library(lmerTest)
library(nlme)
library(car)
library("ggplot2")       
library("lattice")
library("latticeExtra")
library(multcompView)
library(dplyr)
library(plyr)
library(xtable)
library(tidyverse)
library (emmeans)
library("qtl")
library(stringr)
library(data.table)
library(svMisc)
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")
library(Matrix)

filelist.color     <- list.files(path="./Data/rawdata/relevamiento_16_09_2020", pattern = "band.led.sim.",
                       full.names = FALSE)

#filt.list <- "band.led.sim.380_780.txt"

#assuming tab separated values with a header    
df.bands <- bind_rows (lapply(filelist, function(filt.list){
  
  dt <- read.table (file = str_c("./Data/rawdata/relevamiento_16_09_2020/",filt.list), 
                    header = FALSE, sep = "\t",dec = ",")
  #str(dt)
  dt <- dt %>%
        dplyr::rename (units = V1)%>%
        dplyr::rename (value = V2) 
  
  l1 <- dt[1,]
  l2 <- dt[2,]
  
  dt.1 <- dt %>%
          dplyr::filter (units!= "Integration Begin:1107" )%>%
          dplyr::filter (units!= "Integration End:1108") %>%
          dplyr::mutate (l1 = l1[2])%>%
          dplyr::mutate (l2 = l2[2])%>%
          dplyr::mutate (file.origen = filt.list)
  
  
}))


df.leds.1 <- df.leds %>%
             tidyr::separate(pto, c(NA, "ptx"), sep = "_")


unique (df.leds.1$ptx )
df.leds.2 <- df.leds.1 %>%
             tidyr::separate(ptx, c("ptx1" , NA))


colnames(df.leds.2) <- c("lambda", "uW_nm.cm2", "punto")

# Construccion de la matriz de lambdas
d1 <- tibble ( bandwidth = "D1", l1 = 400, l2= 425)
d2 <- tibble ( bandwidth = "D2", l1 = 425, l2= 490)
d3 <- tibble ( bandwidth = "D3", l1 = 490, l2= 560)
d4 <- tibble ( bandwidth = "D4", l1 = 560, l2= 585)
d5 <- tibble ( bandwidth = "D5", l1 = 585, l2= 640)
d6 <- tibble ( bandwidth = "D6", l1 = 640, l2= 700)
dPAR <- tibble ( bandwidth = "PAR", l1 = 400, l2= 700)
dintegr <- tibble ( bandwidth = "total", l1 = 380, l2= 780)

lambdas <- bind_rows (d1, d2, d3, d4, d5, d6, dPAR,dintegr  )

# funcion para la luminaria

list.p1 <- unique (df.leds.2$punto)
#filt.punto <- "A"

df.ppfd.leds <- bind_rows (lapply(list.p1, function (filt.punto){
  
  phot.x <- run_spectrum (df = dist_espectral_uw_raw, lambdas = lambdas, id.survey="sim_sotano", px=filt.punto)
  
  
  print (phot.x)
}))

write_delim (df.ppfd.leds, path="./Data/procdata/df.ppfd.leds.sim.sotano.txt", 
             delim = ",", na = "NA")


df.ppfd.leds$bandwidth

df.ppfd.leds <- df.ppfd.leds %>%
                dplyr::filter (punto != "centro" ) 

promedio.punto <- df.ppfd.leds %>%
  dplyr::group_by(bandwidth)%>%
  dplyr::summarise(mean = mean( umolphotones.s.m2))

class (promedio.punto)
write_delim (promedio.punto, path="./Data/procdata/promedio.punto.sim.sotano.txt", 
             delim = ",", na = "NA")
#unique (df.ppfd.1$punto )

#df.ppfd.2 <- df.ppfd.1 %>%
             #dplyr::mutate (pos = punto) %>%
             dplyr::mutate (pos = fct_recode(pos, "1-1"= "A"))%>%
             dplyr::mutate (pos = fct_recode(pos, "1-2"= "B"))%>%
             dplyr::mutate (pos = fct_recode(pos, "1-3"= "C"))%>%
             dplyr::mutate (pos = fct_recode(pos, "2-1"= "D"))%>%
             dplyr::mutate (pos = fct_recode(pos, "2-2"= "E"))%>%
             dplyr::mutate (pos = fct_recode(pos, "2-3"= "F"))%>%
             dplyr::mutate (pos = fct_recode(pos, "3-1"= "G"))%>%
             dplyr::mutate (pos = fct_recode(pos, "3-2"= "H"))%>%
             dplyr::mutate (pos = fct_recode(pos, "3-3"= "I")) %>%
             tidyr::separate(pos, c("r", "c"), sep ="-")

df.2.PAR <- df.ppfd.2 %>%
       dplyr::filter (bandwidth == "PAR")

# Dummy data
x <- df.2.PAR$r
y <- df.2.PAR$c
data <- expand.grid(X=x, Y=y)
data$Z <- df.2.PAR$umolphotones.s.m2
levelplot(Z ~ X*Y, data=data  ,xlab="X",
          main="")


# Dummy data
x <- df.2.PAR$c
y <-  df.2.PAR$r

x <- LETTERS[1:3]
y <- paste0("var", seq(1,3))
data <- expand.grid(X=x, Y=y)
data$Z <- df.2.PAR$umolphotones.s.m2


library(viridis)

# Heatmap 
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()+ 
  scale_fill_viridis(discrete=FALSE) 
# A =  1-1     
# B =  1-2
# C =  1-3
# D =  2-1
# E=   2-2
# F=   2-3
# G=   3-1
# H =  3-2
# I=   3-3
# 
