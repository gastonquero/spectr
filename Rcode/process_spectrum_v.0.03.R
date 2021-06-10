########################################################
#  Analisis de Datos Luz actinica PAM                  #
#                                                      #
# Gaston Quero - Sebastian Fernandez                   #
#           6-03-2017                                  #
########################################################

getwd()
# seteo el directorio  #### Aca se cambia segun donde tengas los datos 
                          # Esta linea se corre una sola vez en cada caso, y se 
                          # desabilita la linea del directorio de otro usuario
                          # esto lo usamos asi hasta que trabajemos en Git

#setwd("C:/Users/Usuario/OneDrive/Documentos/Grupo_fotobiologia/R_functions/Funcion_process_spectrum")

setwd ("R:/spectr")


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

######### distribucion espectral  ####################################
#########  Ambiente Azul y Rojo  #####################################
#### Actinica 25 lampara nueva  ########################################################

# directorio sebfer
#[dist_espectral]= uW/(cm2.nm)
#dist_espectral_uw_raw <- read.table("./lampara_hpit_AbsoluteIrradiance_0002.txt", header = FALSE, sep = "\t",dec = ",", skip=14)

# Diretorio GQ
dist_espectral_uw_raw <- read.table("./Data/rawdata/lampara_hpit_AbsoluteIrradiance_0002.txt",
                                    header = FALSE, sep = "\t",dec = ",", skip=14)

#View (dist_espectral_uw_raw)
str(dist_espectral_uw_raw)
colnames(dist_espectral_uw_raw) <- c("lambda", "uW_nm.cm2")

dist_espectral_uw_raw <- dist_espectral_uw_raw %>%
                         dplyr::mutate (punto = "X")


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

df  = dist_espectral_uw_raw 
lambdas = lambdas
px = "X"
id.survey="este"

run_spectrum <-  function ( df = NULL, lambdas = NULL, id.survey=NULL, px = NULL) {
   
   id.s <- id.survey
   px <- px
   df.x <- df %>%
           dplyr::filter (punto == px)
   
   dt.1 <- df.x
   dt.lmbd <- lambdas
   
#head(dt.1)
### me quedo solo con el espectro total
dt.lmbd.int <- dt.lmbd %>%
               dplyr::filter (bandwidth == "total" )
    
dt.1.int.total <-  dt.1 %>%
                   dplyr::filter ( lambda >= dt.lmbd.int$l1 ) %>% 
                   dplyr::filter ( lambda <= dt.lmbd.int$l2 ) 

# constantes definidas 
vluz    <- 3e8 #m/s
hplanck <- 6.626070e-34 #J.s
lambda_nm_a_m  <- 1e-9 # nm/m
uW_a_W  <- 1e-6 # uJ/J
Nav <- 6.02e17 # photones/uMol

dt.1.int.phot.total <- dt.1.int.total %>%
                       dplyr::mutate (photones_nm.s.cm2 = uW_nm.cm2*uW_a_W*lambda*lambda_nm_a_m/(hplanck*vluz)) %>%
                       dplyr::mutate (umolphotones_nm.s.cm2 = photones_nm.s.cm2/Nav) %>%
                       dplyr::mutate (umolphotones_nm.s.m2 = umolphotones_nm.s.cm2/1e-4)

dt.lmbd.band <- dt.lmbd %>%
               dplyr::filter (bandwidth != "PAR")%>%
               dplyr::filter (bandwidth != "total")

list.band.1 <- dt.lmbd.band$bandwidth

id.punto <- unique (df.x$punto)

### Este es el grafico de flujo de fotones 
svg (filename= str_c ("./Figures/ppfd_dist.espectral_",id.s,"_", px,".svg"), 
    width=7, 
   height=5, 
  pointsize=12)

plot (umolphotones_nm.s.m2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(dt.lmbd.int$l1,dt.lmbd.int$l2),
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c (id.s, "PPFD", id.punto ,sep="_",collapse = TRUE),
      data= dt.1.int.phot.total)

lapply (list.band.1, function (filt.band) { 
   dt.lmbd.x <- dt.lmbd %>%
                dplyr::filter (bandwidth == filt.band)
   
   abline ( v=dt.lmbd.x$l1, lty=2, col="gray48" )
   abline ( v=dt.lmbd.x$l2 , lty=2, col="gray48" )
   
   text(x = (dt.lmbd.x$l1 + dt.lmbd.x$l2)/2,
         y= max (dt.1.int.phot.total$umolphotones_nm.s.m2 - 1) ,
        label = str_c(dt.lmbd.x$l1 ,"_",dt.lmbd.x$l2),
        cex = 0.9,  col="black")
   })

dev.off ()

plot (umolphotones_nm.s.m2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(dt.lmbd.int$l1,dt.lmbd.int$l2),
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c (id.s, "PPFD", id.punto ,sep="_",collapse = TRUE),
      data= dt.1.int.phot.total)

lapply (list.band.1, function (filt.band) { 
   dt.lmbd.x <- dt.lmbd %>%
      dplyr::filter (bandwidth == filt.band)
   
   abline ( v=dt.lmbd.x$l1, lty=2, col="gray48" )
   abline ( v=dt.lmbd.x$l2 , lty=2, col="gray48" )
   
   text(x = (dt.lmbd.x$l1 + dt.lmbd.x$l2)/2,
        y= max (dt.1.int.phot.total$umolphotones_nm.s.m2 - 1) ,
        label = str_c(dt.lmbd.x$l1 ,"_",dt.lmbd.x$l2),
        cex = 0.9,  col="black")
})


### Este es el grafico de flujo de potencia 
svg (filename= str_c ("./Figures/power_dist.espectral_",id.s,"_", px,".svg"), 
     width=7, 
     height=5, 
     pointsize=12)

plot (uW_nm.cm2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(dt.lmbd.int$l1,dt.lmbd.int$l2),
      col="gray48",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c (id.s, "power", id.punto ,sep="_",collapse = TRUE),
      data= dt.1.int.phot.total)

lapply (list.band.1, function (filt.band) { 
   dt.lmbd.x <- dt.lmbd %>%
      dplyr::filter (bandwidth == filt.band)
   
   abline ( v=dt.lmbd.x$l1, lty=2, col="gray48" )
   abline ( v=dt.lmbd.x$l2 , lty=2, col="gray48" )
   
   text(x = (dt.lmbd.x$l1 + dt.lmbd.x$l2)/2,
        y= max (dt.1.int.phot.total$uW_nm.cm2 - 10) ,
        label = str_c(dt.lmbd.x$l1 ,"_",dt.lmbd.x$l2),
        cex = 0.9,  col="gray48")
})

dev.off ()

plot (uW_nm.cm2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(dt.lmbd.int$l1,dt.lmbd.int$l2),
      col="gray48",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = str_c (id.s, "power", id.punto ,sep="_",collapse = TRUE),
      data= dt.1.int.phot.total)

lapply (list.band.1, function (filt.band) { 
   dt.lmbd.x <- dt.lmbd %>%
      dplyr::filter (bandwidth == filt.band)
   
   abline ( v=dt.lmbd.x$l1, lty=2, col="gray48" )
   abline ( v=dt.lmbd.x$l2 , lty=2, col="gray48" )
   
   text(x = (dt.lmbd.x$l1 + dt.lmbd.x$l2)/2,
        y= max (dt.1.int.phot.total$uW_nm.cm2 - 10) ,
        label = str_c(dt.lmbd.x$l1 ,"_",dt.lmbd.x$l2),
        cex = 0.9,  col="gray48")
})

list.band <- dt.lmbd$bandwidth
#filt.band <- "D1"
df.intervalo <- bind_rows (lapply(list.band, function (filt.band){

   
   dt.lmbd.x <- dt.lmbd %>%
                dplyr::filter (bandwidth == filt.band)
   
   dt.x1.band <- dt.1.int.phot.total  %>%
                        dplyr::filter ( lambda >= dt.lmbd.x$l1 ) %>% 
                        dplyr::filter ( lambda <= dt.lmbd.x$l2 ) %>%
                        dplyr::mutate ( bandwidth = filt.band) %>%
                        dplyr::select ( bandwidth, everything())
   
    dt.x2.band <-  dt.x1.band  %>%
                   dplyr::filter (lambda != min(lambda)) %>%
                   dplyr::select (c(lambda, uW_nm.cm2, umolphotones_nm.s.cm2))

    dt.x3.band <-  dt.x1.band %>%
                           dplyr::filter (lambda != max(lambda)) %>%
                           dplyr::mutate (lmb.2 = dt.x2.band$lambda) %>%
                            dplyr::mutate (uW_nm.cm2.2 = dt.x2.band$uW_nm.cm2)%>%
                           dplyr::mutate (umolphotones_nm.s.cm2.2 = dt.x2.band$umolphotones_nm.s.cm2)%>%
                           dplyr::select (bandwidth, lambda,lmb.2,everything())%>%
                           dplyr::mutate (uW.cm2 = (lmb.2 -lambda) * uW_nm.cm2.2) %>%
                           dplyr::mutate (umolphotones.s.cm2 = (lmb.2 -lambda)*(umolphotones_nm.s.cm2))%>%
                           dplyr::mutate (umolphotones.s.m2 = umolphotones.s.cm2/1e-4)
    
       
    write_delim (dt.x3.band, path = str_c ("./Data/procdata/irradiance_",id.s, "_", filt.band,"_", px,".txt" ), 
                 delim = ",", na = "NA")
    
    
     flujo <-  sum (dt.x3.band$umolphotones.s.m2)
     power <-  sum (dt.x3.band$uW.cm2)
     rad.D <- tibble( dt.lmbd.x, uW.cm2=power, umolphotones.s.m2 =flujo)
     print ( rad.D)
    
}))

df.intervalo <- df.intervalo %>%
                dplyr::mutate (punto = px)%>%
                dplyr::select (punto, everything())


write_delim (df.intervalo, path = str_c ("./Data/procdata/power_ppfd_",id.s,"_", px,".txt"), 
             delim = ",", na = "NA")

return (df.intervalo)
} # aca termina la funcion run_spectrum

phot.x <- run_spectrum (df = dist_espectral_uw_raw, lambdas = lambdas, id.survey="lampara_hpit", px="X" )

#write_delim (phot.x, path= "./Data/procdata/phot.x.txt", delim = ",", na = "NA")



#encuentro indices de lambda = 400 y 700nm
rows_400nm <- which(grepl(400, dist_espectral_uw_raw$lambda))
rows_700nm <- which(grepl(700, dist_espectral_uw_raw$lambda))
index400nm <- rows_400nm[1]
index700nm <- rows_700nm[length(rows_700nm)]

#me quedo solo con el rango par
dist_espectral_uw <- dist_espectral_uw_raw[index400nm:index700nm,]

View (dist_espectral_uw)
str(dist_espectral_uw)


# Energía de un photon de una longitud de onda 
# E = hplanck*vluz/lambda
vluz    <- 3e8 #m/s
hplanck <- 6.626070e-34 #J.s
lambda_nm_a_m  <- 1e-9 # nm/m
uW_a_W  <- 1e-6 # uJ/J
Nav <- 6.02e23 # photones/uMol

E=
dist_espectral_variantes <- dist_espectral_uw %>%
                mutate (photones_nm.s.cm2 = uW_nm.cm2*uW_a_W
                							*lambda*lambda_nm_a_m
                							 /(hplanck* vluz)) %>%
                 mutate (umolphotones_nm.s.cm2 = photones_nm.s.cm2/Nav) 


str(dist_espectral_variantes)    

View(dist_espectral_variantes[index400nm: index700nm,])
 
# estructura de los datos
#svg (filename="./Figures/fig.X2.ambiente.luminico/dist.espectral.ps.sbf3.svg", 
 #    width=7, 
  #   height=5, 
   #  pointsize=12)
par(mfrow=c(2,1))
plot (uW_nm.cm2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(400,750),
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = "uW_nm.cm2",
      data= dist_espectral_variantes)

plot (photones_nm.s.cm2 ~ lambda,
      #ylim = c(0,250),
      xlim=  c(400,750),
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      main = "umolphotones_nm.s.cm2",
      data= dist_espectral_variantes)

