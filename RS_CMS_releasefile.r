################################################################################
###################   red snapper release file for CMS  ########################
###################   M. Karnauskas Jan 8, 2015         ########################
#
#  prepares particle release file for input to Connectivity Modeling System
#  code takes output from statistical RS distribution model 
#  published in Karnauskas et al. 2017 (MCF)
#  
################################################################################


rm(list=ls())

#######################   libraries and functions   ############################
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
library(maps)
if (!"sp" %in% installed.packages()) install.packages("sp", repos='http://cran.us.r-project.org')
library(sp)
if (!"splancs" %in% installed.packages()) install.packages("splancs", repos='http://cran.us.r-project.org')
library(splancs)

source("plotting.r")     # GoM data mapping function                                             
                                  
###########################      import data      ##############################
                                                 
mat <- read.table("KarnauskasMCFoutputs.csv", header=T, sep=",")        #  use Ftotal - total egg production including artificial structures
head(mat)
plotonmap(mat$totalF_Fig7low, mat$longitude, mat$latitude, 15, 0.9, addplus=T); text(-90.5,26.25, "index of fecundity")

dat <- data.frame(cbind(mat$longitude, mat$latitude, mat$totalF_Fig7low))
names(dat) <- c("lon", "lat", "N")

co <- read.table("redSnapperRadialGrid_STATE.xyz", header=F)                    # edited version to align with state boundaries

################################################################################                                                    

################  plot original recruitment habitat grid cells  ################
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, col=j)
  text(m[1,1], m[1,2], m[1,3], cex=0.6)      }
  
points(dat$lon, dat$lat, pch=19, cex=0.5)             # release locations

################  get polygon numbers for release locations  ###################

pts = SpatialPoints(cbind(dat$lon, dat$lat))
m <- co[which(co[,3]==1), ]
m <- rbind(m, m[1,])
oldpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=1)
for (j in 2:max(co$V3))  {
  m <- co[which(co[,3]==j), ]
  m <- rbind(m, m[1,])
  newpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=j)
  oldpol <- c(oldpol, newpol)           }
pol = SpatialPolygons(oldpol)
polynams <- over(pts, pol)
rel <- cbind(polynams, dat)
rel                                     # file coded with polygon numbers

points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/rel$polynams+1)

##############  label release site with closest polygon number  ################
ctrslat <- tapply(co$V2, co$V3, mean)        # find center point of each polygon
ctrslon <- tapply(co$V1, co$V3, mean)        # acts as approximate method for assigning nearest polygon            
points(ctrslon, ctrslat, pch=20, cex=2)

##############  10/13/2017 adjustment for state boundaries      ################
FL <- data.frame(cbind(c(-87.51833333, -80, -80, -87.51833333, -87.51833333), c(24, 24, 32, 32, 32)))
AL <- data.frame(cbind(c(-88.385, -87.51833333, -87.51833333, -88.385, -88.385), c(24, 24, 32, 32, 32)))
MS <- data.frame(cbind(c(-89.16666667, -88.385, -88.385, -89.16666667, -89.16666667), c(24, 24, 32, 32, 32)))
LA <- data.frame(cbind(c(-93.795, -92.88333333, -89.16666667, -89.16666667, -93.85, -93.795), c(29.535, 26.19, 26, 32, 29.74, 29.535)))
TX <- data.frame(cbind(c(-93.85, -93.795, -92.88333333, -100, -100, -93.85), c(29.74, 29.535, 26.19, 25.8, 30, 29.74)))
PT <- cbind(rel$lon, rel$lat)
names(PT) <- c("x","y")
names(TX) <- c("x","y")
names(LA) <- c("x","y")
names(MS) <- c("x","y")
names(AL) <- c("x","y")
names(FL) <- c("x","y")
TXout = inpip(PT, TX, bound=TRUE)
LAout = inpip(PT, LA, bound=TRUE)
MSout = inpip(PT, MS, bound=TRUE)
ALout = inpip(PT, AL, bound=TRUE)
FLout = inpip(PT, FL, bound=TRUE)
points(mat$statlon[TXout], mat$statlat[TXout], pch=15, cex=0.5, col=1)
points(mat$statlon[LAout], mat$statlat[LAout], pch=15, cex=0.5, col=2)
points(mat$statlon[MSout], mat$statlat[MSout], pch=15, cex=0.5, col=3)
points(mat$statlon[ALout], mat$statlat[ALout], pch=15, cex=0.5, col=4)
points(mat$statlon[FLout], mat$statlat[FLout], pch=15, cex=0.5, col=5)

for (i in TXout) {
  if (is.na(rel$polynams[i])) {
      pl <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))
      if (pl > 24) { pl <- 24 }
      rel$polynams[i] <- pl   } }      
for (i in LAout) {
  if (is.na(rel$polynams[i])) {
      pl <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))
      if (pl < 25) { pl <- 25 }
      if (pl > 40) { pl <- 40 }
      rel$polynams[i] <- pl   } }   
for (i in MSout) {
  if (is.na(rel$polynams[i])) {
      pl <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))
      if (pl < 41) { pl <- 41 }
      if (pl > 42) { pl <- 42 }
      rel$polynams[i] <- pl   } }   
for (i in ALout) {
  if (is.na(rel$polynams[i])) {
      pl <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))
      if (pl < 43) { pl <- 43 }
      if (pl > 45) { pl <- 45 }
      rel$polynams[i] <- pl   } }   
for (i in FLout) {
  if (is.na(rel$polynams[i])) {
      pl <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))
      if (pl < 46) { pl <- 46 }      
      rel$polynams[i] <- pl   } }                        
     
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/rel$polynams+1)
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)
     
######################    plot new polygon assignments   #######################
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, lwd=1)
  text(m[1,1], m[1,2], m[1,3], cex=0.6)         }
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])    #  check
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)          #  check
polygon(FL, border=1)
polygon(AL, border=2)
polygon(MS, border=3)
polygon(LA, border=4)
polygon(TX, border=5)

#### adjust depth to make sure spawning releases are not below ocean floor  ####

urla <- "http://coastwatch.pfeg.noaa.gov/erddap/griddap/hycom_gom201D.nc?u[(2010-07-16T00:00:00Z):1:(2010-07-16T00:00:00Z)][(0.0):1:(200)][(20):1:(31)][(-98.0):1:(-81)],v[(2010-07-16T00:00:00Z):1:(2010-07-16T00:00:00Z)][(0.0):1:(200)][(20):1:(31)][(-98.0):1:(-81)]"
#download.file(url=urla, destfile="temp.nc", mode="wb")
nc <- nc_open("temp.nc")        # open netcdf and get temp variable
v1 <- nc$var[[1]]
u <- ncvar_get(nc, v1)
dim(u)
v1 <- nc$var[[2]]
v <- ncvar_get(nc, v1)
dim(v)       
lon <- nc$var[[1]]$dim[[1]]$vals
lat <- nc$var[[1]]$dim[[2]]$vals
dep <- nc$var[[1]]$dim[[3]]$vals
cur <- sqrt(u^2 + v^2)
image(lon, lat, cur[,,1])             # check data inputs

d <- matrix(NA, 426, 307)
for (i in 1:426) {             # convert HYCOM mask 
  for (j in 1:307) {
     d[i,j] <- dep[max(which(!is.na(u[i,j,])))]  }}      
image(d)                       # check assignment
nc_close(nc)

rel$depest <- NA               # determine depth of HYCOM mask at each location
for (i in 1:nrow(rel)) {  rel$depest[i] <- d[which.min(abs(lon - rel$lon[i])), which.min(abs(lat - rel$lat[i]))] }
head(rel)

plot(rel$lon, rel$lat, pch=15, col=gray(1-(rel$depest/max(rel$depest+10)+0.04)))

table(rel$depest)
hist(rel$depest)

length(which(rel$depest<45))
length(which(rel$depest>=45))

rel$spawndep <- rel$depest - 2                     # set spawning 2m above ocean floor
rel$spawndep[which(rel$spawndep > 45)] <- 45       # for >45m, set at 45m

minspawndep <- 15                   #  minimum spawning depth set at 15m based on literature

length(which(rel$depest<minspawndep))
dim(rel)
tapply(rel$N, rel$depest<minspawndep, sum)     # relative spawning biomass above 15m cutoff
rel <- rel[which(rel$depest>=minspawndep),]
dim(rel)

plot(rel$depest, rel$spawndep)                 # check assigned spawning dewpths
table(rel$spawndep < rel$depest)
table(rel$spawndep)
table(rel$depest, rel$spawndep)

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)      # check on map
              
######################### input temporal information ###########################

days <- data.frame(seq(as.Date("2003-01-30"), as.Date("2017-12-31"), by="days"))        #  change ending date if update
names(days) <- "dates"
days$yr <- as.numeric(substr(days$dates, 1, 4))
days$mo <- as.numeric(substr(days$dates, 6, 7))
days$da <- as.numeric(substr(days$dates, 9, 10))

days <- days[which(days$mo > 4 & days$mo < 9),]     # spawning occurs May - Aug
days <- days[seq(1, nrow(days), 3),]                # spawning every 3 days
lis <- days[, 2:4]
dim(lis)

table(lis$yr)
table(lis$yr, lis$mo)
lis <- lis[-which(lis$mo==5 & lis$da==1),]    # take one day out so that there are 10 release days per month

dim(lis)
table(lis$yr, lis$mo)
table(lis$yr)
mean(table(lis$yr))

lis$doy <- NA                                 #  add spawning activity column
dinmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for (i in 1:nrow(lis))  
  {  lis$doy[i] <- (sum(dinmon[1:lis$mo[i]]) + lis$da[i]) / 365  }              # convert to day of year
lis$spawnact <- (lis$doy/0.536)^(0.536/0.024) * exp((0.536-lis$doy)/0.024)      # relationship from Porch et al. 2015 (MCF)
plot(lis$doy, lis$spawnact)                                                     # check relationship
plot(lis$spawnact)
lis <- lis[-c(4)]                                                               # remove DOY column

###############  input spatial site information (from above)  ###################

rel1 <- rel[-c(5)]
head(rel1)
m <- rel1       # 'm' is list of release sites with columns: polygon, lon, lat, number of releases
head(m)

m$N <- m$N / 2.01   # 4.3 # * 3.1   #  18.9              #  / 2.01 works for LoRes
mean(m$N); min(m$N); max(m$N)
which(m$N==0)

prod(nrow(lis), nrow(m))
###################### now, making the release file ############################

mat <- as.data.frame(matrix(data=NA, nrow=nrow(lis)*nrow(m), ncol=9))   # empty matrix to be filled
mat[,1] <- rep(m[,1], nrow(lis))                                        # column 1: release polygon number
mat[,2] <- rep(m[,2], nrow(lis))                                        # column 2: release longitude 
mat[,3] <- rep(m[,3], nrow(lis))                                        # column 3: release latitude 
mat[,4] <- rep(m[,5], nrow(lis))                                        # column 4: release depth
mat[,5] <- rep(m[,4], nrow(lis))                                        # column 5: number of particles per release
mat <- mat[order(mat[,1], mat[,5], mat[,3], mat[,2]), ] # !!!!  CHECK THIS WITH NEW DATA   !!!!                                # resort matrix
# mat
mat[,6] <- rep(lis[,1], nrow(m))                                        # column 6: release year
mat[,7] <- rep(lis[,2], nrow(m))                                        # column 7: release month
mat[,8] <- rep(lis[,3], nrow(m))                                        # column 8: release day
mat[,9] <- 0                                                            # column 9: release hour

mat[,10] <- rep(lis[,4], nrow(m))                                       # column 10: scale by spawning activity!!!!
sum(mat[,5])

mat <- mat[order(mat[,6], mat[,7], mat[,8]), ] # !!!!  REORDER SO DATES ARE TOGETHER   !!!!
head(mat)

x <- mat$V5 * mat$V10
mean(x); min(x); max(x)

mat$V5 <- round(mat$V5 * mat$V10)
sum(mat[,5])
head(mat)

matfin <- mat[-c(10)]
head(matfin)
dim(matfin)
table((matfin$V5 > 0))
mean(matfin$V5); min(matfin$V5); max(matfin$V5)

tapply(x, (matfin$V5 > 0), sum)/sum(x)

matfin <- matfin[which(matfin$V5 >0),]

head(matfin)
dim(matfin)
mean(matfin$V5); min(matfin$V5); max(matfin$V5)
sum(matfin$V5)

getwd()
dim(matfin)
dim(matfin)/8


k <- max(which((nrow(matfin)/33:64 - round(nrow(matfin)/33:64))==0) + 32)    # calculates number of nodes which should be used.
print(paste("Use", k, "nodes."), sep=" ", quote=F)


#mat1 <- mat[which(mat$V6==2014),]; dim(mat1)
#write.table(mat1, file="C:/Users/mkarnauskas/Desktop/red_snapper/Jan2015/RSrel2014.txt", sep="\t", col.names=F, row.names=F)          # WRITE FILE TO TXT

write.table(matfin, file="C:/Users/mkarnauskas/Desktop/red_snapper/Mar2016/RS_FaReg_10perLoss_15mCut.txt", sep="\t", col.names=F, row.names=F)          # WRITE FILE TO TXT


# double check

table(matfin$V6, matfin$V7)
table(matfin$V6)
matplot(table(matfin$V7, matfin$V6), type="l")
diff(table(matfin$V6))

tapply(matfin$V5, list(matfin$V6, matfin$V7), sum)
matplot(tapply(matfin$V5, list(matfin$V7, matfin$V6), sum), type="l")

f <- which(matfin$V6==2003 & matfin$V7 == 8 & matfin$V8 == 5); length(f)
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

f <- which(matfin$V6==2003 & matfin$V7 == 6 & matfin$V8 == 24); length(f)
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)



#################### end release file-making #####################################
