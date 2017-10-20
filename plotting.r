
###########################   MAPPING FUNCTION   ###############################
plotonmap <- function(x, lon, lat, pchnum, cexnum, addplus=F, cutof=NA, startzero=F)  {
    pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 50, 100, 500, 1000)
    if (addplus==TRUE)  {    cutoff <- pos[min(which((quantile(x, probs=0.99)-pos)<0))]
      if (!is.na(cutof))  { cutoff <- cutof }
  x[which(x>cutoff)] <- cutoff  }
  if (startzero==F)  {a <- floor(min(x)) } else { a <- 0 } 
  b <- max((x)-a)*1.03
  pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
#  cols <- rainbow(100, start=0.01, end=0.7)[100:1]; if(is.na(pchnum)) {pchnum=15}
  cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]; if(is.na(pchnum)) {pchnum=15}
  map('usa', fill = 1, interior=F, col = gray(0.95), ylim=c(24.8,30.5), xlim=c(-98,-81))
  points(lon, lat, col=cols[pind], pch=pchnum, cex=cexnum)
  axis(1, at=c(-95, -90, -85), lab=c(expression(paste(95,degree,"W")), expression(paste(90,degree,"W")), expression(paste(85,degree,"W"))))
  axis(2, at=c(26, 28, 30), lab=c(expression(paste(26,degree,"N")), expression(paste(28,degree,"N")), expression(paste(30,degree,"N"))), las=1)
  box()
  xloc <- seq(-95, -86, length.out=100)
  for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(25.5,25.5,25.9,25.9), col=cols[j], border=NA) }
  w <- which.min(abs(((max(x)-min(x))/6) - pos))
  if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
  text(xloc[round((xx-a)/b*100+1)], y=25.25, xx, pos=2)
  if (addplus==TRUE)  {  text(max(xloc[round((xx-a)/b*100+1)]), y=25.25, "+")  }    }
################################################################################

###########################   MAPPING FUNCTION   ###############################
plotonmaphires <- function(xp, lon, lat, pchnum, cexnum, addplus=F, cutof=NA, adjlab = 0)  {

library(maps)

load("bathy.RData")
#x <- x[seq(1,2160,2)]
#y <- y[seq(1,671,2)]
#z <- z[seq(1,2160,2), seq(1,671,2)]

z1 <- z
z1[z1>0] <- 1
z1[z1<0] <- 0

    pos <- c(0.005, 0.01,  0.05, 0.2, 0.5, 1, 2, 4, 5, 10, 50, 100, 200, 500, 1000)
    if (addplus==TRUE)  {    cutoff <- pos[min(which((quantile(xp, probs=0.99)-pos)<0))]
      if (!is.na(cutof))  { cutoff <- cutof }
  xp[which(xp>cutoff)] <- cutoff  }
  a <- floor(min(xp));   b <- max((xp)-a)*1.03
  pind <- round((xp-a)/b*100+1); print(min(pind)); print(max(pind))
  #  cols <- rainbow(100, start=0.01, end=0.7)[100:1]
#  cols <- hcl(seq(1,195,length=100)[100:1], seq(10,300,length=100))
#  cols <- hcl(h=c(seq(300,360,length.out=37), seq(0,100,length.out=63)), c=200, alpha=seq(1, 0.4, length.out=100))[100:1]
  cols <- rainbow(100, start=0.65, end=0.17)[100:1]
#  cols <- heat_hcl(100, h = c(0, 90), c. = c(80, 30), l = c(30, 90), power = c(0.2, 2))[100:1]
# cols <- hcl(h=c(seq(300,360,length.out=37), seq(0,100,length.out=63)), c=200, alpha=seq(1, 0.4, length.out=100))[100:1]

  if(is.na(pchnum)) {pchnum=15} 
  map('usa', fill = 0, interior=F, col = "white", ylim=c(24.8,30.45), xlim=c(-97.8,-81))
  image(x,y,z1, col=c(0, gray(0.5)), ylim=c(24.8,30.45), xlim=c(-97.8,-81), xlab="", ylab="", axes=F, add=T)
  axis(1, at=c(-95, -90, -85), lab=c(expression(paste(95,degree,"W")), expression(paste(90,degree,"W")), expression(paste(85,degree,"W"))))
  axis(2, at=c(26, 28, 30), lab=c(expression(paste(26,degree,"N")), expression(paste(28,degree,"N")), expression(paste(30,degree,"N"))))
  contour(x,y,z, levels=c(-200, -500), col=gray(0.8), add=T, labels=c("", "")) 
  points(lon, lat, col=cols[pind], pch=pchnum, cex=cexnum)
  box()
  text(-86.2, 29, "200m", col=gray(0.6), font=1, cex=0.8, srt = 325)
  text(-86.3, 28.2, "500m", col=gray(0.6), font=1, cex=0.8, srt = 325)
  xloc <- seq(-95, -86, length.out=100)
  for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(25.5,25.5,25.9,25.9), col=cols[j], border=NA) }
  w <- which.min(abs(((max(xp)-min(xp))/6) - pos))
  if(-pos[w]<min(xp)) { xx <- seq(0, max(xp), pos[w]); xx <- xx[xx>min(xp)] } else {  xx <- c(seq(-pos[w], min(xp), -pos[w]), seq(0, max(xp), pos[w])) }
  text(xloc[round((xx-a)/b*100+1)], y=24.65+adjlab, xx, pos=3)
  if (addplus==TRUE)  {  text(max(xloc[round((xx-a)/b*100+1)])+0.55, y=24.65+adjlab, "+", pos=3)  }    }
################################################################################


                                                                                            
