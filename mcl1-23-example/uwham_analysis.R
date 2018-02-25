# cat r*/*.out > data
# R CMD BATCH uwham_analysis.R

#.libPaths("/home/egallicc/R/x86_64-unknown-linux-gnu/3.0")
library("UWHAM")

npot.fcn <- function(e0,ebind, bet, lam) -bet*(e0 + lam*ebind)

uwham.r <- function(label,logQ,ufactormax,ufactormin=1){
  n <- dim(logQ)[1]
  m <- dim(logQ)[2]
  iniz <- array(0,dim=m) 
  uf <- ufactormax
  while(uf >= ufactormin & uf >= 1){
    mask <- seq(1,n,trunc(uf))
    out <- uwham(label=label.cross[mask], logQ=neg.pot[mask,],init=iniz)
    show(uf)
    iniz <- out$ze
    uf <- uf/2
  }
  out$mask <- mask
  out
}

histw <-
function (x, w, xaxis, xmin, xmax, ymax, bar = TRUE, add = FALSE, 
            col = "black", dens = TRUE) 
{
  nbin <- length(xaxis)
  xbin <- cut(x, breaks = xaxis, include.lowest = T, labels = 1:(nbin -  1))
  y <- tapply(w, xbin, sum)
  y[is.na(y)] <- 0
  y <- y/sum(w)
  if (dens) 
    y <- y/(xaxis[-1] - xaxis[-nbin])
  if (!add) {
    plot.new()
    plot.window(xlim = c(xmin, xmax), ylim = c(0, ymax))
    axis(1, pos = 0)
    axis(2, pos = xmin)
  }
  if (bar == 1) {
    rect(xaxis[-nbin], 0, xaxis[-1], y)
  }
  else {
    xval <- as.vector(rbind(xaxis[-nbin], xaxis[-1]))
    yval <- as.vector(rbind(y, y))
    lines(c(min(xmin, xaxis[1]), xval, max(xmax, xaxis[length(xaxis)])), 
          c(0, yval, 0), lty = "11", lwd = 2, col = col)
  }
  invisible()
  list(y = y, breaks = xaxis)
}



data.t <- read.table("data")

temperatures.t <- data.t$V2
lambdas.t <- data.t$V3
binding.energies.t <- data.t$V4
total.energies.t <- data.t$V5
e0.t <- total.energies.t - 0. * lambdas.t * binding.energies.t

lam <- c(0.000001,0.002,0.004,0.008,0.01,0.02,0.04,0.07,0.1,0.17,0.25,0.35,0.5,0.6,0.7,0.8,0.9,1.0)
tempt <- c(300)
bet <- 1.0/(0.001986209*tempt)
mtempt <- length(bet)
mlam <- length(lam)
m <- mlam*mtempt
N <- length(binding.energies.t)


neg.pot <- matrix(0, N,m)
sid <- 1
# note the order of (be,te)
for (be in 1:mlam) {
     for (te in 1:mtempt) {
             neg.pot[,sid] <- npot.fcn(e0=e0.t,ebind=binding.energies.t,bet[te],lam[be])
             sid <- sid + 1
    }
}


# note levels
label.tempt <- factor(temperatures.t, levels=tempt, labels=1:mtempt)
label.lam <- factor(lambdas.t, levels=lam, labels=1:mlam)
label.cross <- (as.numeric(label.lam)-1)*mtempt + as.numeric(label.tempt)
out <- uwham.r(label=label.cross, logQ=neg.pot,ufactormax=4,ufactormin=1)
ze <- matrix(out$ze, nrow=mtempt, ncol=mlam)
-ze/bet
ve <- matrix(out$ve, nrow=mtempt, ncol=mlam)
sqrt(ve)/bet


Vsite <- 4.*pi*(7.0)^3/3.0
dgsite <- -log(Vsite/1668.0)/bet[]
dgbind <- dgsite + (-ze[,mlam]/bet[]) - (-ze[,1]/bet[])
ddgbind <- sqrt(ve[,mlam]+ve[,1])/bet

#average energy
de = mean(binding.energies.t[lambdas.t > 0.99])
dde = sd(binding.energies.t[lambdas.t > 0.99])/sqrt(length(binding.energies.t[lambdas.t > 0.99]))

#DGbind as a function of temperature
dgbind
sink("result.log")
cat("DGb = ", dgbind[1]," +- ",ddgbind[1]," DE = ", de, " +- ",dde,"\n")
sink()

#plot(tempt,dgbind)
#mask <- seq(1,N,1)
#idgomax <-   (18-1)*mtempt + 1
#binding energy distribution at 300K lambda=1
#bins <- seq(-100, 0, 2.0)
#hw <- histw(data.t$V7[mask], w=out$W[,idgomax], xaxis=bins, xmin=-100, xmax=-40, ymax=0.1)
#global rmsd distribution at 300K lambda=1
#bins <- seq(0, 10, 0.5)
#hw <- histw(data.t$V7[mask], w=out$W[,idgomax], xaxis=bins, xmin=0, xmax=10, ymax=1.0, bar=TRUE)

#hw <- histw(data.t$V7[mask], w=out$W[,(18-1)*mtempt + 8], xaxis=bins, xmin=-90, xmax=0, ymax=0.1, bar = TRUE)

