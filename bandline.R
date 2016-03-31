# An implementation of Stephen Few's BandLines
# author: Johan Renaudie (aka plannapus)
# modified from the code written in the context of this StackOverflow question:
# http://stackoverflow.com/questions/14577551/how-to-produce-bandlines-using-ggplot2
#
#
# dat is a data.frame: each column is a new time-serie.
# time is the vector of either numeric or POSIXct dates.
# low.col and high.col are the colors used to identify the outliers.

bandline<-function(dat, time, low.col, high.col){
	par(mfcol=c(ncol(dat), 1))
	for(i in 1:ncol(dat)){
		y <- quantile(dat[,i],c(0.25,0.75))
		z<-boxplot.stats(dat[,i])
		r <- range(dat[,i], na.rm=TRUE)
		ifelse(i==1, par(mar=c(0,3,3,3)), 
					ifelse(i==ncol(dat), par(mar=c(3,3,0,3)), 
										 par(mar=c(0,3,0,3))))
		plot(time, dat[,i], axes=FALSE, bty="n", ylim=r, xaxs="i", type="n")
		rect(time[1],y[1], time[length(time)], r[1], col="grey80", border=NA)
		rect(time[1],y[1], time[length(time)], y[2], col="grey60", border=NA)
		rect(time[1],y[2], time[length(time)], r[2], col="grey40", border=NA)
		abline(h=median(dat[,i]),col="white", lwd=2)
		lines(time, dat[,i])
		zhigh <- zlow <- dat[,i]
		zhigh[zhigh<=z$stats[5]]<-NA
		zlow[zlow>=z$stats[1]]<-NA
		points(time, zlow, bg=low.col, pch=21,cex=2)
		points(time, zhigh, bg=high.col, pch=21, cex=2)
		mtext(colnames(dat)[i], side=4, line=1)
		}
	}