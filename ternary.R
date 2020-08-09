tern2cart <- function(coord){
	x <- coord[1]
	y <- coord[2]
	z <- coord[3]
	tot <- x+y+z
	x <- x/tot
	y <- y/tot
	z <- z/tot
	x1 <- (2*y + z)/(2*(x+y+z))
	y1 <- sqrt(3)*z/(2*(x+y+z))
	c(x1,y1)
	}

ternary.plot <- function(df, ...){
	# df: dataframe containing the xyz values
	# ...: additional arguments to be passed to function points. See ?points for more information
	plot(NA,NA,xlim=c(0,1),ylim=c(0,sqrt(3)/2),asp=1,bty="n",axes=F,xlab="",ylab="")
	n <- sqrt(3)/2
	segments(0,0,0.5,n)
	segments(0.5,n,1,0)
	segments(1,0,0,0)
	text(0,0,labels=colnames(df)[1],pos=1)
	text(1,0,labels=colnames(df)[2],pos=1)
	text(0.5,n,labels=colnames(df)[3],pos=3)
	tern <- t(apply(df,1,tern2cart))
	points(tern, ...)
	}

ternary.grid <- function(col="grey90", lty=2){
	a <- seq(0.9,0.1, by=-0.1)
	b <- rep(0,9)
	c <- seq(0.1,0.9,by=0.1)
	df <- data.frame(x=c(a, b, c, a, c, b),y=c(b, c, a, c, b, a),z=c(c, a, b, b, a, c))
	tern <- t(apply(df,1,tern2cart))
	tern <- cbind(tern[1:27,],tern[28:54,])
	apply(tern,1,function(x){segments(x0=x[1],y0=x[2],x1=x[3],y1=x[4],lty=lty,col=col)})
	}

ternary.contour <- function(df, value, add=FALSE,
	resolution=0.001, breaks, levelcol,
	pch=19, pcol="black", pbg="black", cex=1, ...){
	# df: data.frame with xyz values
	# value: value attributed to each points
	# resolution: resolution of the contour plot (default: one value every 0.001)
	# breaks: vector of levels for the contour plot.
	# levelcol: vector of colors for the contour plot (need to be number of breaks minus 1), default to heat gradient with a transparency effect.
	# pch: point type
	# pcol: point color
	# pbg: point background color
	# cex: point size
	# add: if TRUE, add to previous plot; if FALSE create a new plot.
	# ...: additional parameters to pass to function contour, e. g. lty (for line type), lwd (line width), drawlabels (default to TRUE), labcex (cize of the labels), etc. See ?contour for more info.
	require(akima)
	if(missing(levelcol)){levelcol<-rev(heat.colors(length(breaks)-1, alpha=0.5))}
	if(add==FALSE){
		plot(NA,NA,xlim=c(0,1),ylim=c(0,sqrt(3)/2),asp=1,bty="n",axes=F,xlab="",ylab="")
		n <- sqrt(3)/2
		segments(0,0,0.5,n)
		segments(0.5,n,1,0)
		segments(1,0,0,0)
		text(0,0,labels=colnames(df)[1],pos=1)
		text(1,0,labels=colnames(df)[2],pos=1)
		text(0.5,n,labels=colnames(df)[3],pos=3)
		}
	tern <- t(apply(df,1,tern2cart))
	tern.grid <- interp(tern[,1], tern[,2], z=value, xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution))
	image(tern.grid, breaks=breaks, col=levelcol, add=T)
	contour(tern.grid, levels=breaks, add=T, ...)
	points(tern, pch=pch, col=pcol, cex=cex, bg=pbg)
	}
