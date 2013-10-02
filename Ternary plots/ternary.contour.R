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
		sqrt(3)/2 ->n
		segments(0,0,0.5,n)
		segments(0.5,n,1,0)
		segments(1,0,0,0)
		text(0,0,labels=colnames(df)[1],pos=1)
		text(1,0,labels=colnames(df)[2],pos=1)
		text(0.5,n,labels=colnames(df)[3],pos=3)
		}

	tern2cart <- function(coord){
		coord[1]->x
		coord[2]->y
		coord[3]->z
		x+y+z->tot
		x/tot -> x
		y/tot -> y
		z/tot -> z
		(2*y + z)/(2*(x+y+z)) -> x1
		sqrt(3)*z/(2*(x+y+z)) -> y1
		return(c(x1,y1))
		}
	
	t(apply(df,1,tern2cart)) -> tern

	interp(tern[,1], tern[,2], z=value, xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution)) -> tern.grid

	image(tern.grid, breaks=breaks, col=levelcol, add=T)
	contour(tern.grid, levels=breaks, add=T, ...)
	points(tern, pch=pch, col=pcol, cex=cex, bg=pbg)
	}