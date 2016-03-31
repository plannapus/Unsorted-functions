ternary.plot <- function(df, ...){
	# df: dataframe containing the xyz values
	# ...: additional arguments to be passed to function points. See ?points for more information
	
	plot(NA,NA,xlim=c(0,1),ylim=c(0,sqrt(3)/2),asp=1,bty="n",axes=F,xlab="",ylab="")
	sqrt(3)/2 ->n
	segments(0,0,0.5,n)
	segments(0.5,n,1,0)
	segments(1,0,0,0)
	text(0,0,labels=colnames(df)[1],pos=1)
	text(1,0,labels=colnames(df)[2],pos=1)
	text(0.5,n,labels=colnames(df)[3],pos=3)
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
	points(tern, ...)
	}