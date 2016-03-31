ternary.grid <- function(col="grey90", lty=2){
	a<-seq(0.9,0.1, by=-0.1)
	b<-rep(0,9)
	c<-seq(0.1,0.9,by=0.1)
	df<-data.frame(x=c(a, b, c, a, c, b),y=c(b, c, a, c, b, a),z=c(c, a, b, b, a, c))
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
	cbind(tern[1:27,],tern[28:54,])->tern
	
	apply(tern,1,function(x){segments(x0=x[1],y0=x[2],x1=x[3],y1=x[4],lty=lty,col=col)})
	}