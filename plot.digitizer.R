plot.localizer<-function(file=file.choose()){
	cat("Choose the file to proceed:\n")
	gsub("\\\\ ","\\ ",file)->file
	gsub("^\\ +||| +$", "", file)->file
	if(length(grep("[pP][nN][gG]",file))!=0){
		require(png)
		img <- readPNG(file)
		as.raster(img)->img
		quartz(title=gsub(".[pP][nN][gG]","",file), width=ncol(img)/100, height=nrow(img)/100)	
		par(mar=c(0,0,0,0))
		plot(NA,NA,type="n",xlim=c(1,ncol(img)),ylim=c(1,nrow(img)), asp=1, axes=FALSE,xlab="",ylab="")
		rasterImage(img,1,1,ncol(img),nrow(img))
		}
	if(length(grep("[jJ][pP][eE]?[gG]",file))!=0){
		require(ReadImages)
		img <- read.jpeg(file)
		quartz(title=gsub(".[jJ][pP][eE]?[gG]","",file), width=ncol(img)/100, height=nrow(img)/100)
		par(mar=c(0,0,0,0))
		plot(NA,NA,type="n",xlim=c(1,ncol(img)),ylim=c(1,nrow(img)), asp=1, axes=FALSE,xlab="",ylab="")
		plot(img)
	}
	cat("Click on the lower end of the x-axis\n")
	minx <- locator(1)
	minx$value <- as.numeric(readline("Enter its value: "))
	cat("Click on the upper end of the x-axis\n")
	maxx <- locator(1)
	maxx$value <- as.numeric(readline("Enter its value: "))
	cat("Click on the lower end of the y-axis\n")
	miny <- locator(1)
	miny$value <- as.numeric(readline("Enter its value: "))
	cat("Click on the upper end of the y-axis\n")
	maxy <- locator(1)
	maxy$value <- as.numeric(readline("Enter its value: "))
	cat("Click on any points you want (type ESC when done)\n")
	p <- locator()
	minx$value+(p$x-minx$x)*(maxx$value-minx$value)/(maxx$x-minx$x) -> p$xvalue
	miny$value+(p$y-miny$y)*(maxy$value-miny$value)/(maxy$y-miny$y) -> p$yvalue
	cbind(p$xvalue,p$yvalue) -> res
	colnames(res)<-c("X","Y")
	return(res)
	}