myTransform<-function (x, CRSobj, ...){
	#same as spTransform but does not throw error if points unprojectable: it deletes them instead and cut the lines if needed.
	if(class(res)!="SpatialLinesDataFrame") stop("Only implemented for SpatialLinesDataFrame so far.")
    xSP <- as(x, "SpatialLines")
    xDF <- as(x, "data.frame")
    from_args <- proj4string(xSP)
    to_args <- slot(CRSobj, "projargs")
    input <- slot(xSP, "lines")
    output <- vector(mode = "list", length = length(input))
    for (i in seq_along(input)){
    	y <- input[[i]]
    	ID <- slot(y, "ID")
    	inSP <- slot(y, "Lines")
    	n <- length(inSP)
    	out <- vector(mode = "list", length = n)
    	for (j in 1:n){
    		crds <- slot(inSP[[j]], "coords")
    		m <- nrow(crds)
   			attr(m, "ob_tran") <- 0L
   			res <- .Call("transform", from_args, to_args, m, as.double(crds[,1]), as.double(crds[, 2]), NULL, PACKAGE = "rgdal")
        	crds <- cbind(res[[1]], res[[2]])
        	crds[!is.finite(crds)]<-NA
        	if(all(is.na(crds))){
        		out[[j]]<-NULL
        	}else if(all(!is.na(crds))){
        		out[[j]] <- Line(coords = crds)
        	}else{
        		d <- cumsum(is.na(crds[,1]))
        		d <- d[!is.na(crds[,1])&!is.na(crds[,2])]
        		crds <- crds[!is.na(crds[,1])&!is.na(crds[,2]),]
        		CRDS <- split(as.data.frame(crds),d)
        		sp <- list()
    			for(k in seq_along(CRDS)) sp[[k]] <- Line(coords = as.matrix(CRDS[[k]]))
    			out <- c(out,sp)
    		}
            }
        out <- out[!!sapply(out,length)]
    	if(length(out)) output[[i]] <- Lines(out, ID)

    }
    na <- sapply(output,length)!=0
    output <- output[na]
    xDF <- xDF[na,]
    resSP <- SpatialLines(output, proj4string = CRS(to_args))
    res <- SpatialLinesDataFrame(sl = resSP, data = xDF, match.ID = FALSE)
    res
}