getPangeaList <-
function(DOI){				
	require(XML)
	DOI <- as.character(DOI)
	doi.url <- URLencode(paste("http://doi.pangaea.de/",DOI,sep=""))
	result <- as.data.frame(gsub("^ | $","",do.call(rbind,strsplit(xpathSApply(htmlParse(doi.url),"//ol/li/div/a",xmlValue),":|doi:"))))
	colnames(result)<-c("Source","Name","DOI")
	result
	}