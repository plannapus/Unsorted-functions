#################################
### Get a dataset from Pangea ###
#################################
## Given a dataset DOI, output ##
## a list containing the orig- ##
## inal citation, metadata and ##
## the dataset as a dataframe. ##
#################################

getPangeaData <-
function(DOI){		
	require(XML)
	DOI <- as.character(DOI)
	doi.url <- URLencode(paste("http://doi.pangaea.de/",DOI,sep=""))
	dl <- URLencode(paste("http://doi.pangaea.de/",DOI,"?format=textfile",sep=""))
	p <- readLines(doi.url)
	title <- strsplit(grep("<meta name=\"title\"",p,value=TRUE),"\"")[[1]][4]
	author <- strsplit(grep("<meta name=\"author\"",p,value=TRUE),"\"")[[1]][4]
	date <- as.Date(strsplit(grep("<meta name=\"date\"",p,value=TRUE),"\"")[[1]][4])
	source <- gsub("^In Supplement to: ","",strsplit(grep("<meta name=\"DC.source\"",p,value=TRUE),"\"")[[1]][4])
	p1 <- c(paste(author, " (",format(date,"%Y"),") ", title,sep=""), source, DOI)
	names(p1)<-c("Dataset title","Supplement to:", "DOI:")
	readHTMLTable(p[grep("Parameter",p):tail(grep("/table",p),1)])[[1]]->params
	g <- readLines(dl)
	metaskip <- (1:length(g))[g=="*/"]
	tab <- read.table(dl,sep="\t",skip=metaskip,header=TRUE,check.names=FALSE)
	list(Citation=p1,Parameters=params,Dataset=tab)
	}
