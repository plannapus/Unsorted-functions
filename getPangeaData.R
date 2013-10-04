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
	as.character(DOI)->DOI
	URLencode(paste("http://doi.pangaea.de/",DOI,sep=""))->doi.url
	URLencode(paste("http://doi.pangaea.de/",DOI,"?format=textfile",sep=""))->dl
	readLines(doi.url)->p
	options(warn=-1)
	p[grep("description",p)][1]->p1
	gsub("<meta name=\"description\" content=\"","",p1)->p1
	strsplit(p1,", Supplement to: ")[[1]]->p1
	strsplit(p1[1]," doi:")[[1]][1]->p1[1]
	strsplit(p1[2],", doi:")[[1]][1]->p1[2]
	c(p1,DOI)->p1
	names(p1)<-c("Dataset name:", "Supplement to:", "DOI:")
	readHTMLTable(p[grep("Parameter",p):tail(grep("/table",p),1)])[[1]]->params
	url(dl)->f
	readLines(f)->g
	(1:length(g))[g=="*/"]->metaskip
	strsplit(g[metaskip+1],split="\t",perl=TRUE)[[1]]->cn
	read.table(f,header=FALSE,sep="\t",skip=metaskip+1)->result
	colnames(result)<-cn
	res<-list(Citation=p1,Parameters=params,Dataset=result)
	options(warn=0)
	return(res)
	}
