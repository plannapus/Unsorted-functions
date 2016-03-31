####################################################
#Set of function to retrieve datasets from Pangaea #
####################################################

getPangaeaData<-function(DOI){		
	require(XML)
	DOI <- as.character(DOI)
	if(grepl("^[0-9]",DOI)) DOI <- paste("https://doi.pangaea.de/",DOI,sep="")
	page <- paste(DOI,"?format=html",sep="")
	PAGE <- readLines(page)
	html <- htmlParse(PAGE, encoding="UTF-8")
	title <- xpathSApply(html, "//div[@class='MetaHeaderItem']",xmlValue)[1]
	tables <- readHTMLTable(html)
	params <- tables[[length(tables)-1]]
	data <- tables[[length(tables)]]
	list(Citation=title, Parameters=params, Dataset=data)
	}

getPangaeaMetaData<-function(DOI){ #Grabs site locations in addition to the rest.
	require(XML)
	DOI <- as.character(DOI)
	if(grepl("^[0-9]",DOI)) DOI <- paste("http://doi.pangaea.de/",DOI,sep="")
	page <- paste(DOI,"?format=html",sep="")
	PAGE <- readLines(page)
	html <- htmlParse(PAGE, encoding="UTF-8")
	title <- xpathSApply(html, "//div[@class='MetaHeaderItem']",xmlValue)[1]
	sites <- xpathSApply(html, "//div[@class='MetaHeaderItem geo']/strong",xmlValue)
	latitude <- xpathSApply(html, "//span[@class='latitude']",xmlValue)
	latitude <- as.numeric(latitude[length(latitude)+1-(length(sites):1)])
	longitude <- xpathSApply(html, "//span[@class='longitude']",xmlValue)
	longitude <- as.numeric(longitude[length(longitude)+1-(length(sites):1)])
	tables <- readHTMLTable(html)
	params <- tables[[length(tables)-1]]
	data <- tables[[length(tables)]]
	list(Citation=title, Events = data.frame(Sites=sites, Longitude=longitude, Latitude=latitude), Parameters=params, Dataset=data)
	}
	