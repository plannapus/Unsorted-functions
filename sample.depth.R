sample.depth <- function(samples, intervalsep="/"){
	# Get DSDP/ODP sample depth based on sample name (e. g. 120-747A-1H-1, 45/47cm or 747A-2H-3,32 etc)
	if(length(grep(",",samples))!=length(samples)){
			s <- samples[!grepl(",",samples)]
			samples[!grepl(",",samples)] <- apply(do.call(rbind,strsplit(s," ")), 1,function(x)paste(x[1],x[2], sep=", "))
			}
	if(length(grep(", ",samples))!=length(samples)){
		s <- samples[!grepl(", ",samples)]
		samples[!grepl(", ",samples)] <- apply(do.call(rbind,strsplit(s,",")), 1,
										 	function(x)paste(x[1],x[2], sep=", "))
		}
	S <- do.call(rbind,strsplit(samples,", "))
	mbsf <- rep(NA,length(samples))
	interval_top <- as.integer(do.call(rbind,strsplit(S[,2],split=intervalsep))[,1])
	sect_temp <- strsplit(S[,1],split="-")
	sect <- sapply(sect_temp,function(x)x[length(x)])
	core <- sapply(sect_temp,function(x)x[length(x)-1])
	site <- sapply(sect_temp,function(x)x[length(x)-2])
	leg <- sapply(sect_temp,function(x)ifelse(length(x)>3,x[length(x)-3],""))
	s <- strsplit(site,"")
	is.hole <- lapply(s,`%in%`,LETTERS)
	hole <- lapply(seq_along(s), function(x)s[[x]][is.hole[[x]]])
	site <- sapply(lapply(seq_along(s), function(x)s[[x]][!is.hole[[x]]]),paste,collapse="")
	for(i in seq_along(samples)){if(length(hole[[i]])==0){hole[[i]] <- ""}}
	hole <- unlist(hole)
	C <- strsplit(core,"")
	is.type <- lapply(C,`%in%`,LETTERS)
	type <- sapply(seq_along(C), function(x)C[[x]][is.type[[x]]])
	core <- sapply(lapply(seq_along(C), function(x)C[[x]][!is.type[[x]]]),paste,collapse="")
	s <- strsplit(sect,"")
	is.type2 <- lapply(s,`%in%`,LETTERS)
	sect <- sapply(lapply(seq_along(s), function(x)s[[x]][!is.type2[[x]]]),paste,collapse="",USE.NAMES=F)
	sect[sapply(is.type2,all)]<-sapply(s[sapply(is.type2,all)],paste,collapse="",USE.NAMES=F)
	sect <- sapply(sect,`[`,1)
	
	S <- data.frame(leg=as.character(leg), site=as.character(site), 
					hole=as.character(hole), core=as.character(core), 
					sect=as.character(sect), top=interval_top, stringsAsFactors=FALSE)
	bysite <- split(S, site)
	
	Janus.CorScRequest<-function(site,...){
		require(RCurl)
		require(XML)
		getForm("http://iodp.tamu.edu/janusweb/coring_summaries/coresumm.cgi",site=site,...)->a
		htmlTreeParse(a,useInternalNodes=TRUE)->a1
		xpathApply(a1,path="//pre",xmlValue)[[1]]->a2
		strsplit(a2,split="\n")[[1]]->a3
		a3[a3!=""]->a3
		strsplit(a3,split="\t")->a4
		a5<-array(dim=c(length(a4),10))
		for(i in 1:length(a4)){a5[i,]<-a4[[i]][1:10]}
		gsub(" ","",a5[,1:9])->a5[,1:9]
		a5[-1,]->a5
		as.integer(a5[,1])->leg1
		as.integer(a5[,2])->site1
		as.integer(a5[,4])->core1
		as.character(a5[,3])->hole1
		as.character(a5[,5])->type1
		as.character(a5[,6])->sc1
		as.character(a5[,10])->comment
		as.numeric(a5[,7])->ll
		as.numeric(a5[,8])->cl
		as.numeric(a5[,9])->top
		data.frame(Leg=leg1,Site=site1,H=hole1,Cor=core1,T=type1,Sc=sc1,LL=ll,CL=cl,Top=top,Comment=comment)->res
		return(res)
		}

   for(i in seq_along(bysite)){
   		b <- bysite[[i]]
   		mbsf <- vector(length=nrow(b))
   		cs <- Janus.CorScRequest(site=b$site[1])
   		cs$Sc <- gsub("CC(.*)","CC",cs$Sc)
   		for(j in 1:nrow(b)){
			if(b$hole[j]==""){
				n <- cs[cs$Cor==b$core[j] & cs$Sc==b$sect[j],9]
				}else{n <- cs[cs$H==b$hole[j] & cs$Cor==b$core[j] & cs$Sc==b$sect[j],9]}
   			if(is.na(b$top[j]) & b$sect[j]=="CC"){top <- 0}else{top <- b$top[j]}
   			mbsf[j] <- n + top/100
   			}
   		bysite[[i]]$mbsf <- mbsf
   		if(i==1){cat(i,"site done\r")}else{cat(i,"sites done\r")}
   		}
   	cat("\n")
   	unsplit(bysite,site)
   	}
   			
   		