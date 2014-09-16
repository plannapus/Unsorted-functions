#######################################################################################################################
#                    Subsampling Neptune database (with pacman trimming option) version 3.0                           #
#######################################################################################################################
# author of the implementation: Johan Renaudie.                                                                       #
# authors of the algorithms: Lazarus et al 2012 for pacman, Sanders 1968 for rarefaction, Shinozaki 1963              #
# for UW, Alroy 1996 for OW, Alroy 2000 for O2W and Alroy 2010 for SQS.                                               #
#######################################################################################################################
# id : vector of taxon id                                                                                             #
# age_id : vector of corresponding age values                                                                         #
# sample_id : vector of corresponding samples                                                                         #
# method : can be "all", or any combinations of "SQS","UW","OW","O2W" and "CR" (classical rarefaction)                #
# quota : quota used for each subsampling procedures. If method is "all", the quotas have to be in the order          #
# 	sqs-uw-ow-o2w-cr, otherwise the order given in the method vector. (e.g. if  method=c("CR","OW") then              #
# 	(e.g. if  method=c("CR","OW") then quota=c(100,50) where 100 is the quota for classical rarefaction and           #
#	50 for occurence-weighted by-list subsampling)                                                                    #
# trials : number of trials for each subsampling                                                                      #
# dominant : only useful for SQS subsampling. "include" includes the dominant taxa in each bin, "exclude" excludes it.#
# pac_top, pac_bottom : top and bottom trimming percentage for the pacman procedure                                   #
# bins : time bins. Can be regular (e.g., seq(0,65,by=0.5)) or irregular (stage boundaries)                           #
# write: if TRUE, allows an csv output of the main results to be created (in the default directory).                  #
# rt_method : method used for range-through. Can be either of 'bc' for boundary-crossers, 'rt' for                    #
# 	classic range-through (without singletons) and 'tot' (with singletons).                                           #
# seed : if user wants result to be reproducible, enter a seed for the PRNG, otherwise use current seed of session.   #
#######################################################################################################################

pacnsub <- function(id, age_id, sample_id, 
					method="all", quota=c(.6,10,100,500,100), trials=100, dominant="include", 
					pac_top=0, pac_bottom=0,
					rt_method="rt", rate_method="foote",
					bins=seq(0,65,0.5),
					write=TRUE, seed){
						
	id <- id[!is.na(age_id)]
	sample_id <- sample_id[!is.na(age_id)]
	age_id <- age_id[!is.na(age_id)]
	
	methods <- c("SQS","UW","OW","O2W","CR")
	if("all"%in%method) method <- methods
	quotas <- numeric(5)
	if(length(quota)==length(method)){
		quotas[methods%in%method] <- quota
		}else{
			stop("Please provide a quota for each method (order if all: SQS, UW, OW, O2W and CR)")
			}
	
	if(!missing(seed)) set.seed(seed)
	
	rangethrough<-function(x, rtm){
		if(rtm=="bc"){
			return(((x$bL+x$bt)+(x$Ft+x$bt))/2)
			}else if(rtm=="rt"){
			return(x$bL+x$bt+x$Ft)
			}else if(rtm=="tot"){
			return(x$bL+x$bt+x$Ft+x$FL)
			}else{
			stop("Argument rt_method accept only three methods:\nbc : Standing mean diversity based on boundary-crossers\nrt : Classic range-through excluding singletons\ntot : Classic range-through including singletons")
			}
		}
	
	crossers <- function(sib,nb){
		bL<-bt<-Ft<-FL<-numeric(nb)
		bL[1] <- sum(sib[,1] & rowSums(sib[,2:nb]))
		Ft[nb] <- sum(sib[,nb] & rowSums(sib[,1:(nb-1)]))
		FL[1] <- sum(sib[,1] & !rowSums(sib[,2:nb]))
		FL[nb] <- sum(sib[,nb] & !rowSums(sib[,1:(nb-1)]))
		for(j in 2:(nb-1)){
			bt[j] <- sum(rowSums(as.matrix(sib[,1:(j-1)])) & rowSums(as.matrix(sib[,(j+1):nb])))
			bL[j] <- sum(sib[,j] & !rowSums(as.matrix(sib[,1:(j-1)])) & rowSums(as.matrix(sib[,(j+1):nb])))
			Ft[j] <- sum(sib[,j] & rowSums(as.matrix(sib[,1:(j-1)])) & !rowSums(as.matrix(sib[,(j+1):nb])))
			FL[j] <- sum(sib[,j] & !rowSums(as.matrix(sib[,1:(j-1)])) & !rowSums(as.matrix(sib[,(j+1):nb])))
		}
		return(data.frame(bt,bL,Ft,FL))
		}
		
	rate <- function(sib, nb, rm){
		if(rm=="foote"){
			a <- crossers(sib, nb)
			return(list(E=-log(a$bt/(a$bL+a$bt)),O=-log(a$bt/(a$Ft+a$bt))))
			}else if(rm%in%c("three-timer","3-timer")){
				t2=t3=pt=rep(NA,nb)
				for(i in 1:(nb-1)){
					t2b[i] = colSums(sib[,i+1]&sib[,i])
					if(i>1){
						t2t[i] = colSums(sib[,i]&sib[i-1])
						t3[i] = colSums(sib[,i+1]&sib[,i]&sib[i-1])
						pt[i] = colSums(sib[,i+1]&sib[,i-1]&!sib[,i])
						}
					}
				return(list(E=log(t2b/(t3+pt)),O=log(t2t/(t3+pt))))
				}else if(rm=="gap-filler"){
					t2=t3=pt=g=d=rep(NA,nb)
					for(i in 1:(nb-1)){
						t2b[i] = colSums(sib[,i+1]&sib[,i])
						if(i>1){
							t2t[i] = colSums(sib[,i]&sib[i-1])
							t3[i] = colSums(sib[,i+1]&sib[,i]&sib[i-1])
							pt[i] = colSums(sib[,i+1]&sib[,i-1]&!sib[,i])
							}
						if(i>2) g[i] = colSums(sib[,i+1]&sib[,i-2]&!sib[,i-1])
						if(i<(ncol-2)) d[i] = colSums(sib[,i-1]&sib[,i+2]&!sib[,i+1])
						}
					return(list(E=log((t2b+pt)/(t3+pt+g)),O=log((t2t+pt)/(t3+pt+d))))
					}else{
						stop("Argument rate_method accept only three methods:\nfoote : see Foote (2000)\nthree-timer : see Alroy (2008)\ngap-filler : see Alroy (2014)")
					}
			}
	
	if(pac_top | pac_bottom){
		sp <- unique(id)
		id2 <- age_id2 <- c()
		for(i in sp){
			occ <- sort(age_id[id==i])
			n <- length(occ)
			nb_top <- floor(n*pac_top/100)
			nb_bot <- floor(n*pac_bottom/100)
			id2 <- c(id2, rep(i, n-nb_top-nb_bot))
			age_id2 <- c(age_id2, occ[(nb_top+1):(n-nb_bot)])
			}
		id <- id2
		age_id <- age_id2
		}
		
	sp <- unique(id)
	nsp <- length(sp)
	binned <- cut(age_id,bins)
	res <- as.data.frame(gsub("\\(|\\]","",do.call(rbind,strsplit(as.character(levels(binned)),","))))
	colnames(res)<-c("bin_top", "bin_bottom")
	occbinned <- split(id,binned)
	rawsp <- sapply(occbinned,function(x)length(unique(x)))
	nb <- length(rawsp)
	nocc <- sapply(occbinned, length)
	res$nb_occurrences <- nocc
	if(any(c("UW","OW","O2W")%in%method)){
		collbinned <- split(sample_id,binned)
		ncoll <- sapply(collbinned, length)
		res$nb_collections <- ncoll
		}
	res$nb_species <- rawsp
	if("SQS"%in%method){
		ab <- lapply(occbinned,table)
		singleton <- sapply(ab,function(x)sum(x==1))
		mostfrequent <- lapply(ab,function(x)names(x)[which.max(x)])
		u <- 1-singleton/nocc
		freq <- lapply(ab,function(x)x/sum(x))
		res$good_u <- u
		}
		
	sib <- t(do.call(rbind,lapply(occbinned,function(x){y=rep(0,nsp);y[factor(x,sp)]=1;y})))
	a <- crossers(sib,nb)
	res$rt_raw <- rangethrough(a, rt_method)
	r <- rate(sib,nb,rate_method)
	res$ext_raw <- r$E
	res$orig_raw <- r$O
	
	if("SQS"%in%method){
		q <- quotas[1]
		sp_sqs <- array(NA,dim=c(trials,nb))
		rt_sqs <- ext_sqs <- orig_sqs <- array(NA,dim=c(trials,nb))
		for(t in trials){
			seen <- array(0,dim=c(nsp,nb))
			for(i in seq_len(nb)){
				pool <- occbinned[[i]]
				if(length(pool)){
					sumfreq <- 0
					if(q <= u[i]){
						while(sumfreq<q){
							nib <- length(pool)
							x <- sample(1:nib,1)
							if(!seen[factor(pool[x],sp),i]){
								if(pool[x]!=mostfrequent[i]|dominant=="include"){
									sumfreq <- sumfreq + freq[[i]][names(freq[[i]])==pool[x]]
									seen[factor(pool[x],sp),i] <- 1
									pool <- pool[pool!=pool[x]]
									}
								}
							}
						}
					}
				}
			sqs <- crossers(seen,nb)
			sp_sqs[t,] <- colSums(seen)
			rt_sqs[t,] <- rangethrough(sqs, rt_method)
			r <- rate(seen,nb,rate_method)
			ext_sqs[t,] <- r$E
			orig_sqs[t,] <- r$O
			}
		res$sp_sqs <- colMeans(sp_sqs,na.rm=TRUE)
		res$rt_sqs <- colMeans(rt_sqs,na.rm=TRUE)
		res$ext_sqs <- colMeans(ext_sqs,na.rm=TRUE)
		res$orig_sqs <- colMeans(orig_sqs,na.rm=TRUE)
		}
		
	if("UW"%in%method){
		q <- quotas[2]
		for(t in trials){
			seen <- array(0,dim=c(nsp,nb))
			for(i in seq_len(nb)){
				pool <- collbinned[[i]]
				if(length(pool)){
					collsamp <- 0
					if(q <= ncoll[i]){
						while(collsamp<q){
							nleft <- length(pool)
							x <- sample(1:nleft,1)
							seen[factor(id[sample_id==pool[x]],sp),i] <- 1
							pool <- pool[pool!=pool[x]]
							collsamp <- collsamp + 1
							}
						}
					}
				}
			uw <- crossers(seen)
			sp_uw[t,] <- rowSums(seen)
			rt_uw[t,] <- rangethrough(uw, rt_method)
			r <- rate(seen,nb,rate_method)
			ext_uw[t,] <- r$E
			orig_uw[t,] <- r$O
			}
		res$sp_uw <- colMeans(sp_uw,na.rm=TRUE)
		res$rt_uw <- colMeans(rt_uw,na.rm=TRUE)
		res$ext_uw <- colMeans(ext_uw,na.rm=TRUE)
		res$orig_uw <- colMeans(orig_uw,na.rm=TRUE)
		}
		
	if("OW"%in%method){
		q <- quotas[3]
		for(t in trials){
			seen <- array(0,dim=c(nsp,nb))
			for(i in seq_len(nb)){
				pool <- collbinned[[i]]
				if(length(pool)){
					occsamp <- 0
					if(q <= nocc[i]){
						while(occsamp<q){
							nleft <- length(pool)
							x <- sample(1:nleft,1)
							seen[factor(id[sample_id==pool[x]],sp),i] <- 1
							pool <- pool[pool!=pool[x]]
							occsamp <- occsamp + sum(sample_id==pool[x])
							}
						}
					}
				}
			ow <- crossers(seen)
			sp_ow[t,] <- rowSums(seen)
			rt_ow[t,] <- rangethrough(ow, rt_method)
			r <- rate(seen,nb,rate_method)
			ext_ow[t,] <- r$E
			orig_ow[t,] <- r$O
			}
		res$sp_ow <- colMeans(sp_ow,na.rm=TRUE)
		res$rt_ow <- colMeans(rt_ow,na.rm=TRUE)
		res$ext_ow <- colMeans(ext_ow,na.rm=TRUE)
		res$orig_ow <- colMeans(orig_ow,na.rm=TRUE)
		}
		
	if("O2W"%in%method){
		q <- quotas[4]
		occ2 <- sapply(collbinned,function(x)sum(sapply(x,function(y)sum(sample_id==y)^2)))
		for(t in trials){
			seen <- array(0,dim=c(nsp,nb))
			for(i in seq_len(nb)){
				pool <- collbinned[[i]]
				if(length(pool)){
					occsamp <- 0
					if(q <= occ2[i]){
						while(occsamp<q){
							nleft <- length(pool)
							x <- sample(1:nleft,1)
							seen[factor(id[sample_id==pool[x]],sp),i] <- 1
							pool <- pool[pool!=pool[x]]
							occsamp <- occsamp + sum(sample_id==pool[x])^2
							}
						}
					}
				}
			o2w <- crossers(seen)
			sp_o2w[t,] <- rowSums(seen)
			rt_o2w[t,] <- rangethrough(o2w, rt_method)
			r <- rate(seen,nb,rate_method)
			ext_o2w[t,] <- r$E
			orig_o2w[t,] <- r$O
			}
		res$sp_o2w <- colMeans(sp_o2w,na.rm=TRUE)
		res$rt_o2w <- colMeans(rt_o2w,na.rm=TRUE)
		res$ext_o2w <- colMeans(ext_o2w,na.rm=TRUE)
		res$orig_o2w <- colMeans(orig_o2w,na.rm=TRUE)
		}
		
	if("CR"%in%method){
		q <- quotas[5]
		for(t in trials){
			seen <- array(0,dim=c(nsp,nb))
			for(i in seq_len(nb)){
				if(q <= nocc[i]){
					seen[factor(sample(occbinned[[i]],q),sp),i] <- 1
					}
				}
			cr <- crossers(seen)
			sp_cr[t,] <- rowSums(seen)
			rt_cr[t,] <- rangethrough(cr, rt_method)
			r <- rate(seen,nb,rate_method)
			ext_cr[t,] <- r$E
			orig_cr[t,] <- r$O
			}
		res$sp_cr <- colMeans(sp_cr,na.rm=TRUE)
		res$rt_cr <- colMeans(rt_cr,na.rm=TRUE)
		res$ext_cr <- colMeans(ext_cr,na.rm=TRUE)
		res$orig_cr <- colMeans(orig_cr,na.rm=TRUE)
		}
		
	if(write)write.table(res,file="Results.csv",sep="\t",row.names=FALSE,col.names=TRUE)
	return(res)
	}