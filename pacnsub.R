# Subsampling Neptune database (with pacman trimming option)
# author of the implementation: Johan Renaudie.
# authors of the algorithms: Lazarus et al 2012 for pacman, Sanders 1968 for rarefaction, Shinozaki 1963
# for UW, Alroy 1996 for OW, Alroy 2000 for O2W and Alroy 2010 for SQS.

# Make sure the neptune file contains the following 4 columns: 
# taxon_id (can be an actual id or a species name or whatever but the column must be named taxon_id), 
# sample_age_ma (the age of the sample), hole_id and sample_depth_mbsf (both site and depth_mbsf are 
# used are sample identifier here).
# 
# sub_method can be "all", or any combinations of "SQS","UW","OW","O2W" and "CR" (classical rarefaction)
# sub_quota are the quota used for each subsampling procedures. If sub_method is "all", the quotas have 
# to be in the order sqs-uw-ow-o2w-cr, otherwise the order given in the sub_method vector. (e.g. if 
# sub_method=c("CR","OW") then sub_quota=c(100,50) where 100 is the quota for classical rarefaction and 
# 50 for occurence-weighted by-list subsampling)
# 
# trials is the number of trials for each subsampling
# dominant is only useful for SQS subsampling: "include" includes the dominant taxa in each bin, 
# "exclude" excludes it.
# pac_top and pac_bottom are the top and bottom trimming percentage for the pacman procedure
# bin_length is the length of the bins in Ma: only equal-length bins are supported here.
# age_min and age_max are the temporal boundary of the study (in Ma)
# write=TRUE allows an csv output of the main results to be created (in the default directory).
# rt_method is the method used for range-through. Can be either of 'bc' for boundary-crossers, 'rt' for 
# classic range-through (without singletons) and 'tot' (with singletons).

pacnsub<-function(neptune, write=TRUE,
					sub_method="all", sub_quota=c(0.6,10,100,5000,100), trials=100, dominant="include",
					pac_top=5, pac_bottom=3, rt_method="bc",
					bin_length=0.5, age_min=0, age_max=65){
	
	#Sort out subsampling methods
	c("SQS","UW","OW","O2W","CR") -> methods
	if("all"%in%sub_method){sub_method<-methods}
	quotas <- numeric(5)
	if(length(sub_quota)==length(sub_method)){
		sub_quota -> quotas[methods%in%sub_method]
		}else{stop("Please provide a quota for each method (order if all: SQS, UW, OW, O2W and CR)")}
	
	#Takes care of binning
	bin <- seq(age_min,age_max,by=bin_length)     #set bins
	length(bin)-1 -> vl
	binmid <- numeric(vl)                 		#set the midpoints of each bin (for plotting purposes)
	for(i in 1:vl){binmid[i] <- (bin[i]+bin[i+1])/2}
	
	rangethrough<-function(bL, bt, Ft, FL, rtm){
		if(!rtm%in%c("bc","rt","tot")){stop("Argument rt_method accept only three methods:\nbc : Standing mean diversity based on boundary-crossers\nrt : Classic range-through excluding singletons\ntot : Classic range-through including singletons")}
		if(rtm=="bc"){return(((bL+bt)+(Ft+bt))/2)}
		if(rtm=="rt"){return(bL+bt+Ft)}
		if(rtm=="tot"){return(bL+bt+Ft+FL)}
		}
	
	#Pacman
	pacman<-function(neptune,perc_top,perc_bottom){
		unique(neptune$taxon_id)->species
		sp<-list();noc<-c()
			unique(data.frame(neptune$hole_id,neptune$sample_depth_mbsf))->samples
		outlier<-rep(0,nrow(samples))
		res<-c()
		for(i in 1:length(species)){
			neptune[neptune$taxon_id==species[i],]->sp[[i]]
			sp[[i]][order(sp[[i]]$sample_age_ma),]->sp[[i]]
			noc[i]<-nrow(sp[[i]])
			sp[[i]]->spe
			floor(noc[i]*perc_top/100)->nb_top
			floor(noc[i]*perc_bottom/100)->nb_bottom
			if(nb_top>=1){
				for(k in 1:nb_top){
					samples[,1]==spe$hole_id[k]&samples[,2]==spe$sample_depth_mbsf[k]->out
					outlier[out]<-outlier[out]+1
					}
				spe[-(1:nb_top),]->spe
				}
			if(nb_bottom>=1){
				for(k in ((nrow(spe)-nb_bottom+1):nrow(spe))){
					samples[,1]==spe$hole_id[k]&samples[,2]==spe$sample_depth_mbsf[k]->out
					outlier[out]<-outlier[out]+1
					}
				spe[-((nrow(spe)-nb_bottom+1):nrow(spe)),]->spe}
			rbind(res,spe)->res
			Sys.sleep(.01)      # some loop operations
			cat(i, paste(" of ", length(species)," species\r",sep="")) 
			flush.console()
			}
		cbind(samples,outlier)->outlier
		colnames(neptune)->colnames(res)
		as.data.frame(res)
		colnames(outlier)<-c("hole_id","Depth","Nb Outliers")
		list("After trimming"=res,"Outlier per sample"=outlier)->reslist
		return(reslist)
		}
	if(pac_top!=0 | pac_bottom!=0){
		cat("Pacman trimming :\n")
		pacman(neptune,pac_top,pac_bottom)->PAC
		PAC[[1]]->neptune
		PAC[[2]]->outliers
		}
	else{outliers<-NA}

	#OCCURENCE & COLLECTION BINNING
	
	neptune$taxon_id <- as.factor(neptune$taxon_id)
	nlevels(neptune$taxon_id)->nsp
	
	occ<-numeric(vl)                                  #Nb of Occurences per bin
	rawsp<-numeric(vl)                                #Raw species count per bin
	occbinned<-list()                                 #List of Occurences per bin
	ncoll<-numeric(vl)                                #Nb of Collections per bin
	collbinned<-list()                                #List of Collections per bin
	mostfrequent<-numeric(vl)					#Most abundant species in the bin
	ab<-list()									#Binned abundances
	freq<-list()								 #Binned frequencies
	u<-numeric(vl)								  #Good's u
	single<-numeric(vl)							  #Number of singletons
	for (i in 1:vl){
		occbinned[[i]]<-subset(neptune,bin[i]<sample_age_ma & bin[i+1]>sample_age_ma)
		occ[i]<-nrow(occbinned[[i]])
		rawsp[i]<-length(unique(occbinned[[i]]$taxon_id))
		unique(data.frame(occbinned[[i]]$hole_id,occbinned[[i]]$sample_depth_mbsf))->collbinned[[i]]
		collbinned[[i]][order(collbinned[[i]][,1]),]->collbinned[[i]]
		nrow(collbinned[[i]])->ncoll[i]
		summary(occbinned[[i]]$taxon_id,maxsum=1000)->ab[[i]]
		ab[[i]][ab[[i]]!=0]->ab[[i]]
		single[i]<-sum(ab[[i]]==1)
		u[i]<- 1- single[i]/occ[i]
		if(length(ab[[i]])!=0){
			mf<-names(ab[[i]])[ab[[i]]==max(ab[[i]],na.rm=T)]
			if(length(mf)==1){mostfrequent[i]<-mf}
			}
		freq[[i]]<-ab[[i]]/occ[i]
		}
	
	result<-list()
	result$'Trimmed dataset'<-neptune
	result$'Outliers'<-outliers
	
	if(write){
		write.table(result[[1]],file="Trimmed dataset.csv",sep=",",dec=".")
		write.table(result[[2]],file="Outliers.csv",sep=",",dec=".",row.names=FALSE)
		}
	
	#Raw diversity measurements 
	sib<-array(0,dim=c(nsp,vl))        #Occurence matrix
	neptune[!is.na(neptune$sample_age_ma),]->neptune
	for(j in 1:vl){sib[neptune[neptune$sample_age_ma > bin[j] & neptune$sample_age_ma < bin[j+1],]$taxon_id,j]<-1}
	sib[rowSums(sib)!=0,]->sib
	bL<-bt<-Ft<-FL<-bcRaw<-numeric(vl)
	bL[1] <- sum(sib[,1]==1 & rowSums(sib[,2:vl])>0)
	Ft[vl] <- sum(sib[,vl]==1 & rowSums(sib[,1:(vl-1)])>0)
	FL[1] <- sum(sib[,1]==1 & rowSums(sib[,2:vl])==0)
	FL[vl] <- sum(sib[,vl]==1 & rowSums(sib[,1:(vl-1)])==0)
	for(j in 2:(vl-1)){
		bt[j] <- sum(rowSums(as.matrix(sib[,1:(j-1)]))>0 & rowSums(as.matrix(sib[,(j+1):vl]))>0)
		bL[j] <- sum(sib[,j]==1 & rowSums(as.matrix(sib[,1:(j-1)]))==0 & rowSums(as.matrix(sib[,(j+1):vl]))>0)
		Ft[j] <- sum(sib[,j]==1 & rowSums(as.matrix(sib[,1:(j-1)]))>0 & rowSums(as.matrix(sib[,(j+1):vl]))==0)
		FL[j] <- sum(sib[,j]==1 & rowSums(as.matrix(sib[,1:(j-1)]))==0 & rowSums(as.matrix(sib[,(j+1):vl]))==0)
		}
	bcRaw<-rangethrough(bL,bt,Ft,FL,rt_method)
	rateRaw<- -log(bt/(bL+bt))
	origRaw<- -log(bt/(Ft+bt))
	rbind("Raw Species Richness"=rawsp,"Boundary Crossers"=bcRaw,"Extinction Rate"=rateRaw,"Origination Rate"=origRaw, "Number of Occurences"=occ, "Number of Collections"=ncoll, "Good's u"=u)->result$Raw
	
	if(write){
		rn<-c()
		for(i in 1:vl){rn[i]<-paste(bin[i],bin[i+1],sep=" to ")}
		resmat <- data.frame(`Mid-bin`=binmid, `Species Richness`=rawsp, `Range-Through`=bcRaw, 
							`Number of occurrences`=occ, `Number of collections`=ncoll, `Good's u`=u)
		rn->rownames(resmat)
		}
	
	if("SQS"%in%sub_method){
		quotas[1]->q
		cat("\nSQS :\n")
		SQSspecies<-array(0,dim=c(trials,vl))	#In-bin number of subsampled species
		bcSQS<-array(0,dim=c(trials,vl))		#Boundary crossers after subsampling
		rateSQS<-array(0,dim=c(trials,vl))		#Extinction rate
		origSQS<-array(0,dim=c(trials,vl))		#Origination rate
		for(t in 1:trials){
			seen<-array(0,dim=c(nsp,vl))         #Occurence matrix after subsampling
			for(i in 1:vl){
				if(nrow(occbinned[[i]])!=0){
					pool<-occbinned[[i]]                                           #'Picking jar'
					left<-nrow(pool)                                               #Nb of occurences left to be picked
					sumfreq<-0		                                               #Cumulative frequencies
					if(q<=u[i]){
						while(sumfreq<q){
							x<-floor(runif(1,min=1,max=left+1))
							if(seen[as.integer(pool$taxon_id[x]),i]==0){
								if(pool$taxon_id[x]!=mostfrequent[i] | dominant=="include"){sumfreq<-sumfreq+as.numeric(freq[[i]][names(freq[[i]])==pool$taxon_id[x]])}
								seen[as.integer(pool$taxon_id[x]),i]<-1
								SQSspecies[t,i]<-SQSspecies[t,i]+1
								}
							pool[x,]<-pool[left,]
							left<-left-1
							}
						}
					}
				}
			bL<-bt<-Ft<-FL<-numeric(vl)
			bL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])>0)
			Ft[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])>0)
			FL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])==0)
			FL[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])==0)
			for(j in 2:(vl-1)){
				bt[j] <- sum(rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				bL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				Ft[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
				FL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
				}
			bcSQS[t,]<-rangethrough(bL,bt,Ft,FL,rt_method)
			rateSQS[t,]<- -log(bt/(bL+bt))
			origSQS[t,]<- -log(bt/(Ft+bt))
			Sys.sleep(.01)      # some loop operations
			cat(t, paste(" of ", trials," trials\r",sep="")) 
			flush.console()
			}
		bcSQSM<-SQSspeciesM<-rateSQSM<-origSQSM<-numeric(vl)
		for(i in 1:vl){
			mean(bcSQS[,i],na.rm=T)->bcSQSM[i]
			mean(rateSQS[,i],na.rm=T)->rateSQSM[i]
			mean(origSQS[,i],na.rm=T)->origSQSM[i]
			mean(SQSspecies[,i],na.rm=T)->SQSspeciesM[i]
			}
		rbind("SQS subsampled diversity"=SQSspeciesM,"Boundary Crossers on SQS"=bcSQSM)->sqsmat                        #Output matrix 
		for(i in 1:vl){if(sqsmat[2,i]==0){sqsmat[2,i]<-NA}}
		result$SQS<-list("SQS Subsampling summary"=sqsmat,
						"SQS subsampling value for 100 trials"=SQSspecies, "SQS Boundary crossers in 100 trials"=bcSQS, 
						"Extinction rate"=rateSQSM, "Origination rate"=origSQSM)
		if(write){
			resmat$`SQS diversity`<-SQSspeciesM
			resmat$`SQS range-through`<-bcSQSM
			resmat$`SQS extinction`<-rateSQSM
			resmat$`SQS origination`<-origSQSM
			}
		}
	
	if("UW"%in%sub_method){
		quotas[2]->quota
		cat("\nUW :\n")
		bcUW<-uwSpecies<-rateUW<-origUW<-array(0,dim=c(trials,vl))
		for(t in 1:trials){
			seen<-array(0,dim=c(nsp,vl))
			for(i in 1:vl){
				temp<-collbinned[[i]]
				collsamp<-0                                                     #Nb of Collections already picked
				uwsampled<-c()                                                  #Picked occurences cumulative storage
				ncollLeft<-ncoll[i]
				if(ncoll[i]<quota){uwSpecies[t,i]<-NA}
				else{
					while(collsamp<quota){
						cellNo<-floor(runif(1,min=1,max=ncollLeft+1))           #Pick randomly a sample out of the remaining samples
						neptune[(neptune$hole_id==temp[cellNo,1] & neptune$sample_depth_mbsf==temp[cellNo,2]),]->picked
						rbind(uwsampled,picked)->uwsampled
						seen[picked$taxon_id,i]<-1                            #Fill in the subsampled occurence matrix
						temp[cellNo,]<-temp[ncollLeft,]                         #To avoid picking the same cell again and again
						collsamp<-collsamp+1
						ncollLeft<-ncollLeft-1
					}
					uwSpecies[t,i]<-length(unique(uwsampled$taxon_id))        #Number of taxa picked
				}
			}
			bL<-bt<-Ft<-FL<-numeric(vl)
			bL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])>0)
			Ft[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])>0)
			FL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])==0)
			FL[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])==0)
			for(j in 2:(vl-1)){
				bt[j] <- sum(rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				bL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				Ft[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
				FL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
			}
			bcUW[t,]<-rangethrough(bL,bt,Ft,FL,rt_method)                                             #Boundary Crossers (standing mean diversity)
			rateUW[t,]<- -log(bt/(bL+bt))
			origUW[t,]<- -log(bt/(Ft+bt))
			Sys.sleep(.01)      # some loop operations
			cat(t, paste(" of ", trials," trials\r",sep="")) 
			flush.console()
		}
		uwSpeciesM<-bcUWM<-rateUWM<-origUWM<-numeric(vl)                                    #mean richness & mean bc on t trials 
		for(i in 1:vl){
			mean(uwSpecies[,i],na.rm=T)->uwSpeciesM[i]
			mean(bcUW[,i],na.rm=T)->bcUWM[i]
			mean(rateUW[,i],na.rm=T)->rateUWM[i]
			mean(origUW[,i],na.rm=T)->origUWM[i]
			}
		rbind("UW subsampled diversity"=uwSpeciesM,"Boundary Crossers on UW"=bcUWM)->uwmat                        #Output matrix 
		for(i in 1:vl){if(uwmat[2,i]==0){uwmat[2,i]<-NA}}
		result$UW<-list("UW Subsampling summary"=uwmat,
						"UW subsampling value for 100 trials"=uwSpecies, "UW Boundary crossers in 100 trials"=bcUW, 
						"Extinction rate"=rateUWM, "Origination rate"=origUWM)
		if(write){
			resmat$`UW diversity`<-uwSpeciesM
			resmat$`UW range-through`<-bcUWM
			resmat$`UW extinction`<-rateUWM
			resmat$`UW origination`<-origUWM
			}
		}
		
	if("OW"%in%sub_method){
		quotas[3]->quota
		cat("\nOW :\n")
		bcOW<-owSpecies<-rateOW<-origOW<-array(0,dim=c(trials,vl))
		for(t in 1:trials){
			seen<-array(0,dim=c(nsp,vl))
			for(i in 1:vl){
				temp<-collbinned[[i]]
				occsamp<-0
				owsampled<-c()
				ncollLeft<-ncoll[i]
				if(occ[i]<quota){owSpecies[t,i]<-NA}
				else{
					while(occsamp<quota){
						cellNo<-floor(runif(1,min=1,max=ncollLeft+1))
						neptune[(neptune$hole_id==temp[cellNo,1] & neptune$sample_depth_mbsf==temp[cellNo,2]),]->picked
						rbind(owsampled,picked)->owsampled
						seen[picked$taxon_id,i]<-1
						temp[cellNo,]<-temp[ncollLeft,]
						occsamp<-occsamp+nrow(picked)
						ncollLeft<-ncollLeft-1
					}
					owSpecies[t,i]<-length(unique(owsampled$taxon_id))
				}
			}
			bL<-bt<-Ft<-FL<-numeric(vl)
			bL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])>0)
			Ft[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])>0)
			FL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])==0)
			FL[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])==0)
			for(j in 2:(vl-1)){
				bt[j] <- sum(rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				bL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				Ft[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
				FL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
			}
			bcOW[t,]<-rangethrough(bL,bt,Ft,FL,rt_method)  
			rateOW[t,]<- -log(bt/(bL+bt))
			origOW[t,]<- -log(bt/(Ft+bt))
			Sys.sleep(.01)      # some loop operations
			cat(t, paste(" of ", trials," trials\r",sep="")) 
			flush.console()
			}
		owSpeciesM<-bcOWM<-rateOWM<-origOWM<-numeric(vl)
		for(i in 1:vl){
			mean(owSpecies[,i],na.rm=T)->owSpeciesM[i]
			mean(bcOW[,i],na.rm=T)->bcOWM[i]
			mean(rateOW[,i],na.rm=T)->rateOWM[i]
			mean(origOW[,i],na.rm=T)->origOWM[i]
			}
		rbind("OW subsampled diversity"=owSpeciesM,"Boundary Crossers on OW"=bcOWM)->owmat
		for(i in 1:vl){if(owmat[2,i]==0){owmat[2,i]<-NA}}
		result$OW<-list("OW Subsampling summary"=owmat,
						"OW subsampling value for 100 trials"=owSpecies, "OW Boundary crossers in 100 trials"=bcOW, 
						"Extinction rate"=rateOWM, "Origination rate"=origOWM)
		if(write){
			resmat$`OW diversity`<-owSpeciesM
			resmat$`OW range-through`<-bcOWM
			resmat$`OW extinction`<-rateOWM
			resmat$`OW origination`<-origOWM
			}
	}
	
	if("O2W"%in%sub_method){
		quotas[4]->quota
		cat("\nO2W :\n")
		bc3<-ow2Species<-rateO2W<-origO2W<-array(0,dim=c(trials,vl))
		occ2<-numeric(vl)                   	# Number of occurence-squared per bin
		for(i in 1:vl){
			for(k in 1:nrow(collbinned[[i]])){
				occ2[i]+nrow(subset(neptune,hole_id==collbinned[[i]][k,1] & sample_depth_mbsf==collbinned[[i]][k,2]))^2->occ2[i]
			}
		}
		for(t in 1:trials){
			seen<-array(0,dim=c(nsp,vl))
			for(i in 1:vl){
				temp<-collbinned[[i]]
				spsamp<-0
				ow2sampled<-c()
				ncollLeft<-ncoll[i]
				if(occ2[i]<quota){ow2Species[t,i]<-NA}
				else{
					while(spsamp<quota){
						cellNo<-floor(runif(1,min=1,max=ncollLeft+1))
						subset(neptune,hole_id==temp[cellNo,1] & sample_depth_mbsf==temp[cellNo,2])->selected
						rbind(ow2sampled,selected)->ow2sampled
						seen[selected$taxon_id,i]<-1
						temp[cellNo,]<-temp[ncollLeft,]
						spsamp<-spsamp+nrow(selected)^2
						ncollLeft<-ncollLeft-1
					}
					ow2Species[t,i]<-length(unique(ow2sampled$taxon_id))
				}
			}
			bL<-bt<-Ft<-FL<-numeric(vl)
			bL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])>0)
			Ft[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])>0)
			FL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])==0)
			FL[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])==0)
			for(j in 2:(vl-1)){
				bt[j] <- sum(rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				bL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				Ft[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
				FL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
				}
			bc3[t,]<-rangethrough(bL,bt,Ft,FL,rt_method)  
			rateO2W[t,]<- -log(bt/(bL+bt))
			origO2W[t,]<- -log(bt/(Ft+bt))
			Sys.sleep(.01)      # some loop operations
			cat(t, paste(" of ", trials," trials\r",sep="")) 
			flush.console()
			}
		ow2SpeciesM<-bc3M<-rateO2WM<-origO2WM<-numeric(vl)
		for(i in 1:vl){
			mean(ow2Species[,i],na.rm=T)->ow2SpeciesM[i]
			mean(bc3[,i],na.rm=T)->bc3M[i]
			mean(rateO2W[,i],na.rm=T)->rateO2WM[i]
			mean(origO2W[,i],na.rm=T)->origO2WM[i]
			}
		rbind("OW2 subsampled diversity"=ow2SpeciesM,"Boundary Crossers on OW2"=bc3M)->ow2mat
		for(i in 1:vl){if(ow2mat[2,i]==0){ow2mat[2,i]<-NA}}
		result$O2W<-list("OW2 Subsampling summary"=ow2mat,
						"OW2 subsampling value for 100 trials"=ow2Species, "OW2 Boundary crossers in 100 trials"=bc3, 
						"Extinction rate"=rateO2WM, "Origination rate"=origO2WM)
		if(write){
			resmat$`OW2 diversity`<-ow2SpeciesM
			resmat$`OW2 range-through`<-bc3M
			resmat$`OW2 extinction`<-rateO2WM
			resmat$`OW2 origination`<-origO2WM
			}
		}
		
	if("CR"%in%sub_method){
		quotas[5]->quota
		cat("\nCR :\n")
		rarefied<-bc4<-rateRar<-origRar<-array(dim=c(trials,vl))
		for(t in 1:trials){
			seen<-array(0,dim=c(nsp,vl))
			for(i in 1:vl){
				if(occ[i]<quota){rarefied[t,i]<-NA}
				else{
					occbinned[[i]][sample(nrow(occbinned[[i]]),quota),]->selected
					length(unique(selected$taxon_id))->rarefied[t,i]
					seen[unique(selected$taxon_id),i]<-1
					}
				}
			bL<-bt<-Ft<-FL<-numeric(vl)
			bL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])>0)
			Ft[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])>0)
			FL[1] <- sum(seen[,1]==1 & rowSums(seen[,2:vl])==0)
			FL[vl] <- sum(seen[,vl]==1 & rowSums(seen[,1:(vl-1)])==0)
			for(j in 2:(vl-1)){
				bt[j] <- sum(rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				bL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))>0)
				Ft[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))>0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
				FL[j] <- sum(seen[,j]==1 & rowSums(as.matrix(seen[,1:(j-1)]))==0 & rowSums(as.matrix(seen[,(j+1):vl]))==0)
			}
			bc4[t,]<-rangethrough(bL,bt,Ft,FL,rt_method)
			rateRar[t,]<- -log(bt/(bL+bt))
			origRar[t,]<- -log(bt/(Ft+bt))
			Sys.sleep(.01)      # some loop operations
			cat(t, paste(" of ", trials," trials\r",sep="")) 
			flush.console()
		}
		rarefiedM<-bc4M<-rateRarM<-origRarM<-numeric(vl)
		for(i in 1:vl){
			mean(rarefied[,i],na.rm=T)->rarefiedM[i]
			mean(bc4[,i],na.rm=T)->bc4M[i]
			mean(rateRar[,i],na.rm=T)->rateRarM[i]
			mean(origRar[,i],na.rm=T)->origRarM[i]
			}
		rbind("Rarefied diversity (100 specimens)"=rarefiedM,"Boundary Crossers"=bc4M)->Rarmat
		for(i in 1:vl){if(Rarmat[2,i]==0){Rarmat[2,i]<-NA}}
		result$CR<-list("Rarefaction summary"=Rarmat,
						"Classic Rarefaction value for 100 trials"=rarefied, "Boundary crossers in 100 trials"=bc4, 
						"Extinction rate"=rateRarM, "Origination rate"=origRarM)
		if(write==TRUE){
			resmat$`CR diversity`<-rarefiedM
			resmat$`CR range-through`<-bc4M
			resmat$`CR extinction`<-rateRarM
			resmat$`CR origination`<-origRarM
			}
		}
		
	if(write){write.table(resmat,file="Results.csv",sep=",",dec=".", row.names=F)}
	return(result)
	}