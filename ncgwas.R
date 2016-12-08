#### contact baldassa@email.unc.edu for questions

library(data.table)
library(ncdf4)
library(RcppEigen)
library(parallel)
library(speedglm)

############ start of inputs ####################################################################################

#The objective is to specify the inputs from a bash/tcsh script, to avoid
#Modifying this file 

phenodir <- "/nas02/depts/epi/CVDGeneNas/antoine/ECG_GWAS/WHI/phenotypes/" #phenotype file directory
resdir <- "/proj/epi/CVDGeneNas/antoine/dev_garnetmpi/test_results/" #where to store results
gpath <- "/nas02/depts/epi/Genetic_Data_Center/whi_share/whi_1000g_fh_imp/ncdf-data/" #where the 1KG .nc files live
study <- "GARNET" #the WHI study
outcome <- "jt" #The outcome of interest
form <- ~g+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+region+rr_d+age
pheno <- "ecg_whi_whites_fi.txt" #Name of the .txt phenotype file
type <- "linear" #Right now just "linear" or "logistic"

family <- "" #Optional: family for GLM -- default is binomial. 
link <- "" #Optional: GLM link function -- default is family default (e.g. logit for binomial)
keepchr <- 22 #Chromosomes to restrict to

############# end of inputs ######################################################################################

#Load and pare down data
Epidata <- fread(paste0(phenodir,pheno))
Epidata[,c("g"):=rnorm(nrow(Epidata))]
Epidata <- na.omit(Epidata[,c("id", outcome, all.vars(form)), with = F])
setnames(Epidata, "id", "Common_ID")
setkey(Epidata, Common_ID)

#order ids as in NC file
nc22 <- nc_open(paste0(gpath,study,'-chr',22,'-c.nc'))
ncids <- data.table(Common_ID = ncvar_get(nc22, "Common_ID"), 
                    ncid = seq_len(nc22$dim$Samples$len))
setkey(ncids,Common_ID)
dt_ana <- ncids[Epidata, nomatch = 0]
setkey(dt_ana, ncid)

#Record NCDF indices of participants with phenotype data
nckeep <- dt_ana$ncid

#Create X and y model matrices, for LMs
X <- as.matrix(dt_ana[,all.vars(form),with =FALSE][,int:=1])
y <- dt_ana[,get(outcome)]

#Finalize analytical datasets, for GLMs
dt_ana <- dt_ana[,c(outcome, all.vars(form)),with=FALSE]

#Generate fit function from file inputs. Workhorse functions are
#RcppEigen::fastLmPure for linear models and speedglm::speedglm 
#for generalized linear models
qfit_setup <- function(g, d_type, family = "binomial", link = ""){
  if(d_type == "linear"){
    function(gnow){
      ind <- which(!is.na(gnow))
      X[,g] <- gnow
      tm <- fastLmPure(X[ind,],y[ind])
      list(tm$coefficients[g], tm$se[g])
    }
  }else if(d_type == "glm"){
    if (family == "") family <- "binomial"
    function(gnow){
      ind <- which(!is.na(gnow))
      dt_ana[,g:=gnow]
      #Higher row.chunk will run faster when calculating XX', but will take more
      #memory. 2000 should be fine but change at will. An improvement will be 
      #to deduce XX' from previous iterations using algebra (should work, but
      #will take some code, maybe X2X2' = XX' :/ g :* g2 )
      tm <- speedglm(form,data=dt_ana[ind],
                     family=do.call(family, link),
                     set.default=list(row.chunk=2000))
      as.list(summary(tm)[g,-3])
    }
  }else stop("Specify type as linear or logistic")
}
qfit <- qfit_setup(match("g", all.vars(form))+1, type)

#Split function
splitup <- function(a, n) lapply(split(a[1]:a[2], cut(a[1]:a[2], n)), range)

#Send objects and libraries to worker threads
mpi.bcast.Robj2slave(all = TRUE) 
mpi.bcast.cmd({
  library(data.table); library(ncdf4)
  library(RcppEigen); library(parallel)
  library(speedglm)
})

#The number of workers, -1 because I think the master thread should not count
nworkers <- mpi.comm.size() - 1

#Data for debugging
print(paste("Universe size:", mpi.universe.size()))
print(paste("Comm size:", mpi.comm.size()))

#Start loop over chromosomes
for(i in keepchr){
  
  #misc: get #snps, make output file name, send chromosome # to workers...
  rname <- paste0(resdir,"Chr",i,"_",outcome,"_",study,"_results.csv")
  mpi.bcast.Robj2slave(i)
  nc <- nc_open(paste0(gpath,study,'-chr',i,'-c.nc'))
  nsnp <- nc$dim$SNPs$len
  print(paste("Starting on chromosome:",i,"at:",Sys.time()))
  
  #B
  parts <- splitup(c(1,nsnp), 10)
  for(p in parts){
    #Splitup task into indices for workers
    bits <- splitup(p,nworkers)
    res <- mpi.parLapply(bits, function(k) {
      #Open nc file, get SNP names and create results dataset
      nc <- nc_open(paste0(gpath,study,'-chr',i,'-c.nc'))
      snp_names <- ncvar_get(nc, "SNP_Name",c(1,k[1]), c(-1,k[2]-k[1]+1))
      res_part <- data.table(
        index = as.integer(k[1]:k[2]), 
        snp = snp_names, 
        coded = ncvar_get(nc,"Allele1_Reference", k[1], k[2]-k[1]+1), 
        other = ncvar_get(nc,"Allele2_Reference", k[1], k[2]-k[1]+1),
        caf = as.numeric(NA), b = as.numeric(NA),
        se = as.numeric(NA), p = as.numeric(NA), j = seq_along(k[1]:k[2]))
      
      #Read dosages at relevant indices, and restrict to participants also in phenotype file
      p_aa <- ncvar_get(nc,"Prob_AA", start=c(k[1],1), count=c(k[2]-k[1]+1, -1))[,nckeep]
      p_ab <- ncvar_get(nc,"Prob_AB", start=c(k[1],1), count=c(k[2]-k[1]+1, -1))[,nckeep]
      dos <- p_aa*2 + p_ab
      
      #Add allele frequency, variance and nonomissing N
      res_part[,c("caf", "v", "n") := list(mean(dos[j,]/2, na.rm = TRUE), 
                                           var(dos[j,], na.rm = TRUE),
                                           sum(!is.nan(dos[j,]), na.rm = TRUE)), j]
      
      #Add regression results using qfit over each column of the dosage matrix,
      #wrapped with the data.table by= operator. qfit is defined at the start of 
      #this file 
      
      if(type == "linear"){ res_part[n > 0 & v > 0, c("b", "se") := qfit(dos[j,]), j]
      }else res_part[n > 0 & v > 0, c("b", "se", "p") := qfit(dos[j,]), j]
      
      #Return data.table copy to avoid memory leaks
      copy(res_part)
    })
    #For debugging
    if(!is.data.table(res[[1]])) print(res)
    #warnings()
    
    #Combine worker outputs into single DT
    res <- rbindlist(res)
    
    #Get p-values for linear model (which don't already compute it)
    if(type == "linear") res[,p := 2*(1-pnorm(abs(b/se)))]
    
    #Add chromosome column, and remove j-index column
    res[,chr := i]
    res[,j := NULL]
    
    #Write to file
    if(p[1] == 1){ fwrite(res, rname, sep = ",")
    }else fwrite(res, rname, sep = ",", append = TRUE)
  }
  print(paste("Done with chromosome:",i,"at:",Sys.time()))
}

mpi.close.Rslaves()
mpi.quit()
