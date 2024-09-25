#Function (for R) to estimate allele frequencies from a biallelic multilocus dataset

#input files = geno and locname

#geno is a matrix containing the genotypes: one row per individual; 1st column = individual identifier, then 2 alleles of the same locus given in consecutive cells
#(geno has the same format as the input file of the function GENHET)
#missing data have to be coded as "NA";
#see the example file named "exGENHETgenotinput.txt" provided with this code file
#NOTE: if only one allele of the genotype of a locus is missing, the other one will be included in the estimation of allele freq

#locname is a vector with the names of the different loci (in the same order as in geno)


#copy and paste the commands below into an R window and apply it to your input files


"ALF"=
  function(geno,locname){
    
    nbind=nrow(geno)
    nbloc=(ncol(geno)-1)/2
    
    #creation of the list of alleles
    genov=vector(length=nbind*nbloc*2)
    for (i in 2:ncol(geno)) genov[(nrow(geno)*(i-2)+1):(nrow(geno)*(i-1))]=geno[,i]
    al=sort(na.omit(unique(genov)))
    
    #count of the number of times each allele appears + nb of missing data
    alcount=matrix(nrow=(length(al)+1),ncol=(nbloc+1))
    alcount[,1]=c(al,NA)
    for(j in 1:(nrow(alcount)-1))
      for(k in 1:(ncol(alcount)-1))
        alcount[j,(k+1)]=sum(geno[,(k*2):(k*2+1)]==alcount[j,1],na.rm=T)
    for(l in 2:ncol(alcount))
      alcount[nrow(alcount),l]=(2*nbind-sum(alcount[1:(nrow(alcount)-1),l]))
    
    
    #creation of the table of allele frequencies
    alfreq=matrix(nrow=length(al),ncol=(nbloc+1))
    colnames(alfreq)=c("Allele",locname)
    alfreq[,1]=al
    for(m in (1:nrow(alfreq)))
      for (n in 2:ncol(alfreq)) alfreq[m,n]=alcount[m,n]/(nbind*2-alcount[nrow(alcount),n])
    
    return(alfreq)
    
  }
