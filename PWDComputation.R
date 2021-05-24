

library(gtools)
library(xlsx)
library(gnn)


# important to note, this formulation requires number of samples per observation to be the same! 
# if samples are not the same, need to change, and change averaging procedure (since diff samp sizes)

CpGPWD = function(data,name, Observations, Samples, outputFile)
{
  # Set Variables For Loop
  
  numObs = Observations # number of individuals studied
  
  numSamples= Samples # number of samples per individual 
  
  combinations = combinations(numSamples, 2) # all possible pairwise combinations 
  
  sampleNames = names(data) # names will be used for naming pairwise difference cols
  
  avgMatrix = matrix(0, nrow = dim(data)[1], ncol = numObs) # for easy computation of avg PWD 
  
  ##PWD across samples for each individual, average pwd per cpg site, average PWD per individual##############
  
  for(i in 1:numObs)
  {
    pairs = combinations + (i-1)*numSamples   #generate columns for pairwise differences
    
    pwd = abs(data[pairs[,1]] - data[pairs[,2]])  #compute abs val of pairwise difference
    # will generate Choose(samples, 2) cols
    
    colnames = c()    #initialize, will be used for new col names of pairwise differences
    
    for(j in 1: dim(pairs)[1]) #over all pairs
      colnames = c(colnames,paste("| ",sampleNames[pairs[j,1]], " - ", sampleNames[pairs[j,2]], " |", sep=""))
    #find name of respective pairwise diff col
    
    colnames = c(colnames, paste("Average PWD per CpG Individual" , i) )
    
    average = rowSums(pwd)/dim(pairs)[1]  #col of average PWD over all samples (per individual)
    
    pwd[, dim(pairs)[1]+1] = average      # insert individual average, before merging with dataset
    
    avgMatrix[,i] = average   # save avg pwd in matrix, to compute avg PWD over all individuals
    
    data[colnames] = pwd      #merge the dataset
    
  }
  
  
  data[["AverageAcrossCpG"]] = rowSums(avgMatrix) / numObs    #insert average PWD over all individuals
  ##############################################################################################################
  
  #perform gene-wise analysis, starting here####################################################################
  
  n = dim(data)[1]  #number of CpGs  
  
  geneNames = unique(data$UCSC_RefGene_Name) #get unique gene names
  geneNames = geneNames[which(!is.na(geneNames))]
  
  # will need to automate this part, have gene dataframe name as param
  
  #count the number of cpgs per gene######################################
  genes = rep(" ", n) # to merge with dataframe
  freq = rep(" ", n) #to merge with dataframe
  
  table = as.data.frame(table(data$UCSC_RefGene_Name))
  
  table = table[order(table$Freq, decreasing=TRUE),] #sort by frequency of Cpgs
  
  
  
  genes[1:length(table$Var1)] = as.vector(table$Var1)
  
  freq[1:length(table$Freq)] = as.vector(table$Freq)
  
  genecounts = as.data.frame(cbind(genes, freq))
  
  data = cbind(data, genecounts)
  ########################################################################
  
  #Perform average pwd for each gene across all individuals
  individualGeneAvgPwd = matrix(nrow = n, ncol = numObs+1 )
  
  # generate column names
  columnnames=c()
  for(i in 1: numObs)
  {
    columnnames[i] = paste("Average PWD per Gene Individual", i)
  }
  columnnames=c(columnnames, "Average Over Gene" )
  
  colnames(individualGeneAvgPwd) = columnnames # set the column names
  
  for(i in 1:numObs)
  {
    
    for(j in 1: length(geneNames))
    {
      gene = genes[j] # name of gene for each iteration
      
      colname = paste("Average PWD per CpG Individual" , i) #reference column of individual average pwd for cpg site
      
      CPGs = which(data$UCSC_RefGene_Name == gene) # find cpgs that correspond to gene 
      
      genepwd  = sum(data[[colname]][CPGs])/ length(CPGs) #average pwd over all cpgs corresponding to gene 
      
      individualGeneAvgPwd[j,i] = genepwd 
    }
  }
  if(numObs>1)
  {
    individualGeneAvgPwd[,numObs+1] = rowSums(individualGeneAvgPwd[, 1:numObs])/numObs  # Average Gene PWD
  }else
  {
    individualGeneAvgPwd[,numObs+1] = individualGeneAvgPwd[, 1]
  }
  
  data = cbind(data,individualGeneAvgPwd)
  ########################################################################
  
  # initialize columns
  
  # data$gene= rep(" ", n) # initialize for pwd
  # 
  # data$avgOverGene = rep(" ", n) #initialize 
  # 
  # data$avgOverGene = rep(" ", n)
  # 
  # 
  # ###### Average pwd for each gene over all individuals (using AverageAcrossCpG)
  # for(i in 1: length(geneNames)) #iterate through all unique genes
  # {
  #   gene = geneNames[i] 
  #   
  #   whichCpG = which(data$UCSC_RefGene_Name== gene ) #find which cpgs correspond to gene
  #   
  #   count = length(whichCpG) #how many cpgs correspond to the cpg
  #   
  #   geneAvgPWD = sum(data$AverageAcrossCpG[whichCpG]) / count #find gene average pwd 
  #   
  #   data$gene[i] = gene
  #   
  #   data$avgOverGene[i] = geneAvgPWD
  # }
  # 
  #
  
  
  # Average CpG PWD across each individual
  for(i in 1: numObs)
  {
    colname = paste("Average PWD per CpG Individual" , i)
    
    nonNA = which(!is.na(data[, colname ])) #cannot sum the NAs 
    
    averagePWD = sum(data[nonNA,colname]) /n 
    
    newcolname = paste("Average PWD Individual" , i) 
    
    data[[newcolname]] = c(averagePWD, rep(NA, n-1)) #initialize column
  }
  
  nonNAEntries = which(!is.na(data$AverageAcrossCpG)) #Exclude NAs from sum
  
  averageTissue = sum(data$AverageAcrossCpG[nonNAEntries])/n  #find average pwd across individuals and across CpGs
  
  data[["Average Tissue PWD"]] = c(averageTissue, rep(NA, n-1)) # set column in data frame 
  
  
  sortPWD = data.frame(cbind(genes,freq, data[["Average Over Gene"]])) #create data frame for sorting
  
  names(sortPWD) = c("gene", "freq", "pwd") #names for access to columns of data frames
  
  sortPWD = sortPWD[which(as.numeric(sortPWD$freq) >3) , ] #remove genes with less than 4 cpgs
  
  sortPWD = sortPWD[order(as.numeric(sortPWD$pwd)),] #sort by pwd increasing 
  
  sortPWD[(dim(sortPWD)[1]+1):n, 1:3] = NA #repopulate data frame for proper merge dimensions
  
  names(sortPWD) = c("gene > 3Cpgs ", "freq > 3Cpgs", "pwd > 3Cpgs") #rename for excel sheet 
  
  data= cbind(data, sortPWD)  #merge
  ################################################################################################################
  
  #write.xlsx2(data, outputFile)   #write out to excel file, does not work 
  save(data,file=outputFile)
  rename_rda("data", outputFile, name, outputFile)
}

###### MAIN ####################################################################################

library(readxl)

#Save data locally in R objects


siPWD = read_excel("all 8 SI PWD 10-20.xlsx")
siPWD[24] = c()
save(siPWD,file="siPWD.rda")

allColon = read_excel("all colon x4 9-20 PWD.xlsx")
save(allColon,file="allColon.rda")

uterusGld = read_excel("uterus gld 9-20.xlsx")
save(uterusGld ,file="uterusGld.rda")

EpicArrays= read_excel("list of all EPIC arrays 3-18.xls")
save(EpicArrays,file="EpicArrays.rda")

################################################################################################






load("allColonNew.rda")

CpGPWD(allColon, "allColon", 8, 4, "allColonNew.rda")

load("allColonNew.rda")

write.csv(allColon, "allColonNew.csv")



load("siPWDNew.rda")

CpGPWD(siPWD, "siPWD", 4, 4, "siPWDNew.rda")

load("siPWDNew.rda")

write.csv(siPWD, "siPWDNew.csv")



load("uterusGldNew.rda")

CpGPWD(uterusGld, "uterusGld", 8, 4, "uterusGldNew.rda")

load("uterusGldNew.rda")

write.csv(uterusGld, "uterusGldNew.csv")




