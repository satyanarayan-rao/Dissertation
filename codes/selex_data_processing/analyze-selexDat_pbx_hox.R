args <- commandArgs(trailingOnly=TRUE)
print(args)
workDir <- args[1]
annotation <- args[2]
id <- args[3]
roundAnaly <- args[4]
threshold <- as.integer(args[5])
setoptlen <- as.integer(args[6])

options(java.parameters="-Xmx5000M") #1500
library(SELEX)                                          
selex.config(workingDir=workDir, maxThreadNumber=4)



#####################################################
#Defining samples and creating handles
#
selex.loadAnnotation(annotation)
selex.sampleSummary()
mysamples = selex.sampleSummary()
#..

##create sample handle to be used in other selex functions
#



r0 = selex.sample(seqName=mysamples$seqName[mysamples$round==0], 
                  sampleName=mysamples$sampleName[mysamples$round==0], round=0)
r1 = selex.sample(seqName=mysamples$seqName[mysamples$round==1], 
                  sampleName=mysamples$sampleName[mysamples$round==1], round=1)


##########################################################
## Markov Model
#seqFilter using PCTRL probe:
regex = selex.seqfilter(variableRegionIncludeRegex="NNGAYNNRYNNN") # not using here
kmax.value = selex.kmax(sample=r0, threshold=threshold) 


#mm = selex.mm(sample=r0train, order=NA, crossValidationSample=r0test, Kmax=kmax.value, seqfilter= regex) 
#mm = selex.mm(sample=r0train, order=5, crossValidationSample=r0test, Kmax=kmax.value)
#mm = selex.mm(sample=r0, order=5, seqfilter=regex)
mm = selex.mm(sample=r0, order=5, 
              mmWithLeftFlank=TRUE, 
              Kmax=kmax.value)
#selex.mmSummary()
#..
######################################################
#Optimal motif length
selex.infogain(sample=get(roundAnaly), markovModel=mm, checkBarcode=TRUE)
#selex.infogain(sample=get(roundAnaly) ,markovModel=mm)
igain = selex.infogainSummary()[,1:3]
write.table(igain, paste0(workDir, "/infoGain_", ".txt"), sep="\t")

## CALCULATING FOR OPTIMAL LENGTH
optimalLength=NULL
if (setoptlen !=0){
    optimalLength = setoptlen
}else{ 
    optimalLength = kmax.value
}
table = selex.counts(sample=get(roundAnaly), k=optimalLength, markovModel=mm)
aff = selex.affinities(sample=get(roundAnaly), k=optimalLength, markovModel=mm)
aff = aff[with(aff, order(-Affinity, Probability)), ]
write.table(aff, paste0(workDir, "/results_R0", 
                        roundAnaly , "_", 
                        optimalLength, "mer", "_", 
                        "count_threshold", "_",
                        as.character(threshold),"_",
                        "optimalLength", "_",
                        optimalLength, ".txt"), sep="\t")
