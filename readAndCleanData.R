# Notes For Use
# Should be fairly simple. Read in the clinical data and create a clinicalTable out of it. This data is cleaned by
# removing rows with odd or missing data (detailed below). Another row is created to keep track of who is dead
# or alive and is censored for patients with observations more than 4 years out (this number can be changed)
# Next the geneExpression data is read in. Only cancer columns are kept, and rows that don't meet the threshold
# (4 is chosen though more thought should go into this). Also the data is transformed using asinh. I don't think
# the data has already been transformed but this should be looked into.




GBM.clin.merged <- read.delim("GBM.clin.merged.txt", 
                              header=FALSE, stringsAsFactors=FALSE, row.names = 1)
GBM.clin.merged = t(GBM.clin.merged) # transposing it as it should be.
nobsOriginal = nrow(GBM.clin.merged)
rownames(GBM.clin.merged) = 1:nobsOriginal

# Create clinical table and then convert everything to correct type
clinicalTable = as.data.frame(GBM.clin.merged[, c('patient.bcr_patient_barcode', 'patient.days_to_death', 'patient.days_to_last_followup', 'patient.vital_status')])
clinicalTable[, 'patient.days_to_death'] = as.numeric(as.character(clinicalTable$patient.days_to_death))
clinicalTable[, 'patient.days_to_last_followup'] = as.numeric(as.character(clinicalTable$patient.days_to_last_followup))
clinicalTable[, 'patient.bcr_patient_barcode'] = as.character(clinicalTable$patient.bcr_patient_barcode)
clinicalTable[, 'patient.vital_status'] = as.character(clinicalTable$patient.vital_status)



# Clean up data so things make sense.
# Remove double NA calues in days to death and days since last follow up
# Those who visite after death
# days to death or days to last follow up is zero


# Remove double NA values in days to death and days since last follow up
bothNA = is.na(clinicalTable$patient.days_to_death) & is.na(clinicalTable$patient.days_to_last_followup)
# This vector is all False


# Remove those who visited after death
diedBeforeLastVisit = clinicalTable$patient.days_to_death < clinicalTable$patient.days_to_last_followup
diedBeforeLastVisit[is.na(diedBeforeLastVisit)]= FALSE
# This vector is all False.


# Remove those that have days to death or days to follow up = 0
dayDead_dead_followUp_equal0 = which(clinicalTable$patient.days_to_death == 0 | clinicalTable$patient.days_to_last_followup == 0)
clinicalTable = clinicalTable[-dayDead_dead_followUp_equal0,]


mostRecentEvent = pmax(clinicalTable$patient.days_to_death, clinicalTable$patient.days_to_last_followup, na.rm = TRUE)
clinicalTable = cbind(clinicalTable, mostRecentEvent)


# add one more column that says “alive” if the patient was still alive at the end of observation and “dead” if not.
deadOrAliveCols = apply( cbind(clinicalTable$patient.days_to_death, clinicalTable$patient.days_to_last_followup), 1, which.max)
# ^ Tells you which column to use to determine if person is dead or alive. If days to death is greater than days to follow up
# then they're dead
deadOrAlive = rep('dead', nrow(clinicalTable))
for(i in 1:nrow(clinicalTable)){
  colIdx = deadOrAliveCols[[i]] + 1 # Adding one to make up because these cols correspond to a cbind df
  if(colIdx == 3){
    deadOrAlive[i] = 'alive'
  }
}

clinicalTable = cbind(clinicalTable, deadOrAlive)



cutOff = 4*365
beyondCutoff = which(clinicalTable$mostRecentEvent > cutOff)
clinicalTable$deadOrAlive[beyondCutoff] = 'alive'


# Now dealing with the gene dexpression data

GMBGeneExpOriginal <- read.delim("GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
                         stringsAsFactors=FALSE, header = TRUE, row.names=1)

geneExpressionData = GMBGeneExpOriginal
cancerCols = colnames(geneExpressionData)[grep("01[A-Z].\\d{2}[A-Z]", colnames(geneExpressionData))]
geneExpressionData = geneExpressionData[, cancerCols]
# 01 is cancer and 11 is normal 
geneExpressionData = geneExpressionData[-1,] # remove first row which doesn't have numeric values
# geneExpressionData <- as.data.frame(sapply(geneExpressionData, as.numeric))

# NOTE has the log transform already been taken? Looks like no because there are such high values still
geneExpressionData = data.matrix(geneExpressionData)
geneExpressionData = asinh(geneExpressionData)
# rownames(geneExpressionData) = rownames(GMBGeneExpOriginal[-1,])
# Remove rownames with question marks
geneExpressionData = geneExpressionData[-grep("\\?", rownames(geneExpressionData)),]
hist(geneExpressionData[geneExpressionData != 0])

cutOff = 4
rowMax = apply(geneExpressionData, 1, max)
geneExpressionData = geneExpressionData[rowMax >= cutOff,]
geneExpressionData[geneExpressionData < cutOff] = 0 

genesRemain = nrow(geneExpressionData)
geneExpressionDataAfterCutOff = geneExpressionData
plot(density(geneExpressionDataAfterCutOff[geneExpressionDataAfterCutOff != 0]), main = 'Gene Expression Values')





