print("#######################################################################")
print("Starting subset-allen-to-basal-ganglia.R script...")
sessionInfo()

#Load Allen Brain microarray data and corresponding metadata
print("Loading Allen Brain microarray data...")

#Functions to load data
LoadArrayData <- function (bioRepFolder) {
     arrayDataPath <- paste(
          "../raw_data/allen_data/", bioRepFolder, "/MicroarrayExpression.csv"
          , sep= ""
     )
     read.csv(arrayDataPath, header=FALSE)
}
LoadMetaData <- function (bioRepFolder) {
     metaDataPath <- paste(
          "../raw_data/allen_data/", bioRepFolder, "/SampleAnnot.csv"
          , sep= ""
     )
     read.csv(metaDataPath)
}

#Sub folders with Allen Brain span data, each folder is 1 brain
arrayDataFolderNames <- list(
     "178236545-2015-07-15",
     "178238266-2015-07-15",
     "178238316-2015-07-15",
     "178238359-2015-07-15",
     "178238373-2015-07-15",
     "178238387-2015-07-15"
)

#Load the Allen Brain microarray data
#Columns correspond to samples and samples are listed in the same order in the
#meta and array data
expArrayDataLDF <- lapply(arrayDataFolderNames, LoadArrayData)
metaArrayDataLDF <- lapply(arrayDataFolderNames, LoadMetaData)

# save(array.data, file="../processed.data/all.allen.array.data.rda")
################################################################################

#Subset the data to retain only the substantia nigra samples
print("Subsetting data to select the substantia nigra samples...")

# Areas to subset, identified by unix grep command:
# grep -r --include SampleAnnot* 'accumbens\|\caudate\|putamen\
# |substantia\|subthalamic' 17* | cut -d, -f4,5,6,7 | sort -u | column
areasToSubset <- c("substantia nigra", "nucleus accumbens", "caudate"
                     ,"putamen", "subthalamic")

# Function to subset the meta data
SubsetMetaData <- function (bioRepMetaData, toMatch) {
     subset(bioRepMetaData, grepl(paste(toMatch, collapse="|")
                                     , structure_name))
}

#Function to subset the array expression data
SubsetArrayData <- function (bioRepExpArrayDataLDF, bioRepMetaData
                               , toMatch) {
     subsetArrayData <- bioRepExpArrayDataLDF[ 
          , grepl(paste(toMatch, collapse="|")
                  , bioRepMetaData$structure_name)
          ]
}

#Subset the meta data
metaDataSubsetLDF <- lapply(metaArrayDataLDF
                               , SubsetMetaData
                               , areasToSubset)
#Meta data number of samples and number of features
print("Meta data subset number of samples (Row1) and number of features (Row2):")
sapply(metaDataSubsetLDF, dim)
print("Total number of samples in subset:")
sum(sapply(metaDataSubsetLDF, dim)[1, ])

#Subset the array expression data
arrayDataSubsetDF <- do.call(cbind, mapply(SubsetArrayData
                                           , expArrayDataLDF
                                           , metaArrayDataLDF
                                           , MoreArgs= list(areasToSubset)))
#Add the probe ID numbers as column 1
arrayDataSubsetDF <- cbind(probe=expArrayDataLDF[[1]][,1]
                              , arrayDataSubsetDF)
print("Dimensions of subset expression data, columns are samples:")
print("(row 1 is Probe ID)")
dim(arrayDataSubsetDF)

save(arrayDataSubsetDF, metaDataSubsetLDF,
     file="../processed_data/array_data_subset.rda")

print("End of subset-allen-to-basal-ganglia.R script...")

