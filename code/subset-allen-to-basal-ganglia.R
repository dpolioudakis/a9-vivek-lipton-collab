sessionInfo()

print("Starting subset-allen-to-basal-ganglia.R script...")
################################################################################

#Load Allen Brain microarray data and corresponding metadata
print("Loading Allen Brain microarray data...")

#Functions to load data
load.array.data <- function (bio.rep.folder) {
     array.data.path <- paste(
          "../raw_data/allen_data/", bio.rep.folder, "/MicroarrayExpression.csv"
          , sep= ""
     )
     read.csv(array.data.path, header=FALSE)
}
load.meta.data <- function (bio.rep.folder) {
     meta.data.path <- paste(
          "../raw_data/allen_data/", bio.rep.folder, "/SampleAnnot.csv"
          , sep= ""
     )
     read.csv(meta.data.path)
}

#Sub folders with Allen Brain span data, each folder is 1 brain
array.data.folder.names <- list(
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
exp.array.data <- lapply(array.data.folder.names, load.array.data)
meta.array.data <- lapply(array.data.folder.names, load.meta.data)

# save(array.data, file="../processed.data/all.allen.array.data.rda")
################################################################################

#Subset the data to retain only the substantia nigra samples
print("Subsetting data to select the substantia nigra samples...")

# Function to subset the meta data
subset.sn.meta.data.fxn <- function (bio.rep.meta.data) {
     subset(bio.rep.meta.data, grepl("substantia nigra", structure_name))
}

#Function to subset the array expression data
subset.array.data.fxn <- function (bio.rep.exp.array.data, bio.rep.meta.data) {
     subset.array.data <- bio.rep.exp.array.data[ , grepl(
          "substantia nigra", bio.rep.meta.data$structure_name)]
}

#Subset the meta data
subset.sn.meta.data <- lapply(meta.array.data, subset.sn.meta.data.fxn)
#Meta data substantia nigra number of samples and number of features
print("Meta data substantia nigra number of samples and number of features:")
sapply(subset.sn.meta.data, dim)
#Total number of substantia nigra samples
print("Total number of substantia nigra samples:")
sum(sapply(subset.sn.meta.data, dim)[1, ])

#Subset the array expression data
subset.sn.array.data <- do.call(cbind,
                                mapply(subset.array.data.fxn, exp.array.data, meta.array.data))
#Add the probe ID numbers as column 1
subset.sn.array.data <- cbind(probe=exp.array.data[[1]][,1], subset.sn.array.data)
print("Dimensions of subset expression data, columns are samples:")
print("(row 1 is Probe ID)")
dim(subset.sn.array.data)

save(subset.sn.array.data, subset.sn.meta.data,
     file="../processed_data/subset.sn.array.data.rda")
