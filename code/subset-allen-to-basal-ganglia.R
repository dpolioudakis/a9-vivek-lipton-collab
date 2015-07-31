sessionInfo()

print("Starting subset-allen-to-basal-ganglia.R script...")
################################################################################

#Load Allen Brain microarray data and corresponding metadata
print("Loading Allen Brain microarray data...")

#Functions to load data
Load_array_data <- function (bio.rep.folder) {
     array.data.path <- paste(
          "../raw_data/allen_data/", bio.rep.folder, "/MicroarrayExpression.csv"
          , sep= ""
     )
     read.csv(array.data.path, header=FALSE)
}
Load_meta_data <- function (bio.rep.folder) {
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
exp.array.data <- lapply(array.data.folder.names, Load_array_data)
meta.array.data <- lapply(array.data.folder.names, Load_meta_data)

# save(array.data, file="../processed.data/all.allen.array.data.rda")
################################################################################

#Subset the data to retain only the substantia nigra samples
print("Subsetting data to select the substantia nigra samples...")

# Function to subset the meta data
Subset_meta_data <- function (bio.rep.meta.data, to.match) {
     subset(bio.rep.meta.data, grepl(paste(to.match, collapse="|"), structure_name))
}

# Areas to subset, identified by unix grep command:
     # grep -r --include SampleAnnot* 'accumbens\|\caudate\|putamen\
     # |substantia\|subthalamic' 17* | cut -d, -f4,5,6,7 | sort -u | column
areas.to.subset <- c("substantia nigra", "nucleus accumbens", "caudate"
                    ,"putamen", "subthalamic")

#Function to subset the array expression data
Subset_array_data <- function (bio.rep.exp.array.data, bio.rep.meta.data) {
     subset.array.data <- bio.rep.exp.array.data[ , grepl(
          "substantia nigra", bio.rep.meta.data$structure_name)]
}

#Subset the meta data
meta.data.subset <- lapply(meta.array.data, Subset_meta_data, areas.to.subset)
#Meta data number of samples and number of features
print("Meta data subset number of samples (Row1) and number of features (Row2):")
sapply(meta.data.subset, dim)
print("Total number of samples in subset:")
sum(sapply(meta.data.subset, dim)[1, ])

#Subset the array expression data
array.data.subset <- do.call(cbind, mapply(Subset_array_data
                                           , exp.array.data
                                           , meta.array.data))
#Add the probe ID numbers as column 1
array.data.subset <- cbind(probe=exp.array.data[[1]][,1], array.data.subset)
print("Dimensions of subset expression data, columns are samples:")
print("(row 1 is Probe ID)")
dim(array.data.subset)

save(array.data.subset, meta.data.subset,
     file="../processed_data/array.data.subset.rda")
