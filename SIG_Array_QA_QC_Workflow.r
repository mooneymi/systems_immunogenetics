
## Load libraries and functions for array processing and QA/QC plots
source('array_qa_qc_functions.r')

## Read the annotation spreadsheet into R
annot_dir = '/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/Expression/Bat_Virus_Array'
setwd(annot_dir)
sample_annot = read.xls('BatPlate_Annotation_editedMM.xlsx', header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))
rownames(sample_annot) = sample_annot$ID

## Check annotation dataframe
head(sample_annot[,1:5])

## Set directory where .CEL files are located, and get the list of files
cel_dir = '/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/Expression/Bat_Virus_Array/data'
cel_files = list.celfiles(cel_dir)
## Check file names
cel_files[1:3]

## Create sample names (these must match the annotation file)
sample_names = gsub("_2.CEL", "", cel_files)
## Check sample names
sample_names[1:5]

## Subset sample annotation dataframe
sample_annot = sample_annot[sample_names,]
length(sample_names) == dim(sample_annot)[1]

## Create a phenoData object
phenoData = new("AnnotatedDataFrame", data=sample_annot)
phenoData

## Load the raw expression data
raw.exprs = read.celfiles(file.path(cel_dir, cel_files), pkgname="pd.mogene.2.1.st", 
                          sampleNames=sample_names, phenoData=phenoData)

head(pData(phenoData)[,1:5])

## Create un-normalized ExpressionSet
bgcor.exprs = rma(raw.exprs, normalize=FALSE, target="core")

## Create normalized ExpressionSet
norm.exprs = rma(raw.exprs, normalize=TRUE, target="core")

## Check normalized expression matrix
exprs(norm.exprs)[1:5,1:5]
## Save normalized expression to file (optional)
#save(norm.exprs, file="bat_virus_array_normalized.rda")

describe(make.boxplot)

## Boxplot of un-normalized expression values
make.boxplot(bgcor.exprs, type = 'raw', order.by=c("Mating", "Number"), color.by="Mating", make.pdf=F)

## Boxplot of normalized expression values
make.boxplot(norm.exprs, type = 'norm', order.by=c("Mating", "Sex", "Number"), color.by="Mating", make.pdf=F)

describe(plot.bac.spikes)

## The probe group file from Affy for the array
pg_file = '/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/Expression/MoGene-2_1-st.pgf'

## Create bacterial spike plot
plot.bac.spikes(norm.exprs, pg_file, make.pdf=F)

describe(make.heatmap)

make.heatmap(norm.exprs, base.factors=c('Sex', 'D4_percent'), make.pdf=F)
