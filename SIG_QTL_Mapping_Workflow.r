
## Load R functions and libraries
source('rix_qtl_mapping_functions.r')

library(gdata)

## Read sample annotations (including phenotypes)
annot_dir = '/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/WNV/Cleaned_Data_Releases/r_data_files'

## This loads the 'all_weight' dataframe
load(file.path(annot_dir, 'gale_lund_weight_13-jan-2016_final.rda'))

pheno = all_weight
dim(pheno)

## Select only infected animals with a D10 weight measurement
pheno = pheno[pheno$Timepoint >= 12 & pheno$Virus == 'WNV', ]
dim(pheno)

## Calculate maximum weight loss
w_cols = c('D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'D13', 'D14')
w_cols = paste0(w_cols, '_Percentage')
pheno$d14_percent_wl = apply(pheno[, w_cols], 1, function(x) {x = x[!is.na(x)]; if (length(x) > 0) x[length(x)] else NA})
pheno$d14_percent_wl[is.infinite(pheno$d14_percent_wl)] = NA

library(lattice)
dotplot(reorder(pheno[,'UW_Line'], pheno[,'d14_percent_wl'], mean, na.rm=T) ~ 
        pheno[,'d14_percent_wl'] | pheno[,'Virus'], 
        panel = function(x,y,...) {panel.dotplot(x,y,...); panel.abline(v=0, col.line="red")}, 
        pch=19, ylab='UW Line', xlab="D14 Percent Weight Change")

## Summarize weight loss per line
pheno_per_line = aggregate(pheno[, 'd14_percent_wl'], list(pheno$UW_Line), median, na.rm=T)
colnames(pheno_per_line) = c('UW_Line', 'median_d14_percent_wl')

hist(pheno_per_line$median_d14_percent_wl, breaks=20, main="", xlab="Per-line Median Weight Change at D14")

## View extreme lines
pheno_per_line[with(pheno_per_line, median_d14_percent_wl <= -10 | median_d14_percent_wl >= 5), ]

## Get extreme lines
extreme_lines = pheno_per_line[with(pheno_per_line, median_d14_percent_wl <= -10 | median_d14_percent_wl >= 5), 1]
## Remove lines that recover by D14
extreme_lines = setdiff(extreme_lines, c(19, 55, 76))
## Subset pheno dataframe
pheno = pheno[pheno$UW_Line %in% extreme_lines,]

dotplot(reorder(pheno[,'UW_Line'], pheno[,'d14_percent_wl'], mean, na.rm=T) ~ 
        pheno[,'d14_percent_wl'] | pheno[,'Virus'], 
        panel = function(x,y,...) {panel.dotplot(x,y,...); panel.abline(v=0, col.line="red")}, 
        pch=19, ylab='UW Line', xlab="D14 Percent Weight Change")

## Sort pheno dataframe and set rownames
pheno = pheno[with(pheno, order(Mating, RIX_ID)),]
rownames(pheno) = pheno$ID

## Fix sex column name
colnames(pheno)[which(colnames(pheno) == 'Sex')] = 'sex'
dim(pheno)

## Create covariate dataframe (must include sex)
covar = data.frame(sex = as.numeric(pheno$sex == 'M'))
rownames(covar) = pheno$ID

## Get IDs and Matings for each sample
samples = pheno$ID
matings = unlist(lapply(strsplit(samples, '_'), function(x) {x[1]}))
matings = unique(matings)

## Read strain ID mapping file
mapping_dir = '/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/Mapping'
strain_map = read.xls(file.path(mapping_dir, 'RIXTrueIDs.xlsx'), header=T, as.is=T)
colnames(strain_map) = c('Mating', 'CC_Mating', 'Notes')

head(strain_map)

## Check that matings match the strain mapping file
setdiff(matings, strain_map$Mating)

## Correct the mismatched matings
pheno$Mating[pheno$Mating == '3609x5119'] = '3609x15119'
pheno$Mating[pheno$Mating == '8018x3154'] = '18018x3154'

pheno$ID = paste0(pheno$Mating, '_', pheno$RIX_ID)
rownames(pheno) = pheno$ID
rownames(covar) = pheno$ID

samples = pheno$ID
matings = unlist(lapply(strsplit(samples, '_'), function(x) {x[1]}))
matings = unique(matings)

## Get all file names for CC probability files
cc_prob_files = list.files(mapping_dir, pattern="CC...-.*b38.*\\.csv")

## Map RIX matings to CC matings
cc_matings = sapply(matings, function(x) {strain_map$CC_Mating[strain_map$Mating == x]}, USE.NAMES=F)

## Get unique vector of parental strains
rix_strains = unlist(strsplit(matings, 'x'))
cc_strains = unlist(strsplit(cc_matings, 'x'))

## Check strains
dup_idx = duplicated(cc_strains)
cc_strains = cc_strains[!dup_idx]
names(cc_strains) = rix_strains[!dup_idx]
cc_strains[1:5]

## Iterate through matings and save the 8-state probability file
## for each parental line as a .csv file (RIX IDs will be used for the file names)
rix_dir = '/Users/mooneymi/Documents/MyDocuments/SystemsImmunogenetics/Mapping/RIX'
for (i in 1:length(cc_strains)) {
    ## Get RIX ID
    rix = names(cc_strains[i])
    
    ## Create and save the 8-state probabilities for each parental strain
    cc_file = file.path(mapping_dir, cc_prob_files[grep(paste0(cc_strains[i], '-.*\\.csv'), cc_prob_files)])
    cc_8state = collapse_probs(cc_file)
    
    ## Save Y and M chromosomes separately
    autosome_x_markers = cc_8state$marker[!cc_8state$chromosome %in% c('Y', 'M')]
    y_m_markers = cc_8state$marker[cc_8state$chromosome %in% c('Y', 'M')]
    file_name = file.path(rix_dir, paste0(rix, '.csv'))
    write.table(cc_8state[autosome_x_markers, ], file=file_name, row.names=F, col.names=T, sep=',', quote=F)
    file_name = file.path(rix_dir, paste0(rix, '_Y_M.csv'))
    write.table(cc_8state[y_m_markers, ], file=file_name, row.names=F, col.names=T, sep=',', quote=F)
}

## Create model.probs object needed for DOQTL (an N x 8 x M array; N=samples, M=markers)
## This function will run fastest when samples are ordered by mating
model.probs = make_rix_model_probs(samples, rix_dir)

## Might help free memory
gc()

## Check dimensions and shape of model.probs object
dim(model.probs)
names(dimnames(model.probs))

## Check model.probs object
model.probs[1,,1:5]

## Get vector of all markers
all_markers = dimnames(model.probs)[[3]]

## Check if any markers sum to 0
markers_zero_idx = which(apply(model.probs, 1, colSums) == 0)
markers_zero = rep(all_markers, length(samples))[markers_zero_idx]
markers_zero = markers_zero[!is.na(markers_zero)]
markers_zero = unique(markers_zero)
length(markers_zero)

## Check if any markers don't sum to 1 (Optional)
## WARNING: This can take quite a while
prob_sums = apply(model.probs, 1, colSums)
prob_sums_not_one = sapply(prob_sums, function(x) {!isTRUE(all.equal(x, 1))})
markers_not_one_idx = which(prob_sums_not_one)
markers_not_one = rep(all_markers, length(samples))[markers_not_one_idx]
markers_not_one = markers_not_one[!is.na(markers_not_one)]
markers_not_one = unique(markers_not_one)
length(markers_not_one)

## Save markers with all zeros to file
write.table(data.frame(marker=markers_zero), file='markers_all_zeros.txt', sep='\t', col.names=T, row.names=F, quote=F)

## Remove markers with all zeros
model.probs = model.probs[, , setdiff(all_markers, markers_zero)]
dim(model.probs)

## Make all probabilities non-zero
model.probs[model.probs == 0] = 1e-20

## Check model.probs object
model.probs[1,,1:5]

## Save model.probs object to file (optional)
save(model.probs, file="rix_model_prob_test.rda")

## Create kinship probability matrix
K = kinship.probs(model.probs)

## Might help free memory
gc()

## Check kinship matrix
K[1:5, 1:5]

## First get marker positions
marker_pos = read.csv(file.path(mapping_dir, 'CC001-Uncb38V01.csv'), as.is=T)
marker_pos = marker_pos[,1:3]
marker_pos$position_cM = NA
head(marker_pos)

## Run QTL scan
qtl = scanone(pheno=pheno, pheno.col='d14_percent_wl', probs=model.probs, K=K, addcovar=covar, snps=marker_pos)

perms = scanone.perm(pheno = pheno, pheno.col = 'd14_percent_wl', probs = model.probs, addcovar = covar, snps = marker_pos, nperm = 1000)

thr = quantile(perms, probs = 0.95)

plot(qtl, sig.thr = thr, main = 'd14_percent_wl')

coefplot(qtl, chr = 7)

interval = bayesint(qtl, chr = 7)

ma = assoc.map(pheno = pheno, pheno.col = 'd14_percent_wl', probs = model.probs, K = K, addcovar = covar, 
               snps = marker_pos, chr = interval[1,2], start = interval[1,3], end = interval[3,3])
