library(gdata)
library(abind)
library(DOQTL)


## Help Documentation
describe = function(obj) {
  if ('help' %in% names(attributes(obj))) {
    writeLines(attr(obj, 'help'))
  }
}
attr(describe, 'help') = "
This function prints the contents of the 'help' attribute of any R object. 
It is meant to provide help documentation in the same vein as Docstrings in Python. 
"

geno_states = c("AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "AB", "AC", 
"AD", "AE", "AF", "AG", "AH", "BC", "BD", "BE", "BF", "BG", "BH", 
"CD", "CE", "CF", "CG", "CH", "DE", "DF", "DG", "DH", "EF", "EG", 
"EH", "FG", "FH", "GH")

collapse_probs = function(probs_file, states=geno_states, gbuild='b38') {
  ## Load probabilities
  full_probs = read.csv(probs_file, as.is=T)
  
  ## Founder states
  founder_states = c('AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH')
  
  ## Get heterozygous states for each founder
  het_states = list()
  for (fs in founder_states) {
	letter = strsplit(fs, '')[[1]][1]
	pattern = paste(letter,'[^',letter,']|[^',letter,']',letter,sep='')
    het_states[[fs]] = states[grep(pattern, states)]
  }
  
  ## Initialize the collapsed dataframe with only the homozygous states (AA, BB, etc.)
  coll_probs = full_probs[, 1:11]
  
  ## Collapsed probabilities for each founder is the homozygous probability plus 0.5*probability of all the heterozygous states
  for (fs in founder_states) {
	coll_probs[, fs] = full_probs[, fs] + rowSums(full_probs[, het_states[[fs]]])/2
  }
  
  colnames(coll_probs) = c('marker', 'chromosome', paste0('position_', gbuild), 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
  rownames(coll_probs) = coll_probs$marker
  
  return(coll_probs)
}
attr(collapse_probs, 'help') = "
This function collapses 36-state probabilities to 8 probabilites (one for each founder strain).

Parameters:
probs_file: The filename of a .csv file containing the 36-state probabilities. Reading this file 
   with read.csv() should yield an n x 39 dataframe containing marker info (marker, chromosome, position) 
   for n markers, and the 36-state probabilities.
states: A character vector containing the 36 genotype states (default is the geno_states global variable
   as defined in this script).

Returns:
A n x 11 dataframe, where n is the number of markers and the columns represent the
marker, chromosome, position, and probabilities for each of the 8 founder strains. 
"

make_rix_probs = function(dam_file, sire_file) {
  ## Load probabilities
  dam = read.csv(dam_file, as.is=T)
  sire = read.csv(sire_file, as.is=T)
  
  ## Standardize column names, just in case
  colnames(dam)[1:2] = c('marker', 'chromosome')
  colnames(sire)[1:2] = c('marker', 'chromosome')
  
  ## Check that markers from parental lines match
  if (!identical(dam$marker, sire$marker)) {
    stop("Markers for dam and sire do not match!")
  }
  
  ## Check that genotype states from parental lines match
  dam_states = colnames(dam)[4:length(colnames(dam))]
  sire_states = colnames(sire)[4:length(colnames(sire))]
  if (!identical(dam_states, sire_states)) {
	stop("Genotype states for dam and sire do not match!")
  } else {
	states = dam_states
  }
  
  rownames(dam) = dam$marker
  rownames(sire) = sire$marker
  autosomes = sort(unique(dam$chromosome[!dam$chromosome %in% c('M', 'X', 'Y')]))
  
  ## Create 3D array
  probs = abind(dam[,states], sire[,states], along=3)
  dimnames(probs)[[3]] = c('dam', 'sire')
  names(dimnames(probs)) = c('markers', 'states', 'samples')
  
  ## Calculate marker probabilities for males
  ## Initialize probabilities as 0
  male_df = sire
  male_df[,states] = 0
  
  ## Autosomes
  ## Average the probabilities from dam and sire
  autosome_markers = male_df$marker[male_df$chromosome %in% autosomes]
  if (length(autosome_markers) > 0) {
    male_df[autosome_markers, states] = apply(probs[autosome_markers, states, ], c('markers', 'states'), function(x) {.Internal(mean(x))})
  }
  
  ## Y chromosome
  ## Copy the probabilities from sire
  y_markers = male_df$marker[male_df$chromosome == 'Y']
  if (length(y_markers) > 0) {
    male_df[y_markers, states] = probs[y_markers, states, 'sire']
  }
  
  ## X chromosome
  ## Copy the probabilities from dam
  x_markers = male_df$marker[male_df$chromosome == 'X']
  if (length(x_markers) > 0) {
    male_df[x_markers, states] = probs[x_markers, states, 'dam']
  }
  
  ## M chromosome
  ## Copy the probabilities from dam
  m_markers = male_df$marker[male_df$chromosome == 'M']
  if (length(m_markers) > 0) {
    male_df[m_markers, states] = probs[m_markers, states, 'dam']
  }

  ## Calculate marker probabilities for females
  ## Initialize probabilities as 0
  female_df = dam
  female_df[,states] = 0
  
  ## Autosomes and X chromosome
  ## Average the probabilities from dam and sire
  if (length(c(autosome_markers, x_markers)) > 0) {
    female_df[c(autosome_markers, x_markers), states] = apply(probs[c(autosome_markers, x_markers), states, ], c('markers', 'states'), function(x) {.Internal(mean(x))})
  }

  ## M chromosome
  ## Copy the probabilities from dam
  if (length(m_markers) > 0) {
    female_df[m_markers, states] = probs[m_markers, states, 'dam']
  }
  
  ## Y chromosome
  ## All probabilities = 0
  
  return(list(male=male_df, female=female_df))
}
attr(make_rix_probs, 'help') = "
This function calculates the genotype state probabilities for RIX lines using the probabilities 
from the two parental RI lines.

Parameters:
dam: An n x 3+s dataframe containing marker info (marker, chromosome, position) for n markers,
   and s genotype state probabilities.
sire: An n x 3+s dataframe containing marker info (marker, chromosome, position) for n markers,
   and s genotype state probabilities. 

Returns:
A list with two (male and female) n x 3+s dataframes containing the marker and probability 
information for a RIX (F1-hybrid of the given dam and sire).
"

make_rix_model_probs = function(samples, file_path='.') {
  ## Founder states
  founder_states = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
  
  ## Get probabilities for the first sample
  ## Initialize the 3D array
  mating = unlist(strsplit(samples[1], '_'))[1]
  dam = unlist(strsplit(mating, 'x'))[1]
  sire = unlist(strsplit(mating, 'x'))[2] 
  dam_file = file.path(file_path, paste0(dam, '.csv'))
  sire_file = file.path(file_path, paste0(sire, '.csv'))
  
  rix_probs = make_rix_probs(dam_file, sire_file)$male[, founder_states, ]
  model_prob_array = rix_probs
  
  prev_mating = mating
  ## Iterate through the rest of the samples
  for (s in samples[2:length(samples)]) {
    mating = unlist(strsplit(s, '_'))[1]
    
	## For each mating
	if (mating != prev_mating) {
      dam = unlist(strsplit(mating, 'x'))[1]
      sire = unlist(strsplit(mating, 'x'))[2]
      dam_file = file.path(file_path, paste0(dam, '.csv'))
      sire_file = file.path(file_path, paste0(sire, '.csv'))
	
	  ## Create RIX probabilities
	  rix_probs = make_rix_probs(dam_file, sire_file)$male[, founder_states, ]
	}
    
    model_prob_array = abind(model_prob_array, rix_probs, along=3)
	
	prev_mating = mating
  }
  
  ## Assign names to array
  names(dimnames(model_prob_array)) = c('markers', 'founders', 'samples')
  dimnames(model_prob_array)[['samples']] = samples
  
  ## Reshape the 3D array
  model_prob_array = aperm(model_prob_array, c('samples', 'founders', 'markers'))
  
  return(model_prob_array)
}
attr(make_rix_model_probs, 'help') = "
This function calculates the genotype state probabilities for RIX lines using the probabilities 
from the two parental RI lines.

"

## Read sample annotation

## Get probs for each sample (create RIX probs from RI probs)

## Construct 3d array

## Condense 36 probs to 8 probs

## Get kinship matrix

## Run scan

## Get QTL intervals

## Run association tests in intervals (additive model and dominant model)