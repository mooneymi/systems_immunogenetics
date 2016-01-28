
## Load functions for plotting the flow cytometry data
source('flow_data_plotting_functions.r')

## Load data from an Excel spreadsheet (Warning: this can take a few minutes)
## Note: you may have to change the file paths
#flow_full = read.xls('../Cleaned_Data_Releases/15-Jan-2016/Lund_Flow_Full_11-Jan-2016_final.xlsx', 
#                     header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))
## Load data from an R data file
load('../Cleaned_Data_Releases/r_data_files/lund_flow_full_11-jan-2016_final.rda')

describe(flow_boxplot_data)

## Aggregate the data for boxplots
boxplot_data = flow_boxplot_data(flow_full, c(7,8,9), 'brain', 'treg_T_regs', c('7','12','21','28'))

## Create a list of additional options for the boxplot
opts = list(rm_outliers=F, show_data=F, y_min=0, y_max=60)
## Create the boxplot (the 'cex' parameter controls the size of the x-axis text)
bp = flow_boxplots(c(boxplot_data, opts), cex=0.7)

describe(flow_multiline_plot_data)

lineplot_data = flow_multiline_plot_data(flow_full, c(7,8,9), 'brain', 'treg_T_regs', 1)

## Create a list of additional options for the lineplot
## data_type values: 1 = percentages, 2 = cell counts, 3 = percent ratio, 4 = count ratio
opts2 = list(data_type=1, y_min=0, y_max=60)
## Create a lineplot that compares a single variable across multiple lines
lp = flow_multiline_plots(c(lineplot_data, opts2))

lineplot_data2 = flow_multiline_plot_data(flow_full, c(9), 'brain', c('treg_T_regs', 'tcell_d7_CD8'), 2)

## Create a list of additional options for the lineplot
## data_type values: 1 = percentages, 2 = cell counts, 3 = percent ratio, 4 = count ratio
opts3 = list(data_type=1, y_min=0, y_max=50)
## Create a lineplot that compares multiple variables for a single line
lp2 = flow_multiline_plots(c(lineplot_data2, opts3))

## Load weight, clinical score, and heritability data from the latest data release
## Note: you may have to change the file paths
weights = read.xls('../Cleaned_Data_Releases/15-Jan-2016/Lund_Weight_13-Jan-2016_final.xlsx',
                        header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))
scores = read.xls('../Cleaned_Data_Releases/15-Jan-2016/Lund_Scores_13-Jan-2016_final.xlsx',
                       header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))
heritability = read.xls('../Cleaned_Data_Releases/15-Jan-2016/Lund_Flow_Heritability_11-Jan-2016_final.xlsx',
                        header=T, as.is=T, na.strings=c(""," ", "NA", "#DIV/0!"))

## Set the rownames of the heritability dataframe
rownames(heritability) = heritability$variable

describe(flow_heatmap_data)

heatmap_data = flow_heatmap_data(flow_full, c(7,8,9), 'brain', 
                                 c('treg_T_regs', 'tcell_d7_CD3', 'tcell_d7_CD4', 'tcell_d7_CD8'),
                                 c('7','12','21','28'), heritability)

describe(flow_heatmap_plot)

## Create the heatmap
hm = flow_heatmap_plot(heatmap_data, weights, scores)

## The heatmap without annotations
heatmap_data2 = flow_heatmap_data(flow_full, c(7,8,9), 'brain', 
                                 c('treg_T_regs', 'tcell_d7_CD3', 'tcell_d7_CD4', 'tcell_d7_CD8'),
                                 c('7','12','21','28'), annotations=F)
hm2 = flow_heatmap_plot(heatmap_data2, annotations=F)
