# clear R workspace
rm(list = ls())



library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)
library(plotly)
library(MSnbase)

#read peptides file
header <- read.table(file="/home/compomics/git/Projects/MLMarker_paper/PXD007592_responders_melanoma/QuantifiedPeptides.tsv",sep="\t", nrows=1, header=FALSE, stringsAsFactors = FALSE)
peptidesFile <- read.table(file="/home/compomics/git/Projects/MLMarker_paper/PXD007592_responders_melanoma/QuantifiedPeptides.tsv",sep="\t", skip=1, header=FALSE)
colnames(peptidesFile) <- unlist(header)
peptidesFile <- peptidesFile[!(is.na(peptidesFile$"Protein Groups") | peptidesFile$"Protein Groups"==""), ]
peptidesFile = peptidesFile[order(peptidesFile$"Protein Groups"),]
head(peptidesFile)


peptidesFile <- peptidesFile[,c(
  "Sequence", 
  "Base Sequence", 
  "Protein Groups", 
  "Intensity good_responder_1_1",
  "Intensity good_responder_1_2",
  "Intensity good_responder_2_1",
  "Intensity good_responder_2_2",
  "Intensity good_responder_3_1",
  "Intensity good_responder_3_2",
  "Intensity good_responder_4_1",
  "Intensity good_responder_4_2",
  "Intensity good_responder_5_1",
  "Intensity good_responder_5_2",
  "Intensity poor_responder_1_1",
  "Intensity poor_responder_1_2",
  "Intensity poor_responder_10_1",
  "Intensity poor_responder_10_2",
  "Intensity poor_responder_11_1",
  "Intensity poor_responder_11_2",
  "Intensity poor_responder_12_1",
  "Intensity poor_responder_12_2",
  "Intensity poor_responder_13_1",
  "Intensity poor_responder_13_2",
  "Intensity poor_responder_2_1",
  "Intensity poor_responder_2_2",
  "Intensity poor_responder_3_1",
  "Intensity poor_responder_3_2",
  "Intensity poor_responder_4_1",
  "Intensity poor_responder_4_2",
  "Intensity poor_responder_5_1",
  "Intensity poor_responder_5_2",
  "Intensity poor_responder_6_1",
  "Intensity poor_responder_6_2",
  "Intensity poor_responder_7_1",
  "Intensity poor_responder_7_2",
  "Intensity poor_responder_8_1",
  "Intensity poor_responder_8_2",
  "Intensity poor_responder_9_1",
  "Intensity poor_responder_9_2")
    ]

#get columns with intensity
ecols <- grep(
  "Intensity", 
  names(peptidesFile)
  )
ecols

#store in QFeatures object
pe <- readQFeatures(
  table = peptidesFile,
  fnames = 1,
  ecol = ecols,
  name = "peptideRaw", sep="\t")
pe[["peptideRaw"]]


## Experimental layout

metadata <- read.table(file = "/home/compomics/git/Projects/MLMarker_paper/PXD007592_responders_melanoma/metadata.csv", sep = ",", header = TRUE)
colnames(pe[["peptideRaw"]]) <- metadata$run
#colData(pe)$sample_type <- metadata$condition
colData(pe)$sample_type <- metadata$Predicted
col_data <- colData(pe[["peptideRaw"]])

## Preprocessing <br>
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)  #number of non zero intensities per peptide
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

#plot missingness on peptide level (a lot because of modified peptides)
MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

#log transformation
pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]), legend=FALSE)

#filter for proteins found in at least 12 samples
print(nrow(pe[["peptideLog"]]))
pe <- pe[rowData(pe[["peptideRaw"]])$nNonZero>=(4),,] 
print(nrow(pe[["peptideLog"]])) #number of peptides
colData(pe[["peptideLog"]]) <- colData(pe)

#median normalization
pe <- normalize(pe, i = "peptideLog", method = "center.median", name = "peptideNorm")
limma::plotDensities(assay(pe[["peptideNorm"]]), legend=FALSE)
pe[["peptideNorm"]]
boxplot(assay(pe[["peptideNorm"]]), col = palette()[-1],
       main = "Peptide distribtutions after normalisation", ylab = "intensity")

#multidimensional scaling

limma::plotMDS(assay(pe[["peptideNorm"]]), col = colData(pe)$sample_type, labels = colData(pe)$sample_type)
# distance between two points = leading logFC, leading logFC average of largest absolute logFC between each pari of samples
#aggregate peptide intensities in protein expression values
pe <- aggregateFeatures(pe,
 i = "peptideNorm",
 fcol = "Protein.Groups",
 na.rm = TRUE,
 name = "proteinRobust",
 fun = MsCoreUtils::robustSummary)
 
#multidimensional scaling
limma::plotMDS(assay(pe[["proteinRobust"]]), col = colData(pe)$sample_type)
## Data analysis
### Estimation
### make msqrob models
#encode sampletype as factor
pe <- msqrob(object = pe, i = "proteinRobust", formula = ~sample_type)
getCoef(rowData(pe[["proteinRobust"]])$msqrobModels[[1]])
library(ExploreModelMatrix)
VisualizeDesign(colData(pe),~sample_type)$plotlist[[1]]

#determine the contrast and test your hypothesis
L <- makeContrast(c(
  "sample_typenot_brain=0"
), 
    parameterNames = 
    rowData(pe[["proteinRobust"]])$msqrobModels[[1]] %>% 
    getCoef %>% 
    names
  )
pe <- hypothesisTest(object = pe, i = "proteinRobust", contrast = L, overwrite=TRUE)

results <- as.data.frame(rowData(pe[["proteinRobust"]]))
contrast_names <- c("sample_typenot_brain")

for (contrast in contrast_names) {
  logFC_col <- paste0(contrast, ".logFC")
  pval_col <- paste0(contrast, ".pval")
  adj_pval_col <- paste0(contrast, ".adjPval")


  volcano_data <- results[, c("Protein.Groups", logFC_col, pval_col, adj_pval_col)]
  colnames(volcano_data) <- c("Protein.Groups", "logFC", "pval", "adjPval")
  volcano_data <- na.omit(volcano_data)
  volcano_data$negLogAdjPval <- -log10(volcano_data$adjPval)
  volcano_data$significant <- volcano_data$adjPval < 0.01

  write.csv(volcano_data, file = paste0("/home/compomics/git/Projects/MLMarker_paper/PXD007592_responders_melanoma/prediction_volcano_data_", contrast, ".csv"), row.names = FALSE)

  volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = negLogAdjPval, color = significant)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
    xlim(-7, 5) +
    ylim(-0.2, 5) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = paste("Volcano Plot:", contrast),
      x = "Log Fold-Change",
      y = "-Log10 adjusted P-value"
    )
  ggsave(
    filename = paste0("/home/compomics/git/Projects/MLMarker_paper/Figures_paper/volcano_plot_", contrast, ".png"),
    plot = volcano_plot, width = 8, height = 6, dpi = 300
  )
}

#make a heatmap of volcano_data
#extract significant proteingroups
# Extract significant protein groups
sigNames <- volcano_data %>% filter(sample_typenot_brain.adjPval < 0.01) %>% pull(Protein.Groups)

# Ensure sigNames is a character vector
sigNames <- as.character(sigNames)

# Extract expression values for significant proteins
expression_values <- assay(pe[['proteinRobust']])[rownames(assay(pe[['proteinRobust']])) %in% sigNames, ]

pheatmap(expression_values, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = TRUE, 
         fontsize = 10)