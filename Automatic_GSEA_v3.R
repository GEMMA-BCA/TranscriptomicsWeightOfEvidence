####### Libraries required:
require(optparse)
require(clusterProfiler)
require(tidyverse)
require(tools)
require(readxl)
require(writexl)
require(org.Hs.eg.db)

###### Script options:
rm(list=ls())

option_list = list(
  make_option(c("-p", "--input_file_path"), 
              type = "character", 
              default=NULL,
              help="Working directory where input files are stored", 
              metavar="character"),
  
  make_option(c("-d", "--hallmark_database_file"), 
              type = "character", 
              default="~/Dropbox/Luca/BCA-UniPd/Maci/WOE/Hallmark_database.gmt", 
              help="Name (and full path) of the hallmark database file that contains the Hallmark pathways and the genes that belong to each pathway. It can be a '.gmt' file or an Excel '.xlsx' file with two columns, named 'term' and 'gene' respectively [default= %default]", 
              metavar = "character"),
  
  make_option(c("-c","--hallmark_categories_file"), 
              type="character", 
              default = "~/Dropbox/Luca/BCA-UniPd/Maci/WOE/Hallmark_categories_weights.xlsx", 
              help = "Name (and full path) of the file containing subdivision of Hallmark into categories and the weight of each hallmark [default=%default]", 
              metavar="character"),
  
  make_option(c("-o", "--output_folder"), 
              type="character", 
              default=NULL, 
              help="Output folder name. This option will create a folder with the name that you specify and store the results into that folder. [default= %default]", 
              metavar="character"),
  
  make_option(c("-f", "--input_file"), 
              type="character", 
              default=NULL, 
              help="Name of the file that contains the full Differential Expression results table. It can be multiple file names separated by comma ',' but all files have to be in the same directory", 
              metavar="character")
  
); 

# Example usage of the script:
opt_parser = OptionParser(usage="usage: Rscript %prog --input_file_path Path_to_Input_files --output_folder Results/ --hallmark_database_file Path_to_hallmark_file/Hallmark.gmt --hallmark_categories_file Path_to_hallmark_weights/Hallmark_categories_weights.xlsx --input_file Experiment1.xlsx,Experiment2.xlsx,ExperimentN.xlsx",
                          option_list=option_list);

# Parsing script options:
opt = parse_args(opt_parser)

# Setting working directory:
setwd(opt$input_file_path)

# Define input files:
inputfiles = as.list(strsplit(opt$input_file, ",")[[1]])

###### Analysis:
for (numerofile in c(1:length(inputfiles))){
  
  ## ~~~~~ Import files: ~~~~~ ##
  cat("\n\nProcessing: ", inputfiles[[numerofile]], "\n\n")
  input = inputfiles[[numerofile]]
  
  # Load data:
  if(tools::file_ext(input) == "csv"){
    inputdata=read.csv(input)
  } else {
    inputdata=readxl::read_xlsx(input)
  }
  cat("The imported data contains ", dim(inputdata)[1], " genes\n")
  
  # Output folder:
  out=sub("/$", "", opt$output_folder) #This removes "/" at the end of the string
  dir.create(out)
  output=paste0(out,"/results_",basename(file_path_sans_ext(input)))
  dir.create(output)
  
  
  ## ~~~~~ Process input data and prepare gene lists: ~~~~~ ##
  # Reference organism:
  org="org.Hs.eg.db"
  
  #Convert ids to symbols and fill "NAs":
  dataconv=bitr(inputdata$Geneid, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org)
  inputdata.symbol = merge(inputdata,dataconv, by.x="Geneid",by.y="ENSEMBL",all.x=TRUE)
  inputdata.symbol$SYMBOL[is.na(inputdata.symbol$SYMBOL)] <- as.character(inputdata.symbol$Geneid)[is.na(inputdata.symbol$SYMBOL)]
  
  #Remove Geneid column and get summary of duplicated:
  inputdata.symbol = inputdata.symbol[,-1]
  cat("\n\n Total number of genes: ", dim(inputdata.symbol)[1], "\n Duplicated Symbols:\n")
  summary(duplicated(inputdata.symbol$SYMBOL))
  
  #Filter to keep highest |foldChange|:
  cat("Now removing duplicated IDs/Symbols...\n")
  inputdata.final = inputdata.symbol %>% group_by(SYMBOL) %>% filter(FDR == min(FDR)) %>% filter(abs(logFC) == max(abs(logFC))) %>% ungroup
  inputdata.final = inputdata.final[order(inputdata.final$FDR),]
  inputdata.final = inputdata.final[!duplicated(inputdata.final$SYMBOL),]
  inputdata.final = data.frame(inputdata.final)
  inputdata.final = inputdata.final[,c(6,1,2,3,4,5)]
  cat("\nDuplicated Symbols:\n")
  summary(duplicated(inputdata.final$SYMBOL))
  cat("\nThe dataframe (after removal of duplicated genes) now contains ", dim(inputdata.final)[1], " genes\n")
  
  #Create gene lists:
  geneList_logFC <- inputdata.final[,2] #feature 1: numeric vector
  names(geneList_logFC) <- as.character(inputdata.final[,1]) #feature 2: named vector
  geneList_logFC <- sort(geneList_logFC, decreasing = TRUE) #feature 3: decreasing order
  geneList= geneList_logFC
  
  ## ~~~~~ Load Hallmark data ~~~~~ ##
  #Parameters:
  minsizeGS = 5 #minimum gene set size for inclusion in GSEA analysis
  maxsizeGS= 500 #maximum gene set size for inclusion in GSEA analysis
  pvalCutoffGS=1 #adjusted p-val threshold. Set to 1 to have the results exported for all pathways.
  
  #Create hallmark dataset:
  if(tools::file_ext(opt$hallmark_database_file) == "gmt"){
    hallmark = read.gmt(opt$hallmark_database_file)
  } else if (tools::file_ext(opt$hallmark_database_file) == "xlsx") {
    hallmark = readxl::read_excel(opt$hallmark_database_file)
  } else {
    cat("Make sure the gsea database is either a '.gmt' or a '.xlsx' file")
  }
  hallmark$term = as.factor(hallmark$term)
  
  #Add subcategories:
  if(tools::file_ext(opt$hallmark_categories_file) == "csv"){
    Hallmark_categories <- read_csv(opt$hallmark_categories_file) %>% distinct()
  } else if (tools::file_ext(opt$hallmark_categories_file) == "xlsx") {
    Hallmark_categories <- readxl::read_excel(opt$hallmark_categories_file) %>% distinct()
  } else {
    cat("Make sure the file is a '.xlsx' or '.xsv' file")
  }
  colnames(Hallmark_categories)[1:3] = c("new_name", "process","Weight")
  
  #Create empty dataframe to collect results:
  results = setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalues","rank","leading_edge","core_enrichment"))
  
  
  ## ~~~~~ PERFORM GSEA: ~~~~~ ##
  # Useful for subsequent plot:
  ListaGeni = inputdata.final[,c(1,2)]
  colnames(ListaGeni) <- c("ID","logFC")
  
  # GSEA:
  for (process in levels(as.factor(Hallmark_categories$process))){ #General version for the entire analysis
    GSEA_subset = as.data.frame(Hallmark_categories[Hallmark_categories$process %in% process,]) # We perform GSEA by hallmark's categories
    GSEA_data = as.data.frame(hallmark[as.character(hallmark$term) %in% GSEA_subset$new_name,])
    cat("Analysing: ", process, "\n")
    customGSEA<- GSEA(geneList=geneList, 
                      TERM2GENE = GSEA_data, 
                      minGSSize = minsizeGS,
                      pvalueCutoff = pvalCutoffGS, 
                      maxGSSize = maxsizeGS,
                      verbose = TRUE)
    
    
    #Append results to the empty result dataframe
    results = rbind(results, customGSEA@result)
  }
  
  # Merge everything:
  final_output = merge(Hallmark_categories, results, by.x = "new_name", by.y = "ID")
  write_xlsx(final_output, path=paste0(output,"/results_logFC_GSEA_pval-cutoff",pvalCutoffGS,".xlsx"))
  
}