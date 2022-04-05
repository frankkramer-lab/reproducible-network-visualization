#===============================================================#
#       GENERATION OF MOLECULAR TUMOR BOARD (MTB) REPORT        #
#---------------------------------------------------------------#
#   This script filters SNVs and CNVs using gene-drug public    #
#   databases. Then classifies the variants into levels of      #
#   evidence and finally presents the results in a pdf report   #
#===============================================================#
#   Julia Perera-Bel                                            #
#   Florian J. Auer                                             #
#===============================================================#

require(stringr)
require(xtable)
require(pander)
require(timeSeries)
require(devtools)


#' Generate MTB results for patients
#'
#' @param patients patient ids (must be present as columns in exprs.quant and relevant.genes)
#' @param exprs.quant gene expression levels (HIGH,NORMAL,LOW)
#' @param relevant.genes list of relevant genes for each patient
#' @param exprs.quant.gene.column column name of the gene name in `exprs.quant`
#' @param cancer name of the cancer
#' @param verbose print some information
#'
#' @return list of mtb results (data.frame) for each patient
mtb_analysis = function(patients,
                        exprs.quant,
                        relevant.genes,
                        cancer = "BRCA",
                        exprs.quant.gene.column = "probe",
                        verbose = FALSE){
  
  ######################################
  ## Load additional data and scripts ##
  ######################################
  
  ## Load helper functions
  # source_url("https://raw.githubusercontent.com/jperera-bel/MTB-Report/master/helpers/get_druggable.r")
  # source_url("https://raw.githubusercontent.com/jperera-bel/MTB-Report/master/helpers/get_levels.r")
  ## Local copy
  source("MTB/get_druggable.r")
  source("MTB/get_levels.r")
  
  
  ## Cancer type synonyms between databases
  # synonyms = read.csv(
  #   url("https://raw.githubusercontent.com/jperera-bel/MTB-Report/master/data/cancer_types.csv"), 
  #   header = TRUE,
  #   sep="\t"
  # )
  ## Local copy
  synonyms = read.csv("MTB/cancer_types.csv", 
                      header = TRUE, 
                      sep="\t")
  
  ## Get cancer type synonyms in the databases
  cancerSynonyms = synonyms[synonyms$tcga_cancer==cancer,]
  cancer_GDKD = cancerSynonyms$knowledge
  cancer_CIVIC = unique(unlist(strsplit(cancerSynonyms$civic,",")))
  
  if(verbose){
    cat(paste0('Cancer synonyms for "',cancer,'":\n', collapse = ""))
    cat(paste0(c(cancer_CIVIC,cancer_GDKD), collapse = ", "),"\n")
  }
  
  ## DBs
  # CIVIC = read.csv(
  #   url("https://raw.githubusercontent.com/jperera-bel/MTB-Report/master/data/CIViC.csv"), 
  #   header = TRUE,
  #   sep="\t"
  # )
  # GDKD = read.csv(
  #   url("https://raw.githubusercontent.com/jperera-bel/MTB-Report/master/data/GDKD.csv"), 
  #   header = TRUE,
  #   sep="\t"
  # )
  # TARGET = read.csv(
  #   url("https://raw.githubusercontent.com/jperera-bel/MTB-Report/master/data/TARGET_MERIC.csv"), 
  #   header = TRUE,
  #   sep="\t"
  # )
  ## Local copy
  CIVIC = read.csv("MTB/CIViC.csv", 
                   header = TRUE, 
                   sep="\t")
  GDKD = read.csv("MTB/GDKD.csv", 
                  header = TRUE, 
                  sep="\t")
  TARGET = read.csv("MTB/TARGET_MERIC.csv", 
                    header = TRUE, 
                    sep="\t")
  
  
  ######################
  ## Process Patients ##
  ######################
  
  mtb.results = list()
  for (patient in patients){
    if(verbose){
      cat("Patient: ", patient," ") 
    }
    
    #####################
    ## Read genes list ##
    #####################
    
    ## Assume all are "missense mutation"
    SNV_missense = data.frame(Hugo_Symbol=relevant.genes[[patient]],
                              Variant_Classification="missense mutation",
                              Protein_Change = "")
    
    ## Assume all are "nonsense mutation"
    SNV_nonsense = data.frame(Hugo_Symbol=relevant.genes[[patient]],
                              Variant_Classification="nonsense mutation",
                              Protein_Change = "")
    
    ## Assume all are amplified
    CNV_amplification = data.frame(Hugo_Symbol=relevant.genes[[patient]],
                                   cn_alteration = "amplification")
    
    ## Assume all are deleted
    CNV_deletion = data.frame(Hugo_Symbol=relevant.genes[[patient]],
                              cn_alteration = "deletion")
    
    if( nrow(SNV_missense) == 0 & nrow(SNV_nonsense) == 0 & nrow(CNV_amplification) == 0 & nrow(CNV_deletion) == 0) {
    	warning('The patient has no SNVs nor CNVs. Not worth continuing with the analysis!')
    }
    
    ######################################
    ## Filter SNVs and CNVs by database ##
    ######################################
    
    ## GDKD DB
    druggableGDKD = data.frame()
    druggableGDKD = rbind(
      match_SNV_GDKD(SNV_missense, db = GDKD),
      match_SNV_GDKD(SNV_nonsense, db = GDKD)
    )
    druggableGDKD = rbind(
      druggableGDKD,
      match_CNV_GDKD(CNV_amplification, db = GDKD),
      match_CNV_GDKD(CNV_deletion, db = GDKD)
    )
    druggableGDKD = rbind(
      druggableGDKD,
      match_WT_GDKD(SNV_missense, CNV_amplification, cancer_GDKD, db = GDKD)
    )
    rownames(druggableGDKD) = NULL
    
    ## CIVIC
    druggableCIVIC = data.frame()
    druggableCIVIC = unique(rbind(
      match_SNV_CIVIC(SNV_missense, db = CIVIC),
      match_SNV_CIVIC(SNV_nonsense, db = CIVIC)
    ))
    druggableCIVIC = rbind(
      druggableCIVIC,
      match_CNV_CIVIC(CNV_amplification, db = CIVIC),
      match_CNV_CIVIC(CNV_deletion, db = CIVIC)
    )
    rownames(druggableCIVIC) = NULL
    
    
    ## TARGET DB
    druggableTARGET = data.frame()
    druggableTARGET = rbind(
      match_TARGET_MERIC(SNV_missense, CNV_amplification, db = TARGET), 
      match_TARGET_MERIC(SNV_nonsense, CNV_deletion, db = TARGET)
    )
    
    if( nrow(druggableGDKD) == 0 & nrow(druggableCIVIC) == 0 & nrow(druggableTARGET) == 0) {
      warning('No gene-drug interations were found. Not worth continuing with the analysis!')
    }
    
    
    #########################################################
    ## Classify filtered variants by Levels of Evidence    ##
    #########################################################
    
    ## KNWOLEDGE
    levelsGDKD = get_levels_GDKD(druggableGDKD, cancer_GDKD)
    
    ## CIVIC
    levelsCIVIC = get_levels_CIVIC(druggableCIVIC, cancer_CIVIC)
    
    levels = merge_levels(levelsGDKD, levelsCIVIC)
    
    ## Homogeneize/clean the final table
    table = clean_levels(levels, synonyms, sort_by="genes")
    
    #################################################
    ## Generate report from table and patient data ##
    #################################################
    
    ## Fix format issues
    table$`Pat Var` = sapply(table$`Pat Var` , function(x) gsub("NEW", "", x))
    table$`Pat Var` = sapply(table$`Pat Var` , function(x) gsub("([A-Z]) ([A-Z])", "\\1, \\2", x))
    table$`Pat Var` = sapply(table$`Pat Var` , function(x) gsub(" ", "", x))
    table = unique(table)
    
    #############################
    ## TARGET with other genes ##
    #############################
    
    ## recover genes not in previous table
    druggableTARGET = as.data.frame(
      lapply(druggableTARGET, as.character),
      stringsAsFactors = FALSE
      )
    
    if (nrow(druggableTARGET) != 0){
      new_gene = !(druggableTARGET$Gene %in% table$Gene)
      druggableTARGET = druggableTARGET[new_gene,]
      druggableTARGET = druggableTARGET[,c(1,2,4,3,5)]
      
      if (nrow(druggableTARGET) != 0){
        
        ## clean
        druggableTARGET$Patient_variant = gsub(pattern = "amplification",
                                               replacement = "amp.",
                                               druggableTARGET$Patient_variant)
        druggableTARGET$Patient_variant = gsub(pattern = "deletion",
                                               replacement = "del.",
                                               druggableTARGET$Patient_variant)
        
        ## adjust colnames
        colnames(druggableTARGET) <- c("Gene", "Pat Var", "Known Var", "Predicts", "Drugs")
        
        ## remove duplicated rows
        druggableTARGET = druggableTARGET[!duplicated(druggableTARGET[,-2]),]
        ## keep only predictive (Drugs column full)
        druggableTARGET = druggableTARGET[druggableTARGET$Drugs!="",]
        
        if (nrow(druggableTARGET) != 0){
          ## merge table + target
          result = rbind(table,
                         setNames(data.frame(matrix(NA,nrow = nrow(druggableTARGET),ncol = ncol(table))),
                                  colnames(table))) 
          result[(nrow(table)+1):(nrow(table)+nrow(druggableTARGET)),c(1,2,4,5,6)] = druggableTARGET
        }else{result=table}
      }else{result=table}
    }else{result=table}
    
    ################################
    ## Filter and combine results ##
    ################################
    
    ## add expression
    result$`Patient.Expr` = exprs.quant[match(result$Gene,exprs.quant[,exprs.quant.gene.column]),
                                        patient]
    
    ## match expression with effect observed in evidence
    ## remove normal unless "Known Var=="expression"
    normal = subset(result,
                    subset = Patient.Expr=="NORMAL"&`Known Var`=="expression")
    
    ## match HIGH with GoF
    high = subset(result,
                  subset = Patient.Expr=="HIGH")
    retain = unique(high$Gene[grep("GoF|Amplification|ampl\\.|overexpression",
                                   high$`Known Var`)])
    high = subset(result,
                  subset = Gene %in% retain)
    
    ## match LOW with LoF
    low = subset(result,
                 subset = Patient.Expr=="LOW")
    retain = unique(high$Gene[grep("LoF|deletion|del\\.",
                                   high$`Known Var`)])
    low=subset(low,
               subset = Gene %in% retain)
    
    ## merge normal, high and low
    result = rbind(normal, high, low)
    
    result = result[,c(1,10,3:9)]
    
    if(verbose){
      cat("has ", nrow(result), " results.\n") 
    }
    
    # save
    mtb.results[[patient]] = result
  }
  
  return(mtb.results)
}


