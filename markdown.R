## generate all documentation from RMD files
## html files go to dist/ to be hosted on Github pages
## pdf files go to pdf/
##
## Florian J. Auer
## https://github.com/florian-j-auer
## IT-Infrastructure for Translational Medical Research, University of Augsburg
## https://www.uni-augsburg.de/de/fakultaet/fai/informatik/prof/misit/
## 03/30/2022

## Files to process
# rmd_files = c("NDExEdit.Rmd")
rmd_files = c(
  "Overview.Rmd",
  "SoftwareRequirements.Rmd",
  "Ndexr.Rmd",
  "RCX.Rmd",
  # "RCXvis.Rmd",
  "DataExplorationAndPreparation.Rmd",
  "DataIntegration.Rmd",
  "Visualization.Rmd",
  "MetaRelSubNetVis.Rmd",
  "NDExEdit.Rmd"
)
# exclude = c("navigation.Rmd", "NDExR-and-Cytoscape.Rmd")
# 
# rmd_files = list.files(".", "*.Rmd")
# rmd_files = rmd_files[!rmd_files %in% exclude]


## Generate navigation
source("navigation.R")


## replace the data.frame rendering through DT for html output
pdf <- new.env()
pdf$HTML = FALSE

html <- new.env()
local({
  HTML = TRUE
  print.data.frame = function(x,...){DT::datatable(x, options = list(dom = 't'))}
  
  print.MetaDataAspect = function(x,...){
    cat("Meta-data:\n")
    DT::datatable(x, options = list(dom = 't'))
  }
  
  print.NodeAttributesAspect = function(x,...){
  	cat("Node attributes:\n")
  	DT::datatable(x, options = list(dom = 't'))
  }
  
  print.NetworkAttributesAspect = function(x,...){
  	cat("Network attributes:\n")
  	DT::datatable(x, options = list(dom = 't'))
  }
  
  print.EdgeAttributesAspect = function(x,...){
  	cat("Edge attributes:\n")
  	DT::datatable(x, options = list(dom = 't'))
  }
  
  head_orig = head
  head = function(x,...){
    print(head_orig(x,...))
  }
}, envir = html)

  
for(rmd_file in rmd_files){
  html$HTML_file = paste0("navigation/",rmd_file)
  rmarkdown::render(rmd_file,
                    output_format = "BiocStyle::html_document",
                    output_dir = "docs/",
                    envir = html)
  
  rmarkdown::render(rmd_file,
                    output_format = "BiocStyle::pdf_document",
                    output_dir = "pdf/",
                    envir = pdf)
}

## create the README.md
rmarkdown::render("Overview.Rmd",
									output_file = "README.MD",
									rmarkdown::md_document(variant = "markdown_github"),
									envir = pdf)

