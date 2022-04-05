## generate navigation RMD files in navigation/
##
## Florian J. Auer
## https://github.com/florian-j-auer
## IT-Infrastructure for Translational Medical Research, University of Augsburg
## https://www.uni-augsburg.de/de/fakultaet/fai/informatik/prof/misit/
## 03/30/2022

## Generate navigation
pages = yaml::read_yaml("navigation.yml")

generateRMD = function(pages, indent=0){
  iC = "    "
  results = c()
  for(p in pages){
    results = c(results, paste0(rep(iC, indent),"- [",p$title, "](",p$file,".html) *",p$description,"*\n"))
    if(!is.null(p$subpages)) results = c(results, generateRMD(p$subpages, indent + 1))
  }
  return(results)
}

sink("navigation/navigation.Rmd")
cat(paste(generateRMD(pages)))
sink()

generateLinks = function(pages){
  results = c()
  for(p in pages){
    results = c(results, paste0("[",p$title, "](",p$file,".html)"))
    if(!is.null(p$subpages)) results = c(results, generateLinks(p$subpages))
  }
  return(results)
}

links = generateLinks(pages)
generateFiles = function(pages, total, counter = 1){
  for(p in pages){
    sink(paste0("navigation/",p$file,".Rmd"))
    cat('<div class="pagination">\n')
    if(counter != 1) cat(paste0('<div class="prev">Previous: ',links[counter-1], '</div>\n'))
    cat(paste0('<div class="curr">', p$title, '</div>\n'))
    if(counter != total) cat(paste0('<div class="next">Next: ',links[counter+1], '</div>\n'))
    cat('</div>\n')
    sink()
    counter = counter + 1
    if(!is.null(p$subpages)) counter = generateFiles(p$subpages, total, counter)
  }
  invisible(counter)
}

countFiles = function(pages, counter = 1){
  for(p in pages){
    counter = counter + 1
    if(!is.null(p$subpages)) counter = countFiles(p$subpages, counter)
  }
  return(counter)
}
total = countFiles(pages)

generateFiles(pages, total)
