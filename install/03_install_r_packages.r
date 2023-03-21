

system("mamba install -c  conda-forge r-tarifx r-mixOmics r-devtools r-rmpfr ")
#system("mamba  install -c bioconda r-netcomi r-metacoder bioconductor-dada2 bioconductor-microbiome bioconductor-phyloseq bioconductor-deseq2")


# Install NetCoMi


.cran_packages <- c("shiny", "DT", "tidyverse","svglite", "reshape2", "stringr","dplyr","shinyWidgets", "tibble","png","shinythemes","devtools","remotes","ggsci","ggplot2","purrr","ggpubr","Hmisc","doBy","dplyr","patchwork","stringr","Rmisc","janitor","optparse","here","reader","seqinr","reshape")#"taRifx",

.bioc_packages <- c("BiocManager", "microbiome","phyloseq","mixOmics","dada2","ALDEx2", "DECIPHER","DESeq2","microbiomeMarker")#"biomformat",
.other_packages<-c("metacoder","speedyseq","metagMisc")


chooseCRANmirror(ind=16)#lanzhou

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}



.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if(!.inst[1]) {install.packages("BiocManager")}
  #BiocManager::install('devel')
  install.packages(c("devtools","remotes"))
  devtools::install_github("grunwaldlab/metacoder")#,"metacoder"
  if(any(!.inst[2:length(.inst)])) {
    BiocManager::install(.bioc_packages[!.inst[2:length(.inst)]], ask = F)
    #BiocManager::install('metacoder')
  }
  remotes::install_github("DanielSprockett/reltools")
  remotes::install_github("mikemc/speedyseq")
  remotes::install_github("vmikk/metagMisc")
  remotes::install_github("yiluheihei/microbiomeMarker")
  devtools::install_github("stefpeschel/NetCoMi",  dependencies = c("Depends", "Imports", "LinkingTo"),   repos = c("https://cloud.r-project.org/", BiocManager::repositories()))
  install.packages("metacoder")
}



sapply(c(.cran_packages, .bioc_packages ,.other_packages), require,character.only = TRUE)
