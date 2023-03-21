
# preparatory steps -------------------------------------------------------

# install/load packages -------------------------------------------------


system("apt-get -y install libxml2-dev")
system("conda install -y -c r r-curl")


.cran_packages <- c("shiny", "DT", "tidyverse", "reshape2", "stringr","dplyr", "tibble","png","shinythemes","devtools","remotes","ggsci","ggplot2","purrr","ggpubr","Hmisc","doBy","dplyr",
"taRifx","patchwork","stringr","Rmisc","janitor","optparse","here","reader","seqinr","reshape")#

.bioc_packages <- c("BiocManager", "microbiome","phyloseq","mixOmics","dada2","ALDEx2", "DECIPHER","DESeq2")#"biomformat",
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
  #devtools::install_github("grunwaldlab/metacoder")#,"metacoder"
  if(any(!.inst[2:length(.inst)])) {
    BiocManager::install(.bioc_packages[!.inst[2:length(.inst)]], ask = F)
    #BiocManager::install('metacoder')
  }
  remotes::install_github("DanielSprockett/reltools")
  remotes::install_github("mikemc/speedyseq")
  remotes::install_github("vmikk/metagMisc")
  install.packages("metacoder")
}

system("mamba install -c conda-forge r-tarifx r-mixOmics ")

sapply(c(.cran_packages, .bioc_packages ,.other_packages), require,character.only = TRUE)

system("cp -rf www/database ~")

system("conda info --env ")
