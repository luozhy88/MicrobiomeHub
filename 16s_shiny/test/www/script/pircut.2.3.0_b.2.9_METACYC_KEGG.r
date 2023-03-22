# PICRUSt2
#PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States) is a software for predicting functional abundances based only on marker gene sequences.
#"Function" usually refers to gene families such as KEGG orthologs and Enzyme Classification numbers, but predictions can be made for any arbitrary trait. Similarly, predictions are typically based on 16S rRNA gene sequencing data, but other marker genes can also be used.
#**Note** picrust2 must be run in set_work_path(https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.2.0-beta)).

print("/usr/bin/Rscript /home/lzy/data/script/Pircust/pircut.2.3.0_b.2.9_METACYC_KEGG.r -f /home/lzy/data/metagenome/analysis_project/liudiao_H22_liuyuan_151samples_20200324/tt -s /home/lzy/data/metagenome/analysis_project/liudiao_H22_liuyuan_151samples_20200324/tt/data/FALSE2020_03_20_5_update.rds2020_03_24_2_update.rds
Note: need change activate_local
")

#.libPaths(c("/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library"))
###########Load all packages
library('dada2')
library(ShortRead)
library(Biostrings)
# library(DECIPHER) # only used when switching from naive bayesian to this new approach for taxa assignment
library(phyloseq)
library(ggplot2)
library(microbiome)
library(ape)
library("biomformat")
library(plyr) # only used when phyloseq class is used for beta diversity analysis.
library(tidyr)
library(dplyr)
library(reader)
library("tibble")
library("purrr")
library("reltools")
library(seqinr)
# library("xlsx")
library(rmarkdown)
library(optparse)
library(reshape) # 加载所需的包
###########Load all packages
.libPaths()
#.libPaths("/home/lzy/data/software/anaconda3/envs/pircust2/lib/R/library")
cpus_num <<- 3
################parse_args
option_list <- list(
  make_option(c("-f", "--first"), type = "character", default = FALSE,
              action = "store", help = "This is first! read the work_path"
  ),
  make_option(c("-s", "--second"), type = "character", default = FALSE,
              action = "store", help = "This is second! read the phyloseq_path"
  ), 
              make_option(c("-o", "--out"), type = "character", default = FALSE,
                          action = "store", help = "This is four! out_dir_path"

              # make_option(c("-h", "--help"), type = "logical", default = FALSE,
              #             action = "store_TRUE", help = "This is Help!"
              # )
              ))

  opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))
  print(opt)
  #print(opt$f)
  #
  #
  #
  set_work_path<-opt$f
  phyloseq_path<-opt$s
  out_dir_path<-opt$o
  ################parse_args
  #
  #
  
  ################## Read file
  # set_work_path<-'/home/lzy/data/metagenome/mmu_971/MJ20171226003.yansiyuan.338F_806R.63samples.data'
  # fitGTR_tree <- readRDS(paste0(set_work_path,"/output/tree/tree.rds"))#read all_tree
  # phyloseq_path<-paste0('/home/lzy/data/metagenome/mmu_971/merge_all_phyloseq',"/all_phyloseq_update.rds")
  
  ##set_work_path<-"/home/lzy/data/metagenome/analysis_project/liudiao_H22_liuyuan_151samples_20200324"
  ##phyloseq_path<-paste0(set_work_path,"/data/FALSE2020_03_20_5_update.rds2020_03_24_2_update.rds")
  
  
  setwd(set_work_path)
  all_phyloseq <- readRDS(phyloseq_path) 
  # phy_tree(all_phyloseq)<-fitGTR_tree
  all_phyloseq
  all_phyloseq<-prune_taxa(taxa_sums(all_phyloseq)>0,all_phyloseq)
  #all_phyloseq <- microbiome::transform(all_phyloseq, "compositional")
  all_phyloseq
  ################## Read file
  
  
  #################### 生成中间过程文件
  dir.create(paste0(set_work_path,"/data"))
  dir.create(paste0(set_work_path,"/src"))
  dir.create(paste0(set_work_path,"/analysis"))
  dir.create(paste0(set_work_path,"/output"))
  dir.create(paste0(set_work_path,"/ext"))
  setwd(set_work_path)
  
  all_phyloseq
  phy_otu<-otu_table(all_phyloseq)
  DT::datatable(as.data.frame(head(phy_otu)))
  # view(as.data.frame(head(phy_otu)))
  # write.table(sample_data(all_phyloseq),paste0(set_work_path,'/data/metadata_table.txt'),quote =F,sep ='\t',row.names = FALSE)
  # write_phyloseq(all_phyloseq,type='METADATA',path =paste0(set_work_path,'/data') )
  # saveRDS(all_phyloseq, paste0(set_work_path,"/data/phyloseq.rds"))
  #################### 生成中间过程文件

  
  ## 提取fasta,biom
  # meta_info<- read.table(meta_info,sep="\t",header = T)# if you change the name of ID ,we can use it.
  # meta_info<- read.csv(paste0(set_work_path,'/data/metadata_table.csv'),sep=",", fileEncoding="UTF-16LE")
  # meta_info <- rename(meta_info,c(X = "ID")) 
  
  
  all_phyloseq
  
  # setdiff(rownames(Meta_tab),row.names(as.data.frame(as.matrix(sample_data(all_phyloseq)))))
  # setdiff(rownames(Meta_tab),meta_info$X)
  # phy_<-'101027' %in% row.names(as.data.frame(as.matrix(sample_data(all_phyloseq))))
  
  
  #create otu
  OTU1 = t(as(otu_table(all_phyloseq), "matrix"))
  # view(head(OTU1))
  
  
  #convert the table of all_otu to biom
  # otu<-t(as(otu_table(OTU1,taxa_are_rows = FALSE),"matrix"))
  otu<-as(otu_table(OTU1,taxa_are_rows = FALSE),"matrix")
  #DT::datatable(as.data.frame(otu))
  # view(head(otu))
  otu_biom<-make_biom(data=otu)
  dir.create(paste0(set_work_path,'/output/pircust2/otu_biom'),recursive = TRUE)
  write_biom(otu_biom, paste0(set_work_path,"/output/pircust2/otu_biom/otu_biom.biom"))




  
  # #add tree
  # isnot_build_tree=TRUE
  # if (isnot_build_tree){
  # fitGTR_tree <- readRDS(paste0('/home/lzy/data/metagenome/hsa_meiji_hsa',"/merge_analysis/tree/all_tree.rds"))#read all_tree
  # phy_tree(all_phyloseq)<-fitGTR_tree
  # }
  # all_phyloseq
  
  
  #output tree
  # tree1 = phy_tree(all_phyloseq)#output  tree
  # ape::write.tree(tree1, "tree1.tre")
  
  #output fasta
  df_<-as.data.frame(otu_table(all_phyloseq))
  # %>% colnames() %>% data.frame()
  seqtab<-df_
  
  seqnum <- paste0("Seq", seq(ncol(seqtab)))
  uniqueSeqs <- as.list(colnames(seqtab))
  write.fasta(uniqueSeqs, seqnum,  paste0(set_work_path,"/output/pircust2/uniqueSeqs.fasta"))
  
  seqnum <- as.list(colnames(seqtab))
  write.fasta(uniqueSeqs, seqnum,  paste0(set_work_path,"/output/pircust2/uniqueSeqs.fasta"))
  
  
nf_core=TRUE
if (nf_core==TRUE){
        print("run nf-core output")
	#########nf-core phyloseq to biom
	library("biomformat")
	library("reltools")
	library(seqinr)
	OTU1 = t(as(otu_table(all_phyloseq), "matrix"))
	otu<-as(otu_table(OTU1,taxa_are_rows = TRUE),"matrix") %>% t()
	head(otu)
	otu_biom<-make_biom(data=otu)

	# dir.create(paste0(set_work_path,'/output/pircust2/otu_biom'),recursive = TRUE)
	write_biom(otu_biom, paste0(set_work_path,"/output/pircust2/otu_biom/otu_biom.biom"))
	#########nf-corephyloseq to fasta  
	#output fasta
	phyloseq_new.seq<-refseq(all_phyloseq)  %>% as.list()
	write.fasta(phyloseq_new.seq , names(phyloseq_new.seq),paste0(set_work_path,"/output/pircust2/uniqueSeqs.fasta"))
}




##################### Run picrust2
  ##NOTE:if you can not run picrust2,please source ~/.bashrc or conda activate picrust2.1 picrust2.3 in linux 
str<-sprintf("
#if you can not run picrust2,please source ~/.bashrc or conda activate picrust2.1 picrust2.3 in linux
#source /root/anaconda3/bin/activate pircust2
#conda info --env
#########install#########
#conda create -n pircust2
#conda install -y -c bioconda -c conda-forge picrust2=2.3.0_b #R3.6.1
#conda activate pircust2
#########install#########

#old
  #activate_local=/home/lzy/data/software/anaconda3/bin
  #source $activate_local/activate pircust2
  conda info --env
#new
  #
  #activate_local=/home/zhiyu/data/software/Anaconda3/bin
  activate_local=/home/zhiyu/data/software/Anaconda3/bin
  source $activate_local/activate pircust2
  # source $activate_local/activate pircust2_ceshi
  conda info --env



########################KEGG
echo '####################################KEGG####################################'
#if you can not sun picrust2,please source ~/.bashrc in linux
pwd
mkdir -p output/pircust2/KEGG/picrust2_out_pipeline

#Place reads into reference tree
place_seqs.py -s output/pircust2/uniqueSeqs.fasta -o output/pircust2/KEGG/picrust2_out_pipeline/out.tre -p %s --intermediate output/pircust2/KEGG/picrust2_out_pipeline/intermediate/place_seqs

#Hidden-state prediction of gene families
hsp.py -i 16S -t output/pircust2/KEGG/picrust2_out_pipeline/out.tre -o output/pircust2/KEGG/picrust2_out_pipeline/marker_predicted_and_nsti.tsv.gz -p %s -n
hsp.py -i KO -t output/pircust2/KEGG/picrust2_out_pipeline/out.tre -o output/pircust2/KEGG/picrust2_out_pipeline/EC_predicted.tsv.gz -p %s

#Generate metagenome predictions
# <!-- #Pathway-level inference -->
################KO
##########not strat_out
echo 'KEGG not strat_out'
metagenome_pipeline.py -i output/pircust2/otu_biom/otu_biom.biom -m output/pircust2/KEGG/picrust2_out_pipeline/marker_predicted_and_nsti.tsv.gz -f output/pircust2/KEGG/picrust2_out_pipeline/EC_predicted.tsv.gz -o output/pircust2/KEGG/picrust2_out_pipeline/EC_metagenome_out_unstrat 

# convert_table.py output/pircust2/KEGG/picrust2_out_pipeline/EC_metagenome_out_unstrat/pred_metagenome_unstrat.tsv.gz -c contrib_to_legacy -o output/pircust2/KEGG/picrust2_out_pipeline/EC_metagenome_out_unstrat/pred_metagenome_unstrat.legacy.tsv.gz

pathway_pipeline.py -i output/pircust2/KEGG/picrust2_out_pipeline/EC_metagenome_out_unstrat/pred_metagenome_unstrat.tsv.gz -o output/pircust2/KEGG/picrust2_out_pipeline/KEGG_pathways_unstrat -p %s --no_regroup --map ~/database/KEGG_pathways_to_KO.tsv 

add_descriptions.py -i output/pircust2/KEGG/picrust2_out_pipeline/KEGG_pathways_unstrat/path_abun_unstrat.tsv.gz --custom_map_table ~/database/KEGG_pathways_info.tsv -o output/pircust2/KEGG/picrust2_out_pipeline/KEGG_pathways_unstrat/path_abun_unstrat_descrip.tsv.gz


gunzip output/pircust2/KEGG/picrust2_out_pipeline/KEGG_pathways_unstrat/path_abun_unstrat_descrip.tsv.gz


#python /home/lzy/data/script/python/email_send.py pircust_ok
exit 0

",cpus_num,cpus_num,cpus_num,cpus_num,cpus_num)

  # 2.3.0_b version is not arg (--metagenome_contrib)  in metagenome_pipeline.py
  getwd()
  #system(str)
  # writeLines(str,"/home/lzy/data/script/picrust2.bash")#run the bash in linux;

  #activate_local<-"/root/anaconda3/bin"
  #conda_activate_commond<-paste0(" \n source ",activate_local,"/activate pircust2 \n conda info --env \n")  
  #str<-paste0(str,conda_activate_commond)
  #print(str)

  writeLines(str, paste0(set_work_path,"/src/picrust2.bash"))#run the bash in linux;
  ##################### Run picrust2  conda activate web_test
  sys.commnd<-paste0("/bin/bash -x ", "src/picrust2.bash")
  print(sys.commnd)
  #system(sys.commnd)
  #system(paste0("Rscript ", paste0(set_work_path,"/uu.r") )  )
  
  ##The log of picrust2 will be stored in picrust2.log.
  #read.csv(paste0(set_work_path,'/picrust2.log'))
  #file.copy("/home/lzy/data/script/Pircust/pircut2.9_METACYC_KEGG.r" ,paste0(set_work_path,"/src/pircut2.9_METACYC_KEGG.r" ))
  print('pircust yunxingwan')
  quit()
  
  
  
