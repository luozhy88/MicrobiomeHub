print("/usr/bin/Rscript /home/zhiyu/data/script/merge_metafrom_phyloseq_merge_allphyloseq2.R -s phyloseq  -t meta.table.csv")
#.libPaths("/home/zhiyu/R/x86_64-pc-linux-gnu-library/4.0")
#print(.libPaths())
#message(.libPaths())

# library('dada2')
# library(ShortRead)
# library(Biostrings)
# library(DECIPHER) # only used when switching from naive bayesian to this new approach for taxa assignment
library(phyloseq)
library(ggplot2)
library(microbiome)
library(ape)
library("biomformat")
library(plyr) # only used when phyloseq class is used for beta diversity analysis.
library(tidyr)
library(dplyr)
library("tibble")
library("purrr")
library(seqinr)
library(optparse)
library(metagMisc)
library(here)
# ÃèÊö²ÎÊýµÄ½âÎö·½Ê½
option_list <- list(
  make_option(c("-f", "--first"), type = "character", default = FALSE,
              action = "store", help = "This is first!  phyloseq_work_dir"
  ),
  make_option(c("-s", "--second"), type = "character", default = FALSE,
              action = "store", help = "This is second! read the phyloseq_path"
  ),
  make_option(c("-t", "--third"), type = "character", default = FALSE,
              action = "store", help = "This is third! meta_filename_path"
  ),
  make_option(c("-o", "--out"), type = "character", default = FALSE,
              action = "store", help = "This is four! out_dir_path"
  )
  # make_option(c("-h", "--help"), type = "logical", default = FALSE,
  #             action = "store_TRUE", help = "This is Help!"
  # )
)
# ½âÎö²ÎÊý
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))
print(opt)
#print(opt$f)
#
#
#
phyloseq_workdir_path<-opt$f
phyloseq_path<-opt$s
meta_filename_path<-opt$t
out_dir_path<-opt$o

time_<-format(Sys.Date(), "%Y_%m_%d_%w")
# phyloseq_workdir_path<-'/home/zhiyu/data/metagenome/mmu_different_age/MJ20190627044.yansiyuan.338F_806R.100samples.data'
# phyloseq_path<-'/home/zhiyu/data/metagenome/mmu_different_age/MJ20190627044.yansiyuan.338F_806R.100samples.data/pircust2/data/phyloseq.rds'
# meta_filename_path<-'/home/zhiyu/data/metagenome/mmu_different_age/MJ20190627044.yansiyuan.338F_806R.100samples.data/data/metainfo_chang1.csv'

##################################################
##Êý×Ö¿ªÍ·µÄID»á±¨´í£¬ËùÒÔ¶ÔÊý×Ö¿ªÍ·µÄID¼ÓÉÏ'X'
phyloseq_rename<-function(x){
  vec<-c()
  for (i in 1:length(x)){
    x[i]<-ifelse(grepl("^[[:digit:]]",x[i]),paste0('X',x[i]),x[i])
    vec<-append(vec,x[i])
    # print(x[i])
  }
  return(vec)
}

#ÈÚºÏÄ³¸öÄ¿Â¼ÏÂµÄÈ«²¿phyloseq
if (phyloseq_workdir_path !='FALSE') {
  #phyloseq_workdir_path='/home/zhiyu/data/metagenome/mmu_tumor_16S_analysis_dada2'
  separte_phyloseq<-list.files(path =phyloseq_workdir_path,pattern = ".+hyloseq.rds$",recursive =TRUE ,full.names = TRUE) #
  
  # separte_phyloseq<-list.files(path =phyloseq_workdir_path,pattern = ".+phyloseq/phyloseq.rds",recursive =TRUE ,full.names = TRUE)#
  
  list_phy<-c()
  for (t in 1:length(separte_phyloseq)){
    print(t)
    path_phy=separte_phyloseq[t]
    print(paste0(separte_phyloseq[t],':')) 
    st_ph <- readRDS(path_phy)
    st_ph@phy_tree<-NULL
    st_ph@sam_data<-NULL
    print(st_ph)
    # test<-list(st)
    # list_phy[[separte_phyloseq[t]]] <- st_ph#
    list_phy<-append(list_phy,st_ph)
    # append(list_otu,st)
      
  }
  all_phyloseq <-  do.call(merge_phyloseq,list_phy)#
  phyloseq<-all_phyloseq
  out_file1<-paste0(basename(phyloseq_path),'all_phyloseq_sample.name',time_,'.csv')
  write.csv(sample_names(phyloseq),paste0(out_dir_path,'/',out_file1),quote = F)
  
  sample_names(phyloseq)<-phyloseq_rename(sample_names(phyloseq)) 
  print('do not merge meta for phyloseq:\n')
  phyloseq
}


##################################################
out_file2<-paste0(basename(phyloseq_path),'phyloseq_sample.meta',time_,'.csv') 
print(paste0('shuchulujing:',out_dir_path))
print(paste0('shuchulujing:',out_dir_path))


#¶ÁÈ¡Ä³¸öphyloseq
if (phyloseq_path !='FALSE') {
  phyloseq <- readRDS(phyloseq_path)
  sample_names(phyloseq)<-phyloseq_rename(sample_names(phyloseq))
  phyloseq
}
##################################################
# meta_filename_path='/home/zhiyu/data/metagenome/mmu_tumor_16S_analysis_dada2/id_path.csvmeta_tumor_data_new.csv.csv.csv'
# phyloseq<-FALSE_update

#¼ÓÈëmetaÐÅÏ¢
if (meta_filename_path !='FALSE') {
  #add meta info
  meta_info<-read.csv(meta_filename_path,stringsAsFactors =F,sep=',')
  message(meta_info$SampleID)
  print('ok')
  message(sample_names(phyloseq))
  #meta_info<-read.tsv(meta_filename_path,stringsAsFactors =F,sep='\t')
  meta_info$SampleID<-phyloseq_rename( meta_info$SampleID )
  Meta_tab<-column_to_rownames(meta_info,'SampleID')
  sample_data(phyloseq)<-sample_data(Meta_tab)

  # setdiff( meta_info$SampleID ,rownames(otu_table(all_phyloseq)) )
  # setdiff( rownames(otu_table(all_phyloseq)),meta_info$SampleID)

  # setdiff( rownames(otu_table(all_phyloseqa)),rownames(otu_table(all_phyloseq)) )
  phyloseq
  # phyloseq.to.df<- phyloseq_to_df(phyloseq)%>%data.frame()
  # write.csv(phyloseq.to.df, paste0(dirname(phyloseq_path),'/',paste0(basename(phyloseq_path),time_,"_phyloseq.to.df.csv" )) )
  out_file2<-paste0(basename(phyloseq_path),'phyloseq_sample.meta',time_,'.csv') 
  
  write.csv( sample_data(phyloseq)%>%data.frame() ,paste0(out_dir_path,'/',out_file2),quote = F) 

  
}

print('finnal phyloseq:')
phyloseq
out_file3<-paste0(basename(phyloseq_path),time_,"_phyloseq.names.csv" )
write.csv( rownames(otu_table(phyloseq)) %>% as.data.frame(), paste0(out_dir_path,'/',out_file3)) #save all_phyloseq

out_file4<-paste0(basename(phyloseq_path),"_update.rds" )
saveRDS(phyloseq,paste0(out_dir_path,'/',out_file4)) #save all_phyloseq


# #########nf-core phyloseq to biom
# library("biomformat")
# library("reltools")
# library(seqinr)
# OTU1 = t(as(otu_table(phyloseq), "matrix"))
# otu<-as(otu_table(OTU1,taxa_are_rows = FALSE),"matrix")
# otu_biom<-make_biom(data=otu)

# # dir.create(paste0(set_work_path,'/output/pircust2/otu_biom'),recursive = TRUE)
# write_biom(otu_biom, paste0(out_dir_path,"/otu_biom.biom"))
# #########nf-corephyloseq to fasta  
# #output fasta
# phyloseq_new.seq<-refseq(phyloseq)  %>% as.list()
# write.fasta(phyloseq_new.seq , names(phyloseq_new.seq),paste0(out_dir_path,"/uniqueSeqs.fasta"))
# # df_<-as.data.frame(otu_table(phyloseq_new))




# saveRDS(phyloseq,paste0(dirname(phyloseq_path),'/',paste0(basename(phyloseq_path),time_,"_update.rds" )))#save all_phyloseq





# write.csv(as.data.frame(rownames(otu_table(all_phyloseq))),paste0(phyloseq_workdir_path,'/merge_analysis/all_phyloseq/phyloseq_name.csv' ))
# if (isnot_input_meta){

# }
#############################biom rep.seq to qza  conda activate qiime2-2019-10
# qiime tools import \
#   --input-path otu_biom.biom \
#   --type 'FeatureTable[Frequency]' \
#   --input-format BIOMV100Format   --output-path feature-table.qza

# qiime tools import \
# --input-path uniqueSeqs.fasta \
# --output-path rep.seq.qza \
#  --type "FeatureData[Sequence]"  #--show-importable-formats DNAFASTAFormat
 #############################biom rep.seq to qza
