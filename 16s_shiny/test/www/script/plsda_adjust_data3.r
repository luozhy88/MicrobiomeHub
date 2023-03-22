library(microbiome)
library(phyloseq)
library(mixOmics)
library(metagMisc)
library(tidyverse)
library(dplyr)
library(DT)
library(metacoder)
library(purrr)
library(ggplot2)
library(microbiome)
library(ggplot2)
library(ggpubr)
library('grid')
library(ggsci)
library(png)
`%+%` <- function(a,b) {paste0(a,b)}
# library(ALDEx2)
# library("DESeq2")
color_<-c("black","red","green","blue","maroon","olivedrab2")
mycols <- c(
  "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
  "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
  "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
  "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
  "#000000"
)


aggregate_top_taxa<- function (x, top, level) {
  .Deprecated("aggregate_rare", "The microbiome::aggregate_top_taxa function is deprecated.")
  x <- aggregate_taxa(x, level)
  tops <- top_taxa(x, top)
  tax <- tax_table(x)
  inds <- which(!rownames(tax) %in% tops)
  tax[inds, level] <- "Other"
  tax_table(x) <- tax
  tt <- tax_table(x)[, level]
  tax_table(x) <- tax_table(tt)
  aggregate_taxa(x, level)
}

distribution_project_num<-function(list.meta.path){
    
    meta_num<-list.meta.path %>% purrr::map( function(x) (  read.csv(x,stringsAsFactors = FALSE)$SampleID %>% length()  )  ) %>% data.frame()
    # browser()
    all_meta_id_num <- data.frame( unlist(meta_num) )
    rownames(all_meta_id_num)<-gsub(".meta.csv","",rownames(all_meta_id_num) )
    all_meta_id_num$Project<-rownames(all_meta_id_num)
    
    print(all_meta_id_num)
    c1 <- ggplot(all_meta_id_num, aes(x=Project, y=unlist.meta_num.,fill=Project )) + 
      geom_bar(stat="identity", position=position_dodge() )+  
      scale_fill_manual(values=mycols)+
      #theme_bw()+
      theme_minimal()+
      theme( legend.position = "right",
             axis.title.x = element_blank(),
             axis.title.y = element_blank() ,
             axis.text.x =element_text(angle = 10,size = 20) ,
             axis.text.y =element_text(size = 20), 
             legend.text = element_text(size = 15) ,
             # legend.title = "Project",
             legend.title = element_text(size = 20) ,
             plot.title = element_text(size = 20, face="bold")
      ) +
      coord_flip() +
      ggtitle("The number of SampleID in Meta")
    
    c1
    
  }
  

adjust_data <- function(tax_level,phyloseq,target_2rd_column_){
  
# tax_level<-"Genus"
ps.tax_level<-phyloseq::tax_glom(phyloseq,tax_level)#??
ps.tax_level.df<-phyloseq_to_df(ps.tax_level)

# ps.tax_level.df$sum_<-apply(ps.tax_level.df[,9:length(colnames(ps.tax_level.df))],1,sum)
# ps.tax_level.df$filter<-sapply(ps.tax_level.df$sum,function(x){x/sum(ps.tax_level.df$sum)})
# ps.tax_level.df<-filter(ps.tax_level.df, filter >= 0.00005)

tax_info_<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species","OTU")
pre_tax_level<-tax_info_[which(tax_info_==tax_level)-1]#get the pre tax of tax_level
tax_info_<-tax_info_[ -which(tax_info_==tax_level) ]

ps.tax_level.df_newcolumns<-colnames(ps.tax_level.df)[!colnames(ps.tax_level.df)%in%c(tax_info_)]#filter the columns
 

##copy with the duplicate value
duplicate_loc<-which( duplicated( ps.tax_level.df[,tax_level] )  )#Find the location of the row which is duplicate value.

# ps.tax_level.df[tax_level,duplicate_loc]<-paste(ps.tax_level.df[duplicate_loc,pre_tax_level],ps.tax_level.df[duplicate_loc,tax_level],sep='_')
# ps.tax_level.df$Family[duplicated(ps.tax_level.df$Family)]<-paste(ps.tax_level.df[duplicate_loc,pre_tax_level],ps.tax_level.df[duplicate_loc,tax_level],sep='_')
ps.tax_level.df[duplicate_loc,tax_level]<-paste(ps.tax_level.df[duplicate_loc,pre_tax_level],ps.tax_level.df[duplicate_loc,tax_level],sep='_')#reset the value



ps.tax_level.df.sel<-ps.tax_level.df[,ps.tax_level.df_newcolumns] %>% column_to_rownames(tax_level) %>% t() %>% as.data.frame()#use the new df analysis PLS-DA

# filter ID
# filter_meta<-subset(phyloseq_meta,grepl(filter_str,Class))  #filter ID  KLN  vehicle  	Hepa
filter_meta<-meta(phyloseq)  #all id
filter_meta$SampleID<-rownames(filter_meta)

ps.tax_level.df.sel <- ps.tax_level.df.sel[filter_meta$SampleID, ]%>%as.data.frame()# %>% na.omit()
head(ps.tax_level.df.sel)
# setdiff(filter_meta$SampleID,rownames(ps.tax_level.df.sel))

# filter_meta$SampleID %in% rownames(ps.tax_level.df.sel)


TSS.divide = function(x){
  x/sum(x)}


ps.tax_level.df.sel<-ps.tax_level.df.sel+1
ps.tax_level.df.sel = apply(ps.tax_level.df.sel, 2, TSS.divide)

# Y<-filter_meta$type# zhuyixiiugai
Y<-filter_meta[,target_2rd_column_]# zhuyixiiugai

# Y<-filter_meta$Group
Y<-factor(trimws(Y))
# plsda_result <-plsda(ps.tax_level.df.sel, Y, ncomp = 3)

# ps.tax_level.df.sel[is.na(ps.tax_level.df.sel)] <- 0# fill na 1


# ps.tax_level.df.sel[ps.tax_level.df.sel==0] <- 1# fill 0 1

return(list(ps.tax_level.df.sel,Y))

}



phy_flitered <- function(ps.002.ordered, min_count=5,methods_trans="compositional",prevalence_values=.25){
  # browser()
  print("before:")
  print(ps.002.ordered)
  # min_count<-5
  # methods_trans<-"compositional" #  'Z', 'log10', 'log10p', 'hellinger', 'identity', 'clr'
  # prevalence_values <- .25
  
  filter <- phyloseq::genefilter_sample(ps.002.ordered, filterfun_sample(function(x) x >= min_count))
  phy_fil <- prune_taxa(filter, ps.002.ordered)
  phy_fil
  # methods_trans=NULL
  print(methods_trans)
  if ( !is.null(methods_trans) & !is.null(prevalence_values) ){
  phy_fil <- microbiome::transform(phy_fil, methods_trans)
  phy_fil <- core(phy_fil, detection = 0.0005, prevalence = prevalence_values )
  }
  print("post:")
  print(phy_fil)
  
  return(phy_fil)
}





phy_flitered0 <- function(ps.002.ordered, percent=0.005){
  print("before:")
  print(ps.002.ordered)
  Tax <- tax_table(ps.002.ordered)
  
  # extraction of the metadata
  Metadata <- ps.002.ordered@sam_data
  
  # extract OTU table from phyloseq object
  data.raw <- otu_table(ps.002.ordered) # samples should be in row and variables in column
  
  # offset
  data.offset <- data.raw+1
  
  # filter using the function below
  # function to perform pre-filtering
  low.count.removal = function(
    #data<-data.offset,
    data, # OTU count data frame of size n (sample) x p (OTU)
    percent=0.005 # cutoff chosen
  ){
    print( paste0("low.count.removal :", percent) )
    keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
    data.filter = data[,keep.otu]
    return(list(data.filter = data.filter, keep.otu = keep.otu))
  }
  ##############################################################
  # call the function then apply on the offset data
  result.filter <- low.count.removal(data.offset, percent)
  data.filter <- result.filter$data.filter
  # check the number of variables kept after filtering
  length(result.filter$keep.otu) 
  
  # checking the lib size per sample
  lib.size <- apply(data.filter, 1, sum)
  #### change the figure size according to your sample size
  png(filename = "library.size.png", width = 20, height = 8, units = "in", res = 300)
  barplot(lib.size, las = 2)  # check please to make sure there is no outliers (either too few reads or too many reads)
  dev.off()
  
  #########################################################################
  #########################################################################
  #########################################################################
  #Filter the taxonomy corresponding to the selected OTU
  ############################# you will need this for outputs and further analyses
  Tax.filter <- Tax[unlist(result.filter$keep.otu),]
  #########################################################################
  #########################################################################
  #########################################################################
  
  # check that the OTUs are similar between files
  summary(colnames(data.filter)==rownames(Tax.filter))
  
  
  ps <- phyloseq(tax_table(Tax.filter), otu_table(data.filter, taxa_are_rows = FALSE), sample_data(Metadata))

  print("post:")
  print(ps)
  return(ps)
}


convert_seq_to_asv<-function(phy){
  ## rename sequence to ASV
  print('###########convert_seq_to_asv############')
  library("phyloseq")
  # phy <- phy %>% phyloseq::t()
  print(phy)
  #otu_table(phy) <- otu_table(t( otu_table(phy) %>% data.frame() ),taxa_are_rows = T)
  # print("phyloseq::t ok")
  if(all(rownames(otu_table(phy)) ==rownames(tax_table(phy)))){
    out <- cbind(paste0("ASV",seq(1:nrow(otu_table(phy)))),rownames(tax_table(phy)),tax_table(phy))
    
    write.table(out,file=paste0("out.result/","seq_asv_table.txt"),row.names = F,col.names = F,quote = F)
    taxa_names(phy) <- paste0("ASV",seq(1:nrow(otu_table(phy))))
  }else{
    message("rownames in otu table is not equal to rownames in tax table")
  }
  print('###########END convert_seq_to_asv############')
  return(phy)
}

ordered_phy<-function(phy,target_2rd_column){
  #phy<-phyloseq
  #target_2rd_column<-"group1"
  phy.meta<-meta(phy)
  phy.meta$gp <-phy.meta[,target_2rd_column]
  sample_data(phy)<-phy.meta
  diff_group<-unique(phy.meta$gp)
  linshi<-list()
  
  for (i in diff_group){
    print(i[1])
    print(class(i))
    meta_sub<-subset(phy.meta,grepl(i,gp))
    phy_2<-phy
    sample_data(phy_2)<-meta_sub
  
    #a<-subset_samples(phy, gp == "pre_Tg"  )
    print(phy_2)
    linshi[i]<-phy_2
  }
  #order_phy<-merge_phyloseq(unlist(linshi))
  order_phy <-  do.call(merge_phyloseq,linshi)
  return(order_phy)
}

readpng<-function(filename){

  if (file.exists(filename)){
                uu<-readPNG(filename)
                print(filename)
                grid::grid.raster(uu)
  }else{
        print(paste0("no exist figure:",filename))
        return(NULL)                    
                            }
}

readtable<-function(filename){

  if (file.exists(filename)){
                ff<-read_tsv(filename) %>% data.frame()
                print(filename)
                return(ff)
  }else{
        print(paste0("no exist table:",filename))
        return(NULL)                    
                            }
}


color.gvri <- function (num.vector) 
{
  if (is.factor(num.vector)) 
    num.vector = as.numeric(num.vector)
  if (!is.numeric(num.vector)) {
    stop(paste("num.vector has to be numeric", call. = FALSE))
  }
  ##### lancet color palette
  mixo.col = c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF")
  n = length(num.vector)
  if (isTRUE(num.vector) > length(mixo.col)) {
    stop(paste("We only have a few mix.colors available, n <= ", 
               length(mixo.col)), call. = FALSE)
  }
  if (isTRUE(!is.finite((num.vector))) || (n < 1)) {
    stop("'num.vector' must be an integer vector with positive values.", 
         call. = FALSE)
  }
  return(mixo.col[num.vector])
}

spda_rere_plot<-function(dframe,Tax.filter,n,i){  
  # yy<<-dframe
  # uu<<-Tax.filter
  # dframe<-yy
  # Tax.filter<-uu
  # i<-"Genus"
  # n='n1'
  #dframe=res.splsda.plot_1

  # browser()
  # if(i=="OTU"){Tax.filter$OTU=rownames(Tax.filter)}
  # all_asv<-Tax.filter[, i]
  # all_asv$SampleID<-rownames(all_asv)
  # all_asv$SampleID2<-paste(all_asv[,i],"|",rownames(all_asv))
  # dframe=merge(dframe,all_asv,by="SampleID",all.x = T)
  Tax.filter$SampleID<-rownames(Tax.filter)
  dframe=merge(dframe,Tax.filter,by="SampleID",all.x = T)
  dframe=dframe[order(abs(dframe$importance),decreasing=T),]
  
  
  # browser()
  # dframe %>% select_if(function(importance){abs(importance)>=0.01})
  dframe$importance2<-dframe$importance %>% abs()
  dframe<-dframe %>% dplyr::filter(importance2>=0.01)
  tt<-as.numeric(as.character(dframe$importance) )
  
  
  col <- as.character(dframe$color)
  names(col) <- as.character(dframe$GroupContrib)
  
  
  p <- ggplot(data=dframe,aes(x= reorder(dframe$SampleID, abs(tt) ) ,y= tt,fill = GroupContrib) )+ 
    #scale_y_discrete( limits = c( -1,1 ))+
    scale_fill_manual(values=col) +
    geom_bar(stat="identity",width = 0.5,position=position_dodge())  + 
    coord_flip() +theme_minimal() + 
    #scale_y_discrete( limits = c( -1,1 ))+
    labs(y="importance",x="") + 
    #theme(axis.text.y = element_text(angle = 0, size = 10, face = "bold"))
    theme(strip.text.x = element_text(size = 20, face = "bold"),
	    text = element_text(size = 20),
            axis.text.x = element_text(angle = 0,size=15),
            axis.text.y =element_text(size = 10), 
            axis.title.y = element_blank() ,
            legend.text = element_text(size = 20) ,
            legend.title = element_text(size = 23) ,
	    plot.title = element_text(size = 30, face="bold"))+
    ggtitle( paste0(" comp ",n," ",i) )
   
  ggsave(plot= p, filename = paste0("out.result/splsda/splda_test",".png"), dpi = 300, height = 10, width = 10, units = "in")

  p
  
  return(p)
}


alpha.diversity<-function(tax_level,phyloseq,target_2rd_column){
  
  #filter.meta() 
  #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
  #phyloseq<-readRDS( phyloseq_path   )
  print("##############alpha.diversity#################")
  #phyloseq<-phy.rare
    
    
  meiji<-phyloseq
  #target_2rd_column<-"disease"
  #target_2rd_column<-input$SelectP3_1
  phyloseq::sample_data(meiji)$Group <- meta(meiji)[,target_2rd_column]
  input_levels<-unique(sample_data(meiji)$Group ) %>% as.vector()
  sample_data(meiji)$Group <- factor(sample_data(meiji)$Group, levels = input_levels  )
  cc <- meta(meiji)
  
  #sample_data(meiji)$Diagnosis <- factor(sample_data(meiji)$Diagnosis, levels = c("NC", "MCI", "AD"))
  #color_<-c("black","red","green","blue","maroon","olivedrab2")
  
  
  comb <- combn(input_levels,2)
  #comb <- combn(strsplit("p008_Hepa_Vehicle p023_Hepa_vehicle p027_SP20_vehicle p035_SP20_Vehicle"," ")[[1]],2)
  my_compare <- lapply(1:dim(comb)[2] ,function(i){ return(comb[,i])})
  #my_compare <- list(c("post_Tg","post_Tg_GV971"),c("pre_Tg","post_Tg"),c("pre_Tg","pre_Tg"))#002
  # my_compare <- list(  c( "post_Tg", "post_Tg_GV971_50mpk"),
  #                      c( "post_Tg", "post_Tg_GV971_200mpk"),
  #                      c("pre_Tg","post_Tg"),
  #                      c("pre_Tg","pre_Tg_GV971")    )#004
  
  #'Phylum', 'Order', 'Class', 'Family', 'Genus',    
    if (tax_level=='OTU'){
      print('use OTU')
      pseq <- meiji 
      #%>% microbiome::transform(transform = "compositional")
    }else{
      print('use tax_level')
        pseq <- meiji %>%
        aggregate_taxa(level = tax_level) %>%
        #microbiome::transform(transform = "compositional") %>%
        aggregate_top_taxa(top = 10, tax_level)
    }
    # browser()

    input_levels_<-factor(input_levels)
    names(input_levels_)<-color.gvri(input_levels_)
    
    input_levels_2<-as.vector( names(input_levels_)  )
    names(input_levels_2)<-as.character(input_levels_)
    
  
    alpha <- plot_richness(meiji, x = 'Group',color ='Group', measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "Fisher"),nrow=2 ) + 
      geom_boxplot() +
      geom_jitter() +
      theme_bw() +
      xlab("") +
      #scale_colour_manual(   values = c("post_Tg" = "#E08214", "post_Tg_GV971" = "#E08214","pre_Tg" = "#E08214","pre_Tg_GV971" = "#0099B4FF" )  ) +
      scale_colour_manual(  values =  input_levels_2 ) +
      stat_compare_means(comparisons = my_compare,method = "wilcox.test", label = "p.format") + 
      #stat_compare_means(comparisons = list(input_levels), method = "wilcox.test", label = "p.format")+
      ggplot2::guides(fill = guide_legend(ncol = 1)) +
      ggplot2::guides( fill = guide_legend(tax_level) ) +
      theme(strip.text.x = element_text(size = 20, face = "bold"),
	    text = element_text(size = 20),
            axis.text.x = element_text(angle = 0,size=18),
            axis.text.y =element_text(size = 18), 
            axis.title.y = element_blank() ,
            legend.text = element_text(size = 25) ,
            legend.title = element_text(size = 30) ,
	    plot.title = element_text(size = 30, face="bold"))+
      #geom_text( aes(label= meta(pseq)$SampleID ) )+
      # scale_fill_manual("group")+
      ggtitle("Alpha Diversity Measure")
      #labs(legend = tax_level )
    alpha
    print("ciweizhi alpha:")
    print(getwd())
    dir.create(paste0(getwd(),"/out.result/alpha"))
    ggsave(alpha, filename =paste0(getwd(),"/out.result/alpha/","alpha_",tax_level,"_Shannon_",".png"), dpi = 300, height = 18, width = 20, units = "in")

    
    print("##############END alpha.diversity#################")
    alpha
    
}




beta.diversity<-function(phyloseq,target_2rd_column,pca_method="PCoA"){
  #filter.meta() 
  #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
  #phyloseq<-readRDS( phyloseq_path   )
  print("##############beta.diversity#################")
  meiji<-phyloseq
  #target_2rd_column<-"disease"
  #target_2rd_column<-input$SelectP3_1
  
  library(microbiome)
  library(phyloseq)
  library(ggplot2)
  ##########################
  pseq <- phyloseq
  
  ######################## offset by 1
  
  # pseq.offset <- microbiome::transform(pseq, 'shift', shift=1)
  
  # Convert to compositional data for filtering out the low abundant and low prevalent OTUs
  # pseq.rel <- microbiome::transform(pseq.offset, "compositional")
  
  # Pick core taxa with with the given prevalence and detection limits
  # pseq.core <- core(pseq.rel, detection = .01/100, prevalence = 20/100)
  
  # Use clr transform for the core
  # pseq.core <- microbiome::transform(pseq.core, "clr")
  ####################
  # Ordinate the data
  set.seed(100)
  
  # pca_method <- "PCoA"
  ord <- ordinate(pseq, pca_method , "bray") 
  
  hh<-pseq@sam_data %>% data.frame()
  pseq@sam_data$Group<-hh[,target_2rd_column]
  # browser()
  
  phyloseq::sample_data(pseq)$Group <- meta(pseq)[,target_2rd_column]
  input_levels<-unique(sample_data(pseq)$Group ) %>% as.vector()
  
  input_levels_<-factor(input_levels)
  names(input_levels_)<-color.gvri(input_levels_)
  
  input_levels_2<-as.vector( names(input_levels_)  )
  names(input_levels_2)<-as.character(input_levels_)
  
  
  beta_diversity <- plot_ordination(pseq, ord, color = "Group") +
    geom_point(size = 5) +
    theme_bw() +
    stat_ellipse()+
    theme(strip.text.x = element_text(size = 20, face = "bold"),
	    text = element_text(size = 20),
	    axis.text.x = element_text(angle = 0,size=30),
	    axis.text.y =element_text(size = 30), 
	    axis.title.y = element_blank() ,
	    legend.text = element_text(size = 30) ,
	    legend.title = element_text(size = 35) ,
	    plot.title = element_text(size = 30, face="bold"))+
    scale_colour_manual(  values =  input_levels_2 ) +

    ggtitle("Beta Diversity Measure")
  
  dir.create( paste0(getwd(),"/out.result/beta") ,recursive = T)
  ggsave(beta_diversity, filename =paste0("out.result/beta/","beta_diversity_",".png"), dpi = 300, height = 18, width = 20, units = "in")
  ggsave(beta_diversity, filename =paste0("out.result/beta/","beta_diversity_",".svg"), dpi = 300, height = 18, width = 20, units = "in",device = "svg")
  beta_diversity
  
}


# plot_compositional ------------------------------------------------------
plot_compositional<-function(tax_level,phyloseq,target_2rd_column){
  
  #filter.meta() 
  #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
  #phyloseq<-readRDS( phyloseq_path   )
  print("##############plot_compositional#################")
  meiji<-phyloseq
  #target_2rd_column<-"disease"
  #target_2rd_column<-input$SelectP3_1
  
  library(tidyverse)
  library(phyloseq)
  library(microbiome)
  library(metagMisc)
  
  ### define the color
  mycols <- c(
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
    "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
    "#000000"
  )
  
  phyloseq::sample_data(meiji)$Group <- meta(meiji)[,target_2rd_column]
  sample_data(meiji)$Group <- factor(sample_data(meiji)$Group, levels = unique(sample_data(meiji)$Group ) %>% as.vector()  )
  
  #for (i in c('Phylum', 'Order', 'Class', 'Family', 'Genus')) {
    #tax_level=i
    print(tax_level)
    pseq <- meiji %>%
      aggregate_taxa(level = tax_level) %>%
      microbiome::transform(transform = "compositional") %>%
      aggregate_top_taxa(top = 10, tax_level)
    # browser()

    p.bar.plot_composition <- plot_composition(pseq, average_by = 'Group') +
      scale_fill_manual(values = mycols) +
      guides(fill = guide_legend(ncol = 1)) +
      ylab("Relative Abundance (%)") +
      xlab("") +
      guides(fill = guide_legend(tax_level)) +
      theme_minimal() +
      theme(strip.text.x = element_text(size = 14, face = "bold")) +
      theme(text = element_text(size = 14)) +
      # theme(axis.text.x = element_text(angle = 90, size = 8, face = "bold"))+
    theme(strip.text.x = element_text(size = 20, face = "bold"),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 0,size=30),
          axis.text.y =element_text(size = 30), 
          axis.title.y = element_blank() ,
          legend.text = element_text(size = 30) ,
          legend.title = element_text(size = 35) ,
          plot.title = element_text(size = 30, face="bold"))+
      #scale_x_discrete(labels = abbreviate)+
      scale_x_discrete(labels = label_wrap(10))+
      ggtitle("Compositional analysis")
    
      dir.create( paste0(getwd(),"/out.result/composition") ,recursive = T)
      ggsave(p.bar.plot_composition, filename =paste0("out.result/composition/","composition_.",tax_level,".png"), height = 18, width = 12, units = "in")
  #}
  print("##############END plot_compositional#################")
  p.bar.plot_composition
}

AD_plot_pca_bar<-function(tax_level,phyloseq,target_2rd_column){
  # df_return_g<- adjust_data(tax_level,phyloseq,target_2rd_column)
  # df_new_g<-df_return_g[[1]]
  # Y_g<-df_return_g[2][[1]]
  # tax_level<-"OTU"
  #target_2rd_column<-"disease"
  print("##############PCA#############")
  print("PCA only OTU level ??")
  
  if (tax_level=='OTU'){
    print('use OTU')
    pseq <- phyloseq 
    #%>% microbiome::transform(transform = "compositional")
  }else{
    print('use tax_level')
    pseq <- phyloseq %>%
      aggregate_taxa(level = tax_level) #%>%
    #microbiome::transform(transform = "compositional") %>%
    # aggregate_top_taxa(top = 10, tax_level)
  }
  pseq
  
  

  with_=11
  height_=8
  #paste0("out.result/","pca_",tax_level,".png" )
  # browser()
  dir.create( paste0(getwd(),"/out.result/pca") ,recursive = T)
  png(filename = paste0("out.result/pca/","pca_bar_",tax_level,".png" ), width = with_, height = height_ ,units = "in", res = 300) 
  pca_tab <- t(otu_table(pseq)+1)
  # pca_bar = mixOmics::pca(pca_tab, ncomp = 10, center = T, scale = F,logratio = "CLR")
  pca_bar = mixOmics::pca(pca_tab, ncomp = 10, center = T, scale = F,logratio = "none")
  # par(family="inconsolata")
  plot(pca_bar)
  dev.off()
  
  plot(pca_bar)
}





pca_<-function(tax_level,phyloseq,target_2rd_column,n_componet){
  # df_return_g<- adjust_data(tax_level,phyloseq,target_2rd_column)
  # df_new_g<-df_return_g[[1]]
  # Y_g<-df_return_g[2][[1]]
  # browser()
  # tax_level<-"OTU"
  #target_2rd_column<-"host_disease"
  print("##############PCA#############")
  print("PCA only OTU level ??")
  
  if (tax_level=='OTU'){
    print('use OTU')
    pseq <- phyloseq 
    #%>% microbiome::transform(transform = "compositional")
  }else{
    print('use tax_level')
    pseq <- phyloseq %>%
      aggregate_taxa(level = tax_level) #%>%
    #microbiome::transform(transform = "compositional") %>%
    # aggregate_top_taxa(top = 10, tax_level)
  }
  pseq
  
  
  data.filter<-otu_table(pseq)
  Tax.filter<-tax_table(pseq)
  Metadata<-meta(pseq)
  
  Metadata[,target_2rd_column]<-factor(Metadata[,target_2rd_column])
  Y_g<-Metadata[,target_2rd_column]
  df_new_g<-as.data.frame(data.filter) %>% t()
  
  # pca_g<-pca(df_new_g, ncomp = ifelse(length(colnames(df_new_g))>10,10,length(colnames(df_new_g))), logratio = 'CLR')##pca
  pca_g<-pca(df_new_g, ncomp = ifelse(length(colnames(df_new_g))>10,10,length(colnames(df_new_g))), logratio = 'none')##pca
  plot(pca_g,main='PCA') #main=tax_level
  

  dir.create( paste0(getwd(),"/out.result/pca") ,recursive = T)
  with_=6
  height_=6
  #paste0("out.result/","pca_",tax_level,".png" )
  # browser()
  # n_componet=2
  
  svg(filename = paste0("out.result/pca/","pca_",tax_level,".svg" ), width = 5.5, height = height_ )  
  mixOmics::plotIndiv(pca_g,  comp = c(1,n_componet),
                      ellipse = F,point.lwd = 2,abline = T,
                      ind.names = F, group = Y_g, col.per.group = ggsci::pal_lancet("lanonc")(9)[1:nlevels(Y_g)], 
                      legend = TRUE,legend.position="bottom",title = paste0('PCA 1-',n_componet))#title = tax_level
  dev.off()
  

  png(filename = paste0("out.result/pca/","pca_",tax_level,".png" ), width = 6, height =6, units = "in", res = 300)
  mixOmics::plotIndiv(pca_g,  comp = c(1,n_componet),
                      ellipse = F,point.lwd = 2,abline = T,
                      ind.names = F, group = Y_g, col.per.group = ggsci::pal_lancet("lanonc")(9)[1:nlevels(Y_g)],
                      legend = TRUE,legend.position="bottom",title = paste0('PCA 1-',n_componet))#title = tax_level
  dev.off()
  
  
}














# pca_("Genus",'phyloseq','diease')

splsda.summary <- function(input_pls_obj,outname){
  # plotloading
  
  loadings_result <- c()
  
  for (i in 1:input_pls_obj$ncomp){
    png(filename = paste0(outname,"_",i,".png"),dpi = 300, width = 500,height = 500)
    tab <- plotLoadings(input_pls_obj, comp = i, title = paste('Loadings on comp',i), 
                        contrib = 'max', method = 'mean',
                        legend.col = ggsci::pal_lancet(palette = c("lanonc"), alpha = 1)(2),size.name=0.6)
    tab
    dev.off()
    tab$ncomp <- i
    loadings_result <- rbind(loadings_result,tab)
  }
  return(loadings_result)
}




splsda.replot <- function(input_file,input_levels,input_col){
  # this script will take splsda result file which should be modified manuelly for visualisation purpose then plot a lit bit more beautiful barplot
  
  require(ggplot2)
  require(ggpubr)

  tab <- read.csv(input_file)
  ncomp <- max(tab$ncomp)
  splda.plot.list=list()

  for (i in 1:ncomp){
    comp <- tab[tab$ncomp == i,]
  tabN1 <- apply(comp,1, function(x){

      if(x[7] == input_levels[1]){
        import <- abs(as.numeric(x[9]))
        col <- input_col[1]
      }else{
        import <- -abs(as.numeric(x[9]))
        col <- input_col[2]
      }
      out <- c(x[1],import,col,x[7])
      return(out)
    })
  tabN <- data.frame(t(tabN1),stringsAsFactors = F)
  colnames(tabN) <- c("id","importance","col","Class")
  tabN <- tabN[order(abs(as.numeric(as.character(tabN$importance)))),]
  p <- ggplot(data=tabN,aes(x=tabN$id,y=as.numeric(as.character(tabN$importance)),fill = Class,co=col))+ 
    #scale_y_discrete( limits = c( -1,1 ))+
    geom_bar(stat="identity",width = 0.5,position=position_dodge())  + 
    coord_flip() +theme_minimal() + 
    #scale_y_discrete( limits = c( -1,1 ))+
    labs(y="importance",x="lowest level taxon") + 
    ggtitle(paste0(" comp ",i))+
    theme(axis.text.y = element_text(angle = 0, size = 10, face = "bold"))
  
  ggsave(plot= p, filename = paste0(gsub(".csv","",input_file),"_comp_",i,".png"), dpi = 300, height = 10, width = 10, units = "in")
  print('run jieshu')
  #splda.plot.list[[paste0('t',i)]]=p
  #
  }
  #return(splda.plot.list)
}


pls.da_ <- function(tax_level,phyloseq,target_2rd_column){
	# df_return_g<- adjust_data(tax_level,phyloseq,target_2rd_column)
	# df_new_g<-df_return_g[[1]]
	# Y_g<-df_return_g[2][[1]]
  print('##############PLSDA#############') 
  print(target_2rd_column)
  # browser()
  data.filter<-otu_table(phyloseq)
  Tax.filter<-tax_table(phyloseq)
  Metadata<-meta(phyloseq)

  Metadata[,target_2rd_column]<-factor(Metadata[,target_2rd_column])
  Y_g<-Metadata[,target_2rd_column]
  # df_new_g<-as.data.frame(data.filter)
  df_new_g<-as.data.frame(data.filter) %>% t()



  # gvri.plsda = plsda(X = df_new_g, Y = Y_g, ncomp = nlevels(Y_g)+1, logratio = "CLR")
  gvri.plsda = plsda(X = df_new_g, Y = Y_g, ncomp = nlevels(Y_g)+1, logratio = "none")

  gvri.perf.plsda = perf(gvri.plsda, validation = 'Mfold', folds = 5,
                         progressBar = FALSE, nrepeat = 50)

  #plot(gvri.perf.plsda, overlay = 'measure', sd=TRUE)

  with_=9.5
  height_=9
  #First two components
  #png(filename = paste0("out.result/","pls.da_1-2",tax_level,".png")
  # browser()
  dir.create( paste0(getwd(),"/out.result/plsda") ,recursive = T)
  png(filename = paste0("out.result/plsda/","pls.da_1-2.png"), width = with_, height = height_ ,units = "in", res = 300)  
  #png(filename = "PLSDA_12.png", width = with_, height = 6.5, units = "in", res = 300)
  plotIndiv(gvri.plsda , comp = c(1,2),
            size.xlabel=1.5,size.ylabel=1.5,size.axis=1,size.legend=1.2,
            group = Y_g, ind.names = F, point.lwd = 2, style = "graphics", abline = T,
            ellipse = TRUE, legend = TRUE, title = 'PLSDA comp 1 - 2', col.per.group = ggsci::pal_lancet("lanonc")(9)[1:nlevels(as.factor(Y_g))])
  dev.off()



  #First and the Third components
  png(filename = paste0("out.result/plsda/","pls.da_1-3.png"), width = with_, height = height_ ,units = "in", res = 300)  
  #png(filename = "PLSDA_13.png", width = with_, height = 6.5, units = "in", res = 300)
  plotIndiv(gvri.plsda , comp = c(1,3),
            size.xlabel=1.5,size.ylabel=1.5,size.axis=1,size.legend=1.2,
            group = Y_g,  ind.names = F, point.lwd = 2, style = "graphics", abline = T,
            ellipse = TRUE, legend = TRUE, title = 'PLSDA comp 1 - 3', col.per.group = ggsci::pal_lancet("lanonc")(9)[1:nlevels(as.factor(Y_g))])
  dev.off()

	}
#pls.da_(tax_level,phyloseq,target_2rd_column)





spls.da_new_step1 <- function(tax_level,phyloseq,target_2rd_column,n_componet){
  #tax_level<-"Genus"
  #target_2rd_column<-"Class"
  
  # df_return_g<- adjust_data(tax_level,phyloseq,target_2rd_column)
  # df_new_g<-df_return_g[[1]]
  # Y_g<-df_return_g[2][[1]]
  #phyloseq<-file.phyloseq
  # browser()
  print("##############sPLSDA#############")
  #phyloseq<-convert_seq_to_asv(phyloseq)
  # phy_tmp <- tax_glom(phyloseq,taxrank = "Genus") 
  # tax_level<-"Genus"
  if(tax_level =="OTU"){ 
    print("OTU")
    phy_tmp = phyloseq
    # phy_tmp@tax_table
    # tax_table(phy_tmp)
    }else{
      print(tax_level)
      phy_tmp = phyloseq %>% aggregate_taxa(level = tax_level)
      }
  
  # ifelse(tax_level =="OTU", phy_tmp <- phyloseq , phy_tmp <- phyloseq %>% aggregate_taxa(level = tax_level) )
  # phy_tmp <- phyloseq %>% aggregate_taxa(level = tax_level)
  
  data.filter<-otu_table(phy_tmp)
  Tax.filter<-tax_table(phy_tmp) %>% data.frame()
  Metadata<-meta(phyloseq)
  
  Metadata[,target_2rd_column]<-factor(Metadata[,target_2rd_column])
  Y_g<-Metadata[,target_2rd_column]
  df_new_g<-as.data.frame(data.filter) 
  
  X = as.data.frame(data.filter)%>%t()
  Y = Y_g
  
  trans_method<-"none"#CLR
  with_=12
  height_=6.5
  # first, set a grid of values to test:
  grid.keepX = c(seq(5,150, 5))  
  # if you dont understand what this means, type
  # grid.keepX  # adjust this grid as necessary for your own data
  
  # this chunk takes ~2 min to run
  set.seed(100)  # for reproducible results for this code, remove for your own code
  gvri.tune.splsda = tune.splsda(X = X,
                                 Y = Y,
                                 ncomp = 3,
                                 logratio = trans_method,
                                 test.keepX = grid.keepX,
                                 validation = c('Mfold'),
                                 folds = 5,
                                 dist = 'max.dist', # prediction distance can be chosen according to tune.plsda results
                                 nrepeat = 50,
                                 progressBar = T,cpus=cpus_num)
  # may show some convergence issues for some of the cases, it is ok for tuning
  
  #plot(gvri.tune.splsda)
  
  # optimal number of variables to select on 2 comps:
  select.keepX = gvri.tune.splsda$choice.keepX[1:3]
  select.keepX
  # the sPLS-DA
  gvri.splsda = splsda(X = X,  Y = Y, ncomp = 3, keepX = select.keepX, logratio= trans_method )
  #png(filename = paste0("out.result/","pls.da_1-3",tax_level,".png"), width = with_, height = height_ ,units = "in", res = 300)  
  # browser()
  n_componet
  
  dir.create( paste0(getwd(),"/out.result/splsda") ,recursive = T)
  png(filename = paste0("out.result/splsda/","SPLSDA_1_",n_componet,".png.png"), width = 9, height = 9, units = "in", res = 300)
  plotIndiv(gvri.splsda,
            comp = c(1,n_componet),
            ind.names = F,
            abline = T,
            style = "graphics",
            point.lwd = 2,
            ellipse = TRUE,
            legend = TRUE,
            size.xlabel=1.5,size.ylabel=1.5,size.axis=1.5,size.legend=1.2,
            col.per.group = ggsci::pal_lancet("lanonc")(9)[1:nlevels(as.factor(Y_g ))],
            title = paste0('sPLS-DA comp 1-' , n_componet)  )
  dev.off()
  
  
  set.seed(200)  # for reproducible results for this code, remove for your own code
  
  # gvri.perf.splsda = perf(gvri.splsda, validation = 'Mfold', folds = 5,progressBar = FALSE, nrepeat = 50, dist = 'max.dist')
  # gvri.perf.splsda$error.rate
  # plot(gvri.perf.splsda)
  
  #############################
  
  with_=16
  #tax_group=c('Phylum', 'Order', 'Class', 'Family', 'Genus')
  # tax_group=c(tax_level)
  # for (i in tax_group) {
    # browser()
    # i=tax_level
    ############################
    if(tax_level=="OTU"){Tax.filter$OTU <- rownames(Tax.filter) }
    name.var = Tax.filter[, tax_level]
    # s = paste(Tax.filter[, tax_level], colnames(X), sep = '|')
    ####################comp1#########################
    #paste0("out.result/",i,"_sPLSDA_comp1.png",tax_level,".png")
    n_componet_=n_componet-1
    png(filename = paste0("out.result/splsda/",tax_level,"_sPLSDA_comp",n_componet,".png") , width = with_, height = 6, units = "in", res = 300)
    res.splsda.plot<- plotLoadings(gvri.splsda, comp = n_componet_, method = 'mean', contrib = 'max',
                                     legend.color = ggsci::pal_lancet("lanonc")(9)[1:nlevels(as.factor(Y_g))],
                                     # name.var = name.var,  #comment ?
                                     size.title = 1.0, ndisplay = 50, size.name = .6, size.legend = 1.2, margins=c(120,100)   )
                                     #+ theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))  margins=c(10,5),plot.margin = unit(c(0.8, 0.8, 0.8, 0.5))
    dev.off()
    
    res.splsda.plot$SampleID = rownames(res.splsda.plot)
    #
    p<-NULL
    p<-spda_rere_plot(res.splsda.plot,Tax.filter,paste0('n',n_componet),tax_level)
    p
    ggsave(plot= p, filename = paste0("out.result/splsda/splda_",tax_level,"_plot_n",n_componet,".png"), dpi = 300, height = 10, width = 10, units = "in")
    #
    write.table(res.splsda.plot, file = paste0("out.result/splsda/splda_",tax_level,"_plot_n",n_componet,".xls"),sep = "\t", row.names = F,quote = F)
  # }

  print("##############END sPLSDA#############")
  # return_list<-list(f1=gvri.splsda,f2=Y)
  # return(return_list)
  
  
  print("##############Heatmap#############")
  with_=8
  n_componet
  ##################### Heatmap on outputs of sPLS-DA
  png(filename = paste0("out.result/splsda/","Heatmap_comp1_comp",n_componet,".png"), width = with_, height = 6.5, units = "in", res = 300)
  cim(gvri.splsda, 
      row.sideColors = color.gvri(Y), 
      cluster = "column", 
      comp = c(1,n_componet), ## please make sure you choose the right number of components, or keep all components by unticking this line
      margins = c(5,6)
      #color = viridis::viridis(1000, option = "D") ## change the palette
  )
  dev.off()
  
  
  png(filename = paste0("out.result/splsda/","Heatmap_comp_all.png"), width =with_, height = 6.5, units = "in", res = 300)
  cim(gvri.splsda, 
      row.sideColors = color.gvri(Y), 
      cluster = "column", 
      # comp = c(1,2), ## please make sure you choose the right number of components, or keep all components by unticking this line
      margins = c(5,6)
      #color = viridis::viridis(1000, option = "D") ## change the palette
  )
  dev.off()
  
  print("##############END Heatmap#############")
  
  
  ####### redefine a color palette to make it consistent with the above color palette
} 



# 
# spls.da_new_step2 <- function(gvri.splsda,Y,n_componet){
#   print("##############Heatmap#############")
#   with_=8
#   n_componet
#   ##################### Heatmap on outputs of sPLS-DA
#   png(filename = paste0("out.result/splda/","Heatmap_comp1_comp",n_componet,".png"), width = with_, height = 6.5, units = "in", res = 300)
#   cim(gvri.splsda, 
#       row.sideColors = color.gvri(Y), 
#       cluster = "column", 
#       comp = c(1,n_componet), ## please make sure you choose the right number of components, or keep all components by unticking this line
#       margins = c(5,6)
#       #color = viridis::viridis(1000, option = "D") ## change the palette
#   )
#   dev.off()
#   
#   
#   png(filename = paste0("out.result/splda/","Heatmap_comp_all.png"), width =with_, height = 6.5, units = "in", res = 300)
#   cim(gvri.splsda, 
#       row.sideColors = color.gvri(Y), 
#       cluster = "column", 
#       # comp = c(1,2), ## please make sure you choose the right number of components, or keep all components by unticking this line
#       margins = c(5,6)
#       #color = viridis::viridis(1000, option = "D") ## change the palette
#   )
#   dev.off()
#   
#   print("##############END Heatmap#############")
# }




# Deseq2 ------------------------------------------------------------------

plot_deseq2 <- function(tax_level,phyloseq,target_2rd_column,deseq2_foldchange,deseq2_pvalue,deseq2_padj,Deseq2_height,Deseq2_width){
  # tax_level<-"Genus"
  #tax_level<-"OTU"
  #target_2rd_column<-"Class2"

  #filter.meta()
  #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
  # phyloseq<-readRDS( phyloseq_path   )
  print("##############plot_deseq2#################")
  meiji<-phyloseq
  #target_2rd_column<-"disease"
  #target_2rd_column<-input$SelectP3_1
  library(phyloseq)
  library(tidyverse)
  library(microbiome)
  library(metagMisc)
  library(speedyseq)
  library(ggpubr)




  # browser()
  ###### now you may go for OTU or other levels as input to Deseq2
  ###### most of time, you may need to filter out those very low abundant and prevalence OTUs before going futher.
  ps.sel.core <- phyloseq

  ##  ps.sel.core <- microbiome::core(ps.sel.AD.HC, 0, 30/100)

  ### you may aggregate on specific taxa level such as Genus, Family, etc.
  phyloseq::sample_data(ps.sel.core)$Group <- meta(ps.sel.core)[,target_2rd_column]

  ref_1 <- unique(meta(ps.sel.core)$Group)[1]%>%as.character()
  ref_2 <- unique(meta(ps.sel.core)$Group)[2]%>%as.character()
  ref_=ifelse(length(ref_2)>length(ref_1),ref_2,ref_1)
  print(ref_)
  if (tax_level=='OTU'){
    print('use OTU')
    ps.sel.core=ps.sel.core
  }else{
    print(paste0("use ",tax_level))
    ps.sel.core <- speedyseq::tax_glom(ps.sel.core, tax_level) 
    ps.sel.core <- ps.sel.core %>% aggregate_taxa(level = tax_level)

  }
  #yy<-as.matrix(tax_table(ps.sel.core)) %>% data.frame()
  #ps.sel.core<-if_else(tax_level=='OTU', ps.sel.core, speedyseq::tax_glom(ps.sel.core, tax_level)   )
  #
  ### make sure the reference level setting is right
  phyloseq::sample_data(ps.sel.core)$Group <- as.factor(phyloseq::sample_data(ps.sel.core)$Group)
  phyloseq::sample_data(ps.sel.core)$Group <- relevel(phyloseq::sample_data(ps.sel.core)$Group, ref = ref_)## ? change

  ###
  # ps.sel.core = transform_sample_counts(ps.sel.core, function(x) {1+x} )  #question
  # browser()
  ###### DESeq2 calculation for two-group comparison
  temp <- phyloseq::phyloseq_to_deseq2(ps.sel.core, ~ Group)
  temp <- DESeq2::DESeq(temp)
  temp.results <- DESeq2::results(temp)
  temp.results.df <- as.data.frame(temp.results)

  #### retrieve the taxonomy table
 
  sample_names(ps.sel.core) <-  gsub("-","_",sample_names(ps.sel.core))
  taxa_table <- metagMisc::phyloseq_to_df(ps.sel.core)[, 1:8] %>% tibble::column_to_rownames(., var = "OTU")

  ## get the merged results
  ps.sel.core.DESeq2 <- merge(taxa_table, temp.results.df, by = 0)
  # deseq2_foldchange<-0
  # deseq2_padj<-0.05
  ## select the top ones for downstream plotting, e.g, like the cutoff shown below
  ps.sel.core.DESeq2.top <- subset(ps.sel.core.DESeq2, abs(log2FoldChange) > deseq2_foldchange & padj < deseq2_padj & pvalue < deseq2_pvalue)


  ### define the color
  mycols <- c(
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
    "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
    "#000000"
  )
  dir.create( paste0(getwd(),"/out.result/deseq2") ,recursive = T)
  
  df1 <- ps.sel.core.DESeq2.top
  df1$OTU <- paste0(df1$Family,'|',df1$Genus,'|',df1$Row.names)
  
  df1$tax_columns<-paste0(df1[,tax_level],'|',df1$Row.names)
  write_csv(df1,paste0("out.result/deseq2/","deq2_",tax_level,".csv")  )
  # browser()
#######
df1$order<-df1$log2FoldChange %>%abs()
df1<-df1%>% arrange(dplyr::desc(order))

top_nrow<- ifelse( nrow(df1)<=20,nrow(df1),20 )
AD_plot_deseq2_top_nrow<<-top_nrow
if (tax_level=='OTU'){ df1$unique <-df1$OTU }
# browser()
p1 <- ggpubr::ggbarplot(df1[1:top_nrow,],
                        x = "unique", y = "log2FoldChange", width = 0.5,
                        fill = "red",
                        color = "white",
                        # palette = c("red"), ## if you have 20 different ones, then untick this, otherwise, leave it to default colors
                        sort.val = "asc",
                        sort.by.groups = FALSE,
                        x.text.angle = 90,
                        xlab = "",
                        ylab = "Fold Change (log2)",
                        # legend.title = tax_level,
                        rotate = TRUE,
                        ggtheme = theme_minimal()
)+
  # theme(axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 0,size=10),
        axis.text.y =element_text(size = 15), 
        axis.title.y = element_blank() ,
        # legend.text = element_text(size = 10) ,
        # legend.title = element_text(size = 12) ,
        axis.title.x = element_text(size = 18), 
        plot.title = element_text(size = 18, face="bold"))+
  scale_x_discrete(labels = label_wrap(30))+
  theme(legend.position='none')+
  ggtitle("Deseq2 analysis")
p1
p1 <- ggpubr::ggpar(p1, x.text.angle = 0)

ggsave(p1, filename =paste0("out.result/deseq2/","deq2_",tax_level,".png"), 
       dpi = 300, height = Deseq2_height, width = Deseq2_width, units = "in" )
  print("##############END plot_deseq2#################")
  # return(p1)

}





# spls.da_ <- function(tax_level,phyloseq,target_2rd_column){
# 				# df_return_g<- adjust_data(tax_level,phyloseq,target_2rd_column)
# 				# df_new_g<-df_return_g[[1]]
# 				# Y_g<-df_return_g[2][[1]]
#         print("##############sPLSDA#############")
#         data.filter<-otu_table(phyloseq)
#         Tax.filter<-tax_table(phyloseq)
#         Metadata<-meta(phyloseq)

#         Metadata[,target_2rd_column]<-factor(Metadata[,target_2rd_column])
#         Y_g<-Metadata[,target_2rd_column]
#         df_new_g<-as.data.frame(data.filter)


# 				splsda.tune.genus = tune.splsda(df_new_g,
# 				                          Y = Y_g,
# 				                          ncomp = 2,
# 				                          multilevel = NULL,
# 				                          logratio = 'CLR',
# 				                          test.keepX = c(seq(1,length(df_new_g), 1)),
# 				                          validation = c('Mfold'),
# 				                          folds = 5,
# 				                          dist = 'max.dist',
# 				                          nrepeat = 10,
# 				                          progressBar=FALSE,
# 				                          cpus = 8)

# 				plot(splsda.tune.genus, optimal = TRUE, sd = TRUE)


# 				#########Genus
# 				choice.ncomp.Genus = length(splsda.tune.genus$choice.keepX)
# 				choice.keepX.Genus = splsda.tune.genus$choice.keepX

# 				res.splsda.Genus = splsda(X = df_new_g,
# 				                    Y = Y_g,
# 				                    all.outputs = FALSE,
# 				                    ncomp = choice.ncomp.Genus,
# 				                    keepX = choice.keepX.Genus,
# 				                    logratio= "CLR")

# 				# write.table(df_new_g, file = paste0(work_path,"/output/test6.xls"),sep = "\t", row.names = T,quote = F)
# 				############


# 				##Genus
# 				list.splsda.genus = list()
# 				for(k in 1:choice.ncomp.Genus){
# 				  list.splsda.genus[[k]] = selectVar(res.splsda.Genus, comp = k)$name
# 				}
# 				list.splsda.genus
# 				##Genus
# 				splsda.perf.genus = perf(res.splsda.Genus, validation = 'Mfold', folds = 5,
# 				                   progressBar = FALSE, nrepeat = 100)
# 				plot(splsda.perf.genus)

# 				## VIP taxons
# 				#comp is 1
# 				res.splsda.Genus_plot_n1<-plotLoadings( res.splsda.Genus, comp = 1, method = 'mean', contrib = 'max',title=tax_level,ndisplay=30 )

# 				#dir.create(paste0(work_path,"/output/splsda"))

# 				#write.table(res.splsda.Genus_plot_n1, file = paste0(work_path,"/output/splsda/res.splsda.Genus_plot_n1.xls"),sep = "\t", row.names = T,quote = F)
# 				DT::datatable(res.splsda.Genus_plot_n1)

# 				#comp is 2
# 				#res.splsda.Genus_plot_n2<-plotLoadings(res.splsda.Genus, comp = 2, method = 'mean', contrib = 'max',title=tax_level,ndisplay=30 )
# 				#DT::datatable(res.splsda.Genus_plot_n2)
				
# 				#dir.create(paste0(work_path,"/output/splsda"))

# 				#write.table(res.splsda.Genus_plot_n2, file = paste0(work_path,"/output/splsda/res.splsda.Genus_plot_n2.xls"),sep = "\t", row.names = T,quote = F)

# 				#ggsave(res.splsda.Genus_plot_n1, filename =paste0("out.result/","res.splsda.Genus_plot_n1.",tax_level,".png"), dpi = 300, height = 5, width = 5, units = "in")
				
# 				# png(filename = paste0("out.result/","res.splsda.Genus_plot_n1.",tax_level,".png"))
# 				# #res.splsda.Genus_plot_n1
# 				# plotLoadings( res.splsda.Genus, comp = 1, method = 'mean', contrib = 'max',title=tax_level,ndisplay=30 )
# 			 #    dev.off()


# 			    dir.create("out.result/figures",recursive = T)
# 			    dir.create("out.result/tables",recursive = T)
# 			    res_splsda5 <- splsda.summary(res.splsda.Genus,outname = "out.result/figures/splsda_res_Hepa_1_6_in_situ")
# 				write.csv(res_splsda5,file = "out.result/tables/splsda_res_Hepa_1_6_in_situ.csv")

# 				splsda.replot("out.result/tables/splsda_res_Hepa_1_6_in_situ.csv",c("Hi","MOCK"),mycols[1:2])


# 			 #    png(filename = paste0("out.result/","res.splsda.Genus_plot_n2.",tax_level,".png")  )
# 				# plotLoadings( res.splsda.Genus, comp = 2, method = 'mean', contrib = 'max',title=tax_level,ndisplay=30 )
# 			 #    dev.off()
# 				#res.splsda.Genus_plot_n1
# 				}
# #spls.da_(tax_level,phyloseq,target_2rd_column)

plot_Metacoder0 <- function(phyloseq,target_2rd_column){
  print("##############Metacoder#############")
  #filter.meta() 
  #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
 # phyloseq<-readRDS( phyloseq_path   )
  browser()
  #meiji<-phyloseq
  phyloseq = transform_sample_counts(phyloseq, function(x) {1+x} )  #question
  #target_2rd_column<-"disease"
  #target_2rd_column<-input$SelectP3_1
  
  hmp_data <- parse_phyloseq(phyloseq) 
  sample_data<-sample_data(phyloseq)%>%data.frame() 
  
  sample_data$sample_id<-rownames(sample_data) 
  hmp_data$data$otu_prop <- calc_obs_props(hmp_data,data = "otu_table",cols = sample_data$sample_id) 
  hmp_data$data$tax_prop <- calc_taxon_abund(hmp_data, data = "otu_prop", cols = sample_data$sample_id) 
  
  # Plot data 
  plot_all <- hmp_data %>% 
    mutate_obs("tax_prop", abundance = rowMeans(hmp_data$data$tax_prop[sample_data$sample_id]) ) %>% 
    taxa::filter_taxa(abundance >= 0.001) %>% 
    taxa::filter_taxa(taxon_names != "") %>% # Some taxonomic levels are not named 
    heat_tree(node_size = n_obs, 
              node_size_range = c(0.01, 0.06), 
              node_size_axis_label = "Number of OTUs", 
              node_color = abundance, 
              node_color_axis_label = "Mean proportion of reads", 
              node_label = taxon_names) 
  # ###
  hmp_data$data$diff_table <- compare_groups(hmp_data,
                                             data = "tax_prop",
                                             cols = sample_data$sample_id,
                                             #groups = sample_data$microbiome.time.2 ,
                                             groups = sample_data[,target_2rd_column]
  )
  
  
  hmp_data <- mutate_obs(hmp_data, "diff_table",
                         wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"),
                         log2_median_ratio = ifelse(wilcox_p_value < 0.05 | is.na(wilcox_p_value), log2_median_ratio, 0))
  
  
  # hmp_data %>% 
  #   mutate_obs("tax_prop", abundance = rowMeans(hmp_data$data$tax_prop[sample_data$sample_id])) %>% 
  #   filter_taxa(abundance >= 0.001, reassign_obs = c(diff_table = FALSE)) %>% 
  #   heat_tree_matrix(data = "diff_table", 
  #                    node_size = n_obs, 
  #                    node_size_range = c(0.01, 0.05), 
  #                    node_label = taxon_names, 
  #                    node_color = log2_median_ratio, 
  #                    node_color_range = diverging_palette(), 
  #                    node_color_trans = "linear", 
  #                    node_color_interval = c(-7, 7), 
  #                    edge_color_interval = c(-7, 7), 
  #                    node_size_axis_label = "Number of OTUs", 
  #                    node_color_axis_label = "Log2 ratio median proportions", 
  #                    # initial_layout = "re", 
  #                    # layout = "da", 
  #                    key_size = 0.7, 
  #                    seed = 4) 
  # output_file = result_path("figure_3--hmp_matrix_plot")) 
  dir.create( paste0(getwd(),"/out.result/metacoder") ,recursive = T)
  ggsave(plot_all, filename =paste0("out.result/metacoder/","Metacoder",".png"), dpi = 300, height = 5, width = 5, units = "in")
  return(plot_all)
}


plot_Metacoder <- function(phyloseq,target_2rd_column,min_count=3 ){
  print("##############Metacoder#############")
  #filter.meta() 
  #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
  # phyloseq<-readRDS( phyloseq_path   )
  # browser()
  #meiji<-phyloseq
  
  # phyloseq = transform_sample_counts(phyloseq, function(x) {1+x} )  #question
  #target_2rd_column<-"disease"
  #target_2rd_column<-input$SelectP3_1
  
  hmp_data <- parse_phyloseq(phyloseq) 
  sample_data<-sample_data(phyloseq)%>%data.frame() 
  sample_data$sample_id<-rownames(sample_data) 
  
  #These low-abundance sequences might be the result of sequencing error, so typically we remove any counts/OTUs with less than some number of reads. Lets set all counts with less than 5 reads to zero, overwriting the original table:
  hmp_data$data$otu_table <- zero_low_counts(hmp_data, dataset = "otu_table", min_count = min_count)
  #By setting low abundance counts to zero we might have created OTUs that no longer contain any observations. We can check as follows.
  no_reads <- rowSums(hmp_data$data$otu_table[, sample_data$sample_id]) == 0
  sum(no_reads)
  #We can remove those OTUs and their associated taxa with filter_obs from the taxa package
  hmp_data <- filter_obs(hmp_data, target = "otu_table", ! no_reads, drop_taxa = TRUE)
  print(hmp_data)
  
  #Here we use the function calc_obs_props to divide each sampleâs counts by the total number of counts observed for each sample, resulting in a proportion.
  hmp_data$data$otu_prop <- calc_obs_props(hmp_data,data = "otu_table",cols = sample_data$sample_id) 
  #Currently, we have values for the abundance of each OTU, not each taxon. To get information on the taxa, we can sum the abundance per-taxon and add the results to the taxmap object in a new table:
  hmp_data$data$tax_prop <- calc_taxon_abund(hmp_data, data = "otu_prop", cols = sample_data$sample_id) 
  #We can also easily calculate the number of samples that have reads for each taxon:
  hmp_data$data$tax_occ <- calc_n_samples(hmp_data, "tax_prop", groups = sample_data[,target_2rd_column], cols = sample_data$sample_id)
  
  
  hmp_data$data$diff_table <- compare_groups(hmp_data,
                                        dataset = "tax_prop",
                                        cols = sample_data$sample_id, # What columns of sample data to use
                                        groups =  sample_data[,target_2rd_column] ) # What category each sample is assigned to
  print(hmp_data$data$diff_table)
  
  set.seed(999)
  plot_all1<-heat_tree(hmp_data, 
                       node_label = taxon_names,
                       node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                       node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                       node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
                       node_color_range = c("cyan", "gray", "tan"), # The color palette used
                       node_size_axis_label = "OTU count",
                       node_color_axis_label = "Log 2 ratio of median proportions",
                       layout = "davidson-harel", # The primary layout algorithm
                       initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
  #diff tree
  #Note that we have not taken into account statistical significance when showing this, so lets do that. First, we need to correct for multiple comparisons:
  hmp_data$data$diff_table$wilcox_p_value <- p.adjust(hmp_data$data$diff_table$wilcox_p_value,method = "fdr")
  #If we then look at the distribution of p-values, we can see that none are even close to significant:
  range(hmp_data$data$diff_table$wilcox_p_value, finite = TRUE) 
  #There is no need to graph this, but if there still were some significant differences, we could set any difference that is not significant to zero and repeat the last heat_tree command doing something like:
  hmp_data$data$diff_table$log2_median_ratio[hmp_data$data$diff_table$wilcox_p_value > 0.05] <- 0
  
  plot_all2<-heat_tree(hmp_data, 
                       node_label = taxon_names,
                       node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                       node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                       node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
                       node_color_range = c("cyan", "gray", "tan"), # The color palette used
                       node_size_axis_label = "OTU count",
                       node_color_axis_label = "Log 2 ratio of median proportions",
                       layout = "davidson-harel", # The primary layout algorithm
                       initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
  
  dir.create( paste0(getwd(),"/out.result/metacoder") ,recursive = T)
  ggsave(plot_all1, filename =paste0("out.result/metacoder/","Metacoder1",".png"), dpi = 300, height = 5, width = 8, units = "in")
  ggsave(plot_all2, filename =paste0("out.result/metacoder/","Metacoder2_diff_tree",".png"), dpi = 300, height = 5, width = 8, units = "in")
  return(plot_all1)
}



# Wilcox / Anova analysis
  
library(doBy) #Ê¹ÓÃÆäÖÐµÄ summaryBy() ÒÔ·½±ã°´·Ö×é¼ÆËã¾ùÖµ¡¢ÖÐÎ»Êý
#library(dplyr)
# library(plyr)
library(ggpubr)
library(purrr)
#library(taRifx)
library(patchwork)
library(stringr)
library(ggplot2)
library("Rmisc")
library(janitor)

wilcox_plot_<-function(subset_,two_group_classfication2,Pircust_df,kegg_id_list,p_value_threshold,analysis_method){
  
  color_b1 <- c('#e41a1c','#377eb8','#4daf4a')
  plot_list <- list()
  plot_list2<- list()
  sig_df<-list()
  stat_table_list <- list()
  # browser()
  for (n in 1:nrow(Pircust_df) ) {#nrow(Pircust_df)
    # print(n)
    
    #clean data
    # n=3
    # subset_<-subset_1
    
    Pircust_df_n <- data.frame(  t( Pircust_df[n,colnames(Pircust_df) %in% subset_$sample] ) ) #ÐèÒªÐÞ¸Ä·ÖÎö×é±ðµÄÐÅÏ¢
    #gene_id <- names(Pircust_df_n)[1]
    # Pircust_df_n<-rownames_to_column(Pircust_df_n,var = 'ko_id')
    names(Pircust_df_n)[1] <- 'ko_id'
    
    Pircust_df_n$sample <- rownames(Pircust_df_n) #Pircust_df_n$sample  is SampleID
    Pircust_df_n <- merge(Pircust_df_n, subset_, by = 'sample', all.x = TRUE)
    # Pircust_df_n <- merge(Pircust_df_n, subset_, by.x= 'ko_id' ,by.y='sample', all.x = TRUE)
    
    if (sum(Pircust_df_n$ko_id)==0) next
    Pircust_df_n$group <- factor(Pircust_df_n$group)
    
    ko_id_n <- rownames(Pircust_df[n,])
    Pircust_df_n$Kegg_name <- kegg_id_list[ko_id_n]%>%as.character()
    names(Pircust_df_n)[2] <- 'Value'
    
    #other
    summ <- Pircust_df_n[,c('group','Value')] %>%group_by(group) %>%dplyr::summarize(mean = mean(Value)%>%round(3), median = median(Value)%>%round(3), sd = sd(Value)%>%round(3)  )
    summ$keggid<-ko_id_n
    # summ$kegg_name<-kegg_id_list[ko_id_n]
    
    if (two_group.plot.split){
      Pircust_df_n_mean.df<-Pircust_df_n%>%dplyr::select(group,Value) %>%group_by(group) %>% dplyr::summarize(mean = mean(Value)) %>% data.frame()
      mean_list<-list()
      if (Pircust_df_n_mean.df[,"mean"][1] > Pircust_df_n_mean.df[,"mean"][2]) { 
        Pircust_df_n$NP <- "P"
      }else {
        Pircust_df_n$NP <- "N"
      }
      mean_list
    } 
    
    
    #get significat featrue
    ##get significat table
    compare_means_stat<-compare_means(Value ~ group, data = Pircust_df_n,method = analysis_method)#get p.adj p.format p.signif method
    summ$p.adj<-compare_means_stat$p.adj##get signifcant p.adj
    #print(kegg_id_list[ko_id_n])
    #print( Pircust_df[n,] )
    #print(summ)
    summ_<-merge(summ[1,],summ[2,],by='keggid')
    summ_$p.adj<-summ_$p.adj.x
    summ_<-summ_%>%dplyr::select(-p.adj.x,-p.adj.y)
    
    
    stat_table_list[[ko_id_n]]<-summ_
    
    ###ggbarplot
    if ( nrow( compare_means_stat ) >=1  ){#get signifcant plot
      if ( compare_means_stat$p.adj <p_value_threshold ){
        # print(paste('num:',n,':',ko_id_n))
        # plot_list2[[n]]<-p2
        sig_df[[n]] <- Pircust_df_n
      }
    }
    #get significat featrue
    
  }
  
  
  #plot
  ##ggbarplot overall
  sig_df_overall<-do.call(rbind.data.frame, sig_df)
  
  stat_table_list_zong<-do.call(rbind,stat_table_list) %>% filter(p.adj<p_value_threshold)
  stat_table_list_zong$kegg_name <- purrr::map(stat_table_list_zong$keggid,function(x){kegg_id_list[x]}) %>% unlist()
  
  stat_table_list_zong<-stat_table_list_zong %>% dplyr::select(keggid, kegg_name,everything())#±£³ÖÄ³Ð©±äÁ¿ÔÚÇ°£¬ÆäÓàÔÚºó
  
  # DT::datatable(stat_table_list_zong)
  # print('plot stat_table_list_zong:',length(stat_table_list_zong))
  return(list(sig_df_overall,stat_table_list_zong))
}



##wilcox ¼ìÑéÅú´¦ÀíÊ¾Àý
# work_path<-"/home/zhiyu/data/metagenome/SDQ_F/Pircust-wilco20200330"
# meta_path<-paste0(work_path,"/data/SDQ_F_2020_01_13_1_update.rdsphyloseq_sample.meta2020_01_16_4.csv")
# Pircust_df_name<- paste0(work_path,"/data/path_abun_unstrat_descrip_MJ20200220019-MJ-M-20200220006-SDQ-16samples-rawdata-20200319.tsv")

  
  
pircust_wilcox<-function(work_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2,p_value_threshold=0.05,analysis_method="wilcox.test"){
  #setwd(work_path) 
  #dir.create(paste0(work_path,"/output"))
  # browser()
  print('11111111')
  print(Pircust_df_name)
  Pircust_df <- read.csv(Pircust_df_name, sep = '\t', row.names = 1, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  # Meta_df <- read_tsv(meta_path)
  Meta_df <- read.csv(meta_path, sep = ',', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  two_group.plot.split<<-TRUE
  ###################
  Pircust_df <- Pircust_df  %>%  adorn_percentages("col") # you can convert abundance into percentage.
  
  
  ##ÉèÖÃÍ¨ÓÃKEGGID

  kegg_id_list<-list()
  # kegg_id_list[ko_description$kegg_id]<-ko_description$description
  kegg_id_list[rownames(Pircust_df)]<-Pircust_df$description
  
  ko_description<-Pircust_df%>% dplyr::select("description")
  ko_description$kegg_id<- rownames(Pircust_df)
  
  
  Pircust_df <- Pircust_df[,2:ncol(Pircust_df)]*100
  Meta_df$sample <-Meta_df$SampleID
  

  
  # Wilcox/Anova plot
  #setwd(work_path)
 
  # color_b1 <- c('#e41a1c','#377eb8','#4daf4a')
  # analysis_method<<-"wilcox.test"#anova    wilcox.test
  
  Meta_df$group <- Meta_df[,target_column2]
  # browser()
  subset_1<- subset( Meta_df, grepl( paste0(two_group_classfication1,collapse = '|') ,Meta_df[,target_column1])  )
  print(subset_1)
  #subset_1<- Meta_df
  
  rt<-NULL
  print(analysis_method)
  print(p_value_threshold)
  rt<-wilcox_plot_(subset_1,two_group_classfication2,Pircust_df,kegg_id_list,p_value_threshold,analysis_method)
  
  ## Using ggplot2
  test<-rt[[1]]%>%data.frame()
  
  test$group<-factor(trimws(test$group))#ÖØÖÃfactor
  test$Kegg_name<-factor(trimws(test$Kegg_name))#ÖØÖÃfactor
  test$NP<-factor(trimws(test$NP))#ÖØÖÃfactor
  
  mycols<-c("#00468BFF", "#ED0000FF")
  # mycols <- color.gvri(c("#e6194b","#4363d8"))
  tgc <- summarySE(test, measurevar="Value", groupvars=c("Kegg_name","group","NP"))
  tgc
  # browser()
  tgc$Group<-tgc$group
  c1<-ggplot(tgc, aes(x=Kegg_name, y=Value, fill=Group)) +
    geom_bar(stat="identity", position=position_dodge() )+  
    scale_fill_manual(values=mycols)+
    geom_errorbar(aes(ymin=Value-se, ymax=Value+se),   size=.3, width=.2, position=position_dodge(.9))  +
    coord_flip()+facet_grid(NP ~ ., scales = "free_y", space = "free_y")+
    theme_bw()+theme_minimal()+
    theme( legend.position = "right",
           axis.title.x = element_blank(),
           axis.text.x = element_text(size = 14,color="black"),
           axis.text.y = element_text(size = 14,color="black"), 
           legend.text = element_text(size = 16) ,
           legend.title = element_text(size = 23) ,
           axis.title.y = element_blank()  )
  c1
  
  stat_table_list_zong<-rt[[2]]
  #if (two_group_classfication1 =="*"){two_group_classfication1="all"}
  # ggsave(c1, filename = paste0("out.result/Pircust/output/pircust.",two_group_classfication1,'.',paste0(two_group_classfication2,collapse ='.'),'sig.png' ) ,dpi = 300, height = 12, width = 12, units = "in")
  
  # write.csv(stat_table_list_zong,paste0("out.result/Pircust/output/pircust.",two_group_classfication1,'.',paste0(two_group_classfication2,collapse ='.'),'sig.csv' ),quote=F  )

  ggsave(c1, filename = paste0("out.result/Pircust/output/pircust.",'.',paste0(two_group_classfication2,collapse ='.'),'sig.png' ) ,dpi = 300, height = 12, width = 12, units = "in")
  
  write.csv(stat_table_list_zong,paste0("out.result/Pircust/output/pircust.",'.',paste0(two_group_classfication2,collapse ='.'),'sig.csv' ),quote=F  )
  
  DT::datatable(stat_table_list_zong)
  return(c1)
}











pircust_aldex2<-function(work_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2){
library(ALDEx2)

# target_column2<-'group'
# two_group_classfication2<-c('NC', 'KO4')

    print(Pircust_df_name)
    Pircust_df <- read.csv(Pircust_df_name, sep = '\t', row.names = 1, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    Meta_df <- read.csv(meta_path, sep = ',', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    two_group.plot.split<<-TRUE
    ###################
    Meta_df$sample <-Meta_df$SampleID
    
    ko_description<-Pircust_df%>% dplyr::select("description")
    ko_description$kegg_id<-rownames(Pircust_df)
    kegg_id_list<-list()
    kegg_id_list[ko_description$kegg_id]<-ko_description$description
    
    # Wilcox/Anova plot
    #setwd(work_path)
    color_b1 <- c('#e41a1c','#377eb8','#4daf4a')
    analysis_method<<-"wilcox.test"#anova    wilcox.test
    
    Meta_df$group <- Meta_df[,target_column2]
    #subset_1<- subset(Meta_df, Model %in% two_group_classfication1)
    subset_1<- subset( Meta_df, grepl(two_group_classfication1,Meta_df[,target_column1])  )
    subset_1<- subset(Meta_df, Meta_df[,target_column2] %in% two_group_classfication2)
    print(subset_1)
    #subset_1<- Meta_df

#rownames(Pircust_df)<-Pircust_df_raw$description
Pircust_df <- Pircust_df[,2:ncol(Pircust_df)]

#######################################################################clear data
SampleID<-colnames(Pircust_df)# Pircust ID 
Pircust_df.df<-data.frame(SampleID)# Pircust ID 
Pircust_df.meta<-merge(Pircust_df.df,subset_1,by="SampleID")
conds <- Pircust_df.meta$group# get the vector of group

Pircust_df_screen <-Pircust_df[,Pircust_df.meta$SampleID]
Pircust_df_screen2 <-apply(Pircust_df_screen, 2, as.integer)%>%as.data.frame() #get KEGG table of abuntu
rownames(Pircust_df_screen2)<-rownames(Pircust_df_screen)

#########run aldex###############mainly need Pircust_df_screen2 and conds,conds must have a correct order.########
x <- aldex.clr(Pircust_df_screen2, conds, mc.samples=128, denom="all", verbose=F)

x.kw <- aldex.kw(x)
x.kw$SampleID<-rownames(x.kw)
x.effect <- aldex.effect(x, verbose=FALSE)
x.effect$SampleID<-rownames(x.effect)

result.df<-cbind(x.kw ,x.effect) %>%data.frame()
# df$SampleID<-rownames(df)

df<-result.df %>% dplyr::filter(kw.eBH<0.05) %>% arrange(dplyr::desc(effect))

# df$SampleID<-rownames(df)


    
    
df<-merge(ko_description,df,by.x ="kegg_id",by.y ="SampleID" )

c2<-ggplot(df, aes(x=reorder(description, abs(effect) ), y=effect)) +
  geom_bar(stat="identity", position=position_dodge() ,fill="red")+  
  
  #scale_colour_manual(  c("red" ) )+
  #scale_fill_manual(values=c("#4363d8"))+
  coord_flip()+
  theme_bw()+theme_minimal()+
  theme( legend.position = "right",axis.title.x = element_blank(),axis.text.x =element_text(angle = 0),axis.title.y = element_blank()  )+
    #labs(legend = paste0(two_group_classfication2,collapse ='|') )+
  ggtitle(label = paste0(two_group_classfication2,collapse ='|') )


print(c2)



write.csv(df,paste0("out.result/Pircust/output/pircust.aldex",'.',paste0(two_group_classfication2,collapse ='.'),'sig.csv' ),quote=F  )

ggsave(c2, filename = paste0("out.result/Pircust/output/pircust.aldex",'.',paste0(two_group_classfication2,collapse ='.'),'sig.png' ) ,dpi = 300, height =3, width = 9, units = "in")
DT::datatable(df)
DT::datatable(result.df)
return(c2)

}

# work_path<-"/home/zhiyu/data/luo/16s_shiny/test/out.result/Pircust"
# meta_path<-"/home/zhiyu/data/paper/all.meta/raw/SRP189432.meta.csv"
# Pircust_df_name<- paste0(work_path,"/output/pircust2/KEGG/picrust2_out_pipeline/KEGG_pathways_unstrat/path_abun_unstrat_descrip.tsv")
# target_column1<-'SampleID'
# two_group_classfication1<-c('^male')#LLC-SDQ-1  LLC-YWJ-1 LLC-GV
# target_column2<-'disease'
# two_group_classfication2<-c('Ulcerative Colitis', 'Healthy Control')

# c1<-pircust_wilcox(work_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2)



pircust_deseq2<-function(work_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2){
    # library(ALDEx2)
    # target_column2<-'group'
    # two_group_classfication2<-c('NC', 'KO4')
    print(work_path)
    print(meta_path)
    print(Pircust_df_name)
    print(target_column1)
    print(two_group_classfication1)
    print(target_column2)
    print(two_group_classfication2)

    Pircust_df <- read.csv(Pircust_df_name, sep = '\t', row.names = 1, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    Meta_df <- read.csv(meta_path, sep = ',', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    two_group.plot.split<<-TRUE
    ###################
    Meta_df$sample <-Meta_df$SampleID
    
    ko_description<-Pircust_df%>% dplyr::select("description")
    ko_description$kegg_id<-rownames(Pircust_df)
    kegg_id_list<-list()
    kegg_id_list[ko_description$kegg_id]<-ko_description$description
    
    # Wilcox/Anova plot
    #setwd(work_path)
    color_b1 <- c('#e41a1c','#377eb8','#4daf4a')
    analysis_method<<-"wilcox.test"#anova    wilcox.test
    
    Meta_df$group <- Meta_df[,target_column2]
    #subset_1<- subset(Meta_df, Model %in% two_group_classfication1)
    subset_1<- subset( Meta_df, grepl(two_group_classfication1,Meta_df[,target_column1])  )
    subset_1<- subset(Meta_df, Meta_df[,target_column2] %in% two_group_classfication2)
    print(subset_1)
    #subset_1<- Meta_df

    #rownames(Pircust_df)<-Pircust_df_raw$description
    Pircust_df <- Pircust_df[,2:ncol(Pircust_df)]
    
    #######################################################################clear data
    SampleID<-colnames(Pircust_df)# Pircust ID 
    Pircust_df.df<-data.frame(SampleID)# Pircust ID 
    Pircust_df.meta<-merge(Pircust_df.df,subset_1,by="SampleID")
    conds <- Pircust_df.meta$group# get the vector of group
    
    Pircust_df_screen <-Pircust_df[,Pircust_df.meta$SampleID]
    Pircust_df_screen2 <-apply(Pircust_df_screen, 2, as.integer)%>%as.data.frame() #get KEGG table of abuntu
    rownames(Pircust_df_screen2)<-rownames(Pircust_df_screen)
    
    library("DESeq2")
    pasilla = DESeqDataSetFromMatrix(
      countData = Pircust_df_screen2,
      colData   = Pircust_df.meta,
      design    = ~ group)
    
    pasilla = DESeq(pasilla)
    
    res = results(pasilla)
    res = res[order(res$padj), ] %>% data.frame()
    res$kegg_id <- rownames(res)
    
    df<-merge(ko_description,res,by.x ="kegg_id",by.y ="kegg_id" )
    # browser()
    
    df<-df %>% dplyr::filter(pvalue<0.1 & abs(log2FoldChange) >1 ) %>% arrange(pvalue,log2FoldChange)
    
    c2<-ggplot(df, aes(x=reorder(description, abs(log2FoldChange) ), y=log2FoldChange)) +
      geom_bar(stat="identity", position=position_dodge() ,fill="red")+  
      
      #scale_colour_manual(  c("red" ) )+
      #scale_fill_manual(values=c("#4363d8"))+
      coord_flip()+
      theme_bw()+theme_minimal()+
      theme( legend.position = "right",axis.title.x = element_blank(),axis.text.x =element_text(angle = 0),axis.title.y = element_blank()  )+
        #theme(axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))+
         theme(axis.title.x = element_text(size = 16))+
        theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"))
       
        #labs(legend = paste0(two_group_classfication2,collapse ='|') )+
      ggtitle(label = paste0(two_group_classfication2,collapse ='|') )
    
    
    print(c2)
    return(c2)
    # write.csv(df,paste0(work_path,"/output/bar.",two_group_classfication1,'.',paste0(two_group_classfication2,collapse ='.'),'sig_aldex.csv' ),quote=F  )
    # ggsave(c2, filename = paste0(work_path,"/output/bar.",two_group_classfication1,'.',paste0(two_group_classfication2,collapse ='.'),'sig_aldex.png' ) ,dpi = 300, height =3, width = 9, units = "in")
}


pircust_deseq2_mei<-function(current_dir_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2){
  # browser()
  # library(ALDEx2)
  # library("DESeq2")
  # target_column1<-"host_disease"
  # target_column2<-'host_disease'
  # two_group_classfication1<-c('ulcerative colitis', 'healthy')
  # two_group_classfication2<-c('ulcerative colitis', 'healthy')
  # Pircust_df_name<-"/home/zhiyu/data/luo/16s_shiny/test/out.result/Pircust/output/pircust2/KEGG/picrust2_out_pipeline/EC_metagenome_out_unstrat/pred_metagenome_unstrat.tsv"
  # meta_path<-"/home/zhiyu/data/luo/16s_shiny/test/www/all.meta/filter/filter_meta.tsv"
  # work_path<-"/home/zhiyu/data/luo/16s_shiny/test"
  
  #current_dir_path<-work_path
  work_path<- paste0(current_dir_path,'/out.result/Pircust')

  print(Pircust_df_name)
  Pircust_df <- read.csv(Pircust_df_name, sep = '\t', row.names = 1, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  Meta_df <- read.csv(meta_path, sep = ',', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  two_group.plot.split<<-TRUE
  ###################
  Meta_df$sample <-Meta_df$SampleID
  
  #ko_description<-Pircust_df%>% dplyr::select("description")
  #ko_description$kegg_id<-rownames(Pircust_df)
  #kegg_id_list<-list()
  #kegg_id_list[ko_description$kegg_id]<-ko_description$description
  
  # Wilcox/Anova plot
  #setwd(work_path)
  color_b1 <- c('#e41a1c','#377eb8','#4daf4a')
  analysis_method<<-"wilcox.test"#anova    wilcox.test
  
  Meta_df$group <- Meta_df[,target_column2]
  #subset_1<- subset(Meta_df, Model %in% two_group_classfication1)
  subset_1<- subset( Meta_df, grepl(two_group_classfication1,Meta_df[,target_column1])  )
  subset_1<- subset(Meta_df, Meta_df[,target_column2] %in% two_group_classfication2)
  print(subset_1)
  #subset_1<- Meta_df
  
  #rownames(Pircust_df)<-Pircust_df_raw$description
  Pircust_df <- Pircust_df[,2:ncol(Pircust_df)]
  
  #######################################################################clear data
  SampleID<-colnames(Pircust_df)# Pircust ID 
  Pircust_df.df<-data.frame(SampleID)# Pircust ID 
  Pircust_df.meta<-merge(Pircust_df.df,subset_1,by="SampleID")
  conds <- Pircust_df.meta$group# get the vector of group
  
  Pircust_df_screen <-Pircust_df[,Pircust_df.meta$SampleID]
  Pircust_df_screen2 <-apply(Pircust_df_screen, 2, as.integer)%>%as.data.frame() #get KEGG table of abuntu
  rownames(Pircust_df_screen2)<-rownames(Pircust_df_screen)
  
  
  # browser()
  # library("DESeq2")
  pasilla = DESeqDataSetFromMatrix(
    countData = Pircust_df_screen2,
    colData   = Pircust_df.meta,
    design    = ~ group)
  
  pasilla = DESeq(pasilla)
  
  res = results(pasilla)
  res = res[order(res$padj), ] %>% data.frame()
  res$kegg_id <- rownames(res)
  
  #df<-merge(ko_description,res,by.x ="kegg_id",by.y ="kegg_id" )
  
  df<-res %>% dplyr::filter(pvalue<0.1 & abs(log2FoldChange) >1 ) %>% arrange(pvalue,log2FoldChange)
  df<-df %>% dplyr::select("kegg_id")
  out.deseq2_mei.path<-paste0(work_path,"/output/bar.sig_deseq2_mei.txt")
  print(out.deseq2_mei.path)
  write.csv(df, out.deseq2_mei.path ,quote=F  ,col.names = FALSE,row.names = FALSE)
  return(df)

}

NetCoMi_fc<-function(phyloseq_,target_column2,two_group_classfication2,tax_level){
  library(NetCoMi)
  # browser()
  # group_column<-"disease"
  #group_columns_content<- c("Healthy Control", "Ulcerative Colitis")
  group_column<-target_column2
  group_columns_content<- two_group_classfication2
  # phyloseq_<<-phyloseq
  # phyloseq_<-phyloseq
  #file_meta<<-filter.meta
  
  sample_names(phyloseq_) <-  gsub("-",".",sample_names(phyloseq_))# maybe changed
  phyloseq_2 <- core(phyloseq_, detection = 0.0005, prevalence = 50/100)
  
  otuTable(phyloseq_2)<-otuTable(phyloseq_2) +1
  #sample_data(phyloseq_)<-file_meta 
  otu.mt<-otu_table(phyloseq_2) %>% data.frame() %>% as.matrix() %>% t()
  otu_table(phyloseq_2) <-otu_table(otu.mt,taxa_are_rows=F)
  amgut2.filt.phy<-phyloseq_2
  amgut_split <- metagMisc::phyloseq_sep_variable(amgut2.filt.phy, group_column,drop_zeroes=FALSE)#"SEASONAL_ALLERGIES"
  meta_<-sample_data(amgut2.filt.phy) %>% data.frame()
  amgut_split

  phy1<-amgut_split[[1]]
  print(phy1)
  phy2<-amgut_split[[2]]
  print(phy2)
  
  # net_season <- netConstruct(phy1, phy2, verbose = 2,
  #                            filtTax = "highestVar",filtTaxPar = list(highestVar = 50),measure = "spieceasi",
  #                            measurePar = list(method = "mb", nlambda=20, pulsar.params=list(rep.num=3)),
  #                            normMethod = "none", zeroMethod = "none",sparsMethod = "none", seed = 123456 )
  
  net_season <- netConstruct(phy1, phy2, verbose = 0, 
                              filtTax = "highestVar",
                              filtTaxPar = list(highestVar = 50),
                              measure = "pearson", normMethod = "TSS",
                              sparsMethod = "threshold", thresh = 0.3)
  
  # net_season <- netConstruct(phy1, phy2, verbose = 2,
  #                            filtTax = "pearson",filtTaxPar = list(highestVar = 50),measure = "spieceasi",
  #                            normMethod = "none", zeroMethod = "none",sparsMethod = "none", seed = 123456 )
  # measurePar = list(method = "mb", nlambda=20, pulsar.params=list(rep.num=3)),
  
  props_season <- netAnalyze(net_season, clustMethod = "cluster_fast_greedy")
  labels_t<-paste0(props_season$input$adjaMat1 %>% rownames())
  
  #u can change the lables of nodes
  # otu_tax<-phyloseq_to_df(phyloseq_)
  # rownames(otu_tax)<-otu_tax$OTU
  # otu_tax<-otu_tax[c(labels_t),]
  # otu_tax$tax_otu<-paste0(otu_tax[,tax_level],'|',otu_tax$OTU)

  # browser()
  dir.create( paste0(getwd(),"/out.result/NetCoMi") ,recursive = T)
  png(filename = paste0("out.result/NetCoMi/","NetCoMi_",tax_level,".png") , width =10, height = 7, units = "in", res = 300)
  #netcomi_plot <-plot(props_season, sameLayout = TRUE, layoutGroup = 1,nodeSize = "eigenvector", cexNodes = 1.5, cexLabels = 3,groupNames = group_columns_content,labels=otu_tax$tax_otu)
  netcomi_plot <-plot(props_season, sameLayout = TRUE, layoutGroup = 1,
                      nodeSize = "eigenvector", cexNodes = 1.5, cexLabels = 3,
                      groupNames = group_columns_content)
  dev.off()
  
  # return(netcomi_plot)
  # netcomi_plot <-plot(props_season, sameLayout = TRUE, layoutGroup = 1,
  #                     nodeSize = "eigenvector", cexNodes = 1.5, cexLabels = 3,
  #                     groupNames = group_columns_content,labels=otu_tax$tax_otu)
}





r_lefse_func <- function(group1, group2, phy_core_sub, ref = group1, group = "Groups",outcome_add_info="",level=level,
                         lda_cutoff_values=3,lefse_height=lefse_height,lefse_width=lefse_width){
  library(microbiomeMarker) 
  # browser()
  
  phy_core_sub@sam_data$GROUPS<-phy_core_sub@sam_data[,group] %>% unlist()
  # phy_core_sub2 <- subset_samples(phy_core_sub,(GROUPS == "HC" |GROUPS == "MS" ) )
  # phy_core_sub2 <- subset_samples(phy_core_sub,(GROUPS == group1 |GROUPS == group2 ) )
  phy_core_sub2<-phy_core_sub
  meta_df<-phy_core_sub2@sam_data %>% data.frame()
  
  # level<-"Family"
  # phy_core_sub2 <- phy_core_sub2 %>% aggregate_taxa(level = "Genus")
  ##get marker feature in lefse
  mm <- run_lefse(
    phy_core_sub2, 
    # normalization = 1e6, 
    group = group , 
    taxa_rank=level,
    kw_cutoff = 0.01,
    wilcoxon_cutoff = 0.01,
    lda_cutoff =lda_cutoff_values,
    # multicls_strat = TRUE,
    # summarize= FALSE 
  )
  mm
  
  ## barplot
  #mycol2 <- c("black","red","#EC9EC0","#EBE0D0")
  mycol2<-c("#00468BFF", "#ED0000FF")
  group_vec <- unique(meta_df[,group])
  names(mycol2) <- group_vec
  ysize=15
  lsize=25
  p <-   plot_ef_bar(mm, label_level = 3) + scale_fill_manual(values = mycol2) +
    theme(axis.text.y=element_text(size=ysize),
          axis.text.x=element_text(size=ysize),
          legend.text=element_text(size=lsize),
          axis.title.x=element_text(size=(ysize+2)),
          axis.title.y=element_text(size=(ysize+2)),
          legend.title = element_text(size = lsize) )
  
  dir.create( "out.result/lefse",recursive = TRUE)
  # ggsave(p, filename = work_path %+% "out.result/lefse/" %+% group1 %+% group2 %+% outcome_add_info %+% "_lefse_markder_barplot.svg", device ="svg", height = 14, width =20)
  ggsave(p, filename = "out.result/lefse/"  %+% outcome_add_info %+% "lefse_markder_" %+%  level %+% "_barplot.png", device ="png", height = lefse_height, width =lefse_width,units = "px")
  
  # p
  
}




