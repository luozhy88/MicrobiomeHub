library(shiny)
library(shinythemes)
library(plyr) 
library(dplyr) 
library(DT)
library(metacoder)
library(purrr)
library(ggplot2)
library(microbiome)
library(phyloseq)
library(mixOmics)
library(metagMisc)
library(png)
# install.packages("scales")
library("scales")
# library(future)
library(promises)
`%+%` <- function(a,b) {paste0(a,b)}


print( paste0("###############################################################Start data:",Sys.time()) )

# filename0 <- normalizePath('/home/zhiyu/data/luo/16s_shiny/test/www/zi.png')
# #all.phyloseq.path<-"/home/zhiyu/data/paper/all.phyloseq.upgrade/phyloseq.rds" 
# all.phyloseq.dir<-"/home/zhiyu/data/paper/all.phyloseq" 
# path<-"/home/zhiyu/data/paper/all.meta/raw" 
# filter_meta_path<-"/home/zhiyu/data/paper/all.meta/filter/filter_meta.tsv"
# phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"


print('weizhi!!!!!!!!!!!!!!!!!')
print( getwd() ) 

work_path<- paste0(getwd() ,"/")
#work_path<-"/home/zhiyu/data/luo/16s_shiny/test/"
#filename0 <- normalizePath('www/zi.png')
#all.phyloseq.path<-"/home/zhiyu/data/paper/all.phyloseq.upgrade/phyloseq.rds" 
all.phyloseq.dir<-paste0(work_path,"www/all.phyloseq" )
path<-paste0(work_path,"www/all.meta/raw" )
filter_meta_path<-paste0(work_path,"www/all.meta/filter/filter_meta.csv")
phyloseq_path<-paste0(work_path,"www/phyloseq.upgrade/FALSE_update.rds" )
pre_analysis_meta_path<-paste0(work_path,"www/analysis_file/pre_meta.csv" )
pre_analysis_phyloseq_path<-paste0(work_path,"www/analysis_file/pre_phyloseq.rds" )
pre_analysis_merge_all.phyloseq_path<-paste0(work_path,"www/all.phyloseq/merge/merge_all.phyloseq.rds" )
# conda_site<<-"/home/lzy/data/software/anaconda3_2020"
# conda_site<<-"/home/zhiyu/data/software/Anaconda3/envs/qiime2-2019.10"
cpus_num <<- 3

input.upload.phyloseq.datapath<-NULL
input.upload.meta.datapath<-NULL
dir.create(paste0(work_path,"out.result" ))


###meta path
list.meta.path<-list()
list.meta.path[list.files(path, pattern = "meta", full.names = F)]<-list.files(path, pattern = "meta", full.names = T)
print(list.meta.path)
list.meta.path.basename<-names(list.meta.path)
list.meta.path.basename2<-as.list(list.meta.path.basename)
names(list.meta.path.basename2)<-list.meta.path.basename2
###meta path

paste0(work_path,"www/phyloseq.upgrade" )

system(paste0("rm  -rf ",work_path,"out.result" ))
dir.create(paste0(work_path,"out.result" ),recursive = T)


# server

server <- function(input, output, session) {
  ##################################################################################observe
  rt <- reactiveValues()
  ##P1 observe
  observe({#Select content of column in meta   Select input str(example: ulcerative colitis,healthy )
    lie_name<-input$SelectP1_2
    print(paste0("ddd 第二个参数:",lie_name))
    file.meta<-rt$file_meta
    luo<-file.meta[,lie_name]
    filter_str_lie_name<-unique(luo) 
    print(class(filter_str_lie_name))
    print(paste0("ddd 第3个参数:",filter_str_lie_name ))
    updateSelectizeInput(session, "SelectP1_3",label = paste("Select two classficiations(example: ulcerative colitis,healthy )", length( filter_str_lie_name )),choices = c("All",filter_str_lie_name ) ,options = list(maxItems = 5L) )# page1
  })
  

  
  ##################################################################################reactive
  
  file.all_meta_id_num<- reactive({ distribution_project_num(list.meta.path)  }) #项目情况绘图
  
  filter.meta<- reactive({  print('guolv meta')
    print(input$SelectP1_3)
    print(input$SelectP1_2)
    file.meta<-rt$file_meta
    #print(rt$file_meta)
    #rt$filter_meta <- rt$file_meta %>% filter( grepl(input$SelectP1_3 , input$SelectP1_2 %>% get() ))  
    ifelse (input$SelectP1_3 == "All", rt$filter_meta <- file.meta, rt$filter_meta <- file.meta[file.meta[,input$SelectP1_2] %in% c(input$SelectP1_3), ]  )
    # browser()
    # utils::write.table(rt$filter_meta %>% data.frame(),filter_meta_path,quote = F,row.names = F,sep = "\t")    
    write.csv(rt$filter_meta %>% data.frame(),filter_meta_path,quote = T,row.names = F)    
    rt$filter_meta    }) 
  
  file.phyloseq<- reactive({
    #file.meta.name1<-"/home/zhiyu/data/paper/all.meta/SRP126121.meta.csv"
    system("pwd")
    dir_phyloseq.upgrade<-paste0(work_path,"www/phyloseq.upgrade" )
    system(paste0("rm ",dir_phyloseq.upgrade,"/*"))#筛选后的meta
    
    if (!is.null(input.upload.phyloseq.datapath) | !is.null(input.upload.meta.datapath) ){
      #"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      saveRDS(phyloseq_, phyloseq_path )
      
    }else{
      # browser()
      print(input$SelectP1_1)
      target_phy_path<-paste0(all.phyloseq.dir %+% '/' %+%  str_replace(input$SelectP1_1,".meta.csv",'_phyloseq.rds'))
      # command1<-sprintf("/home/lzy/data/software/anaconda3_2020/envs/pircust2/bin/Rscript  %swww/script/merge_metafrom_phyloseq_merge_allphyloseq2.R -f %s -t %s -o %s",work_path,all.phyloseq.dir,filter_meta_path,dir_phyloseq.upgrade)
      command1<-sprintf(" Rscript  %swww/script/merge_metafrom_phyloseq_merge_allphyloseq2.R -s %s -t %s -o %s",work_path,target_phy_path,filter_meta_path,dir_phyloseq.upgrade)
      # command1<-sprintf(" %s/envs/pircust2/bin/Rscript  %swww/script/merge_metafrom_phyloseq_merge_allphyloseq2.R -s %s -t %s -o %s",conda_site,work_path,target_phy_path,filter_meta_path,dir_phyloseq.upgrade)
      print(command1)
      system( command1 ) 
      system("mv " %+% dir_phyloseq.upgrade %+%'/*_update.rds ' %+% dir_phyloseq.upgrade %+%'/FALSE_update.rds')
      
    }
    
  })
  
  Variable.SelectP4_1_tax <- reactive({ input$SelectP4_1_tax })
  ########################################################################################reactive
  
  ########################################################################################observeEvent  
  observeEvent(input$SelectP1_run1_uploading, { #上传文件
    # 加上进度条显示
    withProgress(message='Program running:', detail="Converting gene IDs", {
      ######处理数据
      input.upload.otu.datapath<-input$upload.otu$datapath
      input.tax.datapath<-input$upload.tax$datapath
      input.upload.meta.datapath<<-input$upload.meta$datapath
      input.upload.phyloseq.datapath<<-input$upload.phyloseq$datapath
      input.upload.sequences.datapath<<-input$upload.sequences$datapath
      # input.upload.otu.datapath<-"/home/zhiyu/data/metabonomics/Sundeqiao_20200609_hong_analysis/SDQ-F-hong-20200609/data/otu.csv"
      # input.tax.datapath<-"/home/zhiyu/data/metabonomics/Sundeqiao_20200609_hong_analysis/SDQ-F-hong-20200609/data/tax.csv"
      # input.upload.meta.datapath<-"/home/zhiyu/data/metabonomics/Sundeqiao_20200609_hong_analysis/SDQ-F-hong-20200609/data/SDQ-F-TILs-Hepa1-6-in-situ-group.csv"
      # browser()
      if (!is.null(input.upload.phyloseq.datapath)){
        print(input.upload.phyloseq.datapath)
        phyloseq_<-readRDS(  input.upload.phyloseq.datapath   )
        sample_names(phyloseq_) <-  gsub("-",".",sample_names(phyloseq_))# maybe changed or comment
        if (!is.null(input.upload.meta.datapath)){
          df.upload.meta <- read.csv(input.upload.meta.datapath,header = TRUE,stringsAsFactors = F,row.names = 1)
          df.upload.meta$SampleID<-rownames(df.upload.meta)
          sample_data(phyloseq_)<-df.upload.meta
        }
        print(phyloseq_)
      }else{
        print(input.upload.otu.datapath)
        df.upload.otu <- read.csv(input.upload.otu.datapath,header = TRUE,stringsAsFactors = F,row.names = 1)
        print(head(df.upload.otu))
        
        print(input.tax.datapath)
        df.upload.tax <- read.csv(input.tax.datapath,header = TRUE,stringsAsFactors = F,row.names = 1)
        print(df.upload.tax)
        
        print(input.upload.meta.datapath)
        df.upload.meta <- read.csv(input.upload.meta.datapath,header = TRUE,stringsAsFactors = F,row.names = 1)
        df.upload.meta$SampleID<-rownames(df.upload.meta)
        print(df.upload.meta)
        
        
        phyloseq_<-phyloseq(otu_table(df.upload.otu,taxa_are_rows = T),tax_table( df.upload.tax%>% as.matrix()) )
        sample_data(phyloseq_)<-df.upload.meta
        
        # browser()
        if (!is.null(input.upload.sequences.datapath)){
          print(input.upload.sequences.datapath)
          library(phyloseq)
          library(Biostrings)
          rep.seqs <- Biostrings::readDNAStringSet(input.upload.sequences.datapath, format = "fasta")
          phyloseq_<-phyloseq(otu_table(df.upload.otu,taxa_are_rows = T),tax_table( df.upload.tax%>% as.matrix()),sample_data(df.upload.meta),refseq(rep.seqs) )
        }
        
        print(phyloseq_)
        # work_path<-"/home/zhiyu/data/luo/16s_shiny/test/"
      }
      incProgress(0.5, detail = paste("half"))
      phyloseq_<<-phyloseq_
      write.csv( df.upload.meta,pre_analysis_meta_path,quote=F  )
      saveRDS(phyloseq_, pre_analysis_phyloseq_path )
      
      
      # file.meta<-read.csv(pre_analysis_meta_path,stringsAsFactors = FALSE)
      file.meta<-df.upload.meta
      rt$file_meta<-file.meta
      
      #Select column of meta
      colname_<-colnames(file.meta) 
      updateSelectizeInput(session, "SelectP1_2",label = paste("Select input label(for example: host_disease)", length( colname_ )),
                           choices = colname_ , options = list(maxItems = 1L))# page1
      # updateSelectizeInput(session, "SelectP1_2",label = paste("Select One Column (example: host_disease)", length( colname_ )),
      #                    choices = colname_  , options = list(maxItems = 1L) )# page3
      rt$file_meta<-file.meta
      ######处理数据
      incProgress(1, detail = paste("Done"))	  
    })
  })
  
  
  observeEvent(input$SelectP1_run2, { #选择数据库中Project的项目ID
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.2, detail = paste("A thousand-li journey is started by taking the first step."))
      ######处理数据
      file_name <- input$SelectP1_1
      print('shifouchenggong')
      print(file_name)
      #file_name <- "SRP189432.meta.csv"
      if (is.null(file_name)){ file_name <- character(0) }
      #file_name<- "SRP126121.meta.csv"
      if ( file_name %in% list.meta.path.basename ){
        print(list.meta.path[file_name])
        file.meta.name1<-list.meta.path[[file_name]]
        print("need know file.meta:")
        print(file.meta.name1)
        system( paste0("cp ",file.meta.name1,' ',pre_analysis_meta_path)  )
        system( paste0("cp ",pre_analysis_merge_all.phyloseq_path, ' ',pre_analysis_phyloseq_path)  )
        file.meta<-read.csv(pre_analysis_meta_path,stringsAsFactors = FALSE)
      }
      #Select column of meta
      colname_<-colnames(file.meta) 
      updateSelectizeInput(session, "SelectP1_2",label = paste("Select input label(for example: host_disease)", length( colname_ )),
                           choices = colname_ , options = list(maxItems = 1L))# page1
      # updateSelectizeInput(session, "SelectP1_2",label = paste("Select One Column (example: host_disease)", length( colname_ )),
      #                    choices = colname_  , options = list(maxItems = 1L) )# page3
      rt$file_meta<-file.meta
      ######处理数据
      incProgress(0.8, detail = paste("A thousand-li journey is started by taking the first step."))
    })
  })
  
  observeEvent(input$SelectP2_run1, { #run 选择的参数：group的列名及列的内容
    # 加上进度条显示
    withProgress(message='Program running:', detail="Converting gene IDs", {
      incProgress(0.5, detail = paste("half"))
      ######处理数据
      #system("rm /home/zhiyu/data/paper/all.meta/filter/*")
      system(paste0("rm ",work_path,"www/all.meta/filter/*" ))
      
      filter.meta() 
      ######处理数据
      incProgress(1, detail = paste("Done"))	  
    })
  })
  
  observeEvent(input$SelectP2_mergephyloseq, { #产生phyloseq
    # 加上进度条显示
    withProgress(message='Program running:', detail="Half of the people who have embarked on a one hundred mile journey may fall by the wayside.", {
      incProgress(0.5, detail = paste("Your time is limited, so don't waste it living someone else's life.…Don't let the noise of others' opinions drown out your own inner voice"))
      ######处理数据
      file.phyloseq() 
      # system(paste0("rm  -rf ",work_path,"out.result" ))
      # dir.create(paste0(work_path,"out.result" ),recursive = T)
      ######处理数据
      incProgress(1, detail = paste("Your time is limited, so don't waste it living someone else's life.…Don't let the noise of others' opinions drown out your own inner voice"))	  
    })
  })
  
  
  
  
  # output$summary <- renderPrint({   summary( rt$file_meta  )  })
  output$meta <- renderDT(  
    rt$file_meta,options = list(columnDefs = list(list(className = 'dt-center', targets = 0: ncol(rt$file_meta))) ) )
  
  ########################################################################################observeEvent  
  
  
  ########################################################################################output
  ###output
  #source("/home/zhiyu/data/script/Shiny/16s/plsda_adjust_data3.r")
  source(paste0(work_path,"www/script/plsda_adjust_data3.r" ))
  
  color_<-c("black","red","green","blue","maroon","olivedrab2")
  #tax_level<-"Genus"
  
  #output$plot <- renderImage( {     list(src = filename0,width=600,height=600) }, deleteFile = FALSE   )
  
  output$plot <- renderPlot( {     #list(src = filename0,width=600,height=600) 
    file.all_meta_id_num()
    
  })
  
  
  output$P1_table <- DT::renderDataTable({     #rt$filter_meta      
    print("kaishi  .....")
    #system("rm /home/zhiyu/data/paper/all.meta/filter/*")
    #filter.meta() 
    rt$filter_meta 
    
  })
  output$AD_plot_table3 <- DT::renderDataTable({     #rt$filter_meta      
    #filter.meta() 
    filter_meta<-subset(rt$filter_meta,grepl(paste(input$SelectP1_3, collapse = "|"),rt$filter_meta[,input$SelectP1_2]))  
    filter_meta
  })
  
  
# alpha

# alpha -------------------------------------------------------------------
  print("############################################ alpha  ###############################33########")
  output$AD_plot_a <- renderPlot({  
    withProgress(message='Program running:', detail="", {
      incProgress(0.2, detail = paste("To become the richest man in the cemetery doesn't matter to me... Said the night to go to bed we have done a great thing... It is important for me."))
      tax_level<-input$SelectP3_2_tax
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      
      percent_<-input$SelectP3_percent
      
      # otu_min_count<-input$otu_min_count#5
      # methods_trans<-input$methods_trans#"compositional"
      # prevalence_values<-input$prevalence_values#0.25
      # browser()
      # phyloseq<-phy_flitered(phyloseq,percent_)
      phyloseq<-phy_flitered(phyloseq, min_count=0,methods_trans=NULL ,prevalence_values=NULL)
      # phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
      
      #phyloseq<-phy_flitered(phyloseq)
      
      print('zheli')
      print(input$SelectP1_2)
      print(input$SelectP3_2_tax)
      incProgress(0.9, detail = paste("To become the richest man in the cemetery doesn't matter to me... Said the night to go to bed we have done a great thing... It is important for me."))
      alpha.diversity(tax_level,phyloseq,input$SelectP1_2)
      
    })
  })
  print("############################################End alpha  ######################################")

# beta --------------------------------------------------------------------
  print("############################################ beta  ##########################################")
  output$AD_plot_b <- renderPlot({  
    tax_level<-input$SelectP3_2_tax
    #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
    phyloseq<-readRDS( phyloseq_path   )
    
    percent_<-input$SelectP3_percent
    # phyloseq<-phy_flitered(phyloseq,percent_)
    phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
    
    #phyloseq<-phy_flitered(phyloseq)
    beta.diversity(phyloseq,input$SelectP1_2,input$beta_method)
    
  })
  print("############################################End beta  #######################################")

# Compositional -----------------------------------------------------------
  print("############################################ compositional  #################################")
  output$AD_plot_compositional <- renderPlot({  
    tax_level<-input$SelectP3_2_tax
    #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
    phyloseq<-readRDS( phyloseq_path   )
    
    percent_<-input$SelectP3_percent
    # phyloseq<-phy_flitered(phyloseq,percent_)
    phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
    
    #phyloseq<-phy_flitered(phyloseq)
    
    plot_compositional(tax_level,phyloseq,input$SelectP1_2)
    
  })
  
  # output$AD_plot_compositional2 <-renderImage({ 
  #   ps_list() #等待上个函数执行后再展示图片
  #   filename <- out.dir %+% "/rarecurve/rarecurve_pv" %+% input$prevalence_values %+% ".png" 
  #   filename =paste0("out.result/composition/","composition_.",tax_level,".png")
  #   width  <- session$clientData$output_image_width
  #   height <- session$clientData$output_image_height
  #   list(src = filename,contentType = "image/png",width = width/2,height = height/2)
  # })
  
  
  print("############################################End compositional  ##############################")

# Deseq2 ------------------------------------------------------------------
  print("############################################deseq2###########################################")
  AD_plot_deseq2_top_nrow<<-1
  output$AD_plot_deseq2 <- renderImage({ 
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.5, detail = paste("The most important thing is to have the courage to follow your heart and intuition, they may already know what you want to be a man"))
      ######处理数据
      #tax_level<-input$SelectP3_2_tax
      tax_level<-input$SelectP4_1_tax
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      percent_<-input$SelectP3_percent
      
      
      # phyloseq<-phy_flitered(phyloseq,percent_)
      phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=NULL,prevalence_values=NULL)
      
      #phyloseq<-phy_flitered(phyloseq)
      phyloseq<-convert_seq_to_asv(phyloseq)
      #phyloseq<-ordered_phy(phyloseq,"Class2")
      phyloseq<-ordered_phy(phyloseq,input$SelectP1_2)
      # p1<-plot_deseq2(tax_level,phyloseq,input$SelectP1_2)
      input$deseq2_foldchange_num
      input$deseq2_padj_num
      
      rt$AD_plot_deseq2<-plot_deseq2(tax_level,phyloseq,
                                     target_2rd_column=input$SelectP1_2,
                                     deseq2_foldchange=input$deseq2_foldchange_num,
                                     deseq2_pvalue=input$deseq2_pvalue,
                                     deseq2_padj=input$deseq2_padj_num,
                                     Deseq2_height=input$Deseq2_height,
                                     Deseq2_width=input$Deseq2_width
                                     )
      # trace(plot_deseq2,edit = T)
      # rt$AD_plot_deseq2<-p1
      filename <- paste0("out.result/deseq2/","deq2_",tax_level,".png" )
      ######处理数据
      # Get width and height of image output
      width  <- session$clientData$output_image_width
      height <- session$clientData$output_image_height
      
      incProgress(1, detail = paste("The most important thing is to have the courage to follow your heart and intuition, they may already know what you want to be a man"))	  
      list(src = filename,contentType = "image/png",width = width/2,height = height/2)
       })
    # rt$AD_plot_deseq2
   
    # },width=1000,height = 300*(AD_plot_deseq2_top_nrow)+100)
  })
  # },width=1000,height = myHeightAlgorithm() ) #
  # output$AD_plot_deseq2_picture <-renderPlot( { 
  #   tax_level<-input$SelectP4_1_tax
  #   filename <- paste0("out.result/deseq2/","deq2_",tax_level,".png" )
  #   readpng(filename)
  # },width=1000,height = 800)
  print("############################################End deseq2#######################################")  

# lefse -------------------------------------------------------------------
  print("############################################ lefse ##########################################") 
  output$AD_plot_lefse_picture <-renderImage({ 
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.5, detail = paste("All for one, one for all"))
      ######处理数据
      #tax_level<-input$SelectP3_2_tax
      tax_level<-input$SelectP4_1_tax
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      percent_<-input$SelectP3_percent
      
      
      # phyloseq<-phy_flitered(phyloseq,percent_)
      phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=NULL,prevalence_values=NULL)
      
      #phyloseq<-phy_flitered(phyloseq)
      phyloseq<-convert_seq_to_asv(phyloseq)
      #phyloseq<-ordered_phy(phyloseq,"Class2")
      phyloseq<-ordered_phy(phyloseq,input$SelectP1_2)
      
      # browser()
      group<-input$SelectP1_2
      phyloseq::sample_data(phyloseq)$Group <- meta(phyloseq)[,group]
      
      group1=unique(meta(phyloseq)$Group)[1]%>%as.character()
      group2=unique(meta(phyloseq)$Group)[2]%>%as.character()
      ref_=ifelse(length(group2)>length(group1),group2,group1)
      outcome_add_info=""
      # tax_level<-input$SelectP3_2_tax
      tax_level<-input$SelectP4_1_tax
      lefse_height<-input$lefse_height
      lefse_width<-input$lefse_width
      rt$AD_plot_lefse<-r_lefse_func(group1,group2,phyloseq,ref_,group,"",level=tax_level,lda_cutoff_values=input$lefse_lda_value,lefse_height=lefse_height,lefse_width=lefse_width)
      
      filename <- paste0("out.result/lefse/lefse_markder_",tax_level,"_barplot.png" )
      ######处理数据
      incProgress(1, detail = paste("All for one, one for all"))	 
      list(src = filename ,contentType = "image/png") 
    
    })

  })
  print("############################################End lefse#######################################") 

# Metacoder ---------------------------------------------------------------
  print("############################################Metacoder#######################################")   
  output$AD_plot_Metacoder <- renderPlot({         
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.5, detail = paste("The most important decisions in life is not what you do, but what you do"))
      ######处理数据
      tax_level<-input$SelectP3_2_tax
      percent_<-input$SelectP3_percent
      #percent_<-0.005
      
      # system(paste0("rm ",work_path,"www/all.meta/filter/*" ))
      
      
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      # phyloseq<-phy_flitered(phyloseq,percent_)
      # phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
      
      plot_MC<-plot_Metacoder(phyloseq,input$SelectP1_2,min_count=input$otu_min_count )
      rt$plot_all_Metacoder<-plot_MC
      ######处理数据
      incProgress(1, detail = paste("The most important decisions in life is not what you do, but what you do"))	  
    })
    rt$plot_all_Metacoder
  })
  print("############################################End Metacoder###################################")   

# PCA ---------------------------------------------------------------------
  print("############################################ PCA ###########################################")  
  output$AD_plot_pca_bar_make <- renderText({  
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.5, detail = paste("One boy is a boy, two boys half a boy, three boys no boy"))
      ######处理数据
      #tax_level<-Variable.SelectP4_1_tax
      tax_level<-input$SelectP4_1_tax
      print(paste0("PCA ",tax_level))
      #tax_level<-input$SelectP3_2_tax
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      
      percent_<-input$SelectP3_percent
      # phyloseq<-phy_flitered(phyloseq,percent_)
      print(input$methods_trans)
      # browser()
      phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
      # phyloseq@otu_table
      #phyloseq<-phy_flitered(phyloseq)
      #setwd("/home/zhiyu/data/luo/16s_shiny/test")
      #pca_("Genus",phyloseq,"diesease")
      rt$AD_plot_pca_bar<-AD_plot_pca_bar(tax_level,phyloseq,input$SelectP1_2)
      ######处理数据
      incProgress(1, detail = paste("One boy is a boy, two boys half a boy, three boys no boy"))	  
    })
    # rt$AD_plot_pca_bar
    
  })
  output$AD_plot_pca_make <- renderText({  
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.5, detail = paste("One boy is a boy, two boys half a boy, three boys no boy"))
      ######处理数据
      #tax_level<-Variable.SelectP4_1_tax
      tax_level<-input$SelectP4_1_tax
      print(paste0("PCA ",tax_level))
      #tax_level<-input$SelectP3_2_tax
      # phyloseq_path<-"/home/lzy/data/luo/16s_shiny/test/www/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      
      percent_<-input$SelectP3_percent
      # phyloseq<-phy_flitered(phyloseq,percent_)
      phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
      
      
      
      rt$plot_pca<-pca_(tax_level,phyloseq,input$SelectP1_2,input$n_componet)
      ######处理数据
      incProgress(1, detail = paste("One boy is a boy, two boys half a boy, three boys no boy"))	  
    })
    # rt$plot_pca
    
  })
  print("############################################End PCA ########################################")  

# PLSDA -------------------------------------------------------------------
  print("############################################ PLSDA  ########################################")  
  output$AD_plot_pls.da <- renderText({  
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.5, detail = paste("Pleasure comes through toil"))
      ######处理数据
      #tax_level<-input$SelectP3_2_tax
      
      tax_level<-input$SelectP4_1_tax
      print(paste0("PLS-DA ",tax_level))
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      
      
      percent_<-input$SelectP3_percent
      # phyloseq<-phy_flitered(phyloseq,percent_)
      phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
      
      #phyloseq<-phy_flitered(phyloseq)
      rt$plot_pls.da<-pls.da_(tax_level,phyloseq,input$SelectP1_2)
      #pls.da_('Genus',phyloseq,'Class')
      ######处理数据
      incProgress(1, detail = paste("Pleasure comes through toil"))	  
    })
  })
  #purrr::map( as.list(diff_group), function(i){ print(i);print(class(i));subset_samples(phyloseq, gp == i[1] )  } ) # 
  #a<-subset_samples(phyloseq, gp== "pre_Tg")
  print("############################################End PLSDA  #####################################")  

# splda -------------------------------------------------------------------
  print("############################################ splda  ########################################")  
  output$AD_plot_spls.da <- renderText({  
    # 加上进度条显示
    withProgress(message='Program running:', detail="No matter how high the mountain is, one can always ascend to its top.", {
      print("########################Start splsda##########################")
      ######处理数据
      #tax_level<-input$SelectP3_2_tax
      tax_level<-input$SelectP4_1_tax
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      phyloseq<-readRDS( phyloseq_path   )
      
      percent_<-input$SelectP3_percent
      #percent_<-0.005
      # phyloseq<-phy_flitered(phyloseq,percent_)
      phyloseq<-phy_flitered(phyloseq, min_count=input$otu_min_count ,methods_trans=input$methods_trans ,prevalence_values=input$prevalence_values)
      
      
      phyloseq<-convert_seq_to_asv(phyloseq)
      #phyloseq<-ordered_phy(phyloseq,"disease")
      
      phyloseq<-ordered_phy(phyloseq,input$SelectP1_2)
      incProgress(0.25, detail = paste("No matter how high the mountain is, one can always ascend to its top"))
      #setwd("/home/zhiyu/data/luo/16s_shiny/test")
      #spls.da_new("Genus",phyloseq,"disease")
      # browser()
      n_componet=input$n_componet
      
      spls.da_new_step1_values<-spls.da_new_step1(tax_level,phyloseq,input$SelectP1_2,n_componet)# return gvri.splsda Y
      incProgress(0.5, detail = paste("No matter how high the mountain is, one can always ascend to its top"))
      # spls.da_new_step2(spls.da_new_step1_values$f1,spls.da_new_step1_values$f2,n_componet)# heatmap
      # incProgress(1, detail = paste("No matter how high the mountain is, one can always ascend to its top"))
      ######处理数据
      print("########################End splsda##########################")
      incProgress(1, detail = paste("Done"))	  
    })
    
  })
  print("############################################End splda  #####################################")  
  
  phyloseq_path<-paste0(work_path,"www/phyloseq.upgrade/FALSE_update.rds" )
  output$AD_plot_pca_bar <-renderPlot( { 
    tax_level<-input$SelectP4_1_tax
    filename <- paste0("out.result/pca/","pca_bar_",tax_level,".png" )
    readpng(filename)
  })
  output$AD_plot_pca <-renderPlot( { 
    tax_level<-input$SelectP4_1_tax
    filename <- paste0("out.result/pca/","pca_",tax_level,".png" )
    readpng(filename)
  })
  output$AD_plot_plda.tu1 <-renderPlot( { 
    filename <- paste0(work_path,"out.result/plsda/","pls.da_1-2.png")
    readpng(filename)
  })
  output$AD_plot_plda.tu2 <-renderPlot( { 
    filename <- paste0(work_path,"out.result/plsda/","pls.da_1-3.png")
    readpng(filename)
  })
  
  output$AD_plot_splda.tu1 <-renderPlot( { 
    n_componet=input$n_componet
    filename <- paste0(work_path,"out.result/splsda/","SPLSDA_1_",n_componet,".png.png")
    readpng(filename)
  })
  output$AD_plot_splda.tu_comp1 <-renderPlot( { 
    #tax_level="Genus"
    tax_level<-input$SelectP4_1_tax
    n_componet=input$n_componet
    filename <- paste0("out.result/splsda/splda_",tax_level,"_plot_n",n_componet,".png")
    readpng(filename)
  })
  output$AD_plot_splda.heatmap1 <-renderPlot( { 
    #tax_level="Genus"
    tax_level<-input$SelectP4_1_tax
    n_componet=input$n_componet
    filename <- paste0("out.result/splsda/","Heatmap_comp1_comp",n_componet,".png")
    readpng(filename)
  })
  output$AD_plot_splda.heatmap3 <-renderPlot( { 
    #tax_level="Genus"
    tax_level<-input$SelectP4_1_tax
    filename <- paste0("out.result/splsda/","Heatmap_all.png")
    readpng(filename)
  })


  print("############################################ NetCoMi  ######################################")  
  output$NetCoMi_make <-renderText( { 
    ##NetCoMi
    # observeEvent(input$runP5_netcomi, { 
    print('#####################Start NetCoMi###########################')
    print(paste0('input$runP5_netcomi:',input$runP5_netcomi))
    # 加上进度条显示
    withProgress(message='Program running:', detail="There are people who will appreciate what I have done but there are also people who will criticize me,ultimately,history will have the final say.", {
      
      ######处理数据
      # browser()
      tax_level<-input$SelectP3_2_tax
      # tax_level<-input$SelectP4_1_tax
      # tax_level<-"Genus"
      #phyloseq_path<-"/home/lzy/data/luo/16s_shiny/test/www/phyloseq.upgrade/FALSE_update.rds"
      phyloseq_temp<-readRDS(  phyloseq_path   )
      # phyloseq_temp <- phyloseq_temp %>% aggregate_taxa(level = tax_level)
      percent_<-input$SelectP3_percent
      #percent_<-0.005
      
      # phyloseq_temp<-phy_flitered(phyloseq_temp,percent_)
      phyloseq_temp<-convert_seq_to_asv(phyloseq_temp)
      
      #phyloseq<-ordered_phy(phyloseq,"disease")
      
      phyloseq_temp<-ordered_phy(phyloseq_temp,input$SelectP1_2)
      
      #setwd("/home/zhiyu/data/luo/16s_shiny/test")
      
      #spls.da_new("Genus",phyloseq,"disease")
      # target_column2<-"host_disease"
      # two_group_classfication2<-c("Healthy Control", "Ulcerative Colitis")
      # two_group_classfication2<-c("pre", "post")
      target_column2<-input$SelectP1_2
      two_group_classfication2<-input$SelectP1_3
      print(phyloseq_temp)
      # print(filter.meta)
      # browser()
      incProgress(0.5, detail = paste("There are people who will appreciate what I have done but there are also people who will criticize me,ultimately,history will have the final say."))
      
      # plan(multisession)
      # availableCores()
      # future::future({        NetCoMi_fc(phyloseq_temp,target_column2,two_group_classfication2,tax_level)      })
      NetCoMi_fc(phyloseq_temp,target_column2,two_group_classfication2,tax_level)
      print('#####################End NetCoMi###########################')
      ######处理数据
      incProgress(1, detail = paste("There are people who will appreciate what I have done but there are also people who will criticize me,ultimately,history will have the final say."))	  
    })
    
    
    # tax_level<-input$SelectP3_2_tax
    # print("The picture of NetCoMi will be showed.")
    # filename <- paste0("out.result/","NetCoMi_",tax_level,".png") 
    # readpng(filename)
  })
  output$NetCoMi_png <-renderPlot( { 
    #tax_level="Genus"
    tax_level<-input$SelectP3_2_tax
    filename <- paste0("out.result/NetCoMi/","NetCoMi_",tax_level,".png") 
    readpng(filename)
  })
  print("############################################End NetCoMi  ###################################")  
  
  print("############################################ pircust  ######################################")  
  observeEvent(input$runP5_1, { 
    print('xiamian............')
    # print(paste0('input$SelectP5_1:',input$SelectP5_1))
    # plan(multisession)
    # 加上进度条显示
    withProgress(message='Program running:', detail=paste("There are people who will appreciate what I have done but there are also people who will criticize me,ultimately,history will have the final say."), {
      incProgress(0.2, detail = paste("There are people who will appreciate what I have done but there are also people who will criticize me,ultimately,history will have the final say.")  )	 
      ######处理数据
      # browser()
      print(paste0('gongzuomulu:',getwd() ))
      ##run pircust
      subDir<-'out.result/Pircust'
      # if (file.exists(subDir)){file.remove(subDir,recursive=T)}
      if ( file.exists(subDir) ){system("rm -rf out.result/Pircust")} # you can remove the last record.
      system('tree out.result/Pircust') 
      system('mkdir -p out.result/Pircust/data')
      #phyloseq_path<-"/home/zhiyu/data/paper/phyloseq.upgrade/FALSE_update.rds"
      # browser()
      # phyloseq_path <- paste0(all.phyloseq.dir,'/',input$SelectP5_1)# you can redefine a new phyloseq_path
      system(paste0('cp ' ,phyloseq_path,' out.result/Pircust/data'))
      system('tree out.result/Pircust') 
      # pircust2.analysis<-sprintf("Rscript /home/zhiyu/data/script/Pircust/pircut.2.3.0_b.2.9_METACYC_KEGG.r -f /home/zhiyu/data/metagenome/analysis_project/liudiao_H22_liuyuan_151samples_20200324/tt -s /home/zhiyu/data/metagenome/analysis_project/liudiao_H22_liuyuan_151samples_20200324/tt/data/FALSE2020_03_20_5_update.rds2020_03_24_2_update.rds")
      # pircust2.analysis<-sprintf("Rscript /home/zhiyu/data/script/Pircust/pircut.2.3.0_b.2.9_METACYC_KEGG.r -f %s/out.result/Pircust -s %s",getwd(),phyloseq_path)
      # browser()
      # pircust2.analysis<-sprintf("%s/envs/pircust2/bin/Rscript %swww/script/pircut.2.3.0_b.2.9_METACYC_KEGG.r -f %s/out.result/Pircust -s %s",conda_site,work_path,getwd(),phyloseq_path)
      pircust2.analysis<-sprintf("Rscript www/script/pircut.2.3.0_b.2.9_METACYC_KEGG.r -f %s/out.result/Pircust -s %s",work_path,phyloseq_path)
      print(pircust2.analysis)
      print(system('pwd'))
      system(pircust2.analysis)
      system('cd ' %+% getwd() %+% '/out.result/Pircust \n ' %+% 'bash -x src/picrust2.bash')
      
      
      # availableCores()
      # print(availableCores())
      incProgress(0.3, detail =paste("There are people who will appreciate what I have done but there are also people who will criticize me,ultimately,history will have the final say.") )
      # browser()
      # output$pircust_value <- renderPrint( yy )
      print('tiaochulaile')
      ######处理数据
      incProgress(0.5, detail = paste("There are people who will appreciate what I have done but there are also people who will criticize me,ultimately,history will have the final say.")  )	  
    })
    print('jindutiao tichu')
  })
  observeEvent(input$runP5_2, { 
    print('xiamian............')
    # browser()
    # 加上进度条显示
    withProgress(message='Program running:', detail="", {
      incProgress(0.5, detail = paste("half"))
      ######处理数据
      print(paste0('gongzuomulu:',getwd() ))
      ##plot pircust
      current_dir_path<-getwd()
      work_path<- paste0(getwd(),'/out.result/Pircust')
      
      meta_path<-filter_meta_path
      # meta_path<-paste0(path,"/",input$SelectP5_2)
      Pircust_df_name<- paste0(work_path,"/output/pircust2/KEGG/picrust2_out_pipeline/KEGG_pathways_unstrat/path_abun_unstrat_descrip.tsv")
      # print('00000000000')
      # print(Pircust_df_name)
      # target_column1<-input$SelectP5_3
      # two_group_classfication1 <- ifelse (input$SelectP5_4 == "All",  "*", paste0('^',input$SelectP5_4 ) )
      # target_column2<-input$SelectP5_5
      # two_group_classfication2<-input$SelectP5_6#c('Ulcerative Colitis', 'Healthy Control')
      
      # work_path<- "/home/zhiyu/data/luo/16s_shiny/test/out.result/Pircust/"
      # meta_path<-"/home/lzy/data/luo/16s_shiny/test/www/all.meta/raw/SRP126121.meta.csv"
      # target_column1<-"disease"
      # two_group_classfication1<-"*"
      # target_column2<-"disease"
      # two_group_classfication2<-c('Ulcerative Colitis', 'Healthy Control')
      
      target_column1<-input$SelectP1_2
      two_group_classfication1<-input$SelectP1_3
      target_column2<-input$SelectP1_2
      two_group_classfication2<-input$SelectP1_3
      
      
      # browser()
      c1<-pircust_wilcox(work_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2,p_value_threshold=input$p_value_threshold,analysis_method=input$analysis_method)
      # c2<-pircust_aldex2(work_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2)
      c3<-pircust_deseq2(work_path,meta_path,Pircust_df_name,target_column1,two_group_classfication1,target_column2,two_group_classfication2)
      #browser()
      #work_path<- "/home/zhiyu/data/luo/16s_shiny/test/out.result/Pircust/"
      #filter_meta_path<-paste0(work_path,"www/all.meta/filter/filter_meta.tsv")
      print(work_path)
      # browser()
      Pircust_df_name_mei<-paste0(work_path,"/output/pircust2/KEGG/picrust2_out_pipeline/EC_metagenome_out_unstrat/pred_metagenome_unstrat.tsv.gz")
      
      # browser()
      c4<-pircust_deseq2_mei(current_dir_path,meta_path,Pircust_df_name_mei,target_column1,two_group_classfication1,target_column2,two_group_classfication2)
      ######处理数据
      incProgress(1, detail = paste("Done"))	  
    })
    rt$AD_plot_Pircust2_percent<-c1
    # rt$AD_plot_Pircust2_aldex<-c2
    rt$AD_plot_Pircust2_deseq2<-c3
    
  })
  output$AD_plot_Pircust2_percent <- renderPlot({  rt$AD_plot_Pircust2_percent   })
  # output$AD_plot_Pircust2_aldex <- renderPlot({  rt$AD_plot_Pircust2_aldex   })
  output$AD_plot_Pircust2_deseq2 <- renderPlot({  rt$AD_plot_Pircust2_deseq2   })
  print("############################################End pircust  ###################################")
  
  print("############################################ Amon  #########################################")
  observeEvent(input$runP5_amon, { 
    print('xiamian............')
    print(paste0('input$runP5_amon:',input$runP5_amon))
    # 加上进度条显示
    withProgress(message='Program running:', detail=paste("As a chinese poem reads we have no fear of the clouds that may block our sights as we are already at the top of the height."), {
      
      ######处理数据
      print(paste0('gongzuomulu:',getwd() ))
      ##run amon
      #browser()
      subDir<-'out.result/Amon'
      if (file.exists(subDir)){file.remove(subDir,recursive=T)}
      #if (file.exists(subDir)){system("rm -rf out.result/Pircust")}
      system('mkdir -p out.result/Amon/data')
      system('mkdir -p out.result/Amon/src')
      # system('mkdir -p out.result/Amon/src')
      incProgress(0.2, detail = paste("As a chinese poem reads we have no fear of the clouds that may block our sights as we are already at the top of the height."))
      conda_site<-"/home/zhiyu/data/software/Anaconda3/bin"
      amon_str<-sprintf(" activate_local=%s \n
              source $activate_local/activate AMON2
              conda info --env
              rm -rf out.result/Amon/output/AMON_result8
              amon.py -o out.result/Amon/output/AMON_result8 -i out.result/Pircust/output/bar.sig_deseq2_mei.txt \n",conda_site)
      writeLines(amon_str, 'out.result/Amon/src/amon.bash'  )#run the bash in linux;
      amon_str<-paste0("/bin/bash -x ", "out.result/Amon/src/amon.bash")
      print(amon_str)
      # plan(multisession)
      # availableCores()
      # future::future({         system(amon_str)      })
      system(amon_str)
      incProgress(0.8, detail = paste("As a chinese poem reads we have no fear of the clouds that may block our sights as we are already at the top of the height."))
      
      filename <- "out.result/Amon/output/AMON_result8/origin_table.tsv"
      rt$amon_table_daixie<-readtable(filename)
      #browser()
      ######处理数据
    })
  })
  output$amon_table_daixie <- DT::renderDataTable({ rt$amon_table_daixie })
  print("############################################End Amon  ######################################")
  
  
  
  ########################################################################################output

  ########################################################################################send the zip of result
  #system( "python /home/zhiyu/data/script/16S/16s_pipline/python_16s_basic_bookdown.py /home/zhiyu/data/paper/analysis/ncbi/IBS/SRP189432  /home/zhiyu/data/paper/data/IBS/ncbi/SRP189432/rawdata  ll FALSE FALSE" )
  #system("cp -r /home/zhiyu/data/paper/analysis/ncbi/IBS/SRP189432/analysis/_book /home/zhiyu/data/luo/16s_shiny/test/www/")
  #

# DownloadData ------------------------------------------------------------
  output$downloadData <- downloadHandler(
    filename = 'files.zip',
    content = function(fname) {
      withProgress(message='Program running:', detail="Still water run deep", {
        incProgress(0.2, detail = paste("Still water run deep"))  
        #tmpdir <- tempdir()
        # browser()
        wk_path1<-getwd()
        wk_path2<-getwd()
        tempdir_ <- paste0(wk_path1,"/out.result")
        print(getwd())
        setwd(tempdir_)
        print(tempdir_)
        
        
        #setwd(tempdir())
        #print(tempdir())
        fs<- dir(tempdir_)
        #fs<- dir(tempdir())
        print(fs)
        print (fs)
        # browser()
        # zip(zipfile=fname, files=fs)
        zip::zip(zipfile=fname, files=fs)
        setwd(wk_path2)
        incProgress(0.3, detail = paste("Still water run deep"))   
      })
    },
    contentType = "application/zip"
  )
  
  

# Downloadexample ---------------------------------------------------------
  output$downloadexample <- downloadHandler(
    filename = 'files_example.zip',
    content = function(fname) {
      file.copy("www/downloadexample.zip", fname)
    }
  )
  
  
  

# Refreshing --------------------------------------------------------------

  ########################################################################################send the zip of result
  print("Initializing")
  
  observeEvent(input$switchtab,{
    aggg_result = -1
    if(aggg_result == -1)
    {
      session$reload()
      return()
      print("session reload not working")
    }
    
    print("Code running this line")
    
    output$code_ran <- renderText("code Ran this line without refreshing")
    
  })
  

# SessionInfo -------------------------------------------------------------
  output$sessionInfo <- renderPrint({
    capture.output(sessionInfo())
    
  })
  
  # session$allowReconnect("force")
  # options(shiny.launch.browser=FALSE)
  # tt <- sessionInfo() %>% as.character()
}



