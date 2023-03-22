library(shiny)
library(shinythemes)
library(dplyr) 
library(DT)

# filename0 <- normalizePath('/home/zhiyu/data/luo/16s_shiny/test/www/zi.png')
# all.phyloseq.path<-"/home/zhiyu/data/paper/all.phyloseq.upgrade/phyloseq.rds"
# path<-"/home/zhiyu/data/paper/all.meta/raw"
# filter_meta_path<-"/home/zhiyu/data/paper/all.meta/filter/filter_meta.tsv"
# all_phyloseq_sigle<-"/home/zhiyu/data/paper/all.phyloseq"

print('weizhi!!!!!!!!!!!!!!!!!')
print( getwd() )
work_path<- paste0(getwd() ,"/")
#work_path<-"/home/zhiyu/data/luo/16s_shiny/test/"

#filename0 <- normalizePath('zi.png')
all.phyloseq.path<-paste0(work_path,"www/all.phyloseq.upgrade/phyloseq.rds")
path<-paste0(work_path,"www/all.meta/raw") 
filter_meta_path<-paste0(work_path,"www/all.meta/filter/filter_meta.tsv")
all_phyloseq_sigle<-paste0(work_path,"www/all.phyloseq")

###meta path
list.meta.path<-list()
list.meta.path[list.files(path, pattern = "meta", full.names = F)]<-list.files(path, pattern = "meta", full.names = T)
list.meta.path.basename<-names(list.meta.path)

list.meta.path.basename2<-as.list(list.meta.path.basename)
names(list.meta.path.basename2)<-list.meta.path.basename2
print("ui :")
print(list.meta.path.basename2)
###meta path
###phyloseq path
list.phyloseq.path<-list()
list.phyloseq.path[list.files(all_phyloseq_sigle, pattern = "phyloseq", full.names = F)]<-list.files(all_phyloseq_sigle, pattern = "phyloseq", full.names = T)
list.phyloseq.path.basename<-names(list.phyloseq.path)
###phyloseq path




# Define UI for app that draws a histogram ----
ui <- fluidPage(
  #theme = shinytheme("united"),
  #shinythemes::themeSelector(),
  #theme = shinytheme("slate"),
  # shinythemes::themeSelector(),
  navbarPage(
    #themeSelector(),
    "Home",
    #theme = shinytheme("united"),
    theme = shinytheme("cerulean"),
    tabPanel('Load Data',
             sidebarLayout(
               
               sidebarPanel(
                 h3(strong("Uploading file:")),
                 fileInput("upload.tax", "Choose Tax File(CSV)",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 fileInput("upload.otu", "Choose Otu File(CSV)",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 fileInput("upload.meta", "Choose Meta File(CSV)",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 fileInput("upload.sequences", "Choose Fasta File(FASTA,Optional)",
                           multiple = FALSE),
                 fileInput("upload.phyloseq", "Choose phyloseq File(Rds,Optional)",
                           multiple = FALSE
                           # ,accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
                 ),
                 # selectizeInput("SelectP1n_2", "Select column of meta",c("SampleID") ,options = list(maxItems = 1L) ),
                 # selectizeInput("SelectP1n_3", "Select content of column in meta",c("Item A", "Item B"),multiple = TRUE,options = list(maxItems = 2L) ),
                 br(),     
                 h6( 'Note: The first column is the row name of the table for otu, tax and meta. \n
                   In addition, if you want to upload phyloseq which hsa otu, tax, meta and fasta. If you can not upload some the fasta file, the pircust2 and amon also are not doing. ' ),
                 h6( 'You also can download a example here. '),
                 tabPanel("downloadexample",downloadButton('downloadexample', 'Download the example')),
                 
                 actionButton("SelectP1_run1_uploading", "Submit"),
                 
                 
                 
                 br(),
                 br(),
                 br(),
                 h3(strong("Database file:")),
                 selectizeInput("SelectP1_1", "Select a project (for example:SRP126121.meta.csv)", list.meta.path.basename2,selected="SRP126121.meta.csv",options = list(maxItems = 1L)  ),
                 # selectizeInput("SelectP1_2", "Select column of meta",c("SampleID") ,options = list(maxItems = 1L) ),
                 # selectizeInput("SelectP1_3", "Select content of column in meta",c("Item A", "Item B"),multiple = TRUE,options = list(maxItems = 2L) ),
                 # br(),
                 actionButton("SelectP1_run2", "Submit"),
                 br(),
                 h6( 'Note: You can choose to upload local data (meta, taxnomy, otu), or upload phyloseq (already merged meta, taxnomy, otu)' ),
                 
                 br()),
               # Main panel for displaying outputs ----
               mainPanel(
                 titlePanel("MicrobiomeHub:Intestinal microbiota Analysis Platform For Diseases"),
                 img(src = "Road_way.png", height = 600, width = 900),
                 br(),
                 # Output: Tabset w/ plot, summary, and table ----
                 tabsetPanel(type = "tabs",
                             tabPanel( "Plot", plotOutput("plot",width = 800,height = '600'),  
                                       #h3("目前所有样本分布数目")
                             )
                             #,tabPanel("Summary", verbatimTextOutput("summary")))
                             ,tabPanel("Meta", DT::dataTableOutput("meta")))
               ))),
    
    tabPanel('Pre-Process', 
             sidebarLayout(
               sidebarPanel(
                 titlePanel("Pre-Process"),
                 selectizeInput("SelectP1_2", "Select column of meta",c("SampleID") ,options = list(maxItems = 1L),selected= "host_disease") ,
                 selectizeInput("SelectP1_3", "Select content of column in meta",c("Item A", "Item B"),multiple = TRUE,options = list(maxItems = 2L) ),
                 br(),
                 actionButton("SelectP2_run1", "Submit"),
                 
                 
                 sliderInput("otu_min_count", "Minimum count of otu",min = 0.0, max = 10, value = 5),
                 selectizeInput("methods_trans", "The methods of transformation",c("compositional",'clr', 'log10', 'log10p', 'hellinger', 'identity', 'Z') ,options = list(maxItems = 1L),selected= "compositional") ,
                 sliderInput("prevalence_values", "Prevalence values",min = 0.0, max = 1, value = 0.25),
                 
                 
                 h3("Create phyloseq"),
                 br(),
                 actionButton("SelectP2_mergephyloseq", "Run phyloseq"),
                 h5("Note:if you chose to uploading file in first section,you can skip the module!"),
                 br()
               ),
               mainPanel(
                 img(src = "16s.png", height = 600, width = 900),
                 br(),
                 br(),
                 tabPanel("Table Filtered", dataTableOutput("P1_table"))  ,
                 #tabsetPanel(type = "tabs",tabPanel("Table2 必选", tableOutput("table2")))
               ))),
    
    tabPanel('Basis analysis', 
             sidebarLayout(
               sidebarPanel(
                 titlePanel("Basis analysis"),
                 #selectizeInput("SelectP3_1", "Select one Column (example: host_disease)",c("SampleID") ,options = list(maxItems = 1L) ),
                 
                 selectizeInput("SelectP3_2_tax", "Select a level of taxonomy" ,c("Phylum","Family", "Order","Class","Genus","OTU"),selected="Genus",options = list(maxItems = 1L))
                 # ,sliderInput("SelectP3_percent", "Threshold_percent",min = 0, max = 0.1, value = 0.005)
                 #percent_<-input$SelectP3_percent
                 #actionButton("run2", "run advanced analysis") ,
                 ,br()
                 ,selectizeInput("beta_method", "Beta method",c("DCA", "CCA", "RDA",  "NMDS", "MDS", "PCoA"),selected="MDS",options = list(maxItems = 1L))
                 
                 #a(em("open report link  请点我 "),href="out_book/gv.16s.total6.html"),
                 ,br(),
                 #sidebarPanel(downloadButton('downloadData', 'Download the results'))
               ),
               mainPanel(
                 tabsetPanel(type = "tabs"
                             #tabPanel("Table3",h5("此表格为目标分类组筛选后tables") , DT::dataTableOutput("AD_plot_table3")), 
                             ,tabPanel("alpha_diversity", plotOutput("AD_plot_a",width = 1200,height = '1000')
                                       ,h3("Alpha diversity refers to the diversity in a specific area or ecosystem. Commonly used ways include chao, shannon, ace, simpson, and coverage. In this function module, you can obtain species diversity by observing various values. And other information, you can also group the samples and use the statistical T test method to detect whether the value between each group is significantly different. \n")
                                       ,h3( "Chao: It is an index that uses the chao1 algorithm to estimate the number of OTUs contained in a sample. chao1 is commonly used in ecology to estimate the total number of species and was first proposed by Chao (1984). \n")
                                       ,h3( " Ace: The index used to estimate the number of OTUs in the community, proposed by Chao, is one of the commonly used indices for estimating the total number of species in ecology, and is different from Chao1's algorithm. \n")
                                       ,h3( " Shannon: Used to estimate one of the microbial diversity indexes in the sample. Shannon and the Simpson diversity index are often used to reflect the alpha diversity index. The larger the Shannon value, the higher the community diversity. \n")
                                       ,h3( " Coverage: refers to the coverage of each sample (clone) library. The higher the value, the higher the probability that the sequence in the sample will be detected, and the lower the probability of not being detected. This index reflects whether the results of this sequencing represent the true conditions of microorganisms in the sample. \nSimpson: It is used to estimate one of the microbial diversity indexes in samples. It was proposed by EdwardHugh Simpson (1949). It is commonly used in ecology to quantitatively describe the biodiversity of an area. \n")
                             )
                             ,tabPanel("Beta_diversity", plotOutput("AD_plot_b",width = 800,height = '700')
                                       ,h3("DCA:Performs detrended correspondence analysis usingdecorana  ")
                                       ,h3("CCA:Performs correspondence analysis, or optionally, constrained correspondence analysis (a.k.a. canonical correspondence analysis), via cca")
                                       ,h3("RDA:Performs redundancy analysis, or optionally principal components analysis, via rda")
                                       ,h3("NMDS:Performs Non-metric MultiDimenstional Scaling of a sample-wise ecological distance matrix onto a user-specified number of axes, k. By default, k=2, but this can be modified as a supplementary argument. This method is ultimately carried out by metaMDS after the appropriate accessions and distance calculations. Because metaMDS includes its own distance calculation wrappers to vegdist, and these provide additional functionality in the form of species scores, ordinate will pass-on the distance argument to metaMDS if it is among the supported vegdist methods. However, all distance methods supported by distance are supported here, including 'unifrac' (the default) and 'DPCoA' ")
                                       ,h3("MDS/PCoA:Performs principal coordinate analysis (also called principle coordinate decomposition, multidimensional scaling (MDS), or classical scaling) of a distance matrix (Gower 1966), including two correction methods for negative eigenvalues.")
                             )
                             ,tabPanel("Composition"
                                       , plotOutput("AD_plot_compositional",width = 630,height = '800')
                                       )
                             ,tabPanel("Metacoder", plotOutput("AD_plot_Metacoder",width = 800,height = '800'))
                             ,tabPanel("Netcomi",br()
                                       # ,actionButton("runP5_netcomi", "run netcomi")
                                       ,textOutput("NetCoMi_make")
                                       ,plotOutput("NetCoMi_png",width = 800,height = '800')
                                       ,h3("Note: The green line represents the positive correlation between the points, and the red line represents the negative correlation between the points.")
                                       ,h3("  The size of a point represents the importance of a node, which includes the number of neighbor nodes of the point, and also depends on the importance of its neighbor nodes.")
                                       
                             )
                             
                 )
               )))
    
    ,tabPanel('Comparative analysis',
              sidebarLayout(
                sidebarPanel(
                  titlePanel("Comparative analysis")
                  ,br()
                  #,sidebarPanel(downloadButton('downloadData', 'Download the results'))
                  ,selectizeInput("SelectP4_1_tax", "Select a level for taxonomy" ,c("Phylum","Family", "Order","Class","Genus","OTU"),selected='Genus',options = list(maxItems = 1L))
                  ,br()
                  # ,selectizeInput("n_componet", "N_componet" ,c(2,3),selected=2,options = list(maxItems = 1L))
                  
                  ,numericInput("n_componet", "N_componet in pca and splda",min = 2, max = 3, value = 2)#1
                  # ,sliderInput("width_value", "The size of picture:width_value",min = 1, max = 1000, value = 1000)#2
                  # ,sliderInput("height_value", "The size of picture:height_value",min = 1, max = 1000, value = 400)#2
                )
                ,mainPanel(
                  tabsetPanel(type = "tabs"
                              #,img(src = "16s.png", height = 600, width = 900)
                              ,tabPanel("Deseq2"
                                        ,h3("  Note:Top 20 will be displayed, but you can also view detailed terms from csv.") 
                                        ,sliderInput("deseq2_foldchange_num", "Deseq2:Foldchange",min = 0.0, max = 10, value = 0)#1
                                        ,sliderInput("deseq2_pvalue", "Deseq2:pvalue",min = 0.001, max = 1, value = 1)#0.05
                                        ,sliderInput("deseq2_padj_num", "Deseq2:Padj",min = 0.00001, max = 1, value = 1)#0.05
                                        ,sliderInput("Deseq2_height", "Deseq2_height",min = 0.0, max = 100, value = 10)
                                        ,sliderInput("Deseq2_width", "Deseq2_width",min = 0.0, max = 100, value = 10)
                                        ,imageOutput("AD_plot_deseq2") 
                                        # ,plotOutput("AD_plot_deseq2_picture",width = 500,height = '500') 
                              )
                              
                              ,tabPanel("Lefse"
                                        # ,textOutput("AD_plot_lefse")
                                        ,sliderInput("lefse_lda_value", "Lefse:LDA value",min = 1, max = 10, value = 2)#2
                                        ,sliderInput("lefse_height", "lefse_height",min = 0.0, max = 5000, value = 4000)
                                        ,sliderInput("lefse_width", "lefse_width",min = 0.0, max = 5000, value = 4000)
                                        ,imageOutput("AD_plot_lefse_picture") 
                              )
                              ,tabPanel("PCA" 
                                        ,textOutput("AD_plot_pca_bar_make")
                                        ,textOutput("AD_plot_pca_make")
                                        ,plotOutput("AD_plot_pca_bar",width = 800,height = '600') 
                                        ,plotOutput("AD_plot_pca",width = 500,height = '500')    )
                              ,tabPanel("Pls.da"
                                        #,plotOutput("AD_plot_pls.da")
                                        ,textOutput("AD_plot_pls.da")
                                        ,h5("Show OTU level")
                                        ,plotOutput("AD_plot_plda.tu1",width = 800,height = '800')
                                        ,plotOutput("AD_plot_plda.tu2",width = 800,height = '800')
                              )
                              ,tabPanel("Spls.da"
                                        #,plotOutput("AD_plot_spls.da",width = 800,height = '800')
                                        ,textOutput("AD_plot_spls.da")
                                        #,img(src =  paste0("out.result/","Genus","_sPLSDA_comp2.png"), height = 600, width = 900)
                                        
                                        #,actionButton("runP4_splda_plot", "run runP4_splda_plot")
                                        #,img(src = "16s.png", height = 600, width = 900)
                                        ,plotOutput("AD_plot_splda.tu1",width = 800,height = '900')
                                        # ,plotOutput("AD_plot_splda.tu2",width = 800,height = '900')
                                        ,plotOutput("AD_plot_splda.tu_comp1",width = 800,height = '800')
                                        # ,plotOutput("AD_plot_splda.tu_comp2",width = 800,height = '800')
                                        # ,plotOutput("AD_plot_splda.tu_comp3",width = 800,height = '800')
                                        
                                        ,plotOutput("AD_plot_splda.heatmap1",width = 800,height = '800')
                                        # ,plotOutput("AD_plot_splda.heatmap2",width = 800,height = '800')
                                        ,plotOutput("AD_plot_splda.heatmap3",width = 800,height = '800')
                              )
                              
                  )
                  
                )))
    
    ,tabPanel('Functional analysis',
              sidebarLayout(
                sidebarPanel(
                  titlePanel("Functional analysis")
                  # ,selectizeInput("SelectP5_2", "Select a project (for example:SRP126121.meta.csv)", names(list.meta.path.basename2),options = list(maxItems = 1L)  )
                  # ,selectizeInput("SelectP5_1", "Select phyloseq name (for example:SRP126121.phyloseq.rds)", list.phyloseq.path.basename,options = list(maxItems = 1L)  )
                  # 
                  # ,selectizeInput("SelectP5_3", "Select column of meta(target 1 column filtered)",c("SampleID") ,options = list(maxItems = 1L) )
                  # ,selectizeInput("SelectP5_4", "Select a classficiation for the column (target 1 column filtered)",c("Item A", "Item B"),options = list(maxItems = 1L) )
                  # 
                  # ,selectizeInput("SelectP5_5", "Select column of meta(target 2 column filtered)",c("SampleID") ,multiple = TRUE,options = list(maxItems = 2L) )
                  # ,selectizeInput("SelectP5_6", "Select a classficiation for the column (target 2 column filtered)",c("Item A", "Item B"),multiple = TRUE,options = list(maxItems = 2L) )
                  # ,br()
                  ,actionButton("runP5_1", "Pircust2")
                  ,br()
                  ,br()
                  ,br()
                  ,selectizeInput("analysis_method", "Select a method" ,c("wilcox.test","t.test"),selected='wilcox.test',options = list(maxItems = 1L))
                  ,sliderInput("p_value_threshold", "p.adj",min = 0.0, max = 0.1, value = 0.05)
                  ,actionButton("runP5_2", "Pircust2 view")
                  
                  #,sidebarPanel(downloadButton('downloadData', 'Download the results'))
                  ,br()
                  ,br()
                  ,br()
                  ,actionButton("runP5_amon", "Amon")
                  #,actionButton('switchtab',"Click this refresh")
                  ,textOutput('code_ran')
                )
                ,mainPanel(
                  tabsetPanel(type = "tabs"
                              #,img(src = "16s.png", height = 600, width = 900)
                              #,tabPanel("AD_plot_deseq2", plotOutput("AD_plot_deseq2",width = 1500,height = '800'))
                              # ,textOutput('pircust_value')
                              ,tabPanel("Pircust2"
                                        # ,textOutput('pircust_value')
                                        , plotOutput("AD_plot_Pircust2_percent",width = 1000,height = '1500')
                                        # , plotOutput("AD_plot_Pircust2_aldex",width = 1500,height = '800')
                              )
                              ,tabPanel("Amon"
                                        ,br()
                                        #,actionButton("runP5_amon", "run Amon")
                                        ,br()
                                        ,br()
                                        ,dataTableOutput("amon_table_daixie")
                              )))
              ))
    ,navbarMenu("More",
                tabPanel("Refresh",actionButton('switchtab',"Click this") )
                ,tabPanel("downloadData",downloadButton('downloadData', 'Download the results'))
                ,tabPanel("R",verbatimTextOutput("sessionInfo")  )      
    )
  ))

