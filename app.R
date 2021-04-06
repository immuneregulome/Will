###load packages
library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(gridExtra)
library(gridGraphics)
library(ggpubr)
library(writexl)

###load data
nk_data <- get(load(file = "data/nk_integrated_labeled.Rdata"))
nk_data_markers <- get(load(file = "data/nk_integrated_markers.Rdata"))
DefaultAssay(nk_data) <- "RNA"

appps1_data <- get(load(file = "data/appps1_lymphocytes_labeled.Rdata"))
appps1_data_markers <- get(load(file = "data/appps1_markers.Rdata"))
DefaultAssay(appps1_data) <- "RNA"

appps1_whole_data <- get(load(file = "data/appps1_whole_labeled.Rdata"))
appps1_whole_markers <- get(load(file = "data/appps1_whole_markers.Rdata"))
DefaultAssay(appps1_whole_data) <- "RNA"

tcellinfil_data <- get(load(file = "data/tcell_infiltration_labeled.Rdata"))
tcellinfil_data_markers <- get(load(file = "data/tcell_infiltration_markers.Rdata"))
DefaultAssay(tcellinfil_data) <- "RNA"

ido_amit_1M_data <- get(load(file = "data/ido_amit_1M_labeled.Rdata"))
ido_amit_1M_markers <- get(load(file = "data/ido_amit_1M_markers.Rdata"))
DefaultAssay(ido_amit_1M_data) <- "RNA"

ido_amit_3M_data <- get(load(file = "data/ido_amit_3M_labeled.Rdata"))
ido_amit_3M_markers <- get(load(file = "data/ido_amit_3M_markers.Rdata"))
DefaultAssay(ido_amit_3M_data) <- "RNA"

ido_amit_6M_data <- get(load(file = "data/ido_amit_6M_labeled.Rdata"))
ido_amit_6M_markers <- get(load(file = "data/ido_amit_6M_markers.Rdata"))
DefaultAssay(ido_amit_6M_data) <- "RNA"



###Set Assays


## User interface

ui <- fluidPage(
  #this establishes the title panel, including the image
  titlePanel(title = div(img(src ="Shiny-solosis.png", height="10%", width="10%"), 
                         "Shih Lab scRNA Repository", br(), h5("Created by William Montgomery"))),
  #this is how the main panel and side panel are organized
  sidebarLayout(
    sidebarPanel(
      #this is the input widget for dataset selection
      selectInput(inputId = "dataset_selec",
                  label = "Choose which Dataset to explore:",
                  choices = list("NK AD Dataset (Zhang, 2020)", 
                                 "APPPS1 Dataset All Cells (Van Hove, 2019)",
                                 "APPPS1 Dataset Lymphocytes (Van Hove, 2019)",
                                 "Aging T Cell Dataset (Dulken, 2019)", 
                                 "DAM Cells 1 MO (Keren-Shaul, 2017)",
                                 "DAM Cells 3 MO (Keren-Shaul, 2017)",
                                 "DAM Cells 6 MO (Keren-Shaul, 2017)"),
                  selected = "APPPS1 Dataset Lymphocytes (Van Hove, 2019)"),
      div(tableOutput("cluster_names"), style = "font-size:60%"),
      width = 3

    ),
    #this is the main panel instructions, organized into tabs
    mainPanel(
      tabsetPanel(
        tabPanel("UMAP",
                 fluidRow(
                   column(6, 
                          checkboxInput("split_u", "Split the graph?"),
                          conditionalPanel(condition = "input.split_u == true",
                                           #display choices to split by
                                           selectInput(inputId = "metadata_split_u", 
                                                       label = "Choose how to split the Seurat data: ", 
                                                       choices = list("Genotype", "Timepoint"))),
                          
                          ),
                   column(6,
                          checkboxInput("group_u", "Group the graph?"),
                          conditionalPanel(condition = "input.group_u == true",
                                           #display choices to group by
                                           selectInput(inputId = "metadata_group_u",
                                                       label = "Choose how to group the Seurat data:",
                                                       choices = list("Genotype", "Timepoint")))
                          )
                 ), 
                 br(),
                 br(),
                 fluidRow(
                   column(4,
                          textInput("save_name_umap",
                                    label = "Enter a file name: ")
                          ),
                   column(4,
                          conditionalPanel(condition = "input.save_name_umap.length > 0",
                                           selectInput("umap_device", 
                                                       label = "Select file type: ",
                                                       choices = list("PNG", "JPEG", "PDF", "TIFF",
                                                                      "BMP", "SVG")))
                          ),
                   column(4,
                          conditionalPanel("input.save_name_umap.length > 0",
                                           br(),
                                           downloadButton("umap_save", label = "Save UMAP")))
                          ),
                 plotOutput("umap")
                 ),
        tabPanel("Feature Plot",
                 fluidRow(
                   column(4, 
                          textInput(inputId = "gene_fp", 
                                    label = "Enter gene(s) of interest here, separated by commas: ")
                          ),
                   column(4,
                          br(),
                          checkboxInput("split_fp", "Split the graph?")
                          ),
                   column(4, 
                          conditionalPanel(condition = "input.split_fp == true",
                                           #display choices to split by
                                           selectInput(inputId = "metadata_split_fp", 
                                                       label = "Choose how to split the Seurat data: ", 
                                                       choices = list("Genotype", "Timepoint")))
                          )
                 ),
                 
                 #ask users if they want to split the graphs
                
                 
                 br(),
                 fluidRow(
                   column(4, 
                          textInput("save_name_fp",
                                    label = "Enter a file name: ")
                          ),
                   column(4, 
                          conditionalPanel(condition = "input.save_name_fp.length > 0",
                                           selectInput("fp_device", 
                                                       label = "Select file type: ",
                                                       choices = list("PNG", "JPEG", "PDF", "TIFF",
                                                                      "BMP", "SVG")))
                          ),
                   column(4, 
                          br(),
                          conditionalPanel(condition = "input.save_name_fp.length > 0",
                                           downloadButton("fp_save", label = "Save Feature Plot"))
                          )
                 ),
                 #plot the actual plot
                 plotOutput("fp")
                 ),
        tabPanel("Violin Plot",
                 fluidRow(
                   column(4,
                          textInput(inputId = "gene_vp", 
                                    label = "Enter gene(s) of interest here, separated by commas: "),
                          ),
                   column(4,
                          checkboxInput("split_vp", "Split the graph?"),
                          conditionalPanel(condition = "input.split_vp == true",
                                           #display choices to split by
                                           br(),
                                           selectInput(inputId = "metadata_split_vp", 
                                                       label = "Choose how to split the Seurat data: ", 
                                                       choices = list("Genotype", "Timepoint")))
                          ),
                   column(4, 
                          checkboxInput("group_vp", "Group the graph?"),
                          conditionalPanel(condition = "input.group_vp == true",
                                           #display choices to group by
                                           br(),
                                           selectInput(inputId = "metadata_group_vp",
                                                       label = "Choose how to group the Seurat data:",
                                                       choices = list("Genotype", "Timepoint")))
                          )
                   
                 ),
                 br(),
                 br(),
                 fluidRow(
                   column(4, 
                          textInput("save_name_vp",
                                   label = "Enter a file name: ")
                          ),
                   column(4,
                          conditionalPanel(condition = "input.save_name_vp.length > 0",
                                           selectInput("vp_device", 
                                                       label = "Select file type: ",
                                                       choices = list("PNG", "JPEG", "PDF", "TIFF",
                                                                      "BMP", "SVG")))
                          ),
                   column(4, 
                          br(),
                          conditionalPanel(condition = "input.save_name_vp.length > 0", 
                                           downloadButton("vp_save", label = "Save Violin Plot")) 
                          )
                 ),
                 plotOutput("vp")
                 ),
        tabPanel("Dot Plot",
                 fluidRow(
                   column(4, 
                          textInput(inputId = "gene_dp", 
                                    label = "Enter gene(s) of interest here, separated by commas: ")
                          ),
                   column(4,
                          checkboxInput("split_dp", "Split the graph?"),
                          br(),
                          conditionalPanel(condition = "input.split_dp == true",
                                           #display choices to split by
                                           selectInput(inputId = "metadata_split_dp", 
                                                       label = "Choose how to split the Seurat data: ", 
                                                       choices = list("Genotype", "Timepoint")))
                          ),
                   column(4, 
                          checkboxInput("group_dp", "Group the graph?"),
                          conditionalPanel(condition = "input.group_dp == true",
                                           #display choices to group by
                                           selectInput(inputId = "metadata_group_dp",
                                                       label = "Choose how to group the Seurat data:",
                                                       choices = list("Genotype", "Timepoint")))
                          )
                 ),
                 br(),
                 br(),
                 fluidRow(
                   column(4,
                          textInput("save_name_dp",
                                    label = "Enter a file name: ")
                          ),
                   column(4, 
                          conditionalPanel(condition = "input.save_name_dp.length == true",
                                           selectInput("vp_device", 
                                                       label = "Select file type: ",
                                                       choices = list("PNG", "JPEG", "PDF", "TIFF",
                                                                      "BMP", "SVG")))
                          ),
                   column(4, 
                          conditionalPanel(condition = "input.save_name_dp.length == true",
                                           downloadButton("dp_save", label = "Save Dot Plot"))
                          )
                 ),
                 plotOutput("dp")
                 ),
        tabPanel("Proportion Graph",
                 fluidRow(
                   column(1),
                   column(8, imageOutput("pic_umap")), 
                   column(3)
                 ),
                 fluidRow(
                   column(5,
                          selectInput(inputId = "proportion_split", 
                                      label = "Choose what variable to find proportions for: ",
                                      choices = NULL)
                          ),
                   column(2),
                   column(5,
                          radioButtons("prop_graph_type", 
                                       label = "Choose style of proportion graph: ", 
                                       choices = list("Side-by-side Bar Graph", "Stacked Bar Graph"),
                                       selected = character(0)))
                   ),
                 br(),
                 br(),
                 fluidRow(
                   column(4, 
                          textInput("save_name_prop", label = "Enter a file name: ")
                          ),
                   column(4, 
                          conditionalPanel(condition = "input.save_name_prop.length == true",
                                           selectInput("prop_device", 
                                                       label = "Select file type: ",
                                                       choices = list("PNG", "JPEG", "PDF", "TIFF",
                                                                      "BMP", "SVG")))
                          ),
                   column(4, 
                           conditionalPanel(condition = "input.save_name_prop.length == true",
                                            br(),
                                            downloadButton("prop_save", label = "Save Proportion Graph"))
                          )
                 ),
                 plotOutput("prop")
                 ),
        tabPanel("Cluster Comparison",
                 fluidRow(
                   column(1),
                   column(8, imageOutput("pic_umap_1")), 
                   column(3)
                 ),
                 fluidRow(
                   column(4, 
                          selectInput("cluster_1", label = h4("Choose the first cluster: "), choices = NULL)
                          ),
                   column(4,
                          selectInput("cluster_2", label = h4("Choose the second cluster: "), choices = NULL)
                          ),
                   column(4,
                          br(),
                          br(),
                          actionButton("calculate", "Find Markers!")
                          )
                 ),
                 br(),
                 #this makes it so that the helptext and save button will only appear after output
                 # has been generated
                 conditionalPanel(condition = "output.cluster_comparison_output",
                                  helpText(h4("Genes high in the first cluster have positive avg_logFC, 
                                              while genes high in the second cluster have negative avg_logFC.")),
                                  br(),
                                  fluidRow(
                                    column(4, 
                                           textInput("cluster_compare_filename", label = "Enter file name: ")
                                           ),
                                    column(4,
                                           selectInput("cluster_compare_device",
                                                       label = "Choose file type: ",
                                                       choices = list("xlsx", "csv"))
                                           
                                           ),
                                    column(4, 
                                           downloadButton("save_cluster_compare", label = "Save this table")
                                           )
                                  )
                                  
                 ),
                 br(),
                 dataTableOutput("cluster_comparison_output")
                 ),
        tabPanel("Cluster Markers",
                 fluidRow(
                   column(1),
                   column(8, imageOutput("pic_umap_2")), 
                   column(3)
                 ),
                 fluidRow(
                   column(4, 
                          selectInput("cluster_select", 
                                      label = "Please select cluster: ",
                                      choices = NULL)
                          ),
                   column(4, 
                          textInput("cluster_marker_filename", label = "Enter file name: "),
                          br(),
                          br(),
                          downloadButton("save_cluster_markers", label = "Save this table")
                          ),
                   column(4, 
                          selectInput("cluster_marker_device",
                                      label = "Choose file type: ",
                                      choices = list("xlsx", "csv"))
                          )
                 ),
                 br(),
                 br(),
                 dataTableOutput("cluster_markers")
                 )
                 
      )
    )
  )
)


##Define server logic
server <- function(input, output, session) {

  #the below allows us to switch between datasets for most of the renderplot functions
  datasetInput <- reactive({
    switch(input$dataset_selec, 
           "NK AD Dataset (Zhang, 2020)" = nk_data, 
           "APPPS1 Dataset All Cells (Van Hove, 2019)" = appps1_whole_data,
           "APPPS1 Dataset Lymphocytes (Van Hove, 2019)" = appps1_data,
           "Aging T Cell Dataset (Dulken, 2019)" = tcellinfil_data,
           "DAM Cells 1 MO (Keren-Shaul, 2017)" = ido_amit_1M_data,
           "DAM Cells 3 MO (Keren-Shaul, 2017)" = ido_amit_3M_data,
           "DAM Cells 6 MO (Keren-Shaul, 2017)" = ido_amit_6M_data)
  })
  
  #the observe expression allows us to use the updateSelectInput() function, so that we can 
  #dynamically select the clusters and metadata labels based on what dataset is loaded
  observe({
    dataset <- switch(input$dataset_selec, 
                      "NK AD Dataset (Zhang, 2020)" = nk_data, 
                      "APPPS1 Dataset All Cells (Van Hove, 2019)" = appps1_whole_data,
                      "APPPS1 Dataset Lymphocytes (Van Hove, 2019)" = appps1_data,
                      "Aging T Cell Dataset (Dulken, 2019)" = tcellinfil_data, 
                      "DAM Cells 1 MO (Keren-Shaul, 2017)" = ido_amit_1M_data,
                      "DAM Cells 3 MO (Keren-Shaul, 2017)" = ido_amit_3M_data,
                      "DAM Cells 6 MO (Keren-Shaul, 2017)" = ido_amit_6M_data)
    #this expression updates the first cluster input
    updateSelectInput(session, 
                      "cluster_1", 
                      choices = c(as.character(sort(unique(dataset@meta.data$Clusters))),
                                  unique(dataset@meta.data$Genotype),
                                  unique(dataset@meta.data$Timepoint)),
                      selected = 0)
    #this expression updates the second cluster input
    updateSelectInput(session,
                      "cluster_2", 
                      choices = c(as.character(sort(unique(dataset@meta.data$Clusters))), 
                                  unique(dataset@meta.data$Genotype),
                                  unique(dataset@meta.data$Timepoint)),
                      selected = 1)
    
    #this expression updates the cluster input on the "Cluster Markers" tab
    updateSelectInput(session,
                      "cluster_select",
                      choices = as.character(sort(unique(dataset@meta.data$Clusters))),
                      selected = 0)
    #this expression updates the options to find proportions for
    updateSelectInput(session,
                      "proportion_split",
                      choices = as.character(colnames(dataset@meta.data)))
  })
  
  
  #this expression needs to be in an eventReactive expression because I want the user to have to
  # click an action button before anything happens. Therefore, I designate the action button at 
  # the beginning, and then define dataset within the expression, then create the desired dataframe
  
  #I also had to create the weird nested flow structure so that I could switch between Idents when
  #necessary 
  cluster_comparison <- eventReactive(input$calculate, {
    dataset <- datasetInput()
    if (input$cluster_1 %in% levels(Idents(dataset))) {
      cluster_comparison <- FindMarkers(dataset, ident.1 = input$cluster_1, ident.2 = input$cluster_2)
    }
    else {
      Idents(dataset) <- "Genotype"
      if (input$cluster_1 %in% levels(Idents(dataset))) {
        cluster_comparison <- FindMarkers(dataset, ident.1 = input$cluster_1, ident.2 = input$cluster_2)
      }
      else {
        Idents(dataset) <- "Timepoint"
        cluster_comparison <- FindMarkers(dataset, ident.1 = input$cluster_1, ident.2 = input$cluster_2)
        }
    }
    
    
    #add the lines below so that gene names are included in the print out
    cluster_comparison$Genes <- rownames(cluster_comparison)
    cluster_comparison <- cluster_comparison[, c(6, 1, 2, 3, 4, 5)]
  })
  
  #the expression below renders the data table generated by the eventReactive function above, and 
  # assigns it to the output object
  output$cluster_comparison_output <- renderDataTable({
    cluster_comparison()
  })
  
  #the expression below saves the output generated by cluster_comparison
  output$save_cluster_compare <- downloadHandler(
    filename = function() {
      paste(input$cluster_compare_filename, input$cluster_compare_device, sep = ".")
    },
    content = function(file) {
      if (input$cluster_compare_device == "csv") {write.csv(cluster_comparison(), file)}
      else {write_xlsx(cluster_comparison(), file)}
    } 
  )

  
  
  
  ###the section below is for the Cluster Markers tab, it switches between pre-loaded marker file
  # based on the sidebarPanel input for dataset choice
  cluster_markers <- reactive({
    markers <- switch(input$dataset_selec, 
                      "NK AD Dataset (Zhang, 2020)" = nk_data_markers, 
                      "APPPS1 Dataset All Cells (Van Hove, 2019)" = appps1_whole_markers,
                      "APPPS1 Dataset Lymphocytes (Van Hove, 2019)" = appps1_data_markers,
                      "Aging T Cell Dataset (Dulken, 2019)" = tcellinfil_data_markers,
                      "DAM Cells 1 MO (Keren-Shaul, 2017)" = ido_amit_1M_markers,
                      "DAM Cells 3 MO (Keren-Shaul, 2017)" = ido_amit_3M_markers,
                      "DAM Cells 6 MO (Keren-Shaul, 2017)" = ido_amit_6M_markers)
    #this subset makes it so that only the desired cluster markers can be seen at one time
    markers <- subset(markers, cluster == input$cluster_select)
  })
  
  output$cluster_markers <- renderDataTable({
    cluster_markers()
  })
  
  output$save_cluster_markers <- downloadHandler(
    filename = function() {
      paste(input$cluster_marker_filename, input$cluster_marker_device, sep = ".")
    },
    content = function(file) {
      if (input$cluster_marker_device == "csv") {write.csv(cluster_markers(), file)}
      else {write_xlsx(cluster_markers(), file)}
    }
  )
  

  output$umap <- renderPlot({
    
    if (input$split_u == TRUE && input$group_u == TRUE) {DimPlot(datasetInput(), label = TRUE, split.by = input$metadata_split_u, group.by = input$metadata_group_u)}
    else if (input$split_u == TRUE && input$group_u == FALSE) {DimPlot(datasetInput(), label = TRUE, split.by = input$metadata_split_u)}
    else if (input$split_u == FALSE && input$group_u == TRUE) {DimPlot(datasetInput(), label = TRUE, group.by = input$metadata_group_u)}
    else {DimPlot(datasetInput(), label = TRUE)}
  })
  
  output$umap_save <- downloadHandler(
    filename = function() {
      paste(input$save_name_umap, tolower(input$umap_device), sep = ".")
    },
    content = function(file) {
      ggsave(file, device = "png")
    }
  )
  
  #create a function for the featureplots so that we can call this when we want to save them
  feat_plot <- reactive({
    fp_genes <- input$gene_fp
    fp_genes <- gsub(" ", "", fp_genes)
    fp_genes <- unlist(strsplit(fp_genes, split = ","))
    
    if (input$split_fp == TRUE) {
      featlist <- list()
      for (i in 1:length(fp_genes)) {
        gene <- fp_genes[i]
        featlist[[i]] <- FeaturePlot(datasetInput(), features = gene, split.by = input$metadata_split_fp)
      }
      plot_grid(plotlist = featlist, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,1))
    }
    else {
      fplots <- FeaturePlot(datasetInput(), features = c(fp_genes), combine = FALSE)
      grid.arrange(grobs = fplots, widths = unit(5.5, units = "in"), heights = unit(c(rep(5, length(fplots))), units = "in"))
    }
  })
  
  #actually render the grid.arrange() object created in feat_plot
  output$fp <- renderPlot({
    
    #this is to prevent an error message from being displayed when a gene hasn't been entered yet
    validate(
      need(input$gene_fp !="", "Please enter a gene.")
    )
    
    feat_plot()

  }, height = reactive({ #this reactive expression allows the height change dynamically with gene inputs
    validate(
      need(input$gene_fp !="", "Please enter a gene.")
    )
    fp_genes <- input$gene_fp
    fp_genes <- gsub(" ", "", fp_genes)
    fp_genes <- unlist(strsplit(fp_genes, split = ","))
    (400*length(fp_genes))
  }), width = reactive({
    validate(
      need(input$gene_fp !="", "Please enter a gene.")
    )
    if (input$split_fp == TRUE) {1000}
    else if (input$split_fp == FALSE) {500}}))
  
  #this function enables us to call this for the ggsave command
  fp_height <- reactive({
    fp_genes <- input$gene_fp
    fp_genes <- gsub(" ", "", fp_genes)
    fp_genes <- unlist(strsplit(fp_genes, split = ","))
    as.numeric(5*length(fp_genes))})
  
  fp_width <- reactive({ 
    if (input$split_fp == TRUE) {11}
    else if (input$split_fp == FALSE) {5.5}})

  output$fp_save <- downloadHandler(
    filename = function() {
      paste(input$save_name_fp, tolower(input$fp_device), sep = ".")
    },
    content = function(file) {
      ggsave(file, feat_plot(), device = tolower(input$fp_device), 
             height = fp_height(), 
             width = fp_width(),
             units = "in")
    }
  )
  
  viol_plot <- reactive({
    vp_genes <- input$gene_vp
    vp_genes <- gsub(" ", "", vp_genes)
    vp_genes <- unlist(strsplit(vp_genes, split = ","))
    
    if (input$split_vp == TRUE && input$group_vp == TRUE) {
      vplots <- VlnPlot(datasetInput(), features = vp_genes, split.by = input$metadata_split_vp, 
              group.by =input$metadata_group_vp, combine = FALSE)
      grid.arrange(grobs = vplots, widths = unit(10, units = "in"), heights = unit(c(rep(5, length(vp_genes))), units = "in"))
      }
    else if (input$split_vp == TRUE && input$group_vp == FALSE) {
      vplots <- VlnPlot(datasetInput(), features = vp_genes, split.by = input$metadata_split_vp, combine = FALSE)
      grid.arrange(grobs = vplots, widths = unit(10, units = "in"), heights = unit(c(rep(5, length(vp_genes))), units = "in"))
      }
    else if (input$split_vp == FALSE && input$group_vp == TRUE) {
      vplots <- VlnPlot(datasetInput(), features = vp_genes, group.by = input$metadata_group_vp, combine = FALSE)
      grid.arrange(grobs = vplots, widths = unit(7, units = "in"), heights = unit(c(rep(5, length(vp_genes))), units = "in"))
      }
    else {
      vplots <- VlnPlot(datasetInput(), features = vp_genes, combine = FALSE)
      grid.arrange(grobs = vplots, widths = unit(7, units = "in"), heights = unit(c(rep(5, length(vp_genes))), units = "in"))
      }
  })
  
  #actually render the grid.arrange() object created in viol_plot
  output$vp <- renderPlot({
    
    #this is to prevent an error message from being displayed when a gene hasn't been entered yet
    validate(
      need(input$gene_vp !="", "Please enter a gene.")
    )
    
    viol_plot()
    
  }, height = reactive({ #this eventReactive expression allows the height change dynamically with gene inputs
    validate(
      need(input$gene_vp !="", "Please enter a gene.")
    )
    vp_genes <- input$gene_vp
    vp_genes <- gsub(" ", "", vp_genes)
    vp_genes <- unlist(strsplit(vp_genes, split = ","))
    (400*length(vp_genes))
  }), width = reactive({
    validate(
      need(input$gene_vp !="", "Please enter a gene.")
    )
    if (input$split_vp == TRUE) {800}
    else if (input$split_vp == FALSE) {500}}))
  
  vp_height <- reactive({
    vp_genes <- input$gene_vp
    vp_genes <- gsub(" ", "", vp_genes)
    vp_genes <- unlist(strsplit(vp_genes, split = ","))
    as.numeric(5*length(vp_genes))})
  
  vp_width <- reactive({ 
    if (input$split_vp == TRUE) {10}
    else if (input$split_vp == FALSE) {7}})
  
  output$vp_save <- downloadHandler(
    filename = function() {
      paste(input$save_name_vp, tolower(input$vp_device), sep = ".")
    },
    content = function(file) {
      ggsave(file, viol_plot(), device = tolower(input$vp_device), 
             height = vp_height(), 
             width = vp_width(),
             units = "in")
    }
  )
  
  output$dp <- renderPlot({
    #this is to prevent an error message from being displayed when a gene hasn't been entered yet
    validate(
      need(input$gene_dp !="", "Please enter a gene.")
    )
    
    #since we want to enable input to be in strings of genes, we need to manipulate the input
    dp_genes <- input$gene_dp
    dp_genes <- gsub(" ", "", dp_genes)
    dp_genes <- unlist(strsplit(dp_genes, split = ","))
    
    if (input$split_dp == TRUE && input$group_dp == TRUE) {DotPlot(datasetInput(), features = c(dp_genes), split.by = input$metadata_split_dp, group.by =input$metadata_group_dp)}
    else if (input$split_dp == TRUE && input$group_dp == FALSE) {DotPlot(datasetInput(), features = c(dp_genes), split.by = input$metadata_split_dp)}
    else if (input$split_dp == FALSE && input$group_dp == TRUE) {DotPlot(datasetInput(), features = c(dp_genes), group.by = input$metadata_group_dp)}
    else {DotPlot(datasetInput(), features = c(dp_genes))}
  })
  
  output$dp_save <- downloadHandler(
    filename = function() {
      paste(input$save_name_dp, tolower(input$dp_device), sep = ".")
    },
    content = function(file) {
      ggsave(file, device = tolower(input$dp_device))
    }
  )
  
  #this section is the part for the proportion graphs, it's a lot longer because these graphs
  #are not covered by wrapper functions, I had to design them myself 
  output$prop <- renderPlot({
    
    #we need to redirect the reactive expression into a variable so that we can subset it
    dataset <- datasetInput()
    
    #now we make the proportion table
    Idents(dataset) <- input$proportion_split
    proportion_table <- prop.table(table(dataset$Clusters, Idents(dataset)), margin = 2)
    proportion_table <- as.data.frame(proportion_table)
    names(proportion_table)[1] <- "Cluster"
    names(proportion_table)[2] <- input$proportion_split
    names(proportion_table)[3] <- "Percentage"
    proportion_table$Percentage <- proportion_table$Percentage * 100
    #We have to use aes_string here because otherwise, the graph won't recognize 
    # input$proportion_split as a variable, it will instead try to fill the graph 
    #using that variable. We could have also used aes(... get(input$proportion_split))
    # but then the legend isn't labeled properly 
    
    #validate prevents an error message from being displayed when the user hasn't chosen a graph
    validate(
      need(input$prop_graph_type !="", "Please select a graph type.")
    )
    
    if (input$prop_graph_type == "Side-by-side Bar Graph") {
      ggplot() +
      geom_bar(aes_string(y = "Percentage", x = "Cluster", fill = input$proportion_split), 
               data = proportion_table, stat = "identity", position = position_dodge(width = 0.7))}
    
    else if (input$prop_graph_type == "Stacked Bar Graph") {
    ggplot(proportion_table, 
           aes_string(x = input$proportion_split, y = "Percentage", fill = "Cluster")) +
      geom_bar(color = "black", position = "fill", stat = "identity", width = 0.3, size = 0.25) + 
      theme_classic() +
      geom_text(aes(label = Cluster), position = position_fill(vjust = 0.5), size = 2.5) +
      coord_flip() + 
      NoLegend() +
      theme(axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10)) + 
      theme(axis.text.x = element_text(vjust = 0.5)) + 
      theme(axis.line = element_line(size = 0.25)) +
      ylab("Proportion")}
  })
  
  output$prop_save <- downloadHandler(
    filename = function() {
      paste(input$save_name_prop, tolower(input$prop_device), sep = ".")
    },
    content = function(file) {
      ggsave(file, device = tolower(input$prop_device))
    }
  )
  
  #have to assign this renderimage to multiple outputs because shiny doesn't let you share output
  #names between multiple output calls
  output$pic_umap <- output$pic_umap_1 <- output$pic_umap_2 <- renderImage({
    filename <- switch(input$dataset_selec,
                    "NK AD Dataset (Zhang, 2020)" = "~/Desktop/Shiny App/Seuratapp/www/NK_umap.jpeg",
                    "APPPS1 Dataset All Cells (Van Hove, 2019)" = "~/Desktop/Shiny App/Seuratapp/www/appps1_whole_umap.jpeg",
                    "APPPS1 Dataset Lymphocytes (Van Hove, 2019)" = "~/Desktop/Shiny App/Seuratapp/www/appps1_lymphocytes.jpeg",
                    "Aging T Cell Dataset (Dulken, 2019)" = "~/Desktop/Shiny App/Seuratapp/www/aging_t_cell_umap.jpeg", 
                    "DAM Cells 1 MO (Keren-Shaul, 2017)" = "~/Desktop/Shiny App/Seuratapp/www/ido_amit_1M_umap.jpeg",
                    "DAM Cells 3 MO (Keren-Shaul, 2017)" = "~/Desktop/Shiny App/Seuratapp/www/ido_amit_3M_umap.jpeg",
                    "DAM Cells 6 MO (Keren-Shaul, 2017)" = "~/Desktop/Shiny App/Seuratapp/www/ido_amit_6M_umap.jpeg")
    list(src = filename, 
         width = 420, 
         height = 300)
    
    
  }, deleteFile = FALSE)
  
  output$cluster_names <- renderTable({
    filename <- switch(input$dataset_selec,
                       "NK AD Dataset (Zhang, 2020)" = "~/Desktop/Shiny App/Seuratapp/data/NK_cluster_names.csv",
                       "APPPS1 Dataset All Cells (Van Hove, 2019)" = "~/Desktop/Shiny App/Seuratapp/data/appps1_whole_cluster_names.csv",
                       "APPPS1 Dataset Lymphocytes (Van Hove, 2019)" = "~/Desktop/Shiny App/Seuratapp/data/appps1_lymphocytes_cluster_names.csv",
                       "Aging T Cell Dataset (Dulken, 2019)" = "~/Desktop/Shiny App/Seuratapp/data/tcell_infiltration_cluster_names.csv",
                       "DAM Cells 1 MO (Keren-Shaul, 2017)" = "~/Desktop/Shiny App/Seuratapp/data/Ido_Amit_M1_cluster_names.csv",
                       "DAM Cells 3 MO (Keren-Shaul, 2017)" = "~/Desktop/Shiny App/Seuratapp/data/Ido_Amit_M3_cluster_names.csv",
                       "DAM Cells 6 MO (Keren-Shaul, 2017)" = "~/Desktop/Shiny App/Seuratapp/data/Ido_Amit_M6_cluster_names.csv")
    read.csv(file = filename, header = TRUE)
  }, striped = TRUE, bordered = TRUE)



}

## Launch App
shinyApp(ui = ui, server = server)


