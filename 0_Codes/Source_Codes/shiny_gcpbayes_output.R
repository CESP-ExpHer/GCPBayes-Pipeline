# 					Shiny App for Visualization of GCPBayes Output
# ================================================================================
# Summary: The code reads an output of GCPBayes (candidate pleiotropic genes)
# and creates different plots based on it
# ================================================================================
# Written first by: Yazdan Asgari
# Modified by: 
# Initial Creation Date: 10/2022
# Edited Date: 12/2022
# ================================================================================
library(shiny)
library(datasets)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(BioCircos)
library(plotly)
library(ggpubr)
library(biomaRt)
library(regioneR)
library(karyoploteR)

ui <- shinyUI(fluidPage(
  titlePanel("GCPBayes Pipeline Plot Interface"),
  tabsetPanel(
    tabPanel("Upload File",
             titlePanel("Uploading Files"),
             sidebarLayout(
               sidebarPanel(
                 fileInput('file1', 'Choose File',
                           accept=c('text/csv', 
                                    'text/comma-separated-values,text/plain', 
                                    '.csv')),
                 
                 # added interface for uploading data from
                 # http://shiny.rstudio.com/gallery/file-upload.html
                 tags$br(),
                 checkboxInput('header', 'Header', TRUE),
                 radioButtons('sep', 'Separator',
                              c(Comma=',',
                                Semicolon=';',
                                Tab='\t',
                                Space=''),
                              ','),
                 radioButtons('quote', 'Quote',
                              c(None='',
                                'Double Quote'='"',
                                'Single Quote'="'"),
                              '"')
                 
               ),
               mainPanel(
                 tableOutput('contents')
               )
             )
    ),
    tabPanel("Scatter Plots",
             pageWithSidebar(
               headerPanel('Different Scatter Plots'),
               sidebarPanel(
                 
                 # "Empty inputs" - they will be updated after the data is uploaded
                 selectInput('xcol', 'X Variable', ""),
                 selectInput('ycol', 'Y Variable', "", selected = "")
                 
               ),
               mainPanel(
                 plotOutput('ScatterPlot'), downloadButton(outputId = "downloadButsc", label = "Download PDF")
               )
             )
    ),
    tabPanel( "Theta Histogram", plotOutput('Thetahist'), downloadButton(outputId = "downloadButth", label = "Download PDF" ) ),
    tabPanel( "Circos Plot", BioCircosOutput('Circoseplot', height = "800px") ),
    tabPanel( "Pie Chart", plotlyOutput('Piechart', height = "800px") ),
    tabPanel( "Box Plot", plotOutput('Boxplot'), downloadButton(outputId = "downloadButbp", label = "Download PDF" ) ),
    tabPanel("Karyotype_chr",
             pageWithSidebar(
               headerPanel('Karyotype Plots in each Chromosome'),
               sidebarPanel(
                 selectInput(inputId = "ChrKP",
                             label = "Select a chromosome",
                             choices = c(seq(1:22)),
                             selected = "1"
                 )
               ),
               mainPanel(
                 plotOutput('Karyotypechrplot'), downloadButton(outputId = "downloadButkchp", label = "Download PDF")
               )
             )
    ),
    tabPanel( "Karyotype_All", plotOutput('Karyotypeplot'), downloadButton(outputId = "downloadButkp", label = "Download PDF" ) )
    
  )
)
)

server <- shinyServer(function(input, output, session) {
  # added "session" because updateSelectInput requires it
  options(shiny.maxRequestSize=30*1024^2)
  
  pleiodata <- reactive({ 
    req(input$file1) ## ?req #  require that the input is available
    
    inFile <- input$file1 
    
    # tested with a following dataset: write.csv(mtcars, "mtcars.csv")
    # and                              write.csv(iris, "iris.csv")
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    
    # Update inputs (you could create an observer with both updateSel...)
    # You can also constraint your choices. If you wanted select only numeric
    # variables you could set "choices = sapply(df, is.numeric)"
    # It depends on what do you want to do later on.
    
    updateSelectInput(session, inputId = 'xcol', label = 'X Variable',
                      choices = names(df), selected = names(df)[4])
    updateSelectInput(session, inputId = 'ycol', label = 'Y Variable',
                      choices = names(df), selected = names(df)[3])
    
    return(df)
  })
  
  output$contents <- renderTable({
    pleiodata()
  })
  
  output$ScatterPlot <- renderPlot({
    x <- pleiodata()[, c(input$xcol, input$ycol)]
    plot(x)
    
    output$downloadButsc <- downloadHandler(
      filename = function() {
        paste("ScatterPlot.pdf")
      },
      content = function(file) {
        pdf(file)
        plot(x)
        dev.off()
      }
    )
  })
  
  output$Thetahist <- renderPlot({
    req(input$file1)
    inFile <- input$file1
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    p1 <- ggplot(df) +
      geom_histogram(aes(x = Theta), boundary = 0, color="black", fill="#00FFFF", binwidth = 0.02)+
      stat_bin( geom="text", colour="black", size=3.5,
                aes(x = Theta, label=..count..), angle = 90, boundary = 0, binwidth = 0.02, hjust = -.05)+
      scale_x_continuous(limits=c(0.5,1))+
      labs(title = "Theta Histogram (0.5 to 1)")
    
    p2 <- ggplot(df) +
      geom_histogram(aes(x=Theta), boundary = 0, color="black", fill="#FFC857", binwidth = 0.05)+
      stat_bin( geom="text", colour="black", size=3.5,
                aes(x = Theta, label=..count..), angle = 90, boundary = 0, binwidth = 0.05, hjust = -.05)+
      labs(title = "Theta Histogram Zoom (Bin = 0.05)")
    
    grid.arrange(p1, p2, ncol=2) 
    
    output$downloadButth <- downloadHandler(
      filename = function() {
        paste("Thetahist.pdf")
      },
      content = function(file) {
        pdf(file)
        grid.arrange(p1, p2, ncol=2)
        dev.off()
      }
    )
    
  })
  
  output$Circoseplot <- renderBioCircos({
    req(input$file1)
    inFile <- input$file1
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.5, maxRadius = 0.8,
                                         borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#FFBBBB")  
    
    BioCircos(tracklist, genomeFillColor = "PuOr",
              chrPad = 0.05, displayGenomeBorder = FALSE, 
              genomeTicksDisplay = FALSE,  genomeLabelTextSize = "9pt", genomeLabelDy = 0)
    
    points_chromosomes = df$chr
    points_coordinates = df$start 
    points_values = df$gene_length
    
    tracklist = BioCircosSNPTrack('myGeneTrack', points_chromosomes, points_coordinates, 
                                  points_values, colors = "Spectral", minRadius = 0.1, maxRadius = 0.45)
    
    
    tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack1", 
                                                     minRadius = 0.1, maxRadius = 0.45,
                                                     borderColors = "#AAAAAA", borderSize = 0.6)  
    
    snps_values = df$snp_number
    
    tracklist = tracklist + BioCircosSNPTrack('mySNPTrack', points_chromosomes, points_coordinates, 
                                              snps_values, colors = "Dark2", minRadius = 0.5, maxRadius = 0.9)
    
    # Background are always placed below other tracks
    tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack2", 
                                                     fillColor = "#FFD700", borderSize = 0.6, maxRadius = 0.9)  
    
    tracklist = tracklist + BioCircosTextTrack("testText", 'Genes/#SNP Distributions', weight = "bold", size = "0.9em",
                                               x = -0.3 , y = - 1.37)
    
    BioCircos(tracklist, genomeFillColor = "PuOr",
              chrPad = 0.05, displayGenomeBorder = FALSE, yChr =  FALSE,
              genomeTicksDisplay = FALSE,  genomeLabelTextSize = 18, genomeLabelDy = 0)
    
  })
  
  output$Piechart <- renderPlotly({
    req(input$file1)
    inFile <- input$file1
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    gene_list <- df %>%
      mutate(chr_character = sub("^", "chr", df$chr ))
    gene_list <- as.data.frame(gene_list)
    gene_df <- data.frame(sort(table(gene_list$chr_character), decreasing = TRUE))
    colnames(gene_df) <- c("Chromosome", "Number")
    
    plot_ly(gene_df, labels = ~Chromosome, values = ~Number, type = 'pie',textposition = 'outside',textinfo = 'label+percent')%>%
      layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    
  })
  
  
  output$Boxplot <- renderPlot({
    req(input$file1)
    inFile <- input$file1
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    df <- data.frame(gene_length=c(df$gene_length),
                     gene_group=rep(c("gene_list"), times=c(length(df$gene_length)) )
    )
    df$gene_group <- factor(df$gene_group , levels=c("gene_list"))
    
    bp <- ggplot(df, aes(x = gene_group, y = gene_length, colour = gene_group)) + 
      geom_boxplot()
    
    jp <- ggplot(df, aes(x = gene_group, y = gene_length, colour = gene_group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15)
    
    vp <- ggplot(df, aes(x = gene_group, y = gene_length, colour = gene_group)) +
      geom_violin() +
      ggpubr::yscale("log10", .format = TRUE)
    
    # ggarrange(bp, jp , vp, 
    #           labels = c("A", "B", "C"),
    #           ncol = 2, nrow = 2)
    grid.arrange(bp, jp, vp, ncol=3)
    
    
    output$downloadButbp <- downloadHandler(
      filename = function() {
        paste("Boxplot.pdf")
      },
      content = function(file) {
        pdf(file)
        grid.arrange(bp, jp, vp, ncol=3)
        dev.off()
      }
    )
    
  })
  
  x_reactive2 <- reactive({
    
    req(input$file1)
    inFile <- input$file1
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    switch(input$ChrKP,  # Look at the value of Chromosome...
           "1" = dplyr::filter(df, chr ==1),
           "2" = dplyr::filter(df, chr ==2), 
           "3" = dplyr::filter(df, chr ==3),
           "4" = dplyr::filter(df, chr ==4),
           "5" = dplyr::filter(df, chr ==5),
           "6" = dplyr::filter(df, chr ==6),
           "7" = dplyr::filter(df, chr ==7),
           "8" = dplyr::filter(df, chr ==8),
           "9" = dplyr::filter(df, chr ==9),
           "10" = dplyr::filter(df, chr ==10),
           "11" = dplyr::filter(df, chr ==11),
           "12" = dplyr::filter(df, chr ==12),
           "13" = dplyr::filter(df, chr ==13),
           "14" = dplyr::filter(df, chr ==14),
           "15" = dplyr::filter(df, chr ==15),
           "16" = dplyr::filter(df, chr ==16),
           "17" = dplyr::filter(df, chr ==17),
           "18" = dplyr::filter(df, chr ==18),
           "19" = dplyr::filter(df, chr ==19),
           "20" = dplyr::filter(df, chr ==20),
           "21" = dplyr::filter(df, chr ==21),
           "22" = dplyr::filter(df, chr ==22),
    )
    
  })
  
  
  output$Karyotypechrplot <- renderPlot({
    
    geneSymbols <- as.character(unlist(x_reactive2()$GeneName))
    
    # reading ensemble data from biomart
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    
    # find ensemble id for our genes
    genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                             filters = 'hgnc_symbol', values =geneSymbols, mart = ensembl))
    
    # checking genes headers
    seqlevelsStyle(genes) <- "UCSC"
    
    kp <- plotKaryotype("hg19", plot.type=2, chromosomes = paste0("chr", input$ChrKP))
    kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "vertical",
                  r1=0.5, cex=0.6, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE ,
                  label.color = "blue")
    output$downloadButkchp <- downloadHandler(
      filename = function() {
        paste("Karyotypechr.pdf")
      },
      content = function(file) {
        pdf(file)
        kp <- plotKaryotype("hg19", plot.type=2, chromosomes = paste0("chr", input$ChrKP))
        kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "vertical",
                      r1=0.5, cex=0.6, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE ,
                      label.color = "blue")
        dev.off()
      }
    )
    
  })   
  
  output$Karyotypeplot <- renderPlot({
    req(input$file1)
    inFile <- input$file1
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    
    geneSymbols <- as.character(unlist(df$GeneName))
    
    # reading ensemble data from biomart
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    
    # find ensemble id for our genes
    genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                             filters = 'hgnc_symbol', values =geneSymbols, mart = ensembl))
    
    # checking genes headers
    seqlevelsStyle(genes) <- "UCSC"
    
    kp <- plotKaryotype("hg19", plot.type=2, chromosomes = c("chr1","chr2", "chr3", "chr4","chr5", "chr6", 
                                                             "chr7","chr8", "chr9", "chr10",
                                                             "chr11", "chr12", "chr13", "chr14", "chr15",
                                                             "chr16", "chr17", "chr18", "chr19",
                                                             "chr20", "chr21", "chr22"))
    kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "vertical",
                  r1=0.5, cex=0.6, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE ,
                  label.color = "blue")
    
    output$downloadButkp <- downloadHandler(
      filename = function() {
        paste("Karyotypeplot.pdf")
      },
      content = function(file) {
        pdf(file)
        kp <- plotKaryotype("hg19", plot.type=2, chromosomes = c("chr1","chr2", "chr3", "chr4","chr5", "chr6", 
                                                                 "chr7","chr8", "chr9", "chr10",
                                                                 "chr11", "chr12", "chr13", "chr14", "chr15",
                                                                 "chr16", "chr17", "chr18", "chr19",
                                                                 "chr20", "chr21", "chr22"))
        kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "vertical",
                      r1=0.5, cex=0.6, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE ,
                      label.color = "blue")
        dev.off()
      }
    )
    
  })   
  
  
})

shinyApp(ui, server)