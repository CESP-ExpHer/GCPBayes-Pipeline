# Load necessary libraries
library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinythemes)

# ================================================================
# Define R parameters
work_dir <- "C:/CESP/AMLAP/4_Manuscripts/GCPBayes-Pipeline/0_code_r_all_in_one/shiny_all_pipeline"
input1 <- "gwas_BCAC_chr5.txt"
input2 <- "gwas_OCAC_chr5.txt"
input_shortname1 <- "BCAC"
input_shortname2 <- "OCAC"
g1_rsid <- "snp"
g1_chr  <- "chr"
g1_pos  <- "bp_hg19"
g1_EA   <- "Effect_A"
g1_nE_A <- "Baseline_A"
g1_beta <- "beta"
g1_se   <- "se"
g1_pval <- "pval"
g1_info <- "info"
g1_EAF  <- "EAF"
g1_MAF  <- "MAF"
g2_rsid <- "snp"
g2_chr  <- "chr"
g2_pos  <- "bp_hg19"
g2_EA   <- "Effect_A"
g2_nE_A <- "Baseline_A"
g2_beta <- "beta"
g2_se   <- "se"
g2_pval <- "pval"
g2_info <- "info"
g2_EAF  <- "EAF"
g2_MAF  <- "MAF"
input_annot_file1 <- "annot_gencode_v38lift37_modified_gene_class.txt"
input_annot_file2 <- "Annot_BCAC_2020_onco_ALL_reformatted_coding.txt"
info_threshold <- 0.8
MAF_threshold <- 0.01
group_clump_threshold <- 700
theta_exploration <- 0.5
toclump <- FALSE
input_ref_path_file <- "C:/Users/Yazdan/Downloads/r_all_in_one/files"
input_ref_name_file <- "EUR"
clump_threshold_r2 <- 0.8
clump_threshold_kb <- 10000
clump_threshold_p <- 0.99
placo_pval_threshold <- 0.05

# ================================================================

# Define a theme
shinytheme <- shinytheme("cerulean")


tags$head(tags$style(HTML('
  .logo {
    background-image: url("logo.png");
    background-repeat: no-repeat;
    background-size: contain;
    height: 50px;
    width: 50px;
    float: right;
    margin-top: 10px;
    margin-right: 10px;
  }
')))

# ================================================================
# Define UI
ui <- fluidPage(theme = shinytheme,

                # # Add a title
                # tags$h1("GCPBayes Pipeline Shiny App"),
                # tags$h5("Created by: Yazdan Asgari", style = "font-style: italic;"),
                # # tags$h5("Author: Yazdan Asgari", style = "color: blue;"),


                # Add a title and logo
                div(
                  style = "display: flex; align-items: center; justify-content: space-between;",
                  tags$h1("GCPBayes Pipeline Shiny App"),
                  tags$img(src = "logo.png", height = "100px", width = "200px")
                ),
                tags$h5("Created by: Yazdan Asgari", style = "font-style: italic;"),

                # Text input widgets for work_dir
                wellPanel(
                  h3("Select Directory PATH for the files"),
                  textInput("work_dir_input", "Work Directory", value = work_dir)
                ),

                tags$hr(),

                # Group related inputs together (File inputs)
                wellPanel(
                  h3("Select GWAS Files Names"),
                  fluidRow(
                    column(width = 6,
                           fileInput("file1_input", "Choose First GWAS file")
                    ),
                    column(width = 6,
                           fileInput("file2_input", "Choose Second GWAS file")
                    )
                  ),
                  # (Shortnames for File inputs)
                  h3("Short_name for your GWAS data"),
                  fluidRow(
                    column(width = 6,
                           textInput("shortname1_input", "Short_name for First GWAS", value = input_shortname1, placeholder = "Enter a short_name for GWAS#1")
                    ),
                    column(width = 6,
                           textInput("shortname2_input", "Short_name for Second GWAS", value = input_shortname2, placeholder = "Enter a short_name for GWAS#1")
                    )
                  )
                ),

                tags$hr(),

                # Group related inputs together (GWAS HEADERS NAMES)
                # GWAS#1
                wellPanel(
                  h3("Definition of GWAS Header File Names - First GWAS"),
                  fluidRow(
                    column(width = 4,
                           textInput("g1_rsid_input", "GWAS#1 rsid", value = g1_rsid, placeholder = "Enter the rsid column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_chr_input", "GWAS#1 chromosome", value = g1_chr, placeholder = "Enter the chromosome column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_pos_input", "GWAS#1 position", value = g1_pos, placeholder = "Enter the position column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_EA_input", "GWAS#1 Ref. Allele", value = g1_EA, placeholder = "Enter the Reference Allele column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_nE_A_input", "GWAS#1 Alt. Allele", value = g1_nE_A, placeholder = "Enter the Alternate Allele column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_beta_input", "GWAS#1 Beta", value = g1_beta, placeholder = "Enter the Beta column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_se_input", "GWAS#1 SE", value = g1_se, placeholder = "Enter the SE column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_pval_input", "GWAS#1 P-value", value = g1_pval, placeholder = "Enter the P-value column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_info_input", "GWAS#1 Info", value = g1_info, placeholder = "Enter the Info column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_EAF_input", "GWAS#1 EAF", value = g1_EAF, placeholder = "Enter the EAF column header name for GWAS#1")
                    ),
                    column(width = 4,
                           textInput("g1_MAF_input", "GWAS#1 MAF", value = g1_MAF, placeholder = "Enter the MAF column header name for GWAS#1")
                    )
                  )
                ),

                tags$hr(),

                # GWAS#2
                wellPanel(
                  h3("Definition of GWAS Header File Names - Second GWAS"),
                  fluidRow(
                    column(width = 4,
                           textInput("g2_rsid_input", "GWAS#2 rsid", value = g2_rsid, placeholder = "Enter the rsid column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_chr_input", "GWAS#2 chromosome", value = g2_chr, placeholder = "Enter the chromosome column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_pos_input", "GWAS#2 position", value = g2_pos, placeholder = "Enter the position column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_EA_input", "GWAS#2 Ref. Allele", value = g2_EA, placeholder = "Enter the Reference Allele column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_nE_A_input", "GWAS#2 Alt. Allele", value = g2_nE_A, placeholder = "Enter the Alternate Allele column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_beta_input", "GWAS#2 Beta", value = g2_beta, placeholder = "Enter the Beta column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_se_input", "GWAS#2 SE", value = g2_se, placeholder = "Enter the SE column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_pval_input", "GWAS#2 P-value", value = g2_pval, placeholder = "Enter the P-value column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_info_input", "GWAS#2 Info", value = g2_info, placeholder = "Enter the Info column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_EAF_input", "GWAS#2 EAF", value = g2_EAF, placeholder = "Enter the EAF column header name for GWAS#2")
                    ),
                    column(width = 4,
                           textInput("g2_MAF_input", "GWAS#2 MAF", value = g2_MAF, placeholder = "Enter the MAF column header name for GWAS#2")
                    )
                  )
                ),


                tags$hr(),

                # Group related inputs together (Annotation File inputs)
                wellPanel(
                  h3("Select Annotation Files Names"),
                  fluidRow(
                    column(width = 6,
                           fileInput("annot_file1_input", "Choose Anotation file")
                    ),
                    column(width = 6,
                           fileInput("annot_file2_input", "Choose GWAS-annotated file")
                    )
                  )
                ),

                tags$hr(),

                # Group related inputs together (Calculation Parameters)
                wellPanel(
                  h3("Calculation Parameters"),
                  fluidRow(
                    column(width = 4,
                           numericInput("info_threshold_input", "Info Threshold", value = info_threshold)
                    ),
                    column(width = 4,
                           numericInput("MAF_threshold_input", "MAF Threshold", value = MAF_threshold)
                    ),
                    column(width = 4,
                           numericInput("group_clump_threshold_input", "Number of SNPs Threshold", value = group_clump_threshold)
                    ),
                    column(width = 4,
                           numericInput("theta_exploration_input", "Pleiotropy Threshold", value = theta_exploration)
                    )
                  )
                ),

                tags$hr(),

                # Checkbox input widget for toclump parameter
                wellPanel(
                  h3("LD Clumping Step"),
                  checkboxInput("toclump_input", "Running LD Clumping", value = toclump),
                  # section for the Conditional Panel
                  uiOutput("condPanel")
                ),

                tags$hr(),

                # Submit button
                wellPanel(
                  h3("Save the Changes"),
                  actionButton("submit_button", "Update Parameters File")
                ),

                tags$hr(),

                # checking a list of packages required for GCPBayes Pipeline
                wellPanel(
                  h3("Checking Required Packages for GCPBayes Pipeline"),
                  actionButton("check_button", "Check Required Packages")
                ),

                tags$hr(),

                # Action button to run the R script
                wellPanel(
                  h3("Running GCPBayes Pipeline"),
                  actionButton("run_button", "Run GCPBayes Pipeline")
                ),

                tags$hr(),

                # Output to display current parameter values
                verbatimTextOutput("current_values"),

                # Add a separator
                tags$hr(),

                # # Add an output to display current parameter values
                # fluidRow(
                #   column(width = 12,
                #          wellPanel(
                #            h3("Current Parameters"),
                #            verbatimTextOutput("current_values")
                #          )
                #   )
                # )

                tags$h6("GitHub Page: https://github.com/CESP-ExpHer/GCPBayes-Pipeline", style = "font-style: italic;")

)

# ================================================================
# Define server
server <- function(input, output, session) {

  # Define reactive values for each parameter
  work_dir_reactive <- reactiveValues(value = work_dir)
  input_reactive1 <- reactiveValues(value = input1)
  input_reactive2 <- reactiveValues(value = input2)
  input_shortname_reactive1 <- reactiveValues(value = input_shortname1)
  input_shortname_reactive2 <- reactiveValues(value = input_shortname2)
  g1_rsid_reactive <- reactiveValues(value = g1_rsid)
  g1_chr_reactive <- reactiveValues(value = g1_chr)
  g1_pos_reactive <- reactiveValues(value = g1_pos)
  g1_EA_reactive <- reactiveValues(value = g1_EA)
  g1_nE_A_reactive <- reactiveValues(value = g1_nE_A)
  g1_beta_reactive <- reactiveValues(value = g1_beta)
  g1_se_reactive <- reactiveValues(value = g1_se)
  g1_pval_reactive <- reactiveValues(value = g1_pval)
  g1_info_reactive <- reactiveValues(value = g1_info)
  g1_EAF_reactive <- reactiveValues(value = g1_EAF)
  g1_MAF_reactive <- reactiveValues(value = g1_MAF)
  g2_rsid_reactive <- reactiveValues(value = g2_rsid)
  g2_chr_reactive <- reactiveValues(value = g2_chr)
  g2_pos_reactive <- reactiveValues(value = g2_pos)
  g2_EA_reactive <- reactiveValues(value = g2_EA)
  g2_nE_A_reactive <- reactiveValues(value = g2_nE_A)
  g2_beta_reactive <- reactiveValues(value = g2_beta)
  g2_se_reactive <- reactiveValues(value = g2_se)
  g2_pval_reactive <- reactiveValues(value = g2_pval)
  g2_info_reactive <- reactiveValues(value = g2_info)
  g2_EAF_reactive <- reactiveValues(value = g2_EAF)
  g2_MAF_reactive <- reactiveValues(value = g2_MAF)
  input_file_annot_reactive <- reactiveValues(value = input_annot_file1)
  input_file_gwas_annot_reactive <- reactiveValues(value = input_annot_file2)
  info_threshold_reactive <- reactiveValues(value = info_threshold)
  MAF_threshold_reactive <- reactiveValues(value = MAF_threshold)
  group_clump_threshold_reactive <- reactiveValues(value = group_clump_threshold)
  theta_exploration_reactive <- reactiveValues(value = theta_exploration)
  toclump_reactive <- reactiveValues(value = toclump)
  input_ref_path_reactive <- reactiveValues(value = input_ref_path_file)
  input_ref_name_reactive <- reactiveValues(value = input_ref_name_file)
  clump_threshold_r2_reactive <- reactiveValues(value = clump_threshold_r2)
  clump_threshold_kb_reactive <- reactiveValues(value = clump_threshold_kb)
  clump_threshold_p_reactive <- reactiveValues(value = clump_threshold_p)
  placo_pval_threshold_reactive <- reactiveValues(value = placo_pval_threshold)

  # for converting Windows Path to a right Path (for "work_dir" and "ref_path_b_files")
  observeEvent(input$work_dir_input, {
    if (grepl("\\\\", input$work_dir_input)) {
      work_dir_reactive$value <- gsub("\\\\", "/", input$work_dir_input)
    } else {
      work_dir_reactive$value <- input$work_dir_input
    }
  })
  observeEvent(input$ref_path_file_input, {
    if (grepl("\\\\", input$ref_path_file_input)) {
      input_ref_path_reactive$value <- gsub("\\\\", "/", input$ref_path_file_input)
    } else {
      input_ref_path_reactive$value <- input$ref_path_file_input
    }
  })

  # GWAS input files
  observeEvent(input$file1_input, {
    input_reactive1$value <- input$file1_input
  })
  observeEvent(input$file2_input, {
    input_reactive2$value <- input$file2_input
  })

  # GWAS input shortnames
  observeEvent(input$shortname1_input, {
    input_shortname_reactive1$value <- input$shortname1_input
  })
  observeEvent(input$shortname2_input, {
    input_shortname_reactive2$value <- input$shortname2_input
  })

  # GWAS input headers names
  observeEvent(input$g1_rsid_input, {
    g1_rsid_reactive$value <- input$g1_rsid_input
  })
  observeEvent(input$g1_chr_input, {
    g1_chr_reactive$value <- input$g1_chr_input
  })
  observeEvent(input$g1_pos_input, {
    g1_pos_reactive$value <- input$g1_pos_input
  })
  observeEvent(input$g1_EA_input, {
    g1_EA_reactive$value <- input$g1_EA_input
  })
  observeEvent(input$g1_nE_A_input, {
    g1_nE_A_reactive$value <- input$g1_nE_A_input
  })
  observeEvent(input$g1_beta_input, {
    g1_beta_reactive$value <- input$g1_beta_input
  })
  observeEvent(input$g1_se_input, {
    g1_se_reactive$value <- input$g1_se_input
  })
  observeEvent(input$g1_pval_input, {
    g1_pval_reactive$value <- input$g1_pval_input
  })
  observeEvent(input$g1_info_input, {
    g1_info_reactive$value <- input$g1_info_input
  })
  observeEvent(input$g1_EAF_input, {
    g1_EAF_reactive$value <- input$g1_EAF_input
  })
  observeEvent(input$g1_MAF_input, {
    g1_MAF_reactive$value <- input$g1_MAF_input
  })
  observeEvent(input$g2_rsid_input, {
    g2_rsid_reactive$value <- input$g2_rsid_input
  })
  observeEvent(input$g2_chr_input, {
    g2_chr_reactive$value <- input$g2_chr_input
  })
  observeEvent(input$g2_pos_input, {
    g2_pos_reactive$value <- input$g2_pos_input
  })
  observeEvent(input$g2_EA_input, {
    g2_EA_reactive$value <- input$g2_EA_input
  })
  observeEvent(input$g2_nE_A_input, {
    g2_nE_A_reactive$value <- input$g2_nE_A_input
  })
  observeEvent(input$g2_beta_input, {
    g2_beta_reactive$value <- input$g2_beta_input
  })
  observeEvent(input$g2_se_input, {
    g2_se_reactive$value <- input$g2_se_input
  })
  observeEvent(input$g2_pval_input, {
    g2_pval_reactive$value <- input$g2_pval_input
  })
  observeEvent(input$g2_info_input, {
    g2_info_reactive$value <- input$g2_info_input
  })
  observeEvent(input$g2_EAF_input, {
    g2_EAF_reactive$value <- input$g2_EAF_input
  })
  observeEvent(input$g2_MAF_input, {
    g2_MAF_reactive$value <- input$g2_MAF_input
  })

  # input for Annotation
  observeEvent(input$annot_file1_input, {
    input_file_annot_reactive$value <- input$annot_file1_input
  })
  observeEvent(input$annot_file2_input, {
    input_file_gwas_annot_reactive$value <- input$annot_file2_input
  })

  # input for Calculation Parameters
  observeEvent(input$info_threshold_input, {
    info_threshold_reactive$value <- input$info_threshold_input
  })
  observeEvent(input$MAF_threshold_input, {
    MAF_threshold_reactive$value <- input$MAF_threshold_input
  })
  observeEvent(input$group_clump_threshold_input, {
    group_clump_threshold_reactive$value <- input$group_clump_threshold_input
  })
  observeEvent(input$theta_exploration_input, {
    theta_exploration_reactive$value <- input$theta_exploration_input
  })

  # Boolean inputs for "LD Clumping"
  observeEvent(input$toclump_input, {
    toclump_reactive$value <- input$toclump_input
  })

  # inputs if "LD Clumping = TRUE"
  # observeEvent(input$ref_path_file_input, {
  #   input_ref_path_reactive$value <- input$ref_path_file_input
  # })
  observeEvent(input$ref_b_file_input, {
    input_ref_name_reactive$value <- input$ref_b_file_input
  })
  observeEvent(input$clump_threshold_r2_input, {
    clump_threshold_r2_reactive$value <- input$clump_threshold_r2_input
  })
  observeEvent(input$clump_threshold_kb_input, {
    clump_threshold_kb_reactive$value <- input$clump_threshold_kb_input
  })
  observeEvent(input$clump_threshold_p_input, {
    clump_threshold_p_reactive$value <- input$clump_threshold_p_input
  })
  observeEvent(input$placo_pval_threshold_input, {
    placo_pval_threshold_reactive$value <- input$placo_pval_threshold_input
  })

  # Render the fluidRow when toclump_input is TRUE
  output$condPanel <- renderUI({
    if (input$toclump_input) {
      fluidRow(
        column(width = 6,
               textInput("ref_path_file_input", "Reference File Directory", value = input_ref_path_file)
        ),
        column(width = 6,
               textInput("ref_b_file_input", "Reference b_files Names", value = input_ref_name_file)
        ),
        column(width = 6,
               numericInput("clump_threshold_r2_input", "LD Clumping r2 Threshold", value = clump_threshold_r2)
        ),
        column(width = 6,
               numericInput("clump_threshold_kb_input", "LD Clumping kb Threshold", value = clump_threshold_kb)
        ),
        column(width = 6,
               numericInput("clump_threshold_p_input", "LD Clumping P-value Threshold", value = clump_threshold_p)
        ),
        column(width = 6,
               numericInput("placo_pval_threshold_input", "Placo P-value Threshold", value = placo_pval_threshold)
        )
      )
    }
  })



  # Define action button
  observeEvent(input$submit_button, {

    # Update R parameter variables with reactive values
    work_dir <- work_dir_reactive$value
    input1 <- input_reactive1$value
    input2 <- input_reactive2$value
    input_shortname1 <- input_shortname_reactive1$value
    input_shortname2 <- input_shortname_reactive2$value
    g1_rsid <- g1_rsid_reactive$value
    g1_chr <- g1_chr_reactive$value
    g1_pos <- g1_pos_reactive$value
    g1_EA <- g1_EA_reactive$value
    g1_nE_A <- g1_nE_A_reactive$value
    g1_beta <- g1_beta_reactive$value
    g1_se <- g1_se_reactive$value
    g1_pval <- g1_pval_reactive$value
    g1_info <- g1_info_reactive$value
    g1_EAF <- g1_EAF_reactive$value
    g1_MAF <- g1_MAF_reactive$value
    g2_rsid <- g2_rsid_reactive$value
    g2_chr <- g2_chr_reactive$value
    g2_pos <- g2_pos_reactive$value
    g2_EA <- g2_EA_reactive$value
    g2_nE_A <- g2_nE_A_reactive$value
    g2_beta <- g2_beta_reactive$value
    g2_se <- g2_se_reactive$value
    g2_pval <- g2_pval_reactive$value
    g2_info <- g2_info_reactive$value
    g2_EAF <- g2_EAF_reactive$value
    g2_MAF <- g2_MAF_reactive$value
    input_annot_file1 <- input_file_annot_reactive$value
    input_annot_file2 <- input_file_gwas_annot_reactive$value
    info_threshold <- info_threshold_reactive$value
    MAF_threshold <- MAF_threshold_reactive$value
    group_clump_threshold <- group_clump_threshold_reactive$value
    theta_exploration <- theta_exploration_reactive$value
    toclump <- toclump_reactive$value
    input_ref_path_file <- input_ref_path_reactive$value
    input_ref_name_file <- input_ref_name_reactive$value
    clump_threshold_r2 <- clump_threshold_r2_reactive$value
    clump_threshold_kb <- clump_threshold_kb_reactive$value
    clump_threshold_p <- clump_threshold_p_reactive$value
    placo_pval_threshold <- placo_pval_threshold_reactive$value

    # Update relevant lines in parameters file
    lines <- readLines("GCPBayes_pipeline_parameters.R")

    lines[grep("^work_dir", lines)] <- paste0("work_dir <- \"", work_dir, "/", "\"")
    lines[grep("^input <-", lines)] <- paste0("input <- c(", "\"", as.character(input1)[1], "\", \"", as.character(input2)[1], "\")")
    lines[grep("^input_shortname <-", lines)] <- paste0("input_shortname <- c(", "\"", as.character(input_shortname1)[1], "\", \"", as.character(input_shortname2)[1], "\")")
    lines[grep("^g1_rsid", lines)] <- paste0("g1_rsid <- \"", g1_rsid, "\"")
    lines[grep("^g1_chr", lines)] <- paste0("g1_chr <- \"", g1_chr, "\"")
    lines[grep("^g1_pos", lines)] <- paste0("g1_pos <- \"", g1_pos, "\"")
    lines[grep("^g1_EA <-", lines)] <- paste0("g1_EA <- \"", g1_EA, "\"")
    lines[grep("^g1_nE_A", lines)] <- paste0("g1_nE_A <- \"", g1_nE_A, "\"")
    lines[grep("^g1_beta", lines)] <- paste0("g1_beta <- \"", g1_beta, "\"")
    lines[grep("^g1_se", lines)] <- paste0("g1_se <- \"", g1_se, "\"")
    lines[grep("^g1_pval", lines)] <- paste0("g1_pval <- \"", g1_pval, "\"")
    lines[grep("^g1_info", lines)] <- paste0("g1_info <- \"", g1_info, "\"")
    lines[grep("^g1_EAF <-", lines)] <- paste0("g1_EAF <- \"", g1_EAF, "\"")
    lines[grep("^g1_MAF", lines)] <- paste0("g1_MAF <- \"", g1_MAF, "\"")
    lines[grep("^g2_rsid", lines)] <- paste0("g2_rsid <- \"", g2_rsid, "\"")
    lines[grep("^g2_chr", lines)] <- paste0("g2_chr <- \"", g2_chr, "\"")
    lines[grep("^g2_pos", lines)] <- paste0("g2_pos <- \"", g2_pos, "\"")
    lines[grep("^g2_EA <-", lines)] <- paste0("g2_EA <- \"", g2_EA, "\"")
    lines[grep("^g2_nE_A", lines)] <- paste0("g2_nE_A <- \"", g2_nE_A, "\"")
    lines[grep("^g2_beta", lines)] <- paste0("g2_beta <- \"", g2_beta, "\"")
    lines[grep("^g2_se", lines)] <- paste0("g2_se <- \"", g2_se, "\"")
    lines[grep("^g2_pval", lines)] <- paste0("g2_pval <- \"", g2_pval, "\"")
    lines[grep("^g2_info", lines)] <- paste0("g2_info <- \"", g2_info, "\"")
    lines[grep("^g2_EAF <-", lines)] <- paste0("g2_EAF <- \"", g2_EAF, "\"")
    lines[grep("^g2_MAF", lines)] <- paste0("g2_MAF <- \"", g2_MAF, "\"")
    lines[grep("^file_annot", lines)] <- paste0("file_annot <- \"", input_annot_file1, "\"")
    lines[grep("^file_gwas_annot", lines)] <- paste0("file_gwas_annot <- \"", input_annot_file2, "\"")
    lines[grep("^info_threshold", lines)] <- paste0("info_threshold <- ", info_threshold)
    lines[grep("^MAF_threshold", lines)] <- paste0("MAF_threshold <- ", MAF_threshold)
    lines[grep("^group_clump_threshold", lines)] <- paste0("group_clump_threshold <- ", group_clump_threshold)
    lines[grep("^theta_exploration", lines)] <- paste0("theta_exploration <- ", theta_exploration)
    lines[grep("^toclump", lines)] <- paste0("toclump <- ", toclump)
    lines[grep("^ref_path_b_files", lines)] <- paste0("ref_path_b_files <- \"", input_ref_path_file, "/", "\"")
    lines[grep("^ref_b_files", lines)] <- paste0("ref_b_files <- \"", input_ref_name_file, "\"")
    lines[grep("^clump_threshold_r2", lines)] <- paste0("clump_threshold_r2 <- ", clump_threshold_r2)
    lines[grep("^clump_threshold_kb", lines)] <- paste0("clump_threshold_kb <- ", clump_threshold_kb)
    lines[grep("^clump_threshold_p", lines)] <- paste0("clump_threshold_p <- ", clump_threshold_p)
    lines[grep("^placo_pval_threshold", lines)] <- paste0("placo_pval_threshold <- ", placo_pval_threshold)

    # Write updated parameters file
    writeLines(lines, "GCPBayes_pipeline_parameters.R")

    # Update current parameter values output
    output$current_values <- renderPrint({
      paste("Work_Directory: ", work_dir,
            "Input_Files: ", as.character(input1)[1], ", ", as.character(input2)[1],
            "Input_Shortname: ", as.character(input_shortname1)[1], ", ", as.character(input_shortname2)[1],
            "GWAS1_rsid: ", g1_rsid,
            "Info_Threshold: ", info_threshold,
            "Running_LD_Clumping: ", ifelse(toclump, "Yes", "No"))
    })

  })


  # Define the action to take when the check button is clicked
  observeEvent(input$check_button, {
    # Call the source() function to run the R script
    source("GCPBayes_pipeline_check_packages.R")
  })


  # Define the action to take when the Run Script button is clicked
  observeEvent(input$run_button, {
    # Call the source() function to run the R script
    source("GCPBayes_pipeline.R")
  })


  # Set maximum file size (to prevent receiving "Maximum upload size exceeded" message)
  options(shiny.maxRequestSize = 10000 * 1024^2)

}

# Run app
shinyApp(ui, server)
