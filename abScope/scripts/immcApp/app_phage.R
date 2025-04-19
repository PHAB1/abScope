# Load libraries
suppressWarnings(library(shiny))
suppressWarnings(library(alakazam))
suppressWarnings(library(shazam))
suppressWarnings(library(dplyr))
suppressWarnings(library(scales))
suppressWarnings(library(ggplot2))
#suppressWarnings(library(dowser))
suppressWarnings(library(igraph))
suppressWarnings(library(tidyr))
suppressWarnings(library("Rtsne"))
suppressWarnings(library('stringdist'))
suppressWarnings(library(plotly))
suppressWarnings(library(ggfortify))
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(rgl))
suppressWarnings(library(DT))
suppressWarnings(library(factoextra))
suppressWarnings(library(tidyverse))
suppressWarnings(library(reshape2))
suppressWarnings(library(htmlwidgets))
suppressWarnings(library(svglite))
suppressWarnings(library(plotrix))
suppressWarnings(library(plot3D))
suppressWarnings(library(MHTdiscrete))

# Function to summarize data, # http://sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  #library(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]])/sqrt(length(x) ))
  }
  #data_sum<-ddply(data, groupnames, .fun=summary_func,
  #                varname)
  #print(data_sum)
  data_sum <- data %>% 
    group_by(subject_id, gene) %>%
    summarize(sd=sd(seq_freq, na.rm = TRUE), se=std.error(seq_freq, na.rm=TRUE), seq_freq=mean(seq_freq, na.rm = TRUE) )

  #data_sum <- rename(data_sum, c("mean" = "seq_freq"))
  return(data_sum)
}


# colnames
#used_cols <- c("locus", "v_call","c_call","d_call","j_call", "fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4","v_score","v_identity","duplicate_count","clone_id","subject_id","cprimer","seq_len","sample_id","clone_size_count","clone_size_freq")
used_cols <- c("junction_aa","locus", "v_call","j_call","duplicate_count","clone_id","subject_id","seq_len","sample_id","clone_size_count","clone_size_freq")

# Functions
df2stat <- function(df, colSamples, colValues, colType, controlPattern, TratPattern, colSubType=FALSE) {
  df <- select(df, colSamples, colValues, colType)
  df$seq_freq <- df$seq_freq*1000000
  df <- df %>% 
    spread(key=colSamples, value=colValues) %>%
    mutate_all(~ifelse(is.na(.), 1, .))

  df$log2FoldChange <- log2((meanByPatter(df,TratPattern) / (meanByPatter(df,controlPattern))))
  
  #print(shapiro.test(c_across(contains(TratPattern))))
  #print(shapiro.test(c_across(contains(controlPattern))))
  df <- df %>%
    rowwise() %>%
    mutate(p_value = wilcox.test(c_across(contains(TratPattern)), c_across(contains(controlPattern)),paired=F)$p.value)
    #mutate(p_value = t.test(c_across(contains(TratPattern)), c_across(contains(controlPattern)))$p.value)
    
  #df$p_adjusted <- p.adjust(df$p_value, method ="BH")
  df$p_adjusted <- Sidak.p.adjust(df$p_value)

  tryCatch({ 
    df <- select(df, colType, log2FoldChange, p_value, p_adjusted, everything())
  }, error = function(e) {
    df <- select(df, colType, log2FoldChange, everything())
  })
}

# Statistical analysis ----
# Sum pattern column                                    
meanByPatter <- function(df, pattern) {                          
  combined_pattern <- paste(pattern, collapse = "|")
  return(rowMeans(df[, grep(combined_pattern, colnames(df))]))
}

# create colnames from column
coll_names <- function(column) {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  name2collor <- data.frame(subject_id=unique(column), collors=col_vector[0:length(unique(column))])
  return(merge(data.frame(subject_id=column), name2collor))
}

options(shiny.maxRequestSize=7000*1024^2)
# Define UI for data upload app ----
ui <- fluidPage(

  # App title ----
  titlePanel(""),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      conditionalPanel(
        'input.dataset == "Table"',
        # Input: Select a file ----
        fileInput("file1", "Choose TSV File",
          multiple = TRUE,
          accept = c("text/csv",
            "text/comma-separated-values,text/plain",
            ".csv", 
            ".tsv")),

        # Horizontal line ----
        tags$hr(),

        # Input: Select number of rows to display ----
        radioButtons("disp", "Display",
          choices = c(Head = "head",
                      All = "all"),
          selected = "head"),
        uiOutput("checkboxes_subject"),
        uiOutput("checkboxes_sample"),
        textAreaInput("junction_aa_input", "Select The Clones To Remove (Separated by space):")
      ),

      conditionalPanel(
        'input.dataset == "Cluster Analysis"',
        # Input: Select subject_id or sample_id ----
        radioButtons("colorByCluster", "Collor by:",
          choices = c(sample_id = "sample_id",
            subject_id = "subject_id"),
            selected = "subject_id"),

        # Input: Simple integer interval ----     
        sliderInput("perplexity", "Perplexity:",
          min = 0, max = 100,
          value = 15)
      ),
      conditionalPanel(
	'input.dataset == "Net Analysis"',
        #radioButtons("NetType", "Network By Type:",
        #  choices = c(sample_id = "sample_id",
        #    subject_id = "subject_id"),
        #    selected = "subject_id"),
		
        # Input: Simple integer interval ----     
        sliderInput("nSample", "Sampling Number:",
          min = 300, max = 1000,
          value = 1000),
	uiOutput("radio_subject_net")
	
      ),
      conditionalPanel(
        'input.dataset == "Diversity"',
        # Input: Select subject_id or sample_id ----
        radioButtons("plotBy", "Plot by:",
          choices = c(sample_id = "sample_id",
              subject_id = "subject_id"),
              selected = "subject_id"),

        # Input: Simple integer interval ----     
        sliderInput("q", "q value:",
                  min = 0, max = 2,
                  value = 0) 
      ),
      conditionalPanel(
        'input.dataset == "Statistics"',
        uiOutput("ctrl"),
        uiOutput("trat"),

        # Input: Select V, D or J ----
        radioButtons("CVDJ_type", "Choose C, V, D or J Type:",
          choices = c(v_call = "v_call",
          j_call = "j_call",
          CDR3 = "junction_aa"),
          selected = "v_call"), 

        # Input: Mode Family or Gene ----
        radioButtons("modeFG", "Choose Family or Gene Mode:",
                choices = c(family = "family",
                    gene = "gene"),
                    selected = "family")
	    ),
      conditionalPanel(
        'input.dataset == "VDJ Usage"',
          # C, V, D or J Plot
          radioButtons("CVDJ_type_usage", "Choose C, V, D or J Type:",
            choices = c(v_call = "v_call",
            j_call = "j_call"),
            selected = "v_call"),

          # Input: Mode Family or Gene ----
          radioButtons("modeFG_usage", "Choose Family or Gene Mode:",
            choices = c(family = "family",
            gene = "gene"),
            selected = "family"),

          # Input: Plot By ----
          radioButtons("plotBy_usage", "Plot by:",
          choices = c(sample_id = "sample_id",
              subject_id = "subject_id"),
              selected = "subject_id"),
              
          textAreaInput("txt_cvdj", "Filter C-V-D or J Type (only one type):")
          
      ),
      conditionalPanel(
        'input.dataset == "DB Search"',
        # Input: Select a file ----
        fileInput("file2", "Choose CSV File",
          multiple = TRUE,
          accept = c("text/csv",
            "text/comma-separated-values,text/plain",
            ".csv", 
            ".tsv")),

        # Horizontal line ----
        tags$hr(),

        # Input: Select number of rows to display ----
        radioButtons("disp2", "Display",
          choices = c(Head = "head",
                      All = "all"),
          selected = "head"),

	uiOutput("cdr3_field")
	),
	
      conditionalPanel(
        'input.dataset == "Ab Length"',
          # Input: Plot By ----
          radioButtons("plotBy_usage_AbLen", "Plot by:",
          choices = c(sample_id = "sample_id",
              subject_id = "subject_id"),
              selected = "subject_id"),

          # Input: Simple integer interval to remove seqs by sample_id or subject_id ----     
          sliderInput("min_group_perCent", "Min percentage c_call (%) per sample_id:",
              min = 0, max = 20,
              value = 3,
	      step = 0.5) 
	  
      )
  ),
  # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(
	      id = 'dataset',
        tabPanel("Table",
          # Output: Data file ----
          DT::dataTableOutput("contents"),

          # Horizontal line ----
          tags$hr(),
        ),
	tabPanel("Cluster Analysis",
          rglwidgetOutput("PCA", width = "512px", height = "512px"),

	  # Horizontal line ----
	  tags$hr(),

	  plotOutput("pca_eig"),

	  # Horizontal line ----
	  tags$hr(),

	  plotlyOutput("tSNE_subject", width = 'auto', height = 'auto'),
	  
	),	
	tabPanel("Net Analysis",
	  plotOutput("igraphPlot"),

	  # Horizontal line ----
	  tags$hr(),
          

          fluidRow(
            column(width = 3,
              numericInput("widthNET", "width (Inches):", value = 8)),
            column(width = 3,
              numericInput("heightNET", "height (Inches):", value = 8)),
            column(width = 3,
              downloadButton("downloadNETPlot", "Download"))
          ),

	),
	tabPanel("Diversity",
          plotOutput("abundance"),

          fluidRow(
            column(width = 3,
              numericInput("widthAbuCur", "width (Inches):", value = 20)),
            column(width = 3,
              numericInput("heightAbuCur", "height (Inches):", value = 10)),
            column(width = 3,
              numericInput("dpiAbuCur", "DPI:", value = 300)),
            column(width = 3,
              selectInput("deviceAbuCur", "File format:", choices = c("png", "jpeg", "pdf", "svg"), selected = "png")),
            column(width = 3,
              downloadButton("downloadAbund", "Download"))
          ),

	  # Horizontal line ----
	  tags$hr(),

	  plotOutput("diversity"),
          
          fluidRow(
            column(width = 3,
              numericInput("widthDivCur", "width (Inches):", value = 20)),
            column(width = 3,
              numericInput("heightDivCur", "height (Inches):", value = 10)),
            column(width = 3,
              numericInput("dpiDivCur", "DPI:", value = 300)),
            column(width = 3,
              selectInput("deviceDivCur", "File format:", choices = c("png", "jpeg", "pdf", "svg"), selected = "png")),
            column(width = 3,
              downloadButton("downloadPlotCurve", "Download"))
          ),

	  # Horizontal line ----
	  tags$hr(),

	  plotOutput("diversity_barplotByQ"),

	  # Horizontal line ----
	  #tags$hr(),

          fluidRow(
            column(width = 3,
              numericInput("width", "width (Inches):", value = 20)),
            column(width = 3,
              numericInput("height", "height (Inches):", value = 10)),
            column(width = 3,
              numericInput("dpi", "DPI:", value = 300)),
            column(width = 3,
              selectInput("device", "File format:", choices = c("png", "jpeg", "pdf", "svg"), selected = "png")),
            column(width = 3,
              downloadButton("downloadPlot", "Download"))
          ),

	  #downloadPlotCurve

	  # Horizontal line ----
	  tags$hr(),

	  plotOutput("diversityByQ"),

          fluidRow(
            column(width = 3,
              numericInput("widthByQ", "width (Inches):", value = 20)),
            column(width = 3,
              numericInput("heightByQ", "height (Inches):", value = 10)),
            column(width = 3,
              numericInput("dpiByQ", "DPI:", value = 300)),
            column(width = 3,
              selectInput("deviceByQ", "File format:", choices = c("png", "jpeg", "pdf", "svg"), selected = "png")),
            column(width = 3,
              downloadButton("downloadPlotByQ", "Download"))
          ),
	),
	tabPanel("Statistics",
    # V, D and J Table, all vs all
    DT::dataTableOutput("vdj_stat"),
          
	  # Horizontal line ----
	  tags$hr(),

    # V, D and J Table, Combinatorial analysis taken two by two
    DT::dataTableOutput("vdj_comb_stat")
	),
	tabPanel("VDJ Usage",
    	  plotlyOutput("vdjBarplot"),


          fluidRow(
            column(width = 3,
              numericInput("widthVDJ", "width (Inches):", value = 20)),
            column(width = 3,
              numericInput("heightVDJ", "height (Inches):", value = 10)),
            column(width = 3,
              numericInput("dpiVDJ", "DPI:", value = 300)),
            column(width = 3,
              selectInput("deviceVDJ", "File format:", choices = c("png", "jpeg", "pdf", "svg"), selected = "png")),
            column(width = 3,
              downloadButton("downloadVDJPlot", "Download"))
          ),

	  # Horizontal line ----
    	  plotlyOutput("vdjSpecBarplot"),
	  tags$hr()
	),
	
	tabPanel("Ab Length",
	  # Horizontal line ----
	  tags$hr(),

	  plotOutput("cdr3_len_plot")
	),

	tabPanel("DB Search",
          # Output: Data file ----
          DT::dataTableOutput("inters_cdr3"),

          # Horizontal line ----
          tags$hr()
	)
      )
    )
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {

  df <- reactive({
    req(input$file1)

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    df <- NULL
    for(f in input$file1$datapath) {
      airr_file <- airr::read_rearrangement(f)
      airr_file$seqorient <- NULL
      df <- tryCatch(
        expr = {
        df <- rbind(df,airr_file)
      },
      error = function(e) {
        df <- airr_file
      }
      )
    }
    df <- df %>%
      filter(!is.na(v_call), !is.na(j_call), !is.na(junction_length)) %>%
      mutate(v_call = sub(",.*","",v_call)) %>%
      mutate(j_call = sub(",.*","",j_call)) %>%
      mutate(v_call = sub("\\*.*","",v_call)) %>%   
      mutate(j_call = sub("\\*.*","",j_call))
  })

  output$contents <- DT::renderDataTable(DT::datatable({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    
    df <- df %>%
        select(all_of(used_cols))
  }))

  # Renderiza os checkbox buttons dinamicamente
  output$checkboxes_subject <- renderUI({
    unique_subjects <- unique(df()$subject_id)
    checkboxGroupInput("subjects", "Choose Subject ids", choices=unique_subjects, selected = unique_subjects)
  })

  # Renderiza os checkbox buttons dinamicamente
  output$checkboxes_sample <- renderUI({
    selected_subjects <- input$subjects
    df <- subset(df(), (subject_id %in% selected_subjects))

    unique_samples <- unique(df$sample_id)
    checkboxGroupInput("samples", "Choose Sample ids", choices=unique_samples, selected = unique_samples)
  })

  # Renderiza os checkbox buttons dinamicamente
  output$radio_subject_net <- renderUI({
    #selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    #selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples, ]
    
    unique_subjects <- unique(df$subject_id)
    radioButtons("subjects_radio", "Choose Subject ids", choices=unique_subjects, selected = unique_subjects[1])
  })

  output$pca_eig <- renderPlot({
      selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
      selected_junctions <- trimws(selected_junctions)
      
      df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]

      # PCA samples plot
      df <- df %>%
        select(sample_id, j_call, v_call) %>% 
        mutate(cjdv_call = paste0(j_call,, v_call)) %>% 
        select(sample_id, cjdv_call) %>% 
        table(.) %>% 
        as.data.frame(.) %>% 
        mutate(Freq = Freq * 1000000 / sum(.$sample_id == sample_id)) %>%
        spread(key="sample_id",value="Freq")

      df <- as.data.frame(df)
      rownames(df) <- df$cjdv_call
      df$cjdv_call <- NULL
      df <- t(df)

      #plot
      df<-df[ , which(apply(df, 2, var) != 0)]

      pca_res <- prcomp(df, scale. = FALSE)
      fviz_eig(pca_res) 
  })

  output$PCA <- renderRglwidget({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]

    df_tags <- df %>%
      select(sample_id,subject_id)
    df_tags <- as.data.frame(unique(df_tags))

    # PCA samples plot
    df <- df %>%
      select(sample_id, j_call, v_call) %>% 
      mutate(cjdv_call = paste0(j_call, v_call)) %>% 
      select(sample_id, cjdv_call) %>% 
      table(.) %>% 
      as.data.frame(.) %>% 
      mutate(Freq = Freq * 1000000 / sum(.$sample_id == sample_id)) %>%
      spread(key="sample_id",value="Freq")

    df <- as.data.frame(df)
    rownames(df) <- df$cjdv_call
    df$cjdv_call <- NULL
    df <- t(df)

    #plot
    df<-df[ , which(apply(df, 2, var) != 0)]
    t <- as.data.frame(df)
    t <- t[order(rownames(t)),]

    pca_res <- prcomp(df, scale. = FALSE)
    
    colorBy <- input$colorByCluster
    colors_names <- coll_names(df_tags$subject_id)
    df_tags <- unique(merge(df_tags,colors_names,by="subject_id"))
    df_tags <- df_tags[order(df_tags$sample_id),]
    eig_perc <- (get_eig(pca_res))$variance.percent %>% round(digits=2)
    
    open3d(useNULL=T)
    plot3d(as.data.frame(pca_res$x[,1:3]),size=8,radius=100,type="s",col=df_tags$collors, box = FALSE, xlab=paste0("PC1 ",eig_perc[1]),ylab=paste0("PC2 ",eig_perc[2]),zlab=paste0("PC3 ",eig_perc[3]))
    if(colorBy == "subject_id") {
      text3d(as.data.frame(pca_res$x[,1:3]),texts=c(df_tags$subject_id),cex= 0.7, pos=3)
    } else {
      text3d(as.data.frame(pca_res$x[,1:3]),texts=c(df_tags$sample_id),cex= 0.7, pos=3)
    }
    rglwidget()
  })
  
  output$tSNE_subject <- renderPlotly({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    # tSNE ----
    # Função para subamostrar 1000 linhas de cada tipo de subject_id
    subsample <- function(data, column, n) {
      data %>%
      group_by(!!sym(column)) %>%
      sample_n(n, replace = TRUE) %>%
      ungroup()
    }

    colorBy <- input$colorByCluster
    bcr.subset <- subsample(df, colorBy, 1000) # 1000 of each

    #bcr.subset <- df %>%
    #  group_by(subject_id) %>%
    #  sample_n(size = 1000, replace = TRUE)

    a <- bcr.subset$junction_aa
    b <- bcr.subset$junction_aa

    # Calculate leveshtein and hamming distance ----
    s.dist.hamming <- stringdistmatrix(a, b, method = "h", useNames = "string")
    s.dist.leveshtein <- stringdistmatrix(a, b, method = "lv", useNames = "string")
    s.dist <- cbind(s.dist.hamming, s.dist.leveshtein)
    s.dist[is.infinite(s.dist)] <- 0

    tsne_results <- Rtsne(as.matrix(s.dist), perplexity=input$perplexity, check_duplicates = FALSE, dims = 2)
    if(colorBy == "subject_id") {
      TsnePlot <<- plot_ly(x = tsne_results$Y[,1], y = tsne_results$Y[,2], type='scatter', color=as.factor(bcr.subset$subject_id), text = as.factor(bcr.subset$junction_aa), mode = "markers")
    } else {
      TsnePlot <<- plot_ly(x = tsne_results$Y[,1], y = tsne_results$Y[,2], type='scatter', color=as.factor(bcr.subset$sample_id), text = as.factor(bcr.subset$junction_aa), mode = "markers")
    }
    TsnePlot
  })


  output$igraphPlot <- renderPlot({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects_radio & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    # tSNE ----
    # Subset the first 1000 rows ----

    # Group by (sample or subject) and sampling by n
    bcr.subset <- df %>%
      group_by(subject_id) %>%
      sample_n(size = input$nSample, replace = FALSE)

    combined_plot <- list()

    a <- bcr.subset$junction_aa
    b <- bcr.subset$junction_aa

    s.dist <- stringdistmatrix(a, b, method = "h", useNames = "string")
    s.dist[s.dist < 2] <- 1
    s.dist[s.dist >= 2] <- 0

    # Create clonality graph 
    network <<- graph_from_adjacency_matrix(s.dist, mode = "undirected", weighted = T, diag = F)
    subNET <<- input$subjects_radio

    plot.igraph(network,vertex.size=1.5, vertex.color = "grey", main=input$subjects_radio, edge.color = "orange",vertex.label=NA,edge.arrow.size=.1,layout=layout.fruchterman.reingold(network, niter=1000))

    #combined_plot <- do.call("subplot", c(combined_plot, nrows = length(combined_plot)))
    #combined_plot
  })

  # Download Net plot
  output$downloadNETPlot <- downloadHandler(

    filename = function() {
            paste0("NETPlot.","pdf")
    },
    content = function(file) {
      # Salvar o gráfico (pdf) com base nos valores de entrada
      pdf(file, width=input$widthNET, height=input$heightNET)
      plot.igraph(network,vertex.size=1.5, vertex.color = "grey", main=subNET, edge.color = "orange",vertex.label=NA,edge.arrow.size=.1,layout=layout.fruchterman.reingold(network, niter=1000))
      dev.off()
    }
  )

  # Immcantation Repertoire Based Analisys # ----
  output$abundance <- renderPlot({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    curve <- estimateAbundance(df, group=input$plotBy, ci=0.95, nboot=100, clone="clone_id")
    #curve <- estimateAbundance(joined_airr, group="sample_id", ci=0.95, nboot=100, clone="clone_id")

    abuCurvePlot <<- plot(curve, legend_title="", main="Abundance Plot") + theme(plot.title=element_text(size = 20), axis.title = element_text(size = 30), axis.text.x=element_text(size=20),axis.text.y=element_text(size=20)) + theme_classic()
    abuCurvePlot
  })

  output$downloadAbund <- downloadHandler(

    filename = function() {
	    paste0("abundCurve.",input$deviceAbuCur)
    },
    content = function(file) {
      # Salvar o gráfico com base nos valores de entrada
      ggsave(file, abuCurvePlot, width = input$widthAbuCur, height = input$heightAbuCur, dpi = input$dpiAbuCur, device = input$deviceAbuCur)
    }
  )

  
  # Diversity Plot ----
  output$diversity <- renderPlot({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    
    # Diversity ----
    abund <- estimateAbundance(df, group=input$plotBy, nboot=100)
    div <- alphaDiversity(abund, max_q=6)
    divCurvePlot <<- plotDiversityCurve(div, legend_title="") + theme(plot.title=element_text(size = 20), axis.title = element_text(size = 30), axis.text.x=element_text(size=20),axis.text.y=element_text(size=20)) + theme_classic()
    divCurvePlot
  })

  output$downloadPlotCurve <- downloadHandler(

    filename = function() {
	    paste0("diversityCurve.",input$deviceDivCur)
    },
    content = function(file) {
      # Salvar o gráfico com base nos valores de entrada
      ggsave(file, divCurvePlot, width = input$widthDivCur, height = input$heightDivCur, dpi = input$dpiDivCur, device = input$deviceDivCur)
    }
  )

  # Diversity by q=n with error ----
  output$diversity_barplotByQ <- renderPlot({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    # Diversity   
    div <- alphaDiversity(df, group=input$plotBy, min_q=0, max_q=2, step_q=1, nboot=100)

    # Define plot values
    # Acessar os valores de diversidade
    df_div <- div@diversity %>%
      filter(q == q) %>%
      mutate(lower = d - d_sd, upper = d + d_sd)

    # Plotar o gráfico de barras com barras de erro
    diversity_plot <<- ggplot(data = df_div, aes_string(x = "q", y = "d", fill = input$plotBy)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08, position = position_dodge(0.9)) +
      labs(x = "q", y = "d", title = "Barplot Diversity") +
      theme_classic() +
      theme(
          axis.text = element_text(size = 16),
	  axis.title = element_text(size = 18),
	  plot.title = element_text(size = 20)
      )

    diversity_plot
  })

  # Diversity by q=n with error ----
  output$diversityByQ <- renderPlot({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    # Diversity   
    div <- alphaDiversity(df, group=input$plotBy, min_q=0, max_q=2, step_q=1, nboot=100)

    plotDivByQ <<- plotDiversityTest(div, input$q, legend_title="") + theme(plot.title=element_text(size = 20), axis.title = element_text(size = 30), axis.text.x=element_text(size=20),axis.text.y=element_text(size=20)) + theme_classic()
    plotDivByQ
  })

  output$downloadPlotByQ <- downloadHandler(

    filename = function() {
	    paste0("diversityBarplotByQ.",input$deviceByQ)
    },
    content = function(file) {
      # Salvar o gráfico com base nos valores de entrada
      ggsave(file, plotDivByQ, width = input$widthByQ, height = input$heightByQ, dpi = input$dpiByQ, device = input$deviceByQ)
    }

  )

  output$downloadPlot <- downloadHandler(

    filename = function() {
	    paste0("diversityBarplot.",input$device)
    },
    content = function(file) {
      # Salvar o gráfico com base nos valores de entrada
      ggsave(file, diversity_plot, width = input$width, height = input$height, dpi = input$dpi, device = input$device)
    }
  )

  # Statistics DataFrames ----
  output$vdj_comb_stat <- DT::renderDataTable(DT::datatable({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    df <- filter(df, subject_id %in% c(input$trat, input$ctrl))

    # combinatory analysis for C-V-D-J calls -> n!/(n-p)!p! = 4!/(4-2)!2! = 6 -> VD, VJ, VC, DJ, DC, JC
    combinatory_list <- list(
      list("j_call","v_call"),
    )

    finalSt_df <- NULL

    tratSamples <- unique(filter(df, subject_id==input$trat)$sample_id)
    controlSamples <- unique(filter(df, subject_id==input$ctrl)$sample_id)
    for(i in 1:length(combinatory_list)) {
      col_up <- unlist(combinatory_list[[i]][1])
      col_down <- unlist(combinatory_list[[i]][2])
      subset_vdj <- countGenes(df,gene=col_up, groups=c("sample_id",col_down),copy="duplicate_count", mode="gene")
  
      teste <- df2stat(subset_vdj, "sample_id", "seq_freq", "gene", controlSamples, tratSamples, col_down)

      colnames(teste)[1] <- "upStream"
      colnames(teste)[4] <- "downStream"
      teste <- select(teste, "upStream", "downStream", everything())
  
      finalSt_df <- tryCatch(
        expr = {
          finalSt_df <- rbind(finalSt_df,teste)
        },
        error = function(e) {
          finalSt_df <- teste
        }
      )
    }
    
    finalSt_df
  }))
  
  # Statistics DataFrames V, D or J ----
  output$vdj_stat <- DT::renderDataTable(DT::datatable({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    df <- filter(df, subject_id %in% c(input$trat, input$ctrl))

    tratSamples <- unique(filter(df, subject_id==input$trat)$sample_id)
    controlSamples <- unique(filter(df, subject_id==input$ctrl)$sample_id)

    family_gene <- countGenes(df, gene=input$CVDJ_type, groups=c("sample_id"), mode=input$modeFG)
    finalSt_df <- df2stat(family_gene, "sample_id", "seq_freq", "gene", controlSamples, tratSamples)

    finalSt_df
  }))

  output$enrich_df <- DT::renderDataTable(DT::datatable({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    df <- filter(df, subject_id %in% c(input$trat, input$ctrl))

  }))
  
  output$ctrl <- renderUI({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]

    subjects <- unique(df$subject_id)
    selectInput(
        'ctrl', 'Control', choices = subjects
    )
  })

  output$trat <- renderUI({
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples, ]
    subjects <- unique(df$subject_id)
    selectInput(
        'trat', 'Tratament', choices = subjects
    )
  })

  output$vdjBarplot <- renderPlotly({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    
    if(input$plotBy_usage=="sample_id") {
    	# Assign sorted levels and subset to all IGHV in subject_id/ sample_id
    	df <- countGenes(df, gene=input$CVDJ_type_usage, groups=input$plotBy_usage, mode=input$modeFG_usage)

    	df$seq_freq <- df$seq_freq*100
    	vdjBarplot_g1 <<- ggplot(df, aes_string(x = "gene", y = "seq_freq", fill = input$plotBy_usage)) +
      	  theme_bw() +
      	  ggtitle("IGH-V-D-J Usage") +
      	  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      	  ylab("Percent of repertoire") +
      	  xlab("") +
      	  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
      	  #scale_fill_brewer(palette = "Set1") +
      	  geom_bar(stat = "identity", position = "dodge", width = 0.7)
      
    	plotly::ggplotly(vdjBarplot_g1)
    }

    else {
    	df <- countGenes(df, gene=input$CVDJ_type_usage, groups=c("sample_id", "subject_id"), mode=input$modeFG_usage)
        
	df <- data_summary(df, varname="seq_freq", 
                    groupnames=c("subject_id","gene"))


    	df$seq_freq <- df$seq_freq*100
	df$sd <- df$sd*100
	df$se <- df$se*100
    	df$subject_id=as.factor(df$subject_id)
	vdjBarplot_g1 <<- ggplot(df, aes(x=gene, y=seq_freq, fill = subject_id)) + 
		geom_bar(stat="identity", color="black",
		    position=position_dodge()) +
		geom_errorbar(aes(ymin=seq_freq, ymax=seq_freq+se), width=.2,
		    position=position_dodge(.9)) + 
		theme_bw() +
                ggtitle("IGH-V-D-J Usage") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                ylab("Percent of repertoire") +
                xlab("") +
                scale_y_continuous(labels = scales::percent_format(scale = 1))
    }

  })

  output$downloadVDJPlot <- downloadHandler(

    filename = function() {
	    paste0("vdjBarplot_g1.",input$deviceVDJ)
    },
    content = function(file) {
      ggsave(file, vdjBarplot_g1, width = input$widthVDJ, height = input$heightVDJ, dpi = input$dpiVDJ, device = input$deviceVDJ)
    }
  )

  output$vdjSpecBarplot <- renderPlotly({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)

    # Selected CVD or J
    selected_cvdj <- strsplit(input$txt_cvdj, " ")[[1]]
    selected_cvdj <- trimws(selected_cvdj)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]

    for(i in selected_cvdj) {
      if(any(grepl(i, df$v_call))) {
        df <- subset(df, grepl(i, v_call))
      } else if(any(grepl(i, df$j_call))) {
        df <- subset(df, grepl(i, j_call))
      }
    }

    #df_count <- countGenes(df, gene="v_call", groups="subject_id", mode="gene")
    gene_count <- countGenes(df, gene = input$CVDJ_type_usage, groups = input$plotBy_usage, mode = input$modeFG_usage)
    gene_count$seq_freq <- gene_count$seq_freq*100

    # Assign sorted levels and subset to downstream
    #gene_count <- gene %>%
    #  mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
    #  filter(getFamily(gene) == "IGHV4")

    # Plot V gene usage in the IGHV4 family by sample
    if(length(selected_cvdj) != 0) {
      vdjSpecBarplot_plot <<- ggplot(gene_count, aes_string(x = "gene", y = "seq_freq", fill = input$plotBy_usage)) +
        theme_bw() +
        ggtitle(paste0("Filtered by ",selected_cvdj)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ylab("Percent of repertoire") +
        xlab("") +
        scale_y_continuous(labels = scales::percent_format(scale = 1)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7)
      plotly::ggplotly(vdjSpecBarplot_plot)
    } else {
      NULL
    }
  })

  # cdr3 len Plot ----
  output$cdr3_len_plot <- renderPlot({
    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)
    
    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ]
    db_props <- aminoAcidProperties(df, seq="junction", trim=TRUE, label="cdr3")
    
    pros_rm <- db_props %>% 
	        count(subject_id,sample_id,v_call) %>%
	        group_by(sample_id) %>%
		mutate(c_freq = n/sum(n)) %>%
		ungroup() %>%
		filter(c_freq > as.numeric(input$min_group_perCent/100))
	            
    pros_rm <- arrange(pros_rm, c_freq)
    db_props <- semi_join(db_props, pros_rm, by=c("sample_id", "v_call"))
    print(pros_rm)

    # The full set of properties are calculated by default
    dplyr::select(db_props[1:3, ], starts_with("cdr3"))

    # Define a ggplot theme for all plots
    tmp_theme <- theme_bw() + theme(legend.position="bottom")

    # Generate plots for all four of the properties
    g1 <- ggplot(db_props, aes(x=v_call, y=cdr3_aa_length)) + tmp_theme +
        ggtitle("CDR3 length") + 
        xlab("Isotype") + ylab("Amino acids") +
        #scale_fill_manual(name="Isotype", values=IG_COLORS) +
        geom_boxplot(aes_string(fill=input$plotBy_usage_AbLen))
    g2 <- ggplot(db_props, aes(x=v_call, y=cdr3_aa_gravy)) + tmp_theme + 
        ggtitle("CDR3 hydrophobicity") + 
        xlab("Isotype") + ylab("GRAVY") +
        #scale_fill_manual(name="Isotype", values=IG_COLORS) +
        geom_boxplot(aes_string(fill=input$plotBy_usage_AbLen))

    # Plot in grid column
    gridPlot(g1, g2, ncol=1)
  })

  output$cdr3_field <- renderUI({
    req(input$file2)
    req(df())

    db <- read.csv(input$file2$datapath, sep=",", header=TRUE)
    db_colnames <- colnames(db) 

    selectInput(
        'cdr3_field', 'CDR3 Colname', choices = db_colnames, selected="CDRH3"
    )
  })

  output$inters_cdr3 <- DT::renderDataTable(DT::datatable({
    req(input$file2)
    req(df())

    selected_junctions <- strsplit(input$junction_aa_input, " ")[[1]]
    selected_junctions <- trimws(selected_junctions)

    df <- df()[ df()$subject_id %in% input$subjects & df()$sample_id %in% input$samples & !(df()$junction_aa %in% selected_junctions), ] %>%
	    group_by(sample_id, v_call, j_call, junction_aa) %>%
	    summarise(dup_count = sum(duplicate_count)) %>%
	    ungroup() %>%
	    arrange(desc(dup_count)) %>%
	    filter(dup_count > 20) %>%
	    as.data.frame(.) %>%
            mutate(junction_aa=substr(junction_aa,2,as.integer(nchar(junction_aa))-1))

    db <- read.csv(input$file2$datapath, sep=",", header=TRUE)
    	
    merged_by_hamm <- data.frame()
    for(i in 1:nrow(db)) {
	print(paste0(i, " de ", nrow(db)))
        for(j in 1:nrow(df)) {
            if ((nchar(db[input$cdr3_field][i,]) != nchar(df["junction_aa"][j,])) | !(grepl(df["v_call"][j,], db["Heavy.V.Gene"][i,])) ) {
               next
           }
           # Calculate hamming distance
           hamm_dist <- sum(strsplit(db[input$cdr3_field][i,], NULL)[[1]] != strsplit(df["junction_aa"][j,], NULL)[[1]])
           if(hamm_dist < 4) {
               new_row <- cbind(db[i, ], df[j, ]) # add new lines where hamming dist < value
               merged_by_hamm <- rbind(merged_by_hamm, new_row)
           }
       }     
    }

    #for (i in 1:nrow(db)) {
    #  for (j in 1:nrow(df)) {
      # Verificar se a distância de Hamming é menor ou igual a 3
    #  if (dist_matrix[i, j] == 1) {
        # Concatenar as linhas de db e df
    #    new_row <- cbind(db[i, ], df[j, ])
        # Adicionar a nova linha ao dataframe resultante
    #    merged_by_hamm <- rbind(merged_by_hamm, new_row)
    #    }
    #  }
    #}
    
    #merged_by_cdr <- merge(db, df, by.x=input$cdr3_field, by.y="junction_aa", all.x = TRUE)
    #merged_by_cdr
    merged_by_hamm
    
  }))
}

# Create Shiny app ----
shinyApp(ui, server)
