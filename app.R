library(shiny)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(DT)

# --- 1. Data Generation & Analysis (Preprocessing) ---
# We run the core DESeq2 logic once outside the UI/Server to initialize data
set.seed(123)
genes <- 1000
samples <- 6
counts <- matrix(rnbinom(genes * samples, mu = 100, size = 1), nrow = genes)
rownames(counts) <- paste0("Gene_", 1:genes)
colnames(counts) <- paste0("Sample_", 1:samples)

metadata <- data.frame(
  condition = factor(c("Control", "Control", "Control", "Disease", "Disease", "Disease")),
  row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res) %>% filter(!is.na(padj))
vsd <- vst(dds, blind = FALSE)

# --- 2. User Interface ---
ui <- fluidPage(
  theme = shinythemes::shinytheme("cosmo"),
  titlePanel("Differential Gene Expression Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Filter Parameters"),
      sliderInput("p_adj", "Adjusted P-value Threshold:", 
                  min = 0.001, max = 0.1, value = 0.05, step = 0.005),
      sliderInput("lfc", "Min Log2 Fold Change:", 
                  min = 0, max = 5, value = 1, step = 0.5),
      hr(),
      h4("Heatmap Settings"),
      numericInput("top_n", "Number of Top Genes:", value = 30, min = 5, max = 100),
      downloadButton("downloadData", "Download Results (.csv)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("PCA Plot", plotOutput("pcaPlot", height = "500px")),
        tabPanel("Volcano Plot", plotOutput("volcanoPlot", height = "500px")),
        tabPanel("Heatmap", plotOutput("heatmapPlot", height = "600px")),
        tabPanel("Results Table", DTOutput("resultsTable"))
      )
    )
  )
)

# --- 3. Server Logic ---
server <- function(input, output) {
  
  # Reactive filtering based on UI inputs
  filtered_res <- reactive({
    res_df %>%
      filter(padj < input$p_adj, abs(log2FoldChange) > input$lfc)
  })
  
  # PCA Plot
  output$pcaPlot <- renderPlot({
    plotPCA(vsd, intgroup = "condition") + 
      theme_minimal() + 
      geom_text(aes(label = name), vjust = 1.5) +
      labs(title = "Principal Component Analysis (PCA)")
  })
  
  # Volcano Plot
  output$volcanoPlot <- renderPlot({
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) + 
      geom_point(aes(color = padj < input$p_adj & abs(log2FoldChange) > input$lfc), alpha = 0.6) +
      scale_color_manual(values = c("black", "red")) +
      geom_hline(yintercept = -log10(input$p_adj), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-input$lfc, input$lfc), linetype = "dashed", color = "blue") +
      theme_minimal() +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
      theme(legend.position = "none")
  })
  
  # Heatmap
  output$heatmapPlot <- renderPlot({
    top_genes_idx <- head(order(res$padj), input$top_n)
    heatmap_data <- assay(vsd)[top_genes_idx, ]
    
    pheatmap(
      heatmap_data,
      scale = "row",
      annotation_col = metadata,
      main = paste("Top", input$top_n, "DEGs Heatmap"),
      color = colorRampPalette(c("blue", "white", "red"))(150)
    )
  })
  
  # Interactive Data Table
  output$resultsTable <- renderDT({
    datatable(filtered_res(), options = list(pageLength = 10)) %>%
      formatRound(columns = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), digits = 4)
  })
  
  # Download Handler
  output$downloadData <- downloadHandler(
    filename = function() { paste("DESeq2_Results_", Sys.Date(), ".csv", sep="") },
    content = function(file) { write.csv(filtered_res(), file) }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)