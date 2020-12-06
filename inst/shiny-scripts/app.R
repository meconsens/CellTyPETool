library(shiny)
source("../../R/bretMarkerEffectOnPathology.R")
load("../../data/bretCellMarkers.rda")
load("../../data/countDf.rda")
load("../../data/metadata.rda")

ui <- fluidPage(
  titlePanel("BRETIGEA Method And the Effect of Marker Genes on Pathology Association"),

  sidebarLayout(
    sidebarPanel(
      helpText("Create a plot showing the effect of each marker genes addition on the cell type proportion calculation's association to the disease pathology."),

      selectInput(inputId = "var",
                  label = "Choose a cell type:",
                  choices = list("Astrocytes",
                                 "Endothelial Cells",
                                 "Microglial Cells",
                                 "Neurons",
                                 "Oligodendrocytes",
                                 "Oligodendrocyte Precursors"),
                  selected = "Astrocytes"),

      sliderInput(inputId = "num",
                  label = "Number of markers",
                  value = 5, min = 1, max = 50),
    ),


    mainPanel(plotOutput("BRET"))
  )
)

server <- function(input, output) {
  output$BRET <- renderPlot({
    data <- switch(input$var,
                   "Astrocytes" = "ast",
                   "Microglial Cells" = "mic",
                   "Endothelial Cells" = "end",
                   "Neurons" = "neu",
                   "Oligodendrocytes" = "oli",
                   "Oligodendrocyte Precurors" = "opc")

    bretMarkerEffectOnPathologyResult <- bretMarkerEffectOnPathology(
      countDf = countDf,
      bretCellMarkers = bretCellMarkers,
      cellName = data,
      metadata = metadata,
      covar = "Covariate",
      pathologyName = "DiseasePhenotypeScore",
      cellTypeNames = unique(bretCellMarkers$cell),
      n= input$num)

    print(bretMarkerEffectOnPathologyResult)

  })
}
shiny::shinyApp(ui = ui, server = server)
# [END]
