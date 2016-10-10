library(shiny)
library(ggplot2)
library(plotly)
# ui.R

js<-"$(function() {
      var $elie = $(document.getElementsByClassName('form-group shiny-input-container'));
rotate(270);
function rotate(degree) {
$elie.css({ WebkitTransform: 'rotate(' + degree + 'deg)'});
$elie.css({ '-moz-transform': 'rotate(' + degree + 'deg)'});
}
});"

# withTags({
#   script(
#     class= "slider",
#     checked = NA,
#     HTML(js))
# })



shinyUI(fluidPage(
  tags$head(tags$style(type="text/css", "
             #loadmessage {
                       position: fixed;
                       top: 0px;
                       left: 0px;
                       width: 100%;
                       padding: 5px 0px 5px 0px;
                       text-align: center;
                       font-weight: bold;
                       font-size: 100%;
                       color: #000000;
                       background-color: #CCFF66;
                       z-index: 105;
                       }
                       "),
            tags$style(HTML("
              input[type=\"number\"] {
                            font-size: 75%;
                            }
                            "))),
  column(3,
         h3("Nanopore State Finder"),
         h6("By: Carlos de Lannoy, v: 0.1.0"),
         wellPanel(
           fileInput("abf",label="",multiple=F),
           fluidRow(
             h6("Sampling frequency"),
            column(8,numericInput("frequencyInput", label = NULL,
                        value = 30)),
            column(1,actionButton("frequencyInputTrigger","adjust"))
           ),
           hr(),
           # FILTER --------------------------------------------------------------
           h4("Filtering"),
           selectInput("selectFilter", label = "",
                       choices = list("no filter" = 1,
                                      "FFT filter" = 2,
                                      "box kernel filter" = 3,
                                      "low pass filter(Gauss.)" = 4), selected = 1),
           fluidRow(
             conditionalPanel("input.selectFilter== 2",
                              column(3, h6("Cut-off: ")),
                              column(5,numericInput("fftCutOff",label=NULL,value=0.10, min = 0, max = 1, step = 0.01, width =80))),
             conditionalPanel("input.selectFilter== 3",
                              column(3, h6("Bandwidth: ")),
                              column(5,numericInput("boxBandwidth",label = NULL,value=0.10, min = 0, max = 1, step = 0.01, width =80))),
             conditionalPanel("input.selectFilter== 4",
                              column(1, h6("W: ")),
                              column(4,numericInput("lowPassW",NULL,value=0.10, min = 0, max = 1, step = 0.01, width =80)),
                              column(1, h6("n: ")),
                              column(4,numericInput("lowPassN",NULL,value=10, min = 1, max = 50, step = 1, width =50))),
             conditionalPanel("input.selectFilter != 1",
                              column(1,actionButton("triggerFilter",label= "apply")))
           ),
           hr(),
           # HMM1 ------------------------------------------------------------------
           h4("Detect open state"),
           checkboxInput("autoTrash1", label = "Auto-trash", value = F),
           selectInput("iresMethod", label = "Ires calculation",
                       choices = list("per open state" = "perEvent",
                                      "periodically (10 min.)" = "periodically",
                                      "overall" = "overall")),
           fluidRow(
             column(5,checkboxInput("toggleHmm1",label="parameters",value=F)),
             column(1,actionButton("triggerHmm1",label= "apply"))
           ),
           conditionalPanel("input.toggleHmm1== 1",
                            br(),
                            fluidRow(
                              column(4,h6("I(pA)")),
                              column(4,h6("sd(pA)")),
                              column(2,h6("Fix value"))
                            ),
                            fluidRow(
                              column(4,numericInput("h1v7",value=NA, label=NULL),
                                       numericInput("h1v9",value=NA,label=NULL)),
                              column(4,numericInput("h1v8",value=NA,label=NULL),
                                       numericInput("h1v10",value=NA,label=NULL)),
                              column(1,checkboxInput("h1c7",value=T,label=NULL),
                                       checkboxInput("h1c9",value=T,label=NULL)),
                              column(1,checkboxInput("h1c8",value=T,label=NULL),
                                       checkboxInput("h1c10",value=T,label=NULL))
                            ),
                            fluidRow(
                              column(8,h6("Transition matrix"))
                              ),
                            fluidRow(
                              column(4,numericInput("h1v3",value=NA,label=NULL),
                                     numericInput("h1v5",value=NA,label=NULL)),
                              column(4,numericInput("h1v4",value=NA,label=NULL),
                                     numericInput("h1v6",value=NA,label=NULL)),
                              column(1,checkboxInput("h1c3",value=T,label=NULL),
                                     checkboxInput("h1c5",value=T,label=NULL)),
                              column(1,checkboxInput("h1c4",value=T,label=NULL),
                                     checkboxInput("h1c6",value=T,label=NULL))
                            )
           ),
           hr(),
           # HMM2 -------------------------------------------------------------------
           h4("Detect blocked states"),
           fluidRow(
             column(4,checkboxInput("autoTrash2", label = "Auto-trash", value = F)),
             column(2,h6("# of states:")),
             column(4,numericInput("nStates",label=NULL,value="0"))
           ),
           fluidRow(
             column(5,checkboxInput("toggleHmm2",label="parameters")),
             column(1,actionButton("triggerHmm2",label= "apply"))
           ),
           conditionalPanel("input.toggleHmm2== 1",
                            br(),
                            fluidRow(
                              column(4,h6("I(pA)")),
                              column(4,h6("sd(pA)")),
                              column(2,h6("Fix value"))
                            ),
                            conditionalPanel("input.nStates==2",
                              fluidRow(
                                column(4,numericInput("h2n2v7",value=NA, label=NULL),
                                       numericInput("h2n2v9",value=NA,label=NULL)),
                                column(4,numericInput("h2n2v8",value=NA,label=NULL),
                                       numericInput("h2n2v10",value=NA,label=NULL)),
                                column(1,checkboxInput("h2n2c7",value=T,label=NULL),
                                       checkboxInput("h2n2c9",value=T,label=NULL)),
                                column(1,checkboxInput("h2n2c8",value=T,label=NULL),
                                       checkboxInput("h2n2c10",value=T,label=NULL))
                              ),
                              fluidRow(
                                column(8,h6("Transition matrix"))
                              ),
                              fluidRow(
                                column(4,numericInput("h2n2v3",value=NA, label=NULL),
                                       numericInput("h2n2v5",value=NA,label=NULL)),
                                column(4,numericInput("h2n2v4",value=NA,label=NULL),
                                       numericInput("h2n2v6",value=NA,label=NULL)),
                                column(1,checkboxInput("h2n2c3",value=T,label=NULL),
                                       checkboxInput("h2n2c5",value=T,label=NULL)),
                                column(1,checkboxInput("h2n2c4",value=T,label=NULL),
                                       checkboxInput("h2n2c6",value=T,label=NULL))
                              )
                            ),
                            conditionalPanel("input.nStates==3",
                                             fluidRow(
                                               column(4,numericInput("h2n3v13",value=NA, label=NULL),
                                                      numericInput("h2n3v15",value=NA, label=NULL),
                                                      numericInput("h2n3v17",value=NA,label=NULL)),
                                               column(4,numericInput("h2n3v14",value=NA,label=NULL),
                                                      numericInput("h2n3v16",value=NA,label=NULL),
                                                      numericInput("h2n3v18",value=NA,label=NULL)),
                                               column(1,checkboxInput("h2n3c13",value=T,label=NULL),
                                                      checkboxInput("h2n3c15",value=T,label=NULL),
                                                      checkboxInput("h2n3c17",value=T,label=NULL)),
                                               column(1,checkboxInput("h2n3c14",value=T,label=NULL),
                                                      checkboxInput("h2n3c16",value=T,label=NULL),
                                                      checkboxInput("h2n3c18",value=T,label=NULL))
                                             ),
                                             fluidRow(
                                               column(8,h6("Transition matrix"))
                                             ),
                                             fluidRow(
                                               column(3,numericInput("h2n3v4",value=NA, label=NULL),
                                                      numericInput("h2n3v7",value=NA, label=NULL),
                                                      numericInput("h2n3v10",value=NA,label=NULL)),
                                               column(3,numericInput("h2n3v5",value=NA,label=NULL),
                                                      numericInput("h2n3v8",value=NA,label=NULL),
                                                      numericInput("h2n3v11",value=NA,label=NULL)),
                                               column(3,numericInput("h2n3v6",value=NA,label=NULL),
                                                      numericInput("h2n3v9",value=NA,label=NULL),
                                                      numericInput("h2n3v12",value=NA,label=NULL)),
                                               column(1,checkboxInput("h2n3c4",value=T,label=NULL),
                                                      checkboxInput("h2n3c7",value=T,label=NULL),
                                                      checkboxInput("h2n3c10",value=T,label=NULL)),
                                               column(1,checkboxInput("h2n3c5",value=T,label=NULL),
                                                      checkboxInput("h2n3c8",value=T,label=NULL),
                                                      checkboxInput("h2n3c11",value=T,label=NULL)),
                                               column(1,checkboxInput("h2n3c6",value=T,label=NULL),
                                                      checkboxInput("h2n3c9",value=T,label=NULL),
                                                      checkboxInput("h2n3c12",value=T,label=NULL))
                                             )
                            )
           ),
           hr(),
           style = "overflow-y:scroll; max-height: 900px")
  ),
  # PLOT-RELATED ---------------------------------------------------------------
  column(9,
         tabsetPanel(
           tabPanel("trace", selectInput("traceSelect", label = NULL, choices = list("current" = "current",
                                                                                     "filtered" = "current.filtered",
                                                                                     "Ires" = "current.res")),
                    plotlyOutput("plotTrace",height = "700px"),
                    conditionalPanel(condition = "input.traceSelect== 'current'",
                                     fluidRow(column(10,sliderInput("vertSlider",label="current",min = -100, max = 100, value = c(-100,100),width ="1200px")),
                                              column(1,br(),br(),checkboxInput("toggleGuides",label = "guides",value=T))),
                                     fluidRow(column(10,sliderInput("horSlider",label="time",min=0, max = 100, value = c(0,100),width ="1200px")),
                                              column(2,actionButton("triggerRestrict","restrict axes"))
                                     )
                    )
                    ),
           tabPanel("scatter-heat", plotlyOutput("plotScatterheat",height = "700px")),
           tabPanel("Histogram", fluidRow(column(6,h3("Event duration")), column(6,h3("Mean residual current"))),
                    plotOutput("plotHist",height = "700px"),
                    fluidRow(
                      column(1, numericInput("bwTime",label="binwidth:", value= NA, step=0.001)),
                      column(5, sliderInput("rangeTime", label="range:", min=NA, max=NA, value=c(NA,NA),width="500px")),
                      column(1, numericInput("bwCurrent",label="binwidth:", value= NA, step=0.001)),
                      column(5, sliderInput("rangeCurrent", label="range:", min=NA, max=NA, value=c(NA,NA),width="500px"))
                    ))
         ),
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div("Loading...",id="loadmessage")
    )
  )
)
)
