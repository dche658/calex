#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinyjs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(irr)


options(shiny.maxRequestSize=30*1024^2)

formatGender <- function(txt) {
  sex <- ifelse(grepl("^[Mm]",txt),"M",
                ifelse(grepl("^[Ff]",txt),"F","U"))
  return(sex)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#formatGender(c("M","Male","m","Male","F","Female","female","I","UNK"))


#' CKD-EPI Estimated GFR
#'
#' Estimated glomerular filtration rate calculated by the CKD-EPI equation.
#' All data should be presented as vectors of values of equal length
#' 
#' @param creatinine Serum creatinine concentration in umol/L
#' @param age Age in years
#' @param gender Sex of the patient as M or F
#' @param black TRUE if African American, otherwise FALSE
#' @param ckdepi.formula Character string indicating "2009" or "2021"
#'
egfr <- function(creatinine, age, gender, black=FALSE, ckdepi.formula="2009") {
  has.gender <- ifelse(!(gender=="M" | gender=="F"), FALSE, TRUE)
  ifelse(!(ckdepi.formula=="2009" | ckdepi.formula=="2021"), stop(simpleError("ckdepi.formula must be '2009' or '2021'")),1)
  
  mu <- ifelse(ckdepi.formula=="2009",141,142)
  kappa <- ifelse(gender=="F", 0.7, 0.9)
  delta <- ifelse(gender=="F" & ckdepi.formula=="2009", 1.018,
                  ifelse(gender=="F" & ckdepi.formula=="2021",1.012, 1.0))
  alpha <- ifelse(gender=="F" & creatinine <= 62 & ckdepi.formula=="2009", -0.329,
                  ifelse(gender=="F" & creatinine > 62 & ckdepi.formula=="2009", -1.209,
                         ifelse(gender=="M" & creatinine <= 80 & ckdepi.formula=="2009", -0.411,
                                ifelse(gender=="M" & creatinine > 80 & ckdepi.formula=="2009", -1.209,
                                       ifelse(gender=="F" & creatinine <= 62 & ckdepi.formula=="2021", -0.241,
                                              ifelse(gender=="F" & creatinine > 62 & ckdepi.formula=="2021", -1.200,
                                                     ifelse(gender=="M" & creatinine <= 80 & ckdepi.formula=="2021", -0.302, -1.200)))))))
  gamma <- ifelse(ckdepi.formula=="2009", 0.9929, 0.9938)
  epsilon <- ifelse(ckdepi.formula=="2009" & black==TRUE, 1.159, 1.0)
  gfr <- ifelse(has.gender==TRUE & age>=18,
                mu * ((creatinine * 0.0113 / kappa)^alpha) * gamma^age * delta * epsilon,
                NA)
  return(gfr)
}


# Define UI for application
ui <- dashboardPage(
  dashboardHeader(title = "Calcium explorer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "intro"),
      menuItem("File Upload", tabName = "upload"),
      menuItem("Column Selection", tabName = "colselect"),
      menuItem("Analysis", tabName = "analysis")
    ) # sidebarMenu
  ),
  dashboardBody(
    useShinyjs(),
    tabItems(
      # First tab content
      tabItem(tabName = "intro",
              fluidRow(
                column(12,
                  p("The purpose of this application is to assist laboratories to 
                      verify the albumin adjusted calcium equation using their
                      own laboratory data."),
                  h4("Upload and select data"),
                  p("Your data should be in the form of a comma separated values (CSV) file
                    with data arrange in columns. Albumin and calcium may be any units. 
                    Gender should be recorded as M, Male, F, or Female. Age must be
                    in years and serum creatinine concentration must be in umol/L in
                    order for the estimated GFR to calculate corretly. The app uses
                    the 2009 CKD-EPI equation."),
                  p("Once the data is loaded a table will display showing the first 
                    few rows of the file."),
                  p("Use the column selection page to indicate
                    which columns should be used for each of the variable"),
                  h4("Analysis"),
                  p("There are two parts to this section. The first uses the supplied
                  data to calculate the regression of calcium on albumin. This can
                  be used to derive the optimal slope for the available data. The
                  filters can be used to select different subsets of data for 
                    analysis. The second part allows the parameters for the albumin
                    adjusted calcium equation to be set to look at the proportion
                    of values below, within, or above the reference interval using
                    the supplied data; as well as assessing the concordance with 
                    ionised calcium if it is available."),
                  style="margin: 4px;"
                )
                
              )
      ),
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Upload Data", width = 12,
                    status = "primary",
                    solidHeader = TRUE,
                    fileInput("file1","Choose CSV file", 
                              accept = "csv",
                              multiple = FALSE),
                    checkboxInput("header", "Header", TRUE))),
              fluidRow(
                box(title = "File data preview", width = 12, solidHeader = TRUE,
                    textOutput("csvInfo"),br(),
                      column(12,tableOutput("table1"),
                             style="height:300px; overflow-y: scroll;overflow-x: scroll;")
                    
                    )
              ),
      ),
      tabItem(tabName = "colselect",
              fluidRow(
                box(title = "Upload Data", width = 12,
                    status = "primary",
                    solidHeader = TRUE,
                    fluidRow(
                    column(3,selectInput("col_ca",
                                         "Calcium column",""),
                           selectInput("col_alb","Albumin column","")),
                    column(3,selectInput("col_ica","Ionised Ca column","")),
                    column(3,selectInput("col_age","Age (years) column",""),
                           selectInput("col_sex","Sex column",""),
                           selectInput("col_crea","Creatinine column","")),
                    ),
                    fluidRow(
                      column(12,h4("Data preview"),
                             tableOutput("table2"))
                    )
                )
              )),
      tabItem(tabName = "analysis",
              fluidRow(
                box(title="Albumin filter", width="4",
                    status = "primary", solidHeader=TRUE,
                    sliderInput("filterAlb", 
                                label = "Albumin", 
                                min=0,max=60, value=c(0,60))),
                box(title="Age filter", width="4",
                    status = "primary", solidHeader=TRUE,
                    sliderInput("filterAge", 
                                label="Age",
                                min=0,max=120, value=c(0,120))),
                box(title="eGFR filter", width="4",
                    status = "primary", solidHeader=TRUE,
                    sliderInput("filterGfr", 
                                label = "eGFR",
                                min=0,max=150, value=c(0,150))),
              ),
              fluidRow(
                box(title="Calcium Reference Interval", width=6,
                    status="primary", solidHeader=TRUE,
                    column(6,numericInput("lrlCal", label="Lower Ref Limit",
                                          value=2.10, step = 0.01)),
                    column(6,numericInput("urlCal", label="Upper Ref Limit",
                                          value=2.60, step = 0.01))),
                box(title="Ionised Ca Reference Interval", width=6,
                    status="primary", solidHeader=TRUE,
                    column(6,numericInput("lrlIca", label="Lower Ref Limit",
                                          value=1.15, step = 0.01)),
                    column(6,numericInput("urlIca", label="Upper Ref Limit",
                                          value=1.30, step = 0.01)))
                
              ),
              fluidRow(
                column(3, numericInput("midAlb", label="Albumin set point",
                                    value="40", step = 1)),
                column(3, p("Albumin concentration to be used in the adjusted 
                            calcium equation."))
              ),
              fluidRow(
                tabBox(id = "tabset1", width = 12,
                       tabPanel("Regression",
                                fluidRow(
                                  column(12, h4("Regression statistics"),
                                         tableOutput("table3"),
                                         textOutput("rse"),
                                         textOutput("mrs"),
                                         plotOutput("scatterplot"),
                                         h4("Albumin adjusted calcium"),
                                         textOutput("equation"),
                                         plotOutput("adjCaHistogram"),
                                         tableOutput("adjCaCatTbl"),
                                         h4("Adjusted calcium and ionised calcium concordance"),
                                         tableOutput("concordance"),
                                         plotOutput("scatterplot2")
                                )
                                
                                )
                       ),
                       tabPanel("Verification",
                                fluidRow(
                                  column(12, h4("Adjusted Calcium Verification"),
                                         p("The adjusted calcium equation is given 
                                           by:"),
                                         p("Adjusted Calcium = Measured Calcium + (Slope x (SetPoint - Measured Albumin))"),
                                         p("where SetPoint is the albumin set point entered above, and slope is given by the
                                           value below."),
                                         numericInput("slope", label = "Slope",
                                                               value = 0.020,
                                                               step = 0.001),
                                         h4("Albumin adjusted calcium"),
                                         plotOutput("adjCaHist2"),
                                         tableOutput("adjCaCatTbl2"),
                                         h4("Adjusted calcium versus ionised calcium"),
                                         fluidRow(
                                           column(6,textOutput("verify.R.sq"))
                                         ),
                                         fluidRow(
                                           box(title="Kappa", width=4,
                                               tableOutput("concordance2")),
                                           box(title="Adjusted Calcium v Ionised Calcium", width=4,
                                               column(12,tableOutput("contingency2"),
                                                      textOutput("agreement2"))),
                                           box(title="Total Calcium v Ionised Calcium", width=4,
                                               column(12,tableOutput("contingency3"),
                                                      textOutput("agreement3")))
                                         ),
                                         plotOutput("concordancePlot2")
                                         )
                                )
                                )
                       )
              )
      )# tabItem
    ) # tabItems
  ) # dashboardBody
) # dashboardPage

# Define server logic 
server <- function(input, output) {
  
  # dat: data frame to hold the imported CSV columns
  global <- reactiveValues(dat = NULL)
  
  # Read data from file
  fileData <- reactive({
    req(input$file1)
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    df <- read.csv(file$datapath, header = input$header)
    csvCols = colnames(df)
    if(length(csvCols)>0) {
      # message("Updating columns")
      updateSelectInput(inputId = "col_ca",
                        choices = c("Select",csvCols),
                        selected = "Select")
      updateSelectInput(inputId = "col_alb",
                        choices = c("Select",csvCols),
                        selected = "Select")
      updateSelectInput(inputId = "col_ica",
                        choices = c("Skip",csvCols),
                        selected = "Skip")
      updateSelectInput(inputId = "col_age",
                        choices = c("Skip",csvCols),
                        selected = "Skip")
      updateSelectInput(inputId = "col_sex",
                        choices = c("Skip",csvCols),
                        selected = "Skip")
      updateSelectInput(inputId = "col_crea",
                        choices = c("Skip",csvCols),
                        selected = "Skip")
    }
    # message("Read file data")
    df
  })
  
  selectedData <- reactive({
    # selected columns
    sh <- c()
    # column names
    ch <- c()
    df <- NULL
    if(input$col_ca != "Select") {
      sh <- c(sh,input$col_ca)
      ch <- c(ch,"calcium")
    }
    if(input$col_alb != "Select") {
      sh <- c(sh,input$col_alb)
      ch <- c(ch,"albumin")
    }
    if(input$col_ica != "Skip") {
      sh <- c(sh,input$col_ica)
      ch <- c(ch,"ica")
    }
    if(input$col_age != "Skip") {
      sh <- c(sh,input$col_age)
      ch <- c(ch,"age")
    }
    if(input$col_sex != "Skip") {
      sh <- c(sh,input$col_sex)
      ch <- c(ch,"sex")
    }
    if(input$col_crea != "Skip") {
      sh <- c(sh,input$col_crea)
      ch <- c(ch,"crea")
    }
    if (length(sh)>1) {
      df <- fileData()[,sh]
      colnames(df) <- ch
      if(input$col_sex != "Skip") {
        df$sex <- formatGender(df$sex)
      }
      if (input$col_age != "Skip" & input$col_sex != "Skip" & input$col_crea != "Skip") {
        df$egfr <- egfr(df$crea, df$age, df$sex)
      }
    }
    df
  })
  
  # Apply filters to the data in dat
  filteredData <- reactive({
    df <- selectedData()
    if (!is.null(df) & 
        "calcium" %in% colnames(df) &
        "albumin" %in% colnames(df)) {
      c <- colnames(df)
      # apply filters
      df <- df[df$albumin > input$filterAlb[1] & df$albumin < input$filterAlb[2],]
      if ("age" %in% c) {
        df <- df[df$age>input$filterAge[1] & df$age<input$filterAge[2],]
      }
      if ("egfr" %in% c) {
        df <- df[df$egfr>input$filterGfr[1] & df$egfr<input$filterGfr[2],]
      }
      # remove rows with null values
      df <- na.omit(df)
    }
    # message("Filtered data")
    df
  })
  
  # Calculate the regression equation
  regression <- reactive({
    r <- NULL
    if (!is.null(filteredData())) {
      if(nrow(filteredData())>1) {
        r <- lm(calcium~albumin,filteredData())
      }
    }
    # message("Calculated regression")
    r
  })
  
  # Calculate the adjusted calcium concentration
  # Chart distribution against the reference intervals and calculate fraction
  # outside the Ca reference interval
  # Determine if adjusted calcium is below, within, or above the calcium 
  # reference interval.
  
  # If ionised calcium data is available
  #  - Determine if ionised calcium is below, within, or above the reference
  #    interval
  #  - Calculate kappa
  #  - Chart distributions against the reference intervals and calculate fraction
  #    outside the iCa reference interval
  
  # Calculate adjusted calcium based on the regression coefficient
  adjustedCalcium <- reactive({
    adjca <- NULL
    if (!is.null(filteredData()) & !is.null(regression())) {
      adjca <- filteredData()$calcium + (regression()$coefficients[2]*(input$midAlb - filteredData()$albumin))
    }
    # message("Calculated adjusted calcium")
    adjca
  })
  
  # Calculate adjusted calcium based on the provided parameters
  adjCa2 <- reactive({
    adjca <- NULL
    if (!is.null(filteredData())) {
      adjca <- filteredData()$calcium + (input$slope*(input$midAlb - filteredData()$albumin))
    }
    # message("Calculated adjusted calcium 2")
    adjca
  })
  
  adjIcaCat <- reactive({
    adjiccat <- NULL
    if(!is.null(filteredData()) & "ica" %in% colnames(filteredData())) {
      adjiccat <- ifelse(filteredData()$ica < input$lrlIca,
                         "below",
                         ifelse(filteredData()$ica > input$urlIca,
                                "above","within"))
      adjiccat <- factor(adjiccat, levels = c("below", "within", "above"))  
    }
    adjiccat
  })
  
  adjCalCategory <- reactive({
    adjcacat <- NULL
    if (!is.null(adjustedCalcium())) {
      adjcacat <- ifelse(adjustedCalcium() < input$lrlCal,
                         "below",
                         ifelse(adjustedCalcium() > input$urlCal,
                                "above","within"))
      adjcacat <- factor(adjcacat, levels = c("below", "within", "above"))
    }
    adjcacat
  })
  
  adjCalCat2 <- reactive({
    adjcacat <- NULL
    if (!is.null(adjCa2())) {
      adjcacat <- ifelse(adjCa2() < input$lrlCal,
                         "below",
                         ifelse(adjCa2() > input$urlCal,
                                "above","within"))
      adjcacat <- factor(adjcacat, levels = c("below", "within", "above"))
    }
    adjcacat
  })
  
  calciumCategory <- reactive({
    cacat <- NULL
    if (!is.null(filteredData())) {
      cacat <- ifelse(filteredData()$calcium < input$lrlCal,
                         "below",
                         ifelse(filteredData()$calcium > input$urlCal,
                                "above","within"))
      cacat <- factor(cacat, levels = c("below", "within", "above"))
    }
    cacat
  })
  
  # Calculate the contingency table for adj calcium and ionised calcium
  # using the user provided adjusted calcium parameters
  ctgcyTbl2 <- reactive({
    df <- NULL
    if(!is.null(adjCalCat2()) & !is.null(adjIcaCat())) {
      df <- table(adjCalCat2(),adjIcaCat())
    }
    df
  })
  
  # Calculate the contingency table for total calcium and ionised calcium
  ctgcyTbl3 <- reactive({
    df <- NULL
    if(!is.null(calciumCategory()) & !is.null(adjIcaCat())) {
      df <- table(calciumCategory(),adjIcaCat())
    }
    df
  })
  
  
  output$table1 <- renderTable({
    global$df <- fileData()
    head(fileData())
  }, bordered = TRUE, width = "400px" )
  
  output$table2 <- renderTable({
    df <- selectedData()
    if(!is.null(df)) {
      return(head(df))
    }
  })
  
  output$csvInfo <- renderText({
    if(!is.null(fileData())) {
      return(paste("Imported file contains",nrow(fileData()),"rows."))
    } else {
      return("")
    }
  })
  
  output$table3 <- renderTable({
    data <- NULL
    if(!is.null(regression())) {
      data <- summary(regression())$coefficients
      rownames(data) <- c("Intercept","Slope")
    }
    data
  }, digits=4, rownames = TRUE)
  
  output$rse <- renderText({
    if(!is.null(regression())) {
      s <- summary(regression())
      return(paste("Residual standard error: ",
                   round(s$sigma, digits = 4),
                   "on",
                   s$df[2],
                   "degrees of freedom"))
    } else {
      return(NULL)
    }
  })
  
  output$mrs <- renderText({
    if(!is.null(regression())) {
      s <- summary(regression())
      return(paste("Multiple R-squared:", 
                   round(s$r.squared, digits = 4)))
    } else {
      return(NULL)
    }
  })
  
  output$scatterplot <- renderPlot({
    plt <- NULL
    if(!is.null(filteredData()) & !is.null(regression())) {
      plt <- ggplot(filteredData())+
        geom_point(aes(albumin,calcium), color="blue",alpha=0.4)+
        geom_abline(slope = regression()$coefficients[2],intercept = regression()$coefficients[1])+
        xlab("Albumin")+
        ylab("Calcium")
    }
    plt
  })
  
  output$equation <- renderText({
    if (!is.null(regression())) {
      paste("Adjusted Calcium = Measured Calcium + ",
            round(regression()$coefficients[2],digits = 4),
            " x (",input$midAlb," - Albumin)")
    }
  })
  
  output$adjCaHistogram <- renderPlot({
    plt <- NULL
    if(!is.null(filteredData())) {
      #global$adjca <- adjustedCalcium()
      if (!is.null(adjustedCalcium())) {
        df <- data.frame(adjca = adjustedCalcium())
        plt <- ggplot(df)+
          geom_histogram(aes(x=adjca), binwidth = 0.02)+
          geom_vline(xintercept = input$lrlCal, color="red")+
          geom_vline(xintercept = input$urlCal, color="red")+
          xlab("Adjusted calcium")  
      }
    }
    plt
  })
  
  output$adjCaHist2 <- renderPlot({
    plt <- NULL
    if(!is.null(filteredData())) {
      if (!is.null(adjCa2())) {
        df <- data.frame(adjca = adjCa2())
        plt <- ggplot(df)+
          geom_histogram(aes(x=adjca), binwidth = 0.02)+
          geom_vline(xintercept = input$lrlCal, color="red")+
          geom_vline(xintercept = input$urlCal, color="red")+
          xlab("Adjusted calcium")  
      }
    }
    plt
  })
  output$adjCaCatTbl <- renderTable({
    data <- NULL
    if(!is.null(filteredData())) {
      if (!is.null(adjCalCategory())) {
        data <- as.data.frame(table(adjCalCategory()))
        data$Percent=round(100*data$Freq/sum(data$Freq),digits = 1)
        colnames(data) <- c("Category","Freq","Percent")
      }
    }
    data
  }, rownames = TRUE)
  
  output$adjCaCatTbl2 <- renderTable({
    data <- NULL
    if(!is.null(filteredData())) {
      if (!is.null(adjCalCat2())) {
        data <- as.data.frame(table(adjCalCat2()))
        data$Percent=round(100*data$Freq/sum(data$Freq),digits = 1)
        colnames(data) <- c("Category","Freq","Percent")
      }
    }
    data
  }, rownames = TRUE)
  
  output$concordance <- renderTable({
    df <- NULL
    if(!is.null(filteredData())) {
      adjca <- adjustedCalcium()
      if (!is.null(adjCalCategory()) & "ica" %in% colnames(filteredData())) {
        adjiccat <- adjIcaCat()
        k <- kappa2(data.frame(adjCalCategory(),adjiccat))
        df <- data.frame(Results=c(as.character(k$subjects),
                                   as.character(k$raters),
                                   as.character(round(k$value,digits = 3)),
                                   as.character(round(k$statistic, digits = 1)),
                                   as.character(k$p.value)))
        rownames(df) <- c("Subject","Raters","Kappa","z","p-value") 
      }
    }
    df
  }, rownames = TRUE)
  
  
  # Compared adjusted calcium calculated by provided parameters to ionised calcium
  
  # Continuous comparison
  output$verify.R.sq <- renderText({
    txt <- ""
    if(!is.null(adjCa2()) & !is.null(filteredData()) & "ica" %in% colnames(filteredData())) {
      txt <- paste("R-squared:",round(cor(adjCa2(),
                                          filteredData()$ica, 
                                          method = "pearson")**2,
                                      digits = 3))
    }
    txt
  })
  
  # Categorical comparison
  output$concordance2 <- renderTable({
    df <- NULL
    if(!is.null(filteredData())) {
      adjca <- adjCa2()
      if (!is.null(adjCalCat2()) & "ica" %in% colnames(filteredData())) {
        adjiccat <- adjIcaCat()
        k <- kappa2(data.frame(adjCalCat2(),adjiccat))
        df <- data.frame(Results=c(as.character(k$subjects),
                                   as.character(k$raters),
                                   as.character(round(k$value,digits = 3)),
                                   as.character(round(k$statistic, digits = 1)),
                                   as.character(k$p.value)))
        rownames(df) <- c("Subject","Raters","Kappa","z","p-value") 
      }
    }
    df
  }, rownames = TRUE)
  
  # Contingency table and percent agreement Adjusted Ca versus Ionised Ca
  output$contingency2 <- renderTable({
    df <- ctgcyTbl2()
    if (!is.null(df)) {
      as.data.frame.matrix(ctgcyTbl2())
    }
  }, rownames = TRUE)
  
  output$agreement2 <- renderText({
    df <- ctgcyTbl2()
    txt <- ""
    if (!is.null(df)) {
      agree <- 0
      for(i in 1:nrow(df)) {
        agree <- agree + df[i,i]
      }
      percentAgree <- 100* agree / nrow(filteredData())
      txt <- paste("Percent agreement:",round(percentAgree,digits = 1),"%")
    }
    txt
  })
  
  # Contingency table and percent agreement Total Ca versus Ionised Ca
  output$contingency3 <- renderTable({
    df <- ctgcyTbl3()
    if (!is.null(df)) {
      as.data.frame.matrix(ctgcyTbl3())  
    }
  }, rownames = TRUE)
  
  output$agreement3 <- renderText({
    df <- ctgcyTbl3()
    txt <- ""
    if (!is.null(df)) {
      agree <- 0
      for(i in 1:nrow(df)) {
        agree <- agree + df[i,i]
      }
      percentAgree <- 100* agree / nrow(filteredData())
      txt <- paste("Percent agreement:",round(percentAgree,digits = 1),"%")
    }
    txt
  })
  
  output$scatterplot2 <- renderPlot({
    plt <- NULL
    if(!is.null(filteredData())) {
      adjca <- adjustedCalcium()
      if (!is.null(adjca) & "ica" %in% colnames(filteredData())) {
        df <- filteredData() %>%
          select(calcium,ica,albumin)
        df <- data.frame(df,adjca) %>% 
          pivot_longer(cols = c(calcium, adjca),
                       names_to = "Analyte",
                       values_to = "calcium")
        df$Analyte <- ifelse(df$Analyte=="calcium",
                             "Total calcium",
                             "Adjusted calcium")
        plt <- ggplot(df)+
          geom_point(aes(ica, calcium, color=albumin))+
          geom_hline(yintercept = input$lrlCal, color="red")+
          geom_hline(yintercept = input$urlCal, color="red")+
          geom_vline(xintercept = input$lrlIca, color="red")+
          geom_vline(xintercept = input$urlIca, color="red")+
          facet_wrap(facets = vars(Analyte))+
          xlab("Ionised calcium (mmol/L)")+
          ylab("Calcium (mmol/L)")
      }
    }
    plt
  })
  
  output$concordancePlot2 <- renderPlot({
    plt <- NULL
    if(!is.null(filteredData())) {
      adjca <- adjCa2()
      if (!is.null(adjca) & "ica" %in% colnames(filteredData())) {
        df <- filteredData() %>%
          select(calcium,ica,albumin)
        df <- data.frame(df,adjca) %>% 
          pivot_longer(cols = c(calcium, adjca),
                       names_to = "Analyte",
                       values_to = "calcium")
        df$Analyte <- ifelse(df$Analyte=="calcium",
                             "Total calcium",
                             "Adjusted calcium")
        plt <- ggplot(df)+
          geom_point(aes(ica, calcium, color=albumin))+
          geom_hline(yintercept = input$lrlCal, color="red")+
          geom_hline(yintercept = input$urlCal, color="red")+
          geom_vline(xintercept = input$lrlIca, color="red")+
          geom_vline(xintercept = input$urlIca, color="red")+
          facet_wrap(facets = vars(Analyte))+
          xlab("Ionised calcium (mmol/L)")+
          ylab("Calcium (mmol/L)")
      }
    }
    plt
  })
  
  observe({
    shinyjs::toggleState("filterAge", input$col_age != "Skip")
  })
  observe({
    shinyjs::toggleState("filterGfr", input$col_crea != "Skip" & 
                           input$col_age != "Skip" &
                           input$col_sex != "Skip")
  })

  observeEvent(input$midAlb, {
    adjustedCalcium()
  })

}

# Run the application 
shinyApp(ui = ui, server = server)
