# library(shiny)
# library(ggplot2)
# library(abf2)
# library(plotly)
# library(depmixS4)

# # Install required packages
# listOfPackages <- c("abf2", "depmixS4","ggplot2","grid", "shiny", "signal")
# newPackages <- listOfPackages[!(listOfPackages %in% installed.packages()[,"Package"])]
# if(length(newPackages)) install.packages(newPackages)

# source("D:/Google Drive/FIG Free Internship Groningen/FIG data/abfData_05092016.R")
options(shiny.maxRequestSize=100*1024^2)
shinyServer(function(input, output, session) {

  # Load a new abfFile
  abfNew = reactive({
    if(is.null(input$abf)){return()}
    if(input$frequencyInputTrigger == 0)
      return(new(Class="AbfData",fileName=input$abf$datapath, samplingFrequency = 30))
    new(Class="AbfData",fileName=input$abf$datapath, samplingFrequency = changeFrequencyTrigger())
  })

  # Restrict axes
  abfRestricted = reactive({
    if(is.null(abfNew())){return()}
    
    abf = abfNew()
    

    if(input$triggerRestrict != 0)
      abf = triggerRestrictFun()
    validTimes = abf@longData$time[!is.na(abf@longData$current)]
    plotExtremes <- round(c(min(validTimes, na.rm = T), max(validTimes, na.rm = T),
                           min(abf@longData$current, na.rm = T), max(abf@longData$current, na.rm = T)), digits=2)

    updateSliderInput(session, "horSlider", min = plotExtremes[1], max = plotExtremes[2], value = plotExtremes[1:2])
    updateSliderInput(session, "vertSlider", min = plotExtremes[3],max = plotExtremes[4],value = plotExtremes[3:4])
    
    return(abf)
  })

  # Apply filter (optional)
  abfFiltered = reactive({
    if(is.null(abfRestricted())){return()}
    
    if(input$triggerFilter==0)
      abf = abfRestricted()
    else{
      abf = triggerFilterFun()
      updateSelectInput(session,"traceSelect",selected = "current.filtered")
    }

    return(abf)
  })

  # apply first HMM at press of button
  abfHmm1 = reactive({
    if(is.null(abfFiltered())){return()}

    if(input$triggerHmm1==0)
      abf = abfFiltered()
    else{
      abf = hmm1TriggerFun()
      # Update params
      params = unlist(unname(getpars(abf@model1)))
      for(i in 3:10){
        parName = paste("h1v",i,sep="")
        updateTextInput(session,parName,value=params[i])
      }
    }
    
    return(abf)
  })
  
  abfHmm2 = reactive({
    if(is.null(abfHmm1())){return()}
    
    if(input$triggerHmm2==0)
      abf = abfHmm1()
    else
      abf = hmm2TriggerFun()
    
    return(abf)
  })

  drawMainPlot = reactive({
    if(is.null(abfHmm2())){return()}
    return(plot(abfHmm2(),input$traceSelect))
  })
  
  output$plotScatterheat = renderPlotly({
    if(is.null(abfHmm2())){return()}
    gg = scatterHeat(abfHmm2())
    ggplotly(gg)
  })
  
  output$plotHist = renderPlot({
    if(is.null(abfHmm2())){return()}
    abf = abfHmm2()
    minDuration = round(min(abf@wideData$duration,na.rm=T),digits=2)
    maxDuration = round(max(abf@wideData$duration,na.rm=T), digits=2)
    minCurrent = round(min(abf@wideData$current.res,na.rm=T), digits=2)
    maxCurrent = round(max(abf@wideData$current.res,na.rm=T), digits=2)
    if(is.na(input$bwTime)){
      updateNumericInput(session, "bwTime", value = (maxDuration - minDuration)/50)
      updateSliderInput(session,"rangeTime", min = minDuration, max=maxDuration, value = c(minDuration,maxDuration))
    }
    if(is.na(input$bwCurrent)){
      updateNumericInput(session, "bwCurrent", value = (maxCurrent - minCurrent)/50)
      updateSliderInput(session,"rangeCurrent", min = minCurrent, max=maxCurrent, value = c(minCurrent,maxCurrent))
    }
    return(hist(abfHmm2(),bwDuration = input$bwTime, bwCurrent = input$bwCurrent,
                rangeDuration = input$rangeTime, rangeCurrent = input$rangeCurrent))
  })
  

  output$plotTrace = renderPlotly({
    if(is.null(drawMainPlot())){return()}
    gg = drawMainPlot()
      if (input$toggleGuides==T){
      gg = gg + geom_hline(yintercept = input$vertSlider[1], colour = "red") + 
                geom_hline(yintercept = input$vertSlider[2], colour = "red") +
                geom_vline(xintercept = input$horSlider[1], colour = "red") + 
                geom_vline(xintercept = input$horSlider[2], colour = "red")
    
    
      # Add guides for HMM1 if applicable
      if(input$toggleHmm1== T & !any(is.na(c(input$h1v7,input$h1v8,input$h1v9,input$h1v10)))){
        gg = gg + geom_hline(yintercept = input$h1v7, colour = "#2ca25f") +
          geom_hline(yintercept = input$h1v7-2*input$h1v8, colour = "#99d8c9") +
          geom_hline(yintercept = input$h1v7+2*input$h1v8, colour = "#99d8c9")
        
        gg = gg + geom_hline(yintercept = input$h1v9, colour = "#3182bd") +
          geom_hline(yintercept = input$h1v9-2*input$h1v10, colour = "#9ecae1") +
          geom_hline(yintercept = input$h1v9+2*input$h1v10, colour = "#9ecae1")
      }
    }
    ggplotly(gg)
  })
  
  
  
  # TRIGGERS ---------------------------------------------------------------
  changeFrequencyTrigger = eventReactive(input$frequencyInputTrigger,{
    input$frequencyInput
  })

  triggerRestrictFun = eventReactive(input$triggerRestrict,{
    abf = abfNew()
    abf = restrictAxesFun(abf, "current", minCurrent = input$vertSlider[1], maxCurrent = input$vertSlider[2],
                                          minTime = input$horSlider[1], maxTime = input$horSlider[2])
    return(abf)
  })
  
  triggerFilterFun = eventReactive(input$triggerFilter,{
    abf = abfRestricted()
    if(input$selectFilter == 2)
      abf = applyFilter(abf, "fft", input$fftCutOff)
    else if(input$selectFilter == 3)
      abf = applyFilter(abf, "boxKernel", input$boxBandwidth)
    else if(input$selectFilter == 4)
      abf = applyFilter(abf, "lowPass", input$lowPassW, input$lowPassN)
    
    return(abf)
  })
  
  hmm1TriggerFun = eventReactive(input$triggerHmm1,{
    abf = abfFiltered()
    
    if(input$toggleHmm1==F)
      abf = detectOpenState(abf, IresMethod= input$iresMethod, autoTrash = input$autoTrash1)
    else{
      params = rep(NA,10)
      for (i in 3:10){
        if(input[[paste("h1c",i,sep="")]]==T)
          params[i] = input[[paste0("h1v",i)]]
      }
      params = as.numeric(params)
      abf = detectOpenState(abf, IresMethod= input$iresMethod, autoTrash = input$autoTrash1, 
                            manualParams = params)
    }
    return(abf)
  })
  
  hmm2TriggerFun = eventReactive(input$triggerHmm2,{
    abf = abfHmm1()
    abf = detectBlockedStates(abf,states=input$nStates)
    
    return(abf)
  })
                            
})
