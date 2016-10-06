# Class definition and initialization -----------------------------------------------------------------------------

#' AbfData Class
#'
#' A class storing information on current data from some ABF file,
#' along with some metadata for this data i.e. original file name,
#' sampling frequency, applied filter and number of desired states.
#'
#' For convenience during analysis, data is saved in "long format"
#' (one row per measurement) and "wide format' (one row per event).
#'
#' @slot fileName Original name of the abf file from which trace was loaded
#' @slot samplingFrequency Sampling frequency to which the trace data was set after it was loaded
#' @slot rawSamplingFrequency Original sampling frequency of trace data in abf-file
#' @slot filter Filter type that was applied to data, or NA if none was applied
#' @slot openPoreState Number denoting the automatically detected open pore state
#' @slot model1 Depmix-model used to discern blocked state from open pore states
#' @slot model2 Depmix-model used to discern different blocked states
#' @slot longData Data frame containing trace data in a long format, i.e. one row per measured time point
#' @slot wideData Data frame containing trace data in a wide format, i.e. one line per event
#'
#' @importClassesFrom depmixS4 depmix.fitted
#' @export
#'
#'
setClass(
  Class="AbfData",
  representation = representation(
    fileName = "character",
    samplingFrequency = "numeric",
    rawSamplingFrequency = "numeric",
    filter = "factor",
    openPoreState = "numeric",
    model1 = "depmix.fitted",
    model2 = "depmix.fitted",
    longData = "data.frame",
    wideData = "data.frame"
  )
)

#' Initialize AbfData object
#'
#' Initialize a new abfData object with given abf file name. Automatically loads the data into
#' the object.
#'
#' @param fileName Path and file name, referrring to the abf file on disk
#' @param samplingFrequency Optionally, sampling frequency in Hz to which trace data should be adjusted
#' @return A new abfData object that stores the trace from the abf-file at the given path, adjusted sampling frequency if applicable, and metadata.
#'
#' @importFrom abf2 abfload
#' @importFrom methods new
#' @export
newAbfData = function(fileName, samplingFrequency){
  return(new(Class="AbfData",fileName=fileName, samplingFrequency=samplingFrequency))
}

#' Initialize AbfData object
#'
#' Initialize a new abfData object with given abf file name. Automatically loads the data into
#' the object. A user friendlier function, "newAbfData", is now in place.
#'
#' @param .Object refers to the object being created
#' @param fileName Path and file name, referrring to the abf file on disk
#' @param samplingFrequency Optionally, sampling frequency in Hz to which trace data should be adjusted
#' @return A new abfData object that stores the trace from the abf-file at the given path, adjusted sampling frequency if applicable, and metadata.
#'
#' @importFrom abf2 abfload
#'
#' @docType methods
#' @rdname initialize-methods
setGeneric(name="initialize", def = function(.Object, fileName, samplingFrequency){standardGeneric("initialize")},
           simpleInheritanceOnly = T)

#' @rdname initialize-methods
#' @aliases initialize,character,character,numeric-method

setMethod(
  f="initialize",
  signature="AbfData",
  definition=function(.Object,fileName,samplingFrequency){
    #require(abf2,quietly =T)
    .Object@fileName = fileName
    .Object@longData = setLongDataFun(fileName)
    .Object@rawSamplingFrequency = round(length(.Object@longData$time)/max(.Object@longData$time))
    if(!missing(samplingFrequency))
      .Object@longData = adjustSamplingFrequencyFun(.Object@longData, samplingFrequency)
    .Object@samplingFrequency= round(length(.Object@longData$time)/max(.Object@longData$time))
    return(.Object)
  }
)


# Load trace data in longData dataframe format
setLongDataFun = function(fileName){
  if (is.na(fileName))
    abfRaw = abfload()
  else
    abfRaw = abfload(fileName)
  longData = data.frame(time= abfRaw$s,current=abfRaw$traces[1,])           # Prepare data frame from raw data
  return(longData)
}

# Adjust sampling frequency in a longData dataframe
adjustSamplingFrequencyFun = function(longData, samplingFrequency){
  stepSize = 1/(samplingFrequency * (longData$time[2]-longData$time[1]))
  longData = longData[seq(1,nrow(longData),stepSize),]
  return(longData)
}

# Noise filtering methods ------------------------------------------------------------------------------------------



#' Restrict axes
#'
#' Restrict current to the range in which data points should be recognized. May solve issues
#' in the convergence of HMM state detection. Any combination of minimum and maximum values is allowed.
#'
#' @param abfData object of AbfData type
#' @param whichSequence Sequence on which axis restriction should be based. Accepted options are "current", "current.filtered" and "current.res"
#' @param minCurrent Minimum current value. points with lower current value are removed.
#' @param maxCurrent Maximum current value. points with higher current value are removed.
#' @param minTime Time in s starting from which datapoints should be retained.
#' @param maxTime Time in s up until which datapoints should be retained.
#'
#' @return Abf-file in which the given restrictions on time and min/max current have been applied
#'
#' @export
#'
#' @docType methods
#' @rdname restrictAxes-methods
setGeneric(name="restrictAxes", def = function(abfData, whichSequence, minCurrent, maxCurrent, minTime, maxTime){standardGeneric("restrictAxes")})
restrictAxesFun = function(abfData, whichSequence, minCurrent, maxCurrent, minTime, maxTime){
  cv = as.vector(unlist(abfData@longData[whichSequence]))
  tv = as.vector(unlist(abfData@longData$time))
  cvSelect = rep(T,length(cv))
  if(!missing(minCurrent))
    cvSelect = cvSelect & cv>minCurrent
  if(!missing(maxCurrent))
    cvSelect = cvSelect & cv<maxCurrent
  if(!missing(minTime))
    cvSelect = cvSelect & tv>minTime
  if(!missing(maxTime))
    cvSelect = cvSelect & tv<maxTime
  cvSelect[is.na(cvSelect)] = F
  abfData@longData[!cvSelect,-1] = NA
  return(abfData)
}

#' @rdname restrictAxes-methods
#' @aliases restrictAxes,AbfData-method

setMethod(f="restrictAxes", signature=c(abfData="AbfData"),
          def=restrictAxesFun)

#' Apply filter
#'
#' Filter white noise from the data, using the given filter type. If state is given, filter only white noise at that level.
#'
#' @param abfData Object of AbfData type, containing data that is to be filtered
#' @param filterType Type of filter that is to be applied. Depending on the given filter type, some additional arguments may be required.
#' Accepted filter types: "fft" (fast fourier transform), "lowPass" (low pass Gaussian), uses signal package,"boxKernel" (moving average kernel).
#' @param cutOff Required for fft-filter. Determines strictness of filtering (0-1, 1= no filter, 0 = flat line)
#' @param w Required for lowPass-filter. Determine  reciprocal of std.dev (so larger is wider window)
#' @param n Required for lowPass-filter. Determines number of coefficients
#' @param bandwidth Required for boxKernel-filter. Determine bandwidth used for moving average (greater = smoother)
#' @param state Optionally, specific state at which filtering should take place
#'
#' @return AbfData-object in which  the column current.filtered, containing filtered current values is added to longData and filter slot is filled with name of used filter, but which in all other ways is identical to the given abfData object.
#'
#' @importFrom signal gausswin filter
#' @export
#'
#' @docType methods
#' @rdname applyFilter-methods
setGeneric(name="applyFilter", def=function(abfData, filterType, cutOff, w, n, bandwidth, state){standardGeneric("applyFilter")})

# Filter white noise from the data, using the given filter type.
applyFilterFun = function(abfData, filterType){
  if(!missing(state)){
    currentVector = abfData@longData$current[!is.na(abfData@longData$state) & abfData@longData$state==state]
  } else{
    currentVector = abfData@longData$current[!is.na(abfData@longData$current)]
    abfData@longData["current.standardized"] = scale(abfData@longData$current)
  }

  if(filterType=="fft" & exists("cutOff"))
    currentVectorFiltered = fftFilterFun(currentVector, cutOff)
  else if(filterType=="lowPass" & exists("w") & exists("n"))
    currentVectorFiltered = lowPassFilterFun(currentVector,w,n)
  else if(filterType=="boxKernel" & exists("bandwidth"))
    currentVectorFiltered = lowPassFilterFun(currentVector, bandwidth)
  else{
    warning("Arguments missing, no filter applied")
    return(currentVector)
  }

  if(!missing(state))
    abfData@longData$current.filtered[!is.na(abfData@longData$state) & abfData@longData$state==state] = currentVectorFiltered
  else{
    abfData@longData["current.filtered"] = NA
    abfData@longData$current.filtered[!is.na(abfData@longData$current)]= currentVectorFiltered
  }
  return(abfData)
}

# home-brew fft filter, to be called from applyFilterFun
fftFilterFun = function(currentVector, cutOff){
  current.ft = fft(currentVector)
  cutOff = length(current.ft) * cutOff*0.5
  current.ft[cutOff:length(current.ft)-cutOff] = 0 + 0i
  current.ft = fft(current.ft, inverse =T)/length(current.ft)
  current.ft = as.vector(Re(current.ft))
  return(current.ft)
}

# Low-pass gaussian filter, to be called from applyFilterFun
lowPassFilterFun = function(currentVector, n, w){
  # require(signal,quietly=T)
  lp = gausswin(n = n , w = w)
  current.filtered = filter(filt = lp, a = rep(1/n, n), x = currentVector)
  return(current.filtered)
}

# Box kernel filter, to be called from applyFilterFun
boxKernelFilterFun = function(currentVector, bandwidth){
  current.filtered = unlist(ksmooth(currentVector, 1:length(currentVector), kernel = "box",bandwidth=bandwidth)[2],use.names=F)
  return(currentVector)
}

#' @rdname applyFilter-methods
#' @aliases applyFilter,AbfData,character-method
setMethod(f="applyFilter", signature= c(abfData="AbfData",filterType="character"),
          def=applyFilterFun)
# HMM methods ------------------------------------------------------------------------------------------------------

#' Detect open pore state
#'
#' Method for discerning open pore state from blocked pore state, using a Hidden Markov Model (HMM). Optionally, provide a vector containing desired values for HMM parameters:
#' - Initial state probabilities (elements 1 and 2)
#' - Transition probabilities, filling the transition matrix row-wise starting top left (elements 3 to 6)
#' - Intercept and standard deviation for state 1 (elements 7 and 8)
#' - Intercept and standard deviation for state 2 (elements 9 and 10)
#' If no parameters are provided, or if some parameters are set to NA, these parameters are filled in automatically using Maximum Likelihood (ML) estimations.
#' After determining open and closed pore states, Ires is automatically caluculated and added to the longData dataframe in the abf-object. By default, Ires Is calculated by division
#' of current through the last open pore state. If open pore states last very short however, this may cause unrealistic fluctuation in the Ires value and one of the other options (
#' blocked pore state averaged per period of 10 minutes or averaged overall) may be preferred.
#' Lastly, the autoTrash option runs a different HMM prior to the one used to discern open and closed pore states and discards the state with highest standard deviation.
#' This may eliminate some extreme peaks deemed artefacts, although restricting the allowed range of currents was found to work better.
#'
#' This method is built around the depmix function from the depmixS4-package. As is stated in the description of this package, ML estimations
#' are made using using the nnet.default routine (nnet package) for initial state and transition model probabilities and a glm for the response model.
#' See the documentation of the DepmixS4-package (\url{https://cran.r-project.org/web/packages/depmixS4/depmixS4.pdf}) for further information.
#'
#' @param abfData an object of class AbfData
#' @param manualParams a vector of length ten, containing required parameters for the HMM
#' @param autoTrash set to TRUE to attempt automatically filtering out extremely high and low values (see description)
#' @param IresMethod Select which method should be used to calculate Ires. Accepted options are "perEvent" (the default), "periodically" (per 10 mins.) and "overall"
#' @param dtmax Optional parameter, defining the maximum time a state may last in s. States lasting longer are removed
#' @param dtmin DO NOT USE, NOT READY YET! Optional parameter, defining the minimum time a state may last in s. If adjacent states are the same, these are joined together.
#' @return an abfData-object, with information on open and closed pore state stored in the longData and wideData dateframes. The model used for determination of those states is stored in
#' the 'model1' variable, as a depmix-object.
#'
#' @importFrom depmixS4 depmix fit getpars setpars posterior
#' @export
#'
#' @docType methods
#' @rdname detectOpenState-methods
setGeneric(name="detectOpenState", def= function(abfData,manualParams, dtmin,dtmax,autoTrash,IresMethod){standardGeneric("detectOpenState")})

# Automatically discern open state and blocked state using HMM or
# use manually provided algorithms to guide the algorithm.
# Todo: consider, is the NoNA step still necessary here?
# Todo: extract function for general HMM application (after writing substate
#       recognition code)
detectOpenStateFun = function(abfData, manualParams, dtmin, dtmax,autoTrash,IresMethod="perEvent"){
  # require(depmixS4,quietly=T)

  # If data has been filtered, fit model on the filtered data. If not, fit on
  # the regular current data.
  if(!is.null(abfData@longData$current.filtered))
    abf.df = data.frame(current = abfData@longData$current.filtered)
  else
    abf.df = data.frame(current = abfData@longData$current)

  # if autoTrash=T, detect trash state and set trash current values to NA
  if(!missing(autoTrash)){
    if(autoTrash==T)
      abf.df = runHmmTrashFun(abf.df)
  }

  # Create and fit HMM model
  abf.dfNoNAs = data.frame(current=abf.df$current[!is.na(abf.df$current)])
  abf.mod = depmix(current~1, family = gaussian(), nstates=2, data=abf.dfNoNAs)
  abfData@model1 = fit(abf.mod,verbose=F)

  # If required levels are provided, refit model using given parameters
  if(!missing("manualParams")){
    params = c(unlist(getpars(abfData@model1)))
    params[!is.na(manualParams)] = manualParams[!is.na(manualParams)]
    abfData@model1 = setpars(abfData@model1, params)
    conpat = as.numeric(is.na(manualParams))
    abfData@model1 = fit(abfData@model1, equal = conpat)
  }

  # Retrieve estimated states, enter them in the abf.df dataframe.
  abf.probs = posterior(abfData@model1)[,1]
  abf.df$current[!is.na(abf.df$current)] = abf.probs

  # Guess which state is the open current state. This should be the state with highest
  # absolute mean. Store this state in the variable openPoreState in the abfData object.
  allPars = unlist(getpars(abfData@model1))
  states = 2 # could have immediately filled in c(7,9) as indices, but may want to generalized code later
  interceptIndices = seq(states+states^2+1, states+states^2+2*states-1, 2)
  interceptValues = allPars[interceptIndices]
  abfData@openPoreState = which(interceptValues==min(interceptValues))


  # Construct wide dataset, filter out min and max event length (if supplied)
  abfData@wideData = setWideDataFun(abfData, abf.df$current, dtmin, dtmax)

  # Add state information to longData
  abfData@longData["state"] = NA
  for(n in 1:nrow(abfData@wideData)){
    abfData@longData$state[abfData@longData$time>=abfData@wideData$start[n]&
                             abfData@longData$time<abfData@wideData$start[n] + abfData@wideData$duration[n]] = abfData@wideData$state[n]
  }

  # Create new column in longData for Ires (current.res), by taking the detected
  # open-pore current and dividing subsequent currents through it up to the next
  # open-pore current.
  abfData@longData["current.res"] = NA
  if(IresMethod == "perEvent"){
    for(n in 1:nrow(abfData@wideData)){
      if(abfData@wideData$state[n] == abfData@openPoreState){
        # Accounts for the cases in which sequence does not start with open pore state. Then stretch
        # before the first open pore state must be treated as well.
        if(n==2){
          si = abfData@longData$time< abfData@wideData$start[n] + abfData@wideData$duration[n]
          abfData@longData$current.res[si] = abfData@longData$current[si] / abfData@wideData$current[n]
        }
        si = abfData@longData$time>=abfData@wideData$start[n]
        abfData@longData$current.res[si] = abfData@longData$current[si] / abfData@wideData$current[n]
      }
    }
  } else if(IresMethod=="overall"){
    meanOpenCurrent = mean(abfData@longData$current[abfData@longData$state==abfData@openPoreState],na.rm=T)
    abfData@longData$current.res = abfData@longData$current/meanOpenCurrent
  } else if(IresMethod=="periodically"){
    periods = seq(0,max(abfData@longData$time),600)
    for(n in periods){
      cp = abfData@longData$time>=n
      openMean = mean(abfData@longData$current[cp & abfData@longData$state==abfData@openPoreState],na.rm=T)
      abfData@longData$current.res[cp] = abfData@longData$current[cp] / openMean
    }
  }

  # Add Ires to wideData datframe
  abfData@wideData["current.res"] = NA
  for (n in 1:nrow(abfData@wideData)){
    longIndex = abfData@longData$time >= abfData@wideData$start[n] & abfData@longData$time < abfData@wideData$start[n] + abfData@wideData$duration[n]
    abfData@wideData$current.res[n] = mean(abfData@longData$current.res[longIndex], na.rm = T)
  }

  return(abfData)
}

#' @rdname detectOpenState-methods
#' @aliases detectOpenState,AbfData-method

setMethod(f="detectOpenState",signature=c(abfData = "AbfData"),
          definition = detectOpenStateFun)

# Detect and filter out data suspected to be faulty blips and voltage changes
# applied to solve blockages. To do so, construct and fit a hidden markov model and remove
# state with highest standard deviation.
# param       abf.df
#             dataframe containing current sequence data only, as given by one of the
#             runHmmFun# functions
# Return      input dataframe, modified such that currents in trash state are now given as NA
# Note        Only intended to be called from other functions!
runHmmTrashFun = function(abf.df){
  # require(depmixS4,quietly=T)
  ns = 4
  abf.mod = depmix(abf.df$current~1, family = gaussian(), nstates=ns, data = abf.df)
  abf.fit = fit(abf.mod,verbose=F)
  abf.post = posterior(abf.fit)[,1]

  sdValues = abs(unlist(getpars(abf.fit))[seq(ns+ns^2+1,ns+ns^2+2*ns-1,2)])
  trashState = which(sdValues==max(sdValues))

  abf.df$current[abf.post==trashState] = NA
  return(abf.df)
}

# Construct wide data dataframe, using the long dataframe and the states per time point as
# estimated by a HMM.
#
# Param     abfData
# Param     hmm.vector
# Param     dtmin
# Param     dtmax
# Return    Current data in wide format
setWideDataFun = function(abfData,hmm.vector, dtmin, dtmax){
  # Pre-allocate memory for 100,000 state-switches (allows overshooting but will be much faster
  # initially)
  N = 100000
  wideData = data.frame(state=rep(NA_integer_,N),
                        start = rep(NA,N),
                        duration=rep(NA,N),
                        current=rep(NA,N),
                        current.filtered = rep(NA,N),
                        sd = rep(NA,N),
                        stringsAsFactors = F)

  # Loop over rows. If state changes: enter new row with state, start time, duration and mean current
  # NOTE: NAs (discarded data) are omited as follows:
  # - at first NA, mark end current event
  # - At subsequent NAs, continue to next.
  # - At start of new event, set new start state (so do not close off set of NA states as event)
  startState = 1
  event = 1
  for (n in 2:length(hmm.vector)){
    if(is.na(hmm.vector[n]) & is.na(hmm.vector[n-1]))
      next
    else if(!is.na(hmm.vector[n]) & is.na(hmm.vector[n-1]))
      startState=n
    else if(is.na(hmm.vector[n]) & !is.na(hmm.vector[n-1]) |
            hmm.vector[n]!=hmm.vector[startState] |
            n==length(hmm.vector)){
      wideData$state[event] = hmm.vector[startState]
      wideData$start[event] = abfData@longData$time[startState]
      wideData$duration[event] = abfData@longData$time[n]-abfData@longData$time[startState]
      wideData$current[event] = mean(abfData@longData$current[startState:n-1],na.rm=T)
      wideData$sd[event] = sd(abfData@longData$current[startState:n-1],na.rm=T)
      #wideData$current.filtered[event] = mean(abfData@longData$current.filtered[startState:n-1],na.rm=T)
      event = event + 1
      if(!is.na(hmm.vector[n]))
        startState = n
    }
  }

  #Remove rows of pre-allocated dataframe that were left unfilled, if any.
  wideData = wideData[!is.na(wideData$state),]

  # If min event duration has been specified, filter out abnorally short events.
  if(!missing(dtmin)){
    wideData = wideData[wideData$duration>dtmin,]
    # # store data temporarily in other dataframe
    # wideStore = wideData
    # wideData = matrix(ncol=ncol(wideData), rep(NA, prod(dim(wideData))))
    # j = 1
    # for (i in 1:nrow(wideData)){
    #   if(wideStore$duration[i]>dtmin){
    #     wideData[j,] = wideStore[i,]
    #     j = j+1
    #   } else{first=F}
    #
    #   if(wideStore$duration[i+1]>dtmin){
    #     wideData[j,"duration"] = sum(wideData[j,"duration"])
    #   }
    # }
    #
    # shortEvents = wideData[wideData$duration>dtmin,"time"]
    #
    # # "Sows" events together if event is removed and adjacent events are of the same state.
    # first = T
    # i = 1
    # for (i in 1:length(wideData)){
    #   if(first==T){
    #     startEvent = i
    #   }
    #   if(!any(i+1==shortEvents)){
    #     if(wideData$state[startEvent]==wideData$state[i]){
    #       wideData$duration = sum(wideData$duration[startEvent:i])
    #       removeEvents = c(removeEvents, i+1)
    #     }
    #     first=T
    #   } else{
    #     first=F
    #   }
    # }
    # wideData = wideData[-removeEvents,]
  }

  # If max event duration has been specified, fiter out abnormally long events.
  # NOTE: unblocked state must be exempt of this filtering; this must be allowed
  # to last as long as it does.
  if(!missing(dtmax)){
    wideData = wideData[wideData$state== abfData@openPoreState | wideData$duration<dtmax,]
  }

  return(wideData)
}


#' Detect blocked states
#'
#' A method for the determination of a given number of blocked sub-states within a blocked current state as detected using the \code{\link{detectOpenState}} method.
#' Optionally, provide a vector containing desired values for HMM parameters in this order:
#' - Initial state probabilities
#' - Transition probabilities, filling the transition matrix row-wise starting top left
#' - Intercept and standard deviation per state
#' If no parameters are provided, or if some parameters are set to NA, these parameters are filled in automatically using Maximum Likelihood (ML) estimations. Furthermore,
#' a "split state" may be defined, within which the HMM will look for the given number of substates. Any previously defined state qualifies as a potential split state.
#' Lastly, the autoTrash option runs a different HMM prior to the one used to discern open and closed pore states and discards the state with highest standard deviation.
#' This may eliminate some extreme peaks deemed artefacts, although restricting the allowed range of currents was found to work better.
#'
#' This method is built around the depmix function from the depmixS4-package. As is stated in the description of this package, ML estimations
#' are made using using the nnet.default routine (nnet package) for initial state and transition model probabilities and a glm for the response model.
#' See the documentation of the DepmixS4-package (\url{https://cran.r-project.org/web/packages/depmixS4/depmixS4.pdf}) for further information.
#'
#' @param abfData An object of class AbfData, in which blocked and open pore state have been discerned using the detectOpenState method
#' @param states Integer denoting the number of blocked sub-states that should be found.
#' @param manualParams A vector containing required parameters for the HMM, in the following order, for n states: prior model (n parameters), transition model (n^2), response model (2n)
#' @param autoTrash set to TRUE to attempt automatically filtering out extremely high and low values (see description)
#' @param splitState Optionally, enter a state within which substates should be recognized.
#' @param dtmax Optional parameter, defining the maximum time a state may last in s. States lasting longer are removed
#' @param dtmin DO NOT USE, NOT READY YET! Optional parameter, defining the minimum time a state may last in s. If adjacent states are the same, these are joined together.
#' @return An abfData-object in which either the blocked state or the given state has been split into the given number of substates.
#'
#' @importFrom depmixS4 depmix fit getpars setpars posterior
#' @export
#'
#' @docType methods
#' @rdname detectBlockedStates-methods
setGeneric(name="detectBlockedStates", def= function(abfData,states,splitState,manualParams,dtmin,dtmax,autoTrash){standardGeneric("detectBlockedStates")})

# Recognize substates in any blocked pore states that are defined in the ABF-file (i.e. when no
# split state has been defined). Takes predicted split state and calles detectBlockedStatesFun, with
# blocked state as splitState.
detectBlockedStatesFun2 = function(abfData,states,manualParams,dtmin,dtmax,autoTrash){
  uniqueStates = unique(abfData@longData$states)
  if(abfData@openPoreState==1)
    blockedState = 2
  else
    blockedState = 1
  abfData@longData$state[abfData@longData$state != abfData@openPoreState] = blockedState
  abfData@wideData$state[abfData@wideData$state != abfData@openPoreState] = blockedState
  return(detectBlockedStatesFun(abfData,states,splitState=blockedState,manualParams,dtmin,dtmax,autoTrash))
}

# Recognize substates in a given blocked pore state.
detectBlockedStatesFun = function(abfData, states, splitState, manualParams, dtmin, dtmax, autoTrash){
  # require(depmixS4,quietly=T)
  if(!is.null(abfData@longData$current.resFiltered))
    abf.df = data.frame(current = abfData@longData$current.resFiltered[!is.na(abfData@longData$state) & abfData@longData$state==splitState])
  else
    abf.df = data.frame(current = abfData@longData$current.res[!is.na(abfData@longData$state) & abfData@longData$state==splitState])


  # if autoTrash=T, detect trash state and set trash current values to NA
  if(!missing(autoTrash))
    if(autoTrash==T)
      abf.df = runHmmTrashFun(abf.df)

  # Construct and fit model
  abf.dfNoNAs = data.frame(current=abf.df$current[!is.na(abf.df$current)])
  abf.mod = depmix(current~1, family = gaussian(), nstates=states, data=abf.dfNoNAs)
  abfData@model2 = fit(abf.mod,verbose=F)

  # If required levels are provided, refit model using given parameters
  if(!missing("manualParams")){
    params = c(unlist(getpars(abf.mod)))
    params[!is.na(manualParams)] = manualParams[!is.na(manualParams)]
    abfData@model2 = setpars(abfData@model2, params)
    conpat = as.numeric(is.na(manualParams))
    abfData@model2 = fit(abfData@model2, equal = conpat)
  }



  # If required levels are provided, refit model using given parameters
  if(!missing("manualParams")){
    params = c(unlist(getpars(abfData@model2)))
    params[!is.na(manualParams)] = manualParams[!is.na(manualParams)]
    abfData@model2 = setpars(abfData@model2, params)
    conpat = as.numeric(is.na(manualParams))
    abfData@model2 = fit(abfData@model2, equal = conpat)
  }

  # Retrieve estimated states, enter them in the abf.df dataframe
  abf.probs = posterior(abfData@model2)[,1]


  #Make sure that values of new states are not equal to states in old sequence
  abf.probs = abf.probs + splitState*10

  # Replace value of split-state with new state values in two steps (i.e. reverse
  # of previous steps)
  abf.df$current[!is.na(abf.df$current)] = abf.probs
  abf.allProbs = abfData@longData$state
  abf.allProbs[!is.na(abf.allProbs) & abf.allProbs==splitState] = abf.df$current
  abf.probs = abf.allProbs

  # Update wide dataset
  abfData@wideData = setWideDataFun(abfData, abf.probs, dtmin, dtmax)
  abfData@longData["state"] = NA

  # Add state information to longData
  for(n in 1:nrow(abfData@wideData)){
    abfData@longData$state[abfData@longData$time>=abfData@wideData$start[n]&
                             abfData@longData$time<abfData@wideData$start[n] + abfData@wideData$duration[n]] = abfData@wideData$state[n]
  }

  # Update Ires to wideData dataframe
  abfData@wideData["current.res"] = NA
  for (n in 1:nrow(abfData@wideData)){
    longIndex = abfData@longData$time >= abfData@wideData$start[n] & abfData@longData$time < abfData@wideData$start[n] + abfData@wideData$duration[n]
    abfData@wideData$current.res[n] = mean(abfData@longData$current.res[longIndex], na.rm = T)
  }

  return(abfData)
}

#' @rdname detectBlockedStates-methods
#' @aliases detectBlockedStates,AbfData,numeric-method

setMethod(f="detectBlockedStates", signature = c(abfData = "AbfData",splitState="numeric"),
          definition = detectBlockedStatesFun)
#' @rdname detectBlockedStates-methods
#' @aliases detectBlockedStates,AbfData-method

setMethod(f="detectBlockedStates", signature = c(abfData = "AbfData"),
          definition = detectBlockedStatesFun2)

# Plotting methods ---------------------------------------------------------------------------------------

#' Draw dwell time and current histograms
#'
#' Overrride of histogram function, automatically produces duration and current histograms of all defined states in a given abfData object.
#'
#' Uses the multiplot function as described in the R cookbook: \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#'
#' @param x abfData-object of which states should be summarized in histograms
#' @param bwDuration Binwidth of histograms on duration. by default, binwidth is set such that entire range of values is divided in 50 bins.
#' @param bwCurrent Binwidth of histograms on Ires. by default, binwidth is set such that entire range of values is divided in 50 bins.
#' @param rangeDuration Vector of two values, denoting min and max duration displayed respectively. By default, entire range of duration values is covered.
#' @param rangeCurrent Vector of two values, denoting min and max Ires displayed respectively. By default, entire range of Ires values is covered.
#' @return A ggplot-object containing the histogram. May be stored, modified using additional ggplot arguments or displayed immediately.
#'
#' @import ggplot2
#' @import grid
#' @export
#'
#' @rdname hist-methods
#' @aliases hist,AbfData-method
setMethod(f="hist",signature="AbfData",
          def=function(x, bwDuration, bwCurrent, rangeDuration, rangeCurrent){

  if(missing(bwDuration))
    bwDuration = (max(x@wideData$duration, na.rm=T) - min(x@wideData$duration,na.rm=T))/50

  ggDuration = ggplot(data=x@wideData,aes(x=duration)) +
    geom_histogram(aes(y=..density..,fill= as.factor(state)),color="grey80",binwidth= bwDuration) +
    facet_grid(state~.) + guides(fill=F) + labs(x="duration (s)")
  if(!missing(rangeDuration))
    ggDuration = ggDuration + coord_cartesian(xlim = rangeDuration)

  if(missing(bwCurrent))
    bwCurrent = (max(x@wideData$current.res, na.rm=T) - min(x@wideData$current.res,na.rm=T))/50

  ggCurrent = ggplot(data=x@wideData,aes(x=current.res)) +
    geom_histogram(aes(y=..density..,fill= as.factor(state)),color="grey80",binwidth= bwCurrent) +
    facet_grid(state~.) + guides(fill=F) + labs(x="Ires (pA)")
  if(!missing(rangeCurrent))
    ggCurrent = ggCurrent + coord_cartesian(xlim = rangeCurrent)

  return(multiplot(ggDuration,ggCurrent,cols=2))
})

# Multiple plot function from the R cookbook, required for override of hist.
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Plot traces
#'
#' Plot sequences stored in an abfData-object and state info (if available).
#'
#' @param x AbfData-object containing the traces (and optionally state information) that are to be plotted.
#' @param y name of the trace that is to be plotted, or a vector of names. Accepted options are "current", "current.filtered",
#' "current.filtered" and "current.res".
#' @return A ggplot-object containing the trace graph. May be stored, modified using additional ggplot arguments or displayed immediately.
#'
#' @import ggplot2
#' @export
#'
#' @rdname plot-methods
#' @aliases plot,AbfData-method
setMethod(f="plot",signature="AbfData",
          def = function(x,y){
  # require(ggplot2,quietly=T)

  ggData = ggplot(data.frame(time=x@longData),aes(x=time))

  # plot any desired combination of current graphs
  if(any(y=="current")){
    height = mean(x@longData$current, na.rm=T)
    ggData = ggData + geom_line(aes(y=current),data=x@longData, na.rm=T)
  }
  if(any(y=="current.standardized")){
    height = mean(x@longData$current.standardized, na.rm=T)
    ggData = ggData + geom_line(aes(y=current.standardized),data=x@longData, na.rm=T)
  }
  if(any(y=="current.filtered")){
    height = mean(x@longData$current.filtered, na.rm=T)
    ggData = ggData + geom_line(aes(y=current.filtered),data=x@longData, na.rm=T)
  }
  if(any(y=="current.res")){
    height = mean(x@longData$current.res, na.rm=T)
    ggData = ggData + geom_line(aes(y=current.res),data=x@longData, na.rm=T)
  }
  if(any(y=="current.resFiltered")){
    height = mean(x@longData$current.resFiltered, na.rm=T)
    ggData = ggData + geom_line(aes(y=current.resFiltered),data=x@longData, na.rm=T)
  }

  tr = range(x@longData$time[!is.na(x@longData$current)], na.rm=T)
  # ggData = ggData + coord_cartesian(xlim = c(tr[1],tr[2]))

  # if state info exists, draw colorbar for states
  if(!is.null(x@wideData$state)){
    ggData = ggData + geom_segment(data=x@wideData,
                                   aes(x=start, xend=start+duration,
                                       y=height,yend=height, color=as.factor(state)),
                                   alpha=0.5,size=25,na.rm=T)
  }
  return(ggData)
})


#' Scatter-heat plot (amplitude (s.d.) vs Ires)
#'
#' Method for making a scatter-heat plot (amplitude vs Ires), using the information stored in an
#' abfData object.
#'
#' @param abfData AbfData object containing the information that is to be displayed in the scatter-heat plot.
#' @return A ggplot-object containing the scatter-heat plot. May be stored, modified using additional ggplot arguments or displayed immediately.
#'
#' @import ggplot2
#' @export
#'
#' @docType methods
#' @rdname scatterHeat-methods
setGeneric(name="scatterHeat", def = function(abfData){standardGeneric("scatterHeat")})
scatterHeatFun = function(abfData){
  # require(ggplot2, quietly=T)
  index = abfData@wideData$state != abfData@openPoreState
  df = data.frame(Ires= abfData@wideData$current.res[index],
                  Amplitude.S.D. = abfData@wideData$sd[index],
                  state = as.factor(abfData@wideData$state[index]))
  gg1 = ggplot(data=df, aes(x=Ires, y=Amplitude.S.D.))
  gg1 = gg1 + stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',na.rm=T) +
    geom_point(aes(shape=state), color="grey50", size=1,na.rm=T) +
    scale_fill_continuous(low="green",high="red") +
    guides(alpha="none",fill="none") +
    theme(panel.grid.major= element_blank(),panel.grid.minor= element_blank())
  return(gg1)
}

#' @rdname scatterHeat-methods
#' @aliases scatterHeat,AbfData-method

setMethod(f="scatterHeat", signature = c(abfData = "AbfData"),
          definition = scatterHeatFun)

# Save/load HMM parameters --------------------------------------------------------

#' Save HMM Models
#'
#' Save the HMM model stored in an abf object, as a file to be re-used for
#' similar sequences, requiring similar analysis.
#'
#' @param abfData The abfData-object from which the models are to be saved.
#' @param fileName Name of the file in which the data should be stored.
#'
#' @export
#'
saveHmmModels = function (abfData, fileName){
  models = list(model1 = abfData@model1, model2 = abfData@model2)
  save(models, file=fileName)
}


#' Load HMM Models
#'
#' Load an existing HMM model to the current data in an object.
#'
#' @note Running loaded existing models not supported yet!
#'
#' @param abfData The abfData object into which the model should be loaded
#' @param fileName The file name (path included if not in home directory), from which HMM models should be loaded.
#' @return AbfData-object with the model1 and model2 variables set to the models in the file denoted by the given file name.
#'
#' @export
#'
loadHmmModels = function(abfData, fileName){
  load(fileName)
  abfData@model1 = models$model1
  abfData@model2 = models$model2
  return(abfData)
}
# Shiny UI ----------------------------------------------------------
#' Shiny State Finder
#'
#' Display a Shiny-based graphical user interface, for a more intuitive application of this package.
#'
#' @importFrom shiny runApp
#' @export
#'
shinyStateFinder = function(){
  appDir = system.file("ShinyStateFinder",package="NanoporeStateFinder")
  shiny::runApp(appDir,display.mode="normal")
}
