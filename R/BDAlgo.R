#' @import inflection
#' @title The Bloom Detecting Algorithm
#' @aliases BDAlgo
#' @description The Bloom Detecting Algorithm enables the detection of  blooms within a time series of species abundance and extracts 22 phenological variables.
#' See details for more information.
#'
#' @usage BDAlgo(data, threshold=c(4,4), perc_of_peak=0.85, nbr_days=300, date_col= 2,
#' station_col=1, Sp = c("Sp1", "Sp2"), min_Log = 2, SP_label = c("Species1", "Species2"),
#' Log=TRUE, PDF=FALSE, saving_path=NULL)
#'
#' @param data, a data frame table containing the species abundance time series with the station and date (dd/mm/yyyy).
#' @param threshold a numeric vector, the abundance threshold from which a detected peak can be considered a bloom.
#' @param perc_of_peak a numeric value between 0 and 1, that creates an upper lower threshold under which the lower point after the peak could be considered the end of the bloom. See the details for more information.
#' @param nbr_days a numeric value corresponding to the maximum of days length of a bloom.
#' @param date_col a numeric value corresponding to the column of the 'date' variable.
#' @param station_col a numeric value corresponding to the column of the 'station' variable.
#' @param Sp a character vector corresponding to the species column in the data to which you wish to apply the function.
#' @param min_Log a numeric value that is used as a lower threshold for the abundance of the data. See details for more information.
#' @param SP_label a character vector corresponding to the species label to use for the output graph.
#' @param Log TRUE or FALSE if you want to log(x+1) transform the abundance data.
#' @param PDF TRUE or FALSE if you wish to save each graph in a single PDF.
#' @param saving_path a character vector used as a directory path to save the output graphs and rda files.
#' @return The BD_Algo function returns:
#'
#'  1. a graph per species and per station.
#'
#'  2. a list containing the species (Sp) and for each station the following data:
#'
#' - \code{smooth_spline}: the results of the smooth spline. See \link[stats:smooth.spline]{smooth.spline} for more information on the return values.
#'
#' - \code{conf_intervall}:  the data of the confidence intervals.
#'
#' - \code{all_date}: the character vector with all the dates used in the smooth spline.
#'
#' - \code{all_bloom}: the phenological data frame of each bloom with the timing variables (DBS, DMF, DMA, DMM, DBE) corrected.
#'
#' - \code{all_bloom_date}: the raw phenological data frame of each bloom.
#'
#' Warnings may occur during the smooth spline applications.
#'
#' @author Stephane Karasiewicz, \email{skaraz.science@gmail.fr}
#' @details
#'
#' The data format required is a simple table with samples as rows and species as columns, with date (dd/mm/yyyy) and station in character. The dates are converted into date format in the algorithm.
#'
#' The Bloom Detecting Algorithm detects the bloom of a species within a time series of abundance according to three conditions. But first, the algorithm locate the high and low points of the curves. For each high point, the closest
#' was considered to be a bloom if:
#'
#'    1. The high points were above the value of \code{threshold} parameter value, which by default is 4, corresponding to the log10 of 10,000 cells/L. The threshold of 10,000 cells/L was used here as the algorithm was created to fit phytoplankton species. See the reference for more information.
#'
#'    2. The low points before and after the high points were inferior to the \code{perc_of_peak} of 0.85 (85\%) of the high point value. In this case, some humps can be merged, as blooms can sometimes be bimodal.
#'
#'    3. The merging of two humps would occur when the value of one of the lowest points did not fit the second condition. The merging of two humps cannot occur if the merging causes the increasing or decreasing phase of the bloom to be greater than \code{nbr_days}, by default 300 days.
#'
#' These three conditions were necessary as they enabled the extraction of the phenological bloom, which in our case corresponded to HABs. The HAB case study helped us define the hump minimum abundance \code{threshold} as well as the amplitude \code{nbr_days} and shape \code{perc_of_peak}.
#'
#' \code{Log} parameter simply transforms the abundance in log(x+1) as it helps with the large variation in the data value. The \code{min_Log} parameter (by default, 2 corresponding to 100 cells/L) was the minimum we fixed for our study.
#'
#' The output graphs are time series of your species abundance (grey dots) for each station with the fitted smooth spline (grey line) and confidence interval (grey shaded area).
#' The colored dots correspond to different timing phenological variables, which need to be in the following order for each bloom detected to confirm the validity of the algorithm fit.
#' Yellow (DBS) > Turquoise (DMF) > Red (DMA) > Purple (DMM) > Pink (DBE)
#'
#' @references Karasiewicz S., and Lefebvre A. (2022). Environmental Impact on Harmful Species Pseudo-nitzschia spp. and Phaeocystis globosa Phenology and Niche. \emph{JMSE} 10(2), 174. \doi{10.3390/jmse10020174}.
#' @examples
#' \donttest{
#' library(BDAlgo)
#' data(Abundance)
#' algo <- BDAlgo(Abundance, threshold=c(4,4), perc_of_peak=0.85, nbr_days=300, date_col= 2,
#' station_col=1, Sp = c("Pseunitz", "Phaeocy"), min_Log = 2,SP_label = c("Pseudo", "Phae"),
#' Log=TRUE,PDF=FALSE, saving_path= NULL)
#'  }
#'
#' @export BDAlgo
#' @rdname BDAlgo
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics axis lines points polygon
#' @importFrom stats quantile predict smooth.spline
#'

BDAlgo <-
  function(data,
           threshold = c(4, 4),
           perc_of_peak = 0.85,
           nbr_days = 300,
           date_col = 2,
           station_col = 1,
           Sp = c("Sp1", "Sp2"),
           min_Log = 2,
           SP_label = c("Species1", "Species2"),
           Log = TRUE,
           PDF = FALSE,
           saving_path = NULL) {
    sp.resampler <- function(x) {
      n <- nrow(x)
      resample.rows <- sample(1:n, size = n, replace = TRUE)
      return(x[resample.rows, ])
    }

    sp.spline.estimator <- function(x, m) {
      # Fit spline to data, with cross-validation to pick lambda
      fit <- smooth.spline(x = x[, 1], y = x[, 2], cv = TRUE)
      # Do the prediction and return the predicted values
      return(predict(fit, x = 1:m)$y) # We only want the predicted values
    }

    sp.spline.cis <- function(x, B, alpha, m) {
      spline.main <- sp.spline.estimator(x, m = m)
      # Draw B boottrap samples, fit the spline to each
      spline.boots <-
        replicate(B, sp.spline.estimator(sp.resampler(x), m = m))
      # Result has m rows and B columns
      cis.lower <-
        2 * spline.main - apply(spline.boots, 1, quantile, probs = 1 - alpha /
                                  2)
      cis.upper <-
        2 * spline.main - apply(spline.boots, 1, quantile, probs = alpha / 2)
      return(list(
        main.curve = spline.main,
        lower.ci = cis.lower,
        upper.ci = cis.upper,
        x = 1:m
      ))
    }

    peaks <- function (y,  y.thresh = 0) {
      pks <- which(diff(sign(diff(y, na.pad = FALSE)),
                        na.pad = FALSE) < 0) + 2
      pks <- cbind(pks, y[pks], rep("max", length(pks)))
      colnames(pks) <- c("x", "y", "Variable")
      pks <- pks[which(pks[, 2] > y.thresh), ]
      return(pks)
    }

    nadir <- function (y, y.thresh = 4) {
      ndr <- which(diff(sign(diff(y, na.pad = FALSE)),
                        na.pad = FALSE) > 0) + 2
      ndr <- cbind(ndr, y[ndr], rep("min", length(ndr)))
      colnames(ndr) <- c("x", "y", "Variable")
      ndr <- ndr[which(ndr[, 2] < y.thresh), ]
      return(ndr)
    }

    oldwd <- getwd()
    on.exit(setwd(oldwd))
    #Saving the results
    BDAlgo_results <- list()
    #Convert the matrix in to dataframe
    datafr <- as.data.frame(data)
    #Create year,month and date vector from the date
    year <- substr(datafr[, date_col], 7, 10)
    #Start of the loop which each jth species selected
    for (j in 1:length(Sp)) {
      #Select the jth species
      value <-  Sp[j]
      #Extract the data for the jth species
      Data <- datafr[, c("station", "date", value)]
      #Extract the dates
      dati <- Data[, 2]
      #find NA
      test <- is.na(as.numeric(Data[, Sp[j]]))
      #Change NA into 0
      Data[test, Sp[j]] <- 0
      #Change the third collumn name
      colnames(Data)[3] <- "value"
      #Order the data in chronological order
      Data <- Data[order(as.Date(dati, format = "%d/%m/%Y")),]
      #Create a year collum
      Data[, "year"] <-
        format(as.Date(Data[, 2], format = "%d/%m/%Y"), format = "%Y")

      if (isTRUE(Log)) {
        #Log(x+1) tranforme the abundance
        Data[, "value"] <- log10(as.numeric(Data[, "value"]) + 1)
        #Set the minimum to 2 (100 cell/L)
        Data[which(Data[, "value"] == 0), "value"] <- min_Log
        #ylabel for graph
        ylabel <-
          paste("Log(x+1)", SP_label[j], "Abundance", sep = " ")
      } else {
        ylabel <- paste(SP_label[j], "Abundance", sep = " ")
      }

      #Extract unique station name
      stat <- unique(Data[, 1])
      #Sort them alphabetically
      stat <- sort(stat)
      #Extract the number of stations
      N <- length(stat)
      #Create an empty list for saving the result
      result <- list()
      #Create a loop for applying the following line for each ith stations
      for (i in 1:N) {
        #Get the station code  of the ith station
        lab_stat <- stat[i]
        #Get the data of the ith statio,
        data <-
          Data[which(Data[, "station"] == lab_stat), c("date", "value", "year")]
        #Create a vector with all the dates
        all_date <- as.character(format(seq.Date(
          from = as.Date(paste("01/01/", min(data[, "year"]), sep = ""), format =
                           "%d/%m/%Y"),
          to = as.Date(paste("31/12/", max(data[, "year"]), sep =
                               ""), format = "%d/%m/%Y"),
          by = "day"
        ), format = "%d/%m/%Y"))
        #Create a vector of the same length than the dates
        all_date_num <- as.integer(1:length(all_date))
        #find the chronological number of the dates
        dy <- all_date_num[all_date %in% data[, "date"]]
        #Make the date vector into Date format
        dati <-
          as.Date(all_date[all_date %in% data[, "date"]], format = "%d/%m/%Y")
        #Calculate the number of days of the period studied
        z = as.numeric(difftime(
          as.Date(paste("31/12/", max(data[, "year"]), sep = ""), format = "%d/%m/%Y")
          ,
          as.Date(paste("01/01/", min(data[, "year"]), sep = ""),
                  format = "%d/%m/%Y")
        ))
        #the number of days of period +1
        Z <- as.numeric(z) + 1
        #Create a new data-frame used for spline
        dat <-
          as.data.frame(cbind(dy, as.numeric(data[, "value"])))
        #Fitting a smooth-spline with cross-validation
        smp <- smooth.spline(x = dat[, 1], y = dat[, 2], cv = TRUE)
        #Calculating the confidence intervalle of the smooth-spline
        sp.cis <- sp.spline.cis(
          x = dat,
          B = 1000,
          alpha = 0.05,
          m = Z
        )
        #Setting all negative values to 0
        sp.cis$lower.ci[which(sp.cis$lower.ci < 0)] <- 0
        #locating the low point in the curve
        nady <- nadir(sp.cis$main.curve, y.thresh = 8)
        #locating the high point of the curves
        peaky <- peaks(sp.cis$main.curve, y.thresh = 0)
        #concatenating the two low and high points
        test <- rbind(nady, peaky)
        #ordering in chronological order the point of high and low
        test <- test[order(as.numeric(test[, 1])),]
        #Creating empty vector for savings the results of the smoothing spline
        res <- c()
        #Creating empty vector for savings the results of date of bloom start
        DBS <- c()
        #Creating empty vector for savings the results of abundance of bloom start
        XO <- c()
        #Creating empty vector for savings the results of date of bloom end
        DBE <- c()
        #Creating empty vector for savings the results of abundance of bloom end
        XE <- c()
        #Creating empty vector for savings the results of date of maximum mortatlity
        DMM <- c()
        #Creating empty vector for savings the results of maximum mortatlity
        MM <- c()
        #Creating empty vector for savings the results of abundance of maximum mortatlity
        XM <- c()
        #Creating empty vector for savings the results of date of maximum fitness
        DMF <- c()
        #Creating empty vector for savings the maximum fitness
        MF <- c()
        #Creating empty vector for savings the results of abundance of maximum fitness
        XF <- c()

        # if the first line of the test data is max, it might be possible that the data
        #prior the max can be missed and therefore recover so start with a low point.
        # if no recovery is possible than remove the max
        if (test[1, 3] == "max") {
          startidate <- smp$x[which(smp$x < as.numeric(test[1, 1]))]
          startival <-
            smp$y[which(smp$x < as.numeric(test[1, 1]))]
          if (isTRUE(as.numeric(test[1, 2]) > startival[length(startival)])) {
            test <- rbind(c(min(startidate), min(startival), "min"), test)
          }
          else{
            test <- test[-1,]
          }
        }
        #delete the last row if it end by a maximum to only keep complete hump
        if (test[nrow(test), 3] == "max") {
          test <- test[-nrow(test),]
        }
        #Extract the number of points in the data -1
        M <- dim(test)[1] - 1
        #Start of loop which will go through each line of the data
        for (k in 1:M) {
          #Extract subdate which are the kth and k+1 lines
          rowi <- test[c(k, k + 1),]
          #if its goes from min to max
          if (rowi[1, 3] == "min" & rowi[2, 3] == "max") {
            #Extract the spline result from min to max
            infl <-
              sp.cis$main.curve[as.numeric(rowi[1, 1]):as.numeric(rowi[2, 1])]
            #Extract the date results from min to max
            dinfl <-
              as.numeric(rowi[1, 1]):as.numeric(rowi[2, 1])
            #Find the maximum different between two points
            MF <- c(MF, max(abs(diff(10 ^ infl))))
            #Extract the date of MF
            dmf <-
              dinfl[which(abs(diff(infl)) == max(abs(diff(infl)))) + 1]
            #Save the DMF
            DMF <- c(DMF, round(dmf))
            #Extract the xf
            xf <-
              infl[which(abs(diff(infl)) == max(abs(diff(infl)))) + 1]
            #Save the xf
            XF <- c(XF, xf)
            #Find the DBS with ese
            dbs <-
              ese((as.numeric(rowi[1, 1]) + 1):as.numeric(rowi[2, 1]),
                  diff(sp.cis$main.curve[as.numeric(rowi[1, 1]):as.numeric(rowi[2, 1])]),
                  0)[3]
            #if ese does not work
            if (is.nan(dbs)) {
              #find dbs with ede
              dbs <-
                ede((as.numeric(rowi[1, 1]) + 1):as.numeric(rowi[2, 1]),
                    diff(sp.cis$main.curve[as.numeric(rowi[1, 1]):as.numeric(rowi[2, 1])]),
                    0)[3]
            }
            #Get the abundance at dbs
            xo <- as.numeric(sp.cis$main.curve[dbs])
            #Save the DBS
            DBS <- c(DBS, round(dbs))
            #SAve the Xo
            XO <- c(XO, xo)
          }
          #if its goes from max to max
          if (rowi[1, 3] == "max" & rowi[2, 3] == "min") {
            #Extract the spline result from max to min
            infl <-
              sp.cis$main.curve[as.numeric(rowi[1, 1]):as.numeric(rowi[2, 1])]
            #Extract the date results from max to min
            dinfl <-
              as.numeric(rowi[1, 1]):as.numeric(rowi[2, 1])
            #Find the maximum different between two points
            MM <- c(MM, max(abs(diff(10 ^ infl))))
            #Extract the date of MM
            dmm <-
              dinfl[which(abs(diff(infl)) == max(abs(diff(infl)))) + 1]
            #Save the dmm
            DMM <- c(DMM, round(dmm))
            #Extract the xf
            xm <-
              infl[which(abs(diff(infl)) == max(abs(diff(infl)))) + 1]
            #Save xm
            XM <- c(XM, xm)
            #Find the DBE with ese
            dbe <-
              ese((as.numeric(rowi[1, 1]) + 1):as.numeric(rowi[2, 1]),
                  diff(sp.cis$main.curve[as.numeric(rowi[1, 1]):as.numeric(rowi[2, 1])]),
                  0)[3]
            #if ese does not work
            if (is.nan(dbe)) {
              #Reduce part of the curve to be checked to after DMM to DBE
              ind <-
                check_curve((as.numeric(dmm) + 1):as.numeric(rowi[2, 1]),
                            diff(sp.cis$main.curve[as.numeric(dmm):as.numeric(rowi[2, 1])]))$index
              #find dbe with ede
              dbe <-
                ede((as.numeric(dmm) + 1):as.numeric(rowi[2, 1]),
                    diff(sp.cis$main.curve[as.numeric(dmm):as.numeric(rowi[2, 1])]),
                    ind)[3]
            }
            #it still does not work
            if (is.nan(dbe)) {
              #find the dbe to after DMM and before DBE
              dbe <-
                ede((as.numeric(dmm) + 1):as.numeric(rowi[2, 1]),
                    diff(sp.cis$main.curve[as.numeric(dmm):as.numeric(rowi[2, 1])]),
                    1 - ind)[3]
            }
            #Extract the xe
            xe <- as.numeric(sp.cis$main.curve[dbe])
            #Save dbe
            DBE <- c(DBE, round(dbe))
            #Save xe
            XE <- c(XE, xe)
          }
        }
        #Calculate difference in abundance between max and min
        testi <-
          cbind(test, c(0, diff(10 ^ as.numeric(test[, 2]))))
        #Detect if there is the 85% difference between Max, min max are not respected
        detec <-
          perc_of_peak * (10 ^ as.numeric(testi[which(testi[, 3] == "max"), 2])) >
          abs(as.numeric(testi[which(testi[, 3] == "max") + 1, 4])) |
          abs(as.numeric(testi[which(testi[, 3] == "max"), 4])) <
          perc_of_peak * (10 ^ as.numeric(testi[which(testi[, 3] ==
                                                        "max"), 2]))
        #The number of detected
        L <- length(detec)
        #Allocate the maximum, its date with the min b4 and after
        maxi <- as.numeric(testi[which(testi[, 3] == "max"), 2])
        datemaxi <-
          as.numeric(testi[which(testi[, 3] == "max"), 1])
        minib4 <-
          as.numeric(testi[which(testi[, 3] == "min"), 2])[-(L + 1)]
        miniaft <-
          as.numeric(testi[which(testi[, 3] == "min"), 2])[-1]

        #Fix skip to 0
        skip <- 0
        #Create an empty vector for saving
        res <- c()
        underthres <- c()
        #Create a loop for each of the potential detected bloom
        for (l in 1:L) {
          #if the max under threshold it is saved for portential checking and move to next bloom
          if (maxi[l] < threshold[j]) {
            underthres <- c(underthres, l)
            next
          }
          #if the difference is higher than 85% than save the bloom
          if (isFALSE(detec[l])) {
            res <-
              rbind(res,
                    c(DBS[l], XO[l], DMF[l], XF[l], MF[l], DMM[l], XM[l], MM[l], DBE[l], XE[l]))
            next
          }
          if (skip %in% l) {
            next
          }
          if (isTRUE(detec[l + 1])) {
            if (isTRUE(detec[l + 2])) {
              datemaximi <- datemaxi[c(l, l + 1, l + 2)]
              keep <- diff(datemaximi)
              if (keep[1] > keep[2]) {
                res <-
                  rbind(res,
                        c(DBS[l], XO[l], DMF[l], XF[l], MF[l], DMM[l], XM[l], MM[l], DBE[l], XE[l]))
                next
              }
              if (keep[1] < keep[2]) {
                skip <- c(l + 1, l + 2)
                res <-
                  rbind(res,
                        c(DBS[l], XO[l], DMF[l], XF[l], MF[l], DMM[l + 1], XM[l + 1], MM[l + 1], DBE[l +
                                                                                                       1], XE[l + 1]))
                res <-
                  rbind(res,
                        c(DBS[l + 2], XO[l + 2], DMF[l + 2], XF[l], MF[l + 2], DMM[l + 2], XM[l +
                                                                                                2], MM[l + 2], DBE[l + 2], XE[l + 2]))
                next
              }
            }
            if (isFALSE(detec[l + 2])) {
              skip <- l + 1
              res <-
                rbind(res,
                      c(DBS[l], XO[l], DMF[l], XF[l], MF[l], DMM[l + 1], XM[l + 1], MM[l + 1], DBE[l +
                                                                                                     1], XE[l + 1]))
              next
            }
            else{
              if (underthres[length(underthres)] == (l - 1)) {
                res <-
                  rbind(res,
                        c(DBS[l - 1], XO[l - 1], DMF[l - 1], XF[l], MF[l - 1], DMM[l], XM[l], MM[l], DBE[l], XE[l]))
              }
              if (datemaxi[l] - datemaxi[l - 1] > nbr_days) {
                res <-
                  rbind(res,
                        c(DBS[l], XO[l], DMF[l], XF[l], MF[l], DMM[l], XM[l], MM[l], DBE[l], XE[l]))
              }
              else{
                last <- dim(res)[1]
                res[last, 5:8] <-
                  c(DMM[l], XM[l], MM[l], DBE[l], XE[l])
              }
            }
          }
        }


        colnames(res) <-
          c("DBS",
            "XO",
            "DMF",
            "XF",
            "MF",
            "DMM",
            "XM",
            "MM",
            "DBE",
            "XE")
        #Create an empty vector for each variables
        DMA <- c()
        MA <- c()
        IL <- c()
        ONS <- c()
        CLI <- c()
        DEC <- c()
        DL <- c()
        HA <- c()
        END <- c()
        BL <- c()
        SI <- c()
        SD <- c()
        #Get the number of lines of the table
        W <- dim(res)[1]
        #Create a
        for (w in 1:W) {
          subi <- smp$x[which(res[w, "DBS"] < smp$x)]
          subi <- subi[which(res[w, "DBE"] > subi)]
          suby <- smp$yin[smp$x %in% subi]
          ma <- max(suby)
          dma <- subi[which(suby == max(suby))]
          if (length(dma) > 1) {
            dma <- dma[1]
          }
          #Calculate all the variables
          il <- as.numeric(dma - res[w, "DBS"])
          ons <- as.numeric(res[w, "DMF"] - res[w, "DBS"])
          cli <- as.numeric(dma - res[w, "DMF"])
          ha <- as.numeric(res[w, "DMM"] - dma)
          dl <- as.numeric(res[w, "DBE"] - dma)
          ha <- as.numeric(res[w, "DMM"] - res[w, "DMF"])
          end <- as.numeric(res[w, "DBE"] - res[w, "DMM"])
          bl <- as.numeric(res[w, "DBE"] - res[w, "DBS"])
          si <-
            log(as.numeric(10 ^ ma) / as.numeric(10 ^ res[w, "XO"])) / as.numeric(il)
          sdi <-
            log(as.numeric(10 ^ res[w, "XE"]) / as.numeric(10 ^ ma)) / as.numeric(dl)
          # Save the variables values
          DMA <- c(DMA, dma)
          MA <- c(MA, ma)
          IL <- c(IL, il)
          ONS <- c(ONS, ons)
          CLI <- c(CLI, cli)
          DEC <- c(DEC, ha)
          DL <- c(DL, dl)
          HA <- c(HA, ha)
          END <- c(END, end)
          BL <- c(BL, bl)
          SI <- c(SI, si)
          SD <- c(SD, sdi)
        }
        #Bind all the results together by collumn
        res <-
          cbind(res[, 1:5], DMA, MA, res[, 6:10], IL, ONS, CLI, DEC, DL, HA, END, BL, SI, SD)
        #>Remove the detected which did not some phases are missing
        res <- res[which(res[, "ONS"] > 0),]
        res <- res[which(res[, "CLI"] > 0),]
        res <- res[which(res[, "DEC"] > 0),]
        res <- res[which(res[, "END"] > 0),]
        #Create a pdf for plotting the results
        if (isTRUE(PDF)) {
          if(!is.null(saving_path)){
            setwd(saving_path)
            }
          pdf(file = paste(
            Sp[j],
            "_",
            paste("smooth-spline", lab_stat, sep = "_"),
            ".pdf",
            sep = ""
          ))
        }
        #Create an empty plotting window
        plot(
          x = all_date_num,
          y = sp.cis$main.curve,
          type = "n",
          xlab = ("Date"),
          ylab = ylabel,
          ylim = c(0, max(sp.cis$upper.ci)),
          main = lab_stat,
          xaxt = "n"
        )
        #Create a vector with the first day of each year
        yeari <-
          paste("01/01/", min(data[, "year"]):max(data[, "year"]), sep = "")
        #Create a vector with the day of the year
        midyeari <-
          paste("30/06/", min(data[, "year"]):max(data[, "year"]), sep = "")
        #Add the x-axis to the graph
        axis(
          side = 1,
          at = c(all_date_num[all_date %in% yeari], max(all_date_num) + 1),
          labels = F
        )
        #Add the tick marks labels of the x-axis
        axis(
          side = 1,
          at = all_date_num[all_date %in% midyeari],
          tick = F,
          labels = min(data[, "year"]):max(data[, "year"]),
          las = 2
        )
        #Add the shaded polygong for the confidence interval
        polygon(
          c(all_date_num, rev(all_date_num)),
          c(sp.cis$lower.ci, rev(sp.cis$main.curve)),
          col = "grey",
          border = NA
        )
        #Add the shaded polygon for the confidence interval
        polygon(
          c(all_date_num, rev(all_date_num)),
          c(sp.cis$main.curve, rev(sp.cis$upper.ci)),
          col = "grey",
          border = NA
        )
        #Add the original data
        points(smp$x , smp$yin, pch = 21, bg = "#737373")
        #Add the spline
        lines(all_date_num , sp.cis$main.curve, lwd = 1.5)
        #add the DMF on the spline
        points(as.numeric(res[, "DMF"]),
               as.numeric(res[, "XF"]),
               pch = 21,
               bg = "turquoise")
        #add the DMM on the spline
        points(as.numeric(res[, "DMM"]),
               as.numeric(res[, "XM"]),
               pch = 21,
               bg = "purple")
        #add the DBS on the spline
        points(as.numeric(res[, "DBS"]),
               as.numeric(res[, "XO"]),
               pch = 21,
               bg = "yellow")
        #add the DBE on the spline
        points(as.numeric(res[, "DBE"]),
               as.numeric(res[, "XE"]),
               pch = 21,
               bg = "pink")
        #add the DMA on the spline
        points(as.numeric(res[, "DMA"]),
               as.numeric(res[, "MA"]),
               pch = 21,
               bg = "red")
        #Close the pdf
        if (isTRUE(PDF)) {
          dev.off()
        }
        #Convert the dates into Julian dates
        julien_date <-
          format(as.Date(all_date, format = "%d/%m/%Y"), "%j")
        #Create a miror image of the results for saving
        all_bloom <- res
        #COnvert the Date variables into julian dates
        all_bloom[, "DBS"] <-
          as.numeric(julien_date[all_bloom[, "DBS"]])
        all_bloom[, "DBE"] <-
          as.numeric(julien_date[all_bloom[, "DBE"]])
        all_bloom[, "DMM"] <-
          as.numeric(julien_date[all_bloom[, "DMM"]])
        all_bloom[, "DMF"] <-
          as.numeric(julien_date[all_bloom[, "DMF"]])
        all_bloom[, "DMA"] <-
          as.numeric(julien_date[all_bloom[, "DMA"]])
        #Create the results list
        result[[i]] <-
          list(
            smooth_spline = smp,
            conf_intervall = sp.cis,
            all_date = all_date,
            all_bloom = all_bloom,
            all_bloom_date = res
          )

      }
      #names the resutls
      names(result) <- stat
      # Save the results
      BDAlgo_results[[Sp[j]]] <- result
    }
    return(BDAlgo_results)
  }
