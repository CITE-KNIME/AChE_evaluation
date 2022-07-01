# AChE_evaluation
NOTE: KNIME workflows used to evaluate the raw data from the AChE assay and the DRC plots for the AChE inhibition; versioned based on the publication in which respective workflows were first used

## KNIME workflows for evaluation of data from acetylcholinesterase (AChE) assay
A KNIME installation is required in order to run the workflows. The workflow can be imported into KNIME by importing each KNIME archive file (.knwf). The workflow has been established and tested with KNIME version 4.2.4. and version 4.4.1. Our workflows have been established and tested with R version 4.0.4 but other versions may work as well. You may need to install various KNIME extentions in order to run the workflows. The following extension is required: KNIME Interactive R Statistics Integration (org.knime.features.r.feature.group).

### Workflow 1: Calculation of % inhibition of AChE for each well treated with a sample/chemical

***Input files required to run the KNIME workflow (Input format is xlsx/xls)***

File 1: Layout of experiment containing the metadata for the experiment
-	One sheet per plate needs to be provided
-	Template YYYY-MM-DD_AChE_layout has to be used, the data given corresponds to the date of the dosing

File 2: Raw data from the reader
-	Absorbance was read in 1 min intervals over a duration of 30 min
-	One sheet per plate needs to be provided
-	Here data is read with a Tecan M1000 Multimode reader
-	If another plate reader is used, then the KNIME workflow may need to be adapted to the format of the raw data given by the reader

***Linear fit of data and calculation of % inhibition of AChE***
-	For each well per plate, a linear fit is performed based on the data from the time series using R.
-	The R code is embedded in a node  ![grafik](https://user-images.githubusercontent.com/108472923/176879980-9d63c93d-a4da-4804-aaa3-75fdf2852b7a.png)
-	R code for linear fit
                            library(plyr)

                            time_in_min <- knime.in$"time_of_measurement_min"
                            OD410 <- knime.in$"absorbance"

                            modlinear <- lm(OD410 ~ time_in_min)
                            df <- summary(modlinear)$coef[2,1:2]
                            dfrenamed <- rename(df, c("Estimate" = "slope", "Std. Error" = "slope_se"))
                            dfrenamed$sample <- knime.in$"sample"
                            dfrenamed$sample_category <- knime.in$"sample_category"
                            dfrenamed$cell_line <- knime.in$"Cell line"
                            dfrenamed$concentration <- knime.in$"concentration"
                            dfrenamed$well <- knime.in$"Well Position"
                            dfrenamed$plate_number <- knime.in$"plate_number"
                            dfrenamed$Exposure_in_h <- knime.in$"Exposure duration"
                            dfrenamed$date_time_measurement <- knime.in$"Date_Time"

                            knime.out <- as.data.frame(dfrenamed)

-	Slope values were obtained by the linear fit to determine enzyme velocity. Subsequently, the slope values were used to calculate the % Inhibition of AChE as follows: <br>
<p align="center">
<img src ="https://user-images.githubusercontent.com/108472923/176880272-79497435-2744-494a-bc82-5a485a3a9db1.png" width=40% align="center">
</p>

***Output of the workflow***<br>
The result are two Excel files which are named in a standardized manner and saved automatically in the same folder as the raw data.

`Output 1 = YYYY-MM-DD_AChE_metadata_rawdata.xlsx (Note: the date is the same date as for the raw data)`<br>

<dl>
  <dt>Columns</dt>
  <dd>Sample – Sample name<br>
  Sample category – SC = …, TC = …, NC w cells = negative control with cells, NC wo cells = negative control without cells <br>
  column, row and Well position – based on standard nomenclature <br>
  concentration – concentration dosed in bioassay in M for chemicals and REF for environmental samples <br>
  Exposure duration – time of incubation of cells with sample <br>
  Plate number <br>
  Cell line – name of cell line <br>
  time_of_measurement_min – time of time series measurement of absorbance <br>
  absorbance – result of absorbance measurement <br>
  date & time – of absorbance measurement </dd>
  </dd>
</dl>

`Output 2 = YYYY-MM-DD_AChE_%_inhibition`<br>

<dl>
  <dt>Columns</dt>
  <dd>Sample – Sample name <br>
	Sample category – SC = …, TC = …, NC w cells = negative control with cells, NC wo cells = negative control without cells <br>
  Cell line – name of cell line <br>
  concentration – concentration dosed in bioassay in M for chemicals and REF for environmental samples <br>
	well  – based on standard nomenclature <br>
  Plate number <br>
	Exposure_in_h – time of incubation of cells with sample <br>
	Date_time_measurement - of absorbance measurement <br>
	Inhibition_% - of AChE
  </dd>
</dl>


### Workflow 2: Fitting dose response curves of replicates based on results for % inhibition of AChE

***Input files required to run the KNIME workflow (Input format is xlsx/xls)***<br>

Output files per replicate from workflow 1 (`YYYY-MM-DD_AChE_%_inhibition.xlsx`) need to be saved into one folder. This folder needs to be selected in the "List files"-Node.

***Log-logistic fit of data of % inhibition of AChE***
-	Log-logistic fit of data of % inhibition of AChE is performed using R with the bottom and top fixed to 0 and 100, respectively. In addition, the possible lowest value for the EC50 was fixed to 0 (see lowerl=c(NA,0)).
-	The R code is embedded in a node ![grafik](https://user-images.githubusercontent.com/108472923/176885436-2f0baf8b-818e-4fb8-8ce3-a3f8ceeb03e9.png)

        library("drc")

        LL_4P_parameters=NULL

        #Fitting 4-parameter log-logistic model
        LL_4P_effect<- drm(knime.in$"Inhibition_%"~knime.in$"concentration" , fct = LL.4(names = c("slope", "lower_limit", 
        "upper_limit", "EC50"), fixed = c(NA, 0, 100, NA)), lowerl=c(NA,0))

        #calculate Akaike information criterion (AIC), which is an estimator of the relative quality of statistical models 
        #for a given set of data
        AIC_effect <- AIC(LL_4P_effect)

        #extract model coefficients (slope, top_effect ...)
        LL_4P_coef_effect <-t(as.data.frame(coef(LL_4P_effect)))

        #calculation of ECX std errors
        EC_results <- ED(LL_4P_effect, c(10,20,50))
        EC_results <- as.matrix(unlist(EC_results))

        #calculation and extraction of EC50 std error
        #EC_results <- ED(LL_4P_effect, c(50))
        #EC50_error <- as.matrix(EC_results["e:1:50",2])
        #colnames(EC50_error)=c("EC50_error")

        #Returns the EC10 & EC20
        EC10 <- as.matrix(EC_results["e:1:10",1])
        colnames(EC10)= c("EC10")
        EC20 <- as.matrix(EC_results["e:1:20",1])
        colnames(EC20)= c("EC20")

        #Returns the ECX error
        EC10_error <- as.matrix(EC_results["e:1:10",2])
        colnames(EC10_error)=c("EC10_error")
        EC20_error <- as.matrix(EC_results["e:1:20",2])
        colnames(EC20_error)=c("EC20_error")
        EC50_error <- as.matrix(EC_results["e:1:50",2])
        colnames(EC50_error)=c("EC50_error")

        #summary of model required for extraction of stats
        LL_4P_sum_effect<-summary(LL_4P_effect)
        LL_4P_sum_effect<-as.matrix(unlist(LL_4P_sum_effect))

        #Returns the residual degrees-of-freedom
        degrees_of_freedom <- as.matrix(LL_4P_sum_effect["df.residual",1])
        colnames(degrees_of_freedom)=c("degrees_of_freedom")

        #Returns residual sum of squares
        rse_effect<-as.matrix(LL_4P_sum_effect["rseMat1",1])
        colnames(rse_effect)=c("Res_sum_squares")

        # Confidence interval for a single parameter (e=x50)
        conf<-confint(LL_4P_effect, "EC50", level = 0.9) 

        knime.out <- as.data.frame(cbind(LL_4P_coef_effect, EC50_error, EC10, EC10_error, EC20, EC20_error, conf, 
        AIC_effect, rse_effect, degrees_of_freedom))
 
-	In addition, the result was plotted using R as well. ![grafik](https://user-images.githubusercontent.com/108472923/176885789-cf34f937-e5f8-4645-bf51-ea541248d874.png)

library("plotrix")
library("drc")

        #Model: y = e0 + eMax (x^h / (ed50^h +x^h))

        # %EFFECT - probe and conc column selection account for potential missing values in any column
        LL_4P_effect<- drm(knime.in$"Inhibition_%"~knime.in$"concentration" , fct = LL.4(fixed = c(NA, 0, 100, NA)), 
        lowerl=c(NA,0))

        sample_name <-as.matrix(knime.in$"sample")
        sample_name <-sample_name[1,]

        ### Prepare for showing zero values on x-axis
        ## add extra space to right margin of plot within frame
        par(mar=c(4, 5.5, 4, 4) + 0.1)

        max<-max(knime.in$"concentration")
        plot(LL_4P_effect, type="all", main=paste(sample_name,"\n"),xlab = "",ylab="",ylim=c(-20,110), axes=FALSE, col ="red")
        # axes = FALSE, thus add labels to all axes with axis and mtext
        axis(2, las=1) # las=1 makes horizontal labels
        mtext ("Inhibition in %", side=2, line=4)
        axis(1, col = "black") 
        mtext ("c [M]", side=1, line=2.5)
        par(new=TRUE)
        plot(LL_4P_effect, type="confidence", xlab = "", ylab = "", ylim=c(-20,110), axes=FALSE, col = "red")
 
-	Modelling of the dose response curves may fail for negative controls or samples showing no activity. Thus, an error catch was introduced to deal with this case.

***Output of the workflow***<br>
The results are Excel files which are named in a standardized manner and saved automatically in the same folder as the input data.











