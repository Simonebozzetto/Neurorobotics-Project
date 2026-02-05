# **ASSIGNMENT 1**

###### **Simone Bozzetto, Davide Gasparini, Samuele Pinello** 



## **MAIN FILES:**



1. preparation\_grand\_average.m

&nbsp;	preparation and preprocessing of all the data for each subject for the grand average analysis

2\. grand\_average\_analysis.m

&nbsp;	temporal visualization and topoplot for each subject and for the average of subjects

3\. preparation\_all\_data.m

 	preparation and preprocessing of all the data for each subject for the BMI decoding

4\. cal\_feature\_identification\_extraction.m

&nbsp;	Calibration phase: only the offline runs already processed. Computation of the features, selection of the most discriminant features and creatation of a classifier based on those 	features + metrics for calibration results

5\. evaluation\_phase.m

&nbsp;	Evaluation phase: only the online runs already processed. Computation of the features and extraction of those already selected during the calibration phase. Evaluation of the 	classifier created during the calibration phase. Implementation and application of 2 evidence accumulation frameworks on the posterior probabilities + metrics for evaluation 	results

6\. BMI\_metrics.m

&nbsp;	plots and results for each subject and average result of the BMI decoding part	





## **FUNCTIONS (devided based on the main files):**



preparation\_grand\_average.m

* label\_vectors\_1.m (extraction of labels required, in samples)
* preprocessing\_1.m (preprocessing required for calculation of ERD)
* conc\_to\_trial.m (concatenation over trials of the EEG signal and cut based on the length of the trials)



grand\_average\_analysis.m

* concat\_avg\_z\_ERD.m (concatenation and averaging over trials of the ERD of each subject)
* avg\_subj\_ERD.m (concatenation and averaging over subjects)



preparation\_all\_data.m

* preprocessing.m (preprocessing required for calculation of PSD)
* PSDcalculation.m (PSD calculation)
* proc\_spectrogram.m (used for the PSD calculation)
* proc\_pos2win.m (used for the PSD calculation)



cal\_feature\_identification\_extraction.m \& evaluation\_phase.m

* label\_vectors.m (extraction of labels required, in windows)
* fisher.m (Fisher score computation)



## **HOW TO RUN ALL THE SCRIPTS**

* all the GDF files of all the subjects have to be in the same folder, the same folder in which there are all the mains and functions
* to run the scripts in order, just follow the order of the main files above



* in particular here below there are all the features that has to be selected during the cal\_feature\_identification\_extraction.m for each subject:

&nbsp;		Subject 1: \[18 9; 20 9; 22 9]

&nbsp;		Subject 2: \[14 11; 14 7]

&nbsp;		Subject 3: \[14 11; 14 7]

&nbsp;		Subject 4: \[12 7; 12 11; 10 7; 10 11]

&nbsp;		Subject 5: \[14 7; 12 7; 12 6]

&nbsp;		Subject 6: \[12 7; 14 7; 12 8; 14 8]

&nbsp;		Subject 7: \[12 11; 10 11]

&nbsp;		Subject 8: \[12 8; 12 12; 12 16; 14 16]

&nbsp;  in particular the feature are pairs of \[frequency,channels] and has to be selected in this matrix way of the dimension (N\_features x 2)

&nbsp;  just copy and paste the \[...] when required



## **CONTRIBUTES:**



* Simone Bozzetto: 

&nbsp;	Matlab code of preprocessing, BMI decoding

&nbsp;	Methods and results discussion

&nbsp;	Report writing

* Davide Gasparini: 

&nbsp;	Matlab code of preprocessing and grand average analysis

&nbsp;	Methods and results discussion

&nbsp;	Report writing

* Samuele Pinello: 

&nbsp;	Methods discussion

&nbsp;	Report writing





