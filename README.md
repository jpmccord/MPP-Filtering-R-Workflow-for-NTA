# MPP-Filtering-R-Workflow-for-NTA

R Processing script for filtering and flagging output from Agilent Mass Profiler Professional (MPP) based on experimental validations parameters.

DataProcessingWrapper loads config file and necessary functions from other scripts to process.

config file contains variables detailing MPP input files, processing variables, and output settings. 

SampleInput_XXX.R scripts are intended to reformat known output formats into the format expected by the processing workflow. Expecdted columns are:
FeatureNumber
ExpID
IonMode
ModeDetail
SampleDetail
Rep
volume
Mass
RetentionTime
Compound
Score
Frequency

