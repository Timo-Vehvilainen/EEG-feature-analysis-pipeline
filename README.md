# EEG-feature-analysis-pipeline

This pipeline can be used to analyze changes in computational features of EEG data during drug trials. The computational features provided include:

    - Power Spectral Density (PSD)
    - cross Power Spectral Density (cPSD)
    - Activation Synchrony Index (ASI)
    - Nestedness Coefficient (NC)
    - weighted Phase Lag Index (wPLI)
    - amplitude-integrated EEG (aEEG):
        - mean
        - interquartile range 
    - range-EEG (rEEG):
        - mean
        - interquartile range
        - lower 5th percentile
    - Multifractal Detrended Fluctuation Analysis (MFDFA):
        - height
        - width
        - peak
        - tail
    - Suppression Curve (SC)
    
the pipeline computes and displays the following statistical representations for the aforementioned features:
    
    - absolute time trend
    - time trend relative to a specified baseline
    - group comparison relative to a specified baseline
    - correlation testing relative to a specified baseline
    
# Input
    
The pipeline takes in EEG data in the form of .e-files (Nicolet) and uses the package available in https://github.com/ieeg-portal/Nicolet-Reader to read them in. It also can take in information about artifacts in the data, in the format of .mat-files (more on their formatting below).

The .e-files need to include an annotation somewhere indicating the starting time of the drug administration, as this annotation is used for pinpointing the moment in the data. This annotation needs be searchably similar in all the files, e.g. by including the name of the drug in the annotation ("fentanyl"). Unmodified, this pipeline is only able to handle 4-channel EEG data designed for neonates. It was specifically designed for use with the F3, F4, P3 and P4 channels used in neonate studies. By default it will calculate the Left (F3-P3), Right (F4-P4), Frontal (F3-F4) and Parietal (P3-P4) bipolar channel montages, and perform feature calculations on them aswell.

# Usage

To use this pipeline, follow the instructions provided in the comments of the feature_analysis_pipeline.m -file, and run it through section by section. Commentation is provided in the subfunctions, so that moderate handling of matlab should suffice for modifying the pipeline to custom needs. This pipeline was built and first used with MatLab version 2019a.

# Output

After running the final part of the pipeline succesfully, the following three files are saved for each figure window that is produced:
    
    -.fig-file of the full figure
    -.pdf-image of the full figure
    -.mat file, which includes a cell array of tables. Each cell in the cell array visually corresponds to a subplot in the figure of the same name. Eg. The data table for the subplot in the upper left-hand corner is in cell with index (1,1).

# Artifact data formatting

The automatic reading of artifact data is only implemented for monopolar channels F3, F4, P3 and P4. Information about artifacts in the data should be provided in .mat-files, as cell arrays. The cell arrays need to follow the following format:

Each row in the cell array will represent a single artifact period in the data, and it will be described in the first 3 columns in the following manner:

1st column:

    Text annotation starting with the string "Artifact". 
    If nothing else is provided, the artifact is assumed to affect all channels. 
    Individual channels may be specified by continuing the string with a dash surrounded by spaces " - ",  
    and following this by naming the individual channels affected, separated by commas. 
    An example annotation: "Artifact - F3,P3".
    
2nd column:

    The starting time of the artifact, in "HH:mm:ss" format, for example "00:26:38"
    
3rd column:

    The duration of the artifact, in "mm:ss" format (the amount of minutes may exceed 60 if needed). For example "01:05" 
