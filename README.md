# PPGFeat - Extended  
This repository contains an extended version of the [PPGFeat](https://github.com/saadsur/PPGFeat) application.  
The app takes raw unfiltered PPG data and tries the find the best cycle in each segment. From each cycle, the magnitude and tome are taken from 15 fiducial points which are used to calculate a set of 145 different features.

## Additions from PPGFeat  
- Allows loading of PPG data from different file formats (.csv, .mat, .txt, .xlsx)
- Accepts data collected with any sampling frequency and allows resampling of data to other frequencies
- Can segment larger PPG entries to smaller segments
- Calculates signal quality index values for both whole segments and individual cycles using skewness and template matching
- Calculates an a set of 145 biomarkers in addition to the 15 fiducial points

## Read more
[Expanding PPGFeat for Multi-Dataset PPG Analysis, Comprehensive Biomarker Extraction and Machine learning for Cardiovascular Risk Prediction](https://mdh.diva-portal.org/smash/record.jsf?pid=diva2%3A2040268&dswid=7506)

## Acknowledgements  
This project builds on the PPGFeat project by saadsur (Saad Abdullah): https://github.com/saadsur/PPGFeat  

This project uses the Natural-Order Filename Sort function by Stephen23 (Stephen Cobeldick): https://se.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort
