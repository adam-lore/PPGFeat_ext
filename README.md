# PPGFeat - Extended  
This repository contains an extended version of the [PPGFeat](https://github.com/saadsur/PPGFeat) application.  
The app takes raw unfiltered PPG data and tries the find the best cycle in each segment. From each cycle, the magnitude and tome are taken from 15 fiducial points which are used to calculate a set of 145 different features.  

![Image of PPGFeat_ext app in use](/media/PPGFeat_example.png)

## Additions from PPGFeat  
- Allows loading of PPG data from different file formats (.csv, .mat, .txt, .xlsx, WFDB)
- Accepts data collected with any sampling frequency and allows resampling of data to other frequencies
- Can segment larger PPG entries to smaller segments
- Option to invert loaded data
- Calculates signal quality index values for both whole segments and individual cycles using skewness and template matching
- Calculates an a set of 145 biomarkers in addition to the 15 fiducial points

## Usage Notes
1. Enter information about the data to be loaded which includes: sampling frequency (not needed when loading WFDB data), whether to resample the data and to what frequency if so, low and high-pass values for the filter, whether to segment the entries and how long each segment should be and finally if the data needs to be inverted.
2. Data can then be loaded either by selecting a file or by selecting a folder and then choosing what types of files should be loaded from the folder. (Both .hea and .dat files are needed for WFDB)
3. The toolbox will try to find find the highest quality cycle to use. The cycle used can be changed by manually editing Min1 and/or Min2.
4. Pressing the Plot button will cause the toolbox to calculate all the fiducial points and the features for the current cycle. The values can be manually updated by pressing the Update button which will save the fiducial values from the edit fields and recalculate the features.
5. The entries can either be processed one by one using the Next button or the toolbox can process every entry automatically using the Process every entry button.

[Video of PPGFeat - Extended in use](https://www.youtube.com/watch?v=HmEc5A7E4Vw)

## Dependencies
MATLAB toolboxes:
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

Needed to load WFDB data:
- [WFDB MATLAB Toolbox](https://physionet.org/content/wfdb-matlab/0.10.0/)

## Read more
[Expanding PPGFeat for Multi-Dataset PPG Analysis, Comprehensive Biomarker Extraction and Machine learning for Cardiovascular Risk Prediction](https://mdh.diva-portal.org/smash/record.jsf?pid=diva2%3A2040268&dswid=7506)  

## Acknowledgements  
This project builds on the PPGFeat project by saadsur (Saad Abdullah): https://github.com/saadsur/PPGFeat  

This project uses the Natural-Order Filename Sort function by Stephen23 (Stephen Cobeldick): https://se.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort

This project utilizes the Waveform Database Software Package (WFDB) for MATLAB and Octave created by Ikaro Silva, Benjamin Moody, George Moody: https://physionet.org/content/wfdb-matlab/0.10.0/


