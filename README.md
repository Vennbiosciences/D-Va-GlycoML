# D-Va: Glycopeptide Fragment Prediction
[![License: CC BY-ND 4.0](https://img.shields.io/badge/License-CC_BY--ND_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nd/4.0/)

D-Va project for glycopeptide fragment prediction comprises two parts: spectral library generation and deep learning.

## Spectral library
This pipeline constructs annotated glycopeptide spectral libraries from MGF files and gLabel 
results (from [pGlyco3](https://www.nature.com/articles/s41592-021-01306-0)). The spectral library format is in NIST MSP format. In order to analyze the MSP files, we need to parse the data
and convert to CSV file, which could be feed into the deep learning model.

Following [pGlyco3](https://www.nature.com/articles/s41592-021-01306-0) nomenclature, we use b- and y- ion annotations for the peptide backbone, as well as Y-ions (and not B-ions) for the intact glycopeptide. Specifically, in order to generate intact glycopeptide fragmentation pathway, we need do the deconvolution for different charges. After sorting Y-ions with the total number of monosaccharides, we use *de novo* sequencing to 
determine the structure's linear sequence.

## Deep learning model
Finally, we select the top ten candidates with the maximum summary of the intensities, then attach them to the peptide sequence.  After one-hot encoding the linearized sequence, we feed them into the deep learning model.

The deep learning model is a bidirectional LSTM(Long Short-Term Memory) with multi-head attention. 
This deep neural network can predict both peptide and glycan fragmentation patterns.


## Folders and files:

* deep learning: deep learning models for the training and prediction
  * biLSTM: basic biLSTM model and multi head attention model
  * models: construct the deep learning network based on biLSTM
  * run_model: the model could be used for the prediction after training
  * sample_predict: input the csv file as the testing
  * train_model: input the data to train the model
  * trainer: pipeline for handling the input data
* input_files: some example files for msp, and csv, etc.
* output_files: some example csv files for the prediction
* saved_models: save the models after training
* spectral_library: pipeline for constructing spectral library, parse data, de novo sequencing, etc.
  * calculate_fragment_mz: calculate the fragment mz for b-ions, y-ions, and Y-ions
  * construct_all_library: construct the spectral library which contains all the information: annotated one or not
  * construct_annotated_library: construct the spectral library which only contains annotated ions
  * de_novo_sequencing: determine the intact glycopeptide structure without using glycan database
  * parse_msp_data: read msp files, and use deno sequecning, then generate the csv file for the deep learning model
  * stack_queue: basic data structure for analyzing the glycans
  
## Contributors:

* Zhewei Liang (zliang@venn.bio)
* Mingqi Liu (mliu@venn.bio)
* Richard "DJ" Shipman (richard.shipman@venn.bio)
* Norton Kitagawa (norton.kitagawa@venn.bio) - corresponding author
