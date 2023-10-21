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

## Installation and development
This project uses Git Large File Storage (`git-lfs`) for versioning of the larger Pytorch `*.pth` files.  To gain access, you will need to install `git-lfs` on your system as follows:
```
sudo apt-get install git-lfs  # Ubuntu/Debian Linux
brew install git-lfs          # macOS with Homebrew
port install git-lfs          # macOS with MacPorts
```

Python dependencies are managed by `poetry`.  Install dependencies and enter the virtual environment as follows
```
python -m pip install --user poetry
python -m poetry install
python -m poetry shell
```

## Spectral library building scripts (in order of operation):
### `construct_annotated_library.py`

This tool reads Mascot Generic Format (.mgf) files, with companion gLabel (.txt) results files, then fuses them together into a NIST MSP format file containing annotated peak lists.

<details>
<summary>Usage:</summary>

```bash
python ./src/spectral_library/construct_annotated_library.py --IT {N-NlinkedGlycoPeptide or O-OlinkedGlycoPeptide} --IP /path/to/mgf_and_gLabel_files
```

#### Parameters:

`--InputType` or `-IT`:  specify either of `N` for N-linked or `O` for O-linked glycopeptides<br>
`--InputPath` or `-IP`:  path to `.mgf` and `-glabel.txt` files.<br>
NOTE:  for this function to work properly, you must pair the `.mgf` and `-glabel.txt` files as follows:<br>
If the spectral file is named `filename.msp`, then the corresponding gLabel file needs to be named `filename.mgf-glabel.txt`
</details>

### `parse_msp_to_csv.py`

This tool performs multiple functions:
1. perform *de novo* sequencing on `.msp` spectra (defined in `src/spectral_library/de_novo_sequencing.py`)
2. perform one-hot encoding of peptide sequence to an array of 20 amino acid codes, with the following added special codes:
   - X: represents I and L
   - J: glycosylated N
   - Z: zero-padding for consistent embedding length
3. create a linearized glycan sequence representation (as described in publication)

<details>
<summary>Usage:</summary>

```bash
python parse_msp_to_csv.py -IT=N -TN=10 -IP=/path/to/spectral/msp_files
```

#### Parameters:

`--InputType` or `-IT`:  specify either of `N` for N-linked or `O` for O-linked glycopeptides<br>
`--TopNumber` or `-TN`:  specify "top N" scoring *de novo* sequencing matches<br>
`--InputPath` or `-IP`:  path to `.msp` files.<br>

This tool assumes a strict naming to subfolders for the data:
1. for N-linked files:
   - .msp files should be placed under `N-GP-MSP`
   - *de novo* sequencing files will be written to `N-GP-DENOVO-MSP-TOP-{N}`, where `N` is the "top N" value specified in the `-TN` parameter
   - .csv output files will be written to `N-GP-CSV-TOP-{N}`, where `N` is the "top N" value specified in the `-TN` parameter
2. for O-linked files:
   - all folders will be prefixed by `O-` instead of `N-`
</details>

### `parse_csv_to_pkl.py`
This tool loads the `.csv` file containing one-hot encoded glycopeptide sequences into a Numpy array, then pickles it out to a `.pkl` file.
This prepares the processed data for hand-off to the ML model.

<details>
<summary>Usage:</summary>

```bash
python parse_csv_to_pkl.py -IT=N -TN=10 -IP=/path/to/csv_files
```

#### Parameters:

`--InputType` or `-IT`:  specify either of `N` for N-linked or `O` for O-linked glycopeptides<br>
`--TopNumber` or `-TN`:  specify "top N" scoring *de novo* sequencing matches (NOTE: this helps locate files from the previous script stage)<br>
`--InputPath` or `-IP`:  path to `.msp` files.<br>

</details>

## Model training and testing scripts:
### `train_model_*.py`
These scripts perform the training stage for ML models.  They are designed to train over a user-specified number of epochs.
Data is loaded from a folder containing annotated spectral `.pkl` files produced in the spectral library building stage.

These files are named for the specific training datasets used in model training in the publication.

<details>
<summary>Usage:</summary>

```bash
python train_model_batch.py -TA=/data/Training/N-GP-PKL -TE=/data/Testing/N-GP-PKL -SM=/data/saved_models
```

#### Parameters:

`--TrainingPath` or `-TA`:  Training data `.pkl` path, i.e. `/data/Training/N-GP-PKL`<br>
`--TestingPath` or `-TE`:  Testing data `.pkl` path, i.e. `/data/Testing/N-GP-PKL`<br>
`--SavedModelPath` or `-SM`:  Saved model path, i.e. `/data/saved_models`<br>

NOTE: some of the specific `train_model_*.py` scripts may have slightly different parameters.  Please inspect each script to determine the arguments needed.

</details>

### `predict_different_models.py` and `predict_different_models_batch.py`

These scripts perform model training and testing, whilst reporting real-time statistics to the console and summary statistics post-training.

<details>
<summary>Usage:</summary>

```bash
python predict_different_models.py -MP=/data/saved_models/model_path -IT=N -TN=10 -TP=/data/testing
```

#### Parameters:

`--ModelPath` or `-MP`: Model path, i.e. `/data/saved_models/model_path`<br>
`--InputType` or `-IT`:  specify either of `N` for N-linked or `O` for O-linked glycopeptides<br>
`--TopNumber` or `-TN`:  specify "top N" scoring *de novo* sequencing matches (NOTE: this helps locate files from the previous script stage)<br>
`--TestingPath` or `-TP`: Input and Output testing path, i.e. `/data/testing`<br>

</details>


## Results reporting scripts:
### `parse_msp_to_stat.py`

This tool reads a collection of `.msp` files and outputs summary statistics, such as:
- number of files, spectra, spectra for deep learning, samples
- total number of peptides, charge states, glycan charge states
- intensity statistics

Statistics are written back to the same containing folder as the `.msp` files, in two formats:  `.txt` and `.csv`.  The same statistics are also output to the console.

<details>
<summary>Usage:</summary>

```bash
python parse_msp_to_stat.py -IT=N -TN=10 -IP=/path/to/spectral/msp_files
```

#### Parameters:

`--InputType` or `-IT`:  specify either of `N` for N-linked or `O` for O-linked glycopeptides<br>
`--TopNumber` or `-TN`:  specify "top N" scoring *de novo* sequencing matches<br>
`--InputPath` or `-IP`:  path to `.msp` files.<br>

</details>

## Contributors:

* Zhewei Liang (zliang@venn.bio)
* Mingqi Liu (mliu@venn.bio)
* Richard "DJ" Shipman (richard.shipman@venn.bio)
* Norton Kitagawa (norton.kitagawa@venn.bio) - corresponding author
