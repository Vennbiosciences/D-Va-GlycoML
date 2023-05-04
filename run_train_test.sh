#!/bin/bash
set -eu
DIR=$(cd "$(dirname "$0")"; pwd)

if [[ $# -lt 3 ]]; then
  echo "Usage: $(basename "$0") TRAIN_DIR TEST_DIR MODEL_DIR"
  exit 1
fi

TRAIN_DIR=$1
TEST_DIR=$2
MODEL_DIR=$3

if [[ -d ${MODEL_DIR} ]]; then
  echo "Model already exists"
  exit 2
fi

mkdir -p "${MODEL_DIR}"

python ./src/spectral_library/construct_annotated_library.py -IT=N -IP="${TRAIN_DIR}"
python ./src/spectral_library/parse_msp_to_csv.py -IT=N -TN=1 -IP="${TRAIN_DIR}"
python ./src/spectral_library/parse_csv_to_pkl.py -IT=N -TN=1 -IP="${TRAIN_DIR}"
python ./src/spectral_library/parse_msp_to_stat.py -IT=N -TN=1 -IP="${TRAIN_DIR}"
python -u ./src/deep_learning/train_model_cyno.py \
  -TA="${TRAIN_DIR}/N-GP-PKL-TOP-1" \
  -SM="${MODEL_DIR}"
python -u ./src/deep_learning/predict_different_models_background.py \
  -MP="${MODEL_DIR}" \
  -IT=N \
  -TN=1 \
  -TP="${TEST_DIR}"
