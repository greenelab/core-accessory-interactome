#!/bin/bash
conda env create --file environment.yml
conda activate core_acc
pip install -e .