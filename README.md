# Neural encoding of acoustic and semantic features during speech and music perception: Matlab to Python code translation

## Table of contents
1. [Project Description](https://github.com/neurogima/bhack_td/blob/main/README.md#project-description)
2. [Communication Channels](https://github.com/neurogima/bhack_td/blob/main/README.md#communication-channels)
3. [Goals for BrainHack Global](https://github.com/neurogima/bhack_td#goals-for-brainhack-global)
4. [Environment preparation](https://github.com/neurogima/bhack_td#environment-preparation)

## Project Description

In everyday life, humans are particularly attuned to listening to two
particular types of sound: speech and music. We apply a novel analysis method
to shed light on how the brain is almost effortlessly able to use acoustic
features to assign meaning to sounds. To do so, we use an original
cross-validated Representational Similarity Analyses (RSA) approach implemented
in Matlab to estimate the similarity between acoustic or semantic features of
an auditory stream (speech, music) and neural activity (here intracranial EEG
recordings decomposed into frequency bands).

## Communication channels

https://mattermost.brainhack.org/brainhack/channels/brainhack_marseille_2022_speech_music_representation

## Goals for Brainhack Global

The main goal of this project is to translate the Matlab code into Python:

- [ ] task 0: identify python libraries that can speed up the code translation effort
- [ ] task 1: translate preliminaries, data massaging
- [ ] task 2: translate the temporal folding on the neural signal in Python
- [ ] task 3: translate distance-computation metrics
  - [ ] sub-task 3.1: translate `BLG_CosDistND.m`
  - [ ] sub-task 3.1=2: translate `BLG_EucDistND.m`
- [ ] task 4: translate GLM computation and cross-validation stats
  - [ ] sub-task 4.1: translate `BLG_GLM_ND.m`
  - [ ] sub-task 4.2: implement CV stats  
- [ ] task 5: implement native plotting using the Python libraries (e.g. Seaborn, matplotlib, etc)
- [ ] task unassigned: build up Documentation pages that can be rendered as "read the doc" style

## Environment preparation

1. Install [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html)
   - If you do not have conda installed in your operating system, install [mambaforge](https://github.com/conda-forge/miniforge/releases)
   a wrapper around conda and a community project of the conda-forge
   community.
   - If you already have conda installed:
   `conda install --channel=conda-forge --name=base mamba`

2. Download the environment file with a starting set of packages [bh22_environment.yml](https://raw.githubusercontent.com/neurogima/bhack_td/main/bh22_environment.yml)
3. Install the environment for the BrainHack, using the file:
   `mamba env create -f bh22_environment.yml`
