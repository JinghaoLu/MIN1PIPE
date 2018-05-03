# MIN1PIPE
[A **MIN**iscope **1**-photon-based Calcium Imaging Signal Extraction **PIPE**line](https://www.biorxiv.org/content/early/2018/04/30/311548). 

MIN1PIPE is a fully automatic, Matlab-based toolbox, solving the full range problems in 1-photon calcium imaging (**data enhancement** + **movement morrection** + **signal extraction**) in one package. It requires minimal parameter-tuning and integrates the semi-auto options. The inidividual modules can also be easily adapted for 2-photon imaging data processing.

## Contents
1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Usage](#usage)
    1. [Demo](#demo)
    2. [Custom data](#custom-data)

## Introduction
MIN1PIPE contains the following three steps:
- **`Neural Enhancing`**: remove spatial noise and then adaptively remove non-neural background in the field of view in a framewise manner.
- **`Movement Correction`**: remove field of view movement with a hierarchical designed nonrigid movement correction module (integrating KLT Tracker and LogDemons deformation registration method), which is free of assumption about movement type and amplitude.
- **`Neural Signal Extraction`**: identify neuronal ROIs and corresponding calcium traces with minimal false positive rates (incorporating GMM, LSTM as true neuron selector and modified CNMF as spatiotemporal calcium signal identifier)

**`Semi-auto options`**: We also provide semi-auto options, including
  - Automated manual seeds selection module: for users who want to manually select seeds of neuron ROIs that will result in **ZERO** false positives
  - Post-process exclusion of "bad" neural components
  
## Dependencies
This Matlab implementation has the following dependencies (included under [`utilities`](./utilities)):
- [CVX](http://cvxr.com/cvx/)

## Usage
- Download/clone the git repository of the codes
- Open Matlab and set the MIN1PIPE folder as the current folder in Matlab
- Run [`min1pipe.m`](./min1pipe.m), the code automatically sets the package to the path, and processes the data the user specifies. For argument options, please see [`min1pipe.m`](./min1pipe.m)

### Demo
As a demo, we demonstrate the use of 1-photon calcium imaging video recorded with UCLA miniscope:
```
demo_min1pipe.m
```
The same code can also be adapted to custom scripts for the processing.

### Custom data
To use the code on a custom dataset, you need to put your dataset under the foler [demo](./demo) as "*.tif"
  
## References
Please cite [the MIN1PIPE paper](https://www.biorxiv.org/content/early/2018/04/30/311548) if it helps your research.
```
@article{lu2018MIN1PIPE,
  title={MIN1PIPE: A Miniscope 1-photon-based Calcium Imaging Signal Extraction Pipeline.},
  author={Lu, J., Li, C., Singh-Alvarado, J., Zhou, Z., Fröhlich, F., Mooney, R., & Wang, F.},
  journal={bioRxiv, 311548},
  year={2018}
}
```
or the related [arXiv version](https://arxiv.org/abs/1704.00793):
```
@article{lu2017seeds,
  title={Seeds Cleansing CNMF for Spatiotemporal Neural Signals Extraction of Miniscope Imaging Data},
  author={Lu, Jinghao and Li, Chunyuan and Wang, Fan},
  journal={arXiv preprint arXiv:1704.00793},
  year={2017}
}
```

## Questions?
Please email to `min1pipe2018@gmail.com` for additional questions.
