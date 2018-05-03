# MIN1PIPE
A **MIN**iscope **1**-photon-based Calcium Imaging Signal Extraction **PIPE**line. 

This fully automatic, Matlab-based toolbox solves the full range problems of 1-photon calcium imaging (**data enhancement** + **movement morrection** + **signal extraction**) in one package, with minimal parameter tuning and semi-auto options integrated. The inidividual modules can also be easily adapted for two-photon imaging data processing.

## Usage
- Download/clone the git repository of the codes
- Open Matlab and set the MIN1PIPE folder as the current folder in Matlab
- Run `min1pipe.m`, the code automatically sets the package to the path, and processes the data the user specifies.

## Demo
Users can run `demo_min1pipe.m` for a demo of 1-photon calcium imaging video recorded with UCLA miniscope. The same code can also be adapted to custom scripts for the processing.

## Brief Intro
MIN1PIPE contains essentially three steps: neural enhancing, movement correction, and neural signal extraction.
- Neural Enhancing: remove spatial noise and then adaptively remove non-neural background in the field of view in a framewise manner.
- Movement Correction: remove field of view movement with a hierarchical designed nonrigid movement correction module (integrating KLT Tracker and LogDemons deformation registration method), which is free of assumption about movement type and amplitude.
- Neural Signal Extraction: identify neuronal ROIs and corresponding calcium traces with minimal false positive rates (incorporating GMM, LSTM as true neuron selector and modified CNMF as spatiotemporal calcium signal identifier)
- Semi-auto options: we also provide semi-auto options for certain sections including
  - Automated manual seeds selection module: for users who want to manually select seeds of neuron ROIs that will result in **ZERO** false positives
  - Post-process exclusion of "bad" neural components
  
  
## References
Please cite this paper if MIN1PIPE help your research:
```
**Lu, J.**, Li, C., Singh-Alvarado, J., Zhou, Z., Fröhlich, F., Mooney, R., & Wang, F. (2018). 
MIN1PIPE: A Miniscope 1-photon-based Calcium Imaging Signal Extraction Pipeline. 
*bioRxiv*, 311548.
```
or related reference:

**Lu, J.**, Li, C., & Wang, F. (2017). Seeds Cleansing CNMF for Spatiotemporal Neural Signals Extraction of Miniscope Imaging Data. *arXiv preprint* arXiv:1704.00793.

## Questions?
Please email to `min1pipe2018@gmail.com` provided additional questions.
