# iCBIHT
Implementation of journal paper entitled '1-Bit direction of arrival estimation via improved complex-valued binary iterative hard thresholding' for DOA estimation [Digital Signal Processing]

**Notation:** Codes for comparison methods are not available. Please remove them from `main.m` manually.

## Citation
```
@article{WANG2022103265,
title = {1-Bit direction of arrival estimation via improved complex-valued binary iterative hard thresholding},
journal = {Digital Signal Processing},
volume = {120},
pages = {103265},
year = {2022},
issn = {1051-2004},
doi = {https://doi.org/10.1016/j.dsp.2021.103265},
url = {https://www.sciencedirect.com/science/article/pii/S1051200421003043},
author = {Pengyu Wang and Huichao Yang and Zhongfu Ye},
keywords = {DOA estimation, Sparse representation, 1-bit quantization, Sensor array},
abstract = {Aiming to estimate direction-of-arrival (DOA) using 1-bit quantized observation of sensor arrays, an improved complex-valued binary iterative hard thresholding (iCBIHT) algorithm is proposed in this research. In this work, an error function of signal reconstruction is defined. The signals are estimated by gradient descending. A refined hard threshold function and a novel stopping criterion are designed to judge the convergence. A bases-updating strategy is introduced to solve off-grid DOAs, which improves the accuracy of DOA estimation. A backtracking strategy is utilized to refine the estimated signals processed by the hard thresholding function as well as promote convergence. Simulations and analyses show that the proposed algorithm is superior to 1-bit multiple signal classification (MUSIC) and complex-valued binary iterative hard thresholding (CBIHT) with few snapshots and low signal-to-noise ratio (SNR).}
}
```
