# TRACMIP_pa_notebooks
Jupyter notebooks used to analyze polar amplification in the [TRACMIP](https://sites.google.com/site/tracmip/) multi-model ensemble.

This code is for a paper by Rick Russotto and Michela Biasutti published in Geophysical Research Letters in 2020: https://doi.org/10.1029/2019GL086771

Notebooks are now permanently archived on Zenodo, accessible here: https://doi.org/10.5281/zenodo.3711669



### Notebooks used to make each figure: 
Figure 1: EBM_Gregory_analyze.ipynb

Figure 2, S2: AnalyzeGregoryFeedbacks_v2.ipynb

Figure 3, S3, S8: EBM_Gregory_analyze_part3.ipynb

Figure 4, S4, S5, S9: EBM_BigScatter_TropicsPoles.ipynb

Figure S5, S6: SupplementFigures.ipynb




### Other notebooks and scripts that the above depend on: 
#### Notebooks that ran EBM experiments: 
EBM_Gregory_run_noG.ipynb

EBM_LocalRemote_Gregory_run_noG.ipynb

EBM_noQ_run.ipynb

#### Notebooks that ran Gregory regressions:
GregoryTRACMIP.ipynb

GregoryIndividualFeedbacks.ipynb

#### Notebooks that calculated monthly TOA radiative changes for regressions:
MonthlyAPRP.ipynb

MonthlyKernels.ipynb

MonthlyKernels_CloudLW.ipynb

#### Miscellaneous:

voigtColors.py

APRPX.py--adapted from https://github.com/rdrussotto/pyaprp

