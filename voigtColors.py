#Color definitions to match Voigt et al. paper
#(Import this to use in any script)

import numpy as np

voigtColors = {'AM2': np.array([255,204,153])/255,
               'CAM3': np.array([128,128,128])/255,
               'CAM4': np.array([148,255,181])/255,
               'CNRM-AM6-DIA-v2': np.array([0,51,128])/255, 
               'CaltechGray': np.array([255,164,5])/255, 
               'ECHAM-6.1': np.array([0,117,220])/255,
               'ECHAM-6.3': np.array([153,63,0])/255,
               'GISS': np.array([157,204,0])/255,
               'IPSL-CM5A': np.array([76,0,92])/255, 
               'MIROC5': np.array([43,206,72])/255, 
               'MPAS': np.array([143,124,0])/255, 
               'MetUM-GA6-CTL': np.array([25,25,25])/255, 
               'MetUM-GA6-ENT': np.array([0,92,49])/255, 
               'NorESM2': np.array([194,0,136])/255}