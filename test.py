import pandas as pd
import numpy as np
from utils import read_excel, graphTool

df = read_excel('NRW~1\MeasData(19-12-23 17-10-06)_CNT.xls')
df.columns
gT = graphTool(figshape=(1,1), figsize= (5,5))


freq = np.asarray(df.index)
S11 = np.exp(df[('LOG MAG [dB]', 'S11')].to_numpy())
S11_i  =
S21 = np.exp(df[('LOG MAG [dB]', 'S21')].to_numpy())

gT.plot(ax_num=(0,0), x = freq, y = S11_logMagnitude, c='r', label ='S11')
gT.plot(ax_num=(0,0), x = freq, y = S21_logMagnitude, c='b', label ='S21')
**2 
S11**2 