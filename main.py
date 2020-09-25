import pandas as pd
import numpy as np
from preprocessing import get_S_parameter
from utils import read_excel, graphTool
from core import NRWsolver

df = read_excel('NRW~1\MeasData(19-12-23 17-05-23)_air.xls')

freq = np.asarray(df.index)

S11_DB = df[('LOG MAG [dB]', 'S11')]
S11_phase = df[(( 'PHASE [deg]', 'S11'))]
S11 = get_S_parameter(S11_DB, S11_phase)

S12_DB = df[('LOG MAG [dB]', 'S12')]
S12_phase = df[(( 'PHASE [deg]', 'S12'))]
S12 = get_S_parameter(S12_DB, S12_phase)

solver = NRWsolver(freq, S11, S12)
solver.solve()

real_mu = np.real(solver.mu)
imag_mu = np.imag(solver.mu)
real_epsilon = np.real(solver.epsilon)
imag_epsilon = np.imag(solver.epsilon)

gT = graphTool((2,1), (10,10))
gT.plot((0,0), x=freq, y=real_mu, c='r', label='real_mu')
gT.plot((0,0), x=freq, y=imag_mu, c='b', label='imag_mu')
gT.plot((1,0), x=freq, y=real_epsilon, c='r', label = 'real_eps')
gT.plot((1,0), x=freq, y=imag_epsilon, c='b', label = 'imag_eps')

gT.show()