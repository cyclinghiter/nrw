import pandas as pd
import numpy as np
from utils import read_excel, graphTool
from preprocessing import get_S_parameter

class NRWsolver():
    def __init__(self, freq, S11, S21, cutoff=0.02285*2, d=0.00969):
        c = 299792458 # m/s
        self.cutoff = cutoff
        self.d = d
        self.freq = freq
        self.lambda_0 = c / (self.freq * 10**6)
        self.S11 = S11
        self.S21 = S21


    def _get_reflection_coeff(self):
        reflection_coeff_1 = self.K + np.sqrt(self.K**2 - 1)
        reflection_coeff_2 = self.K - np.sqrt(self.K**2 - 1)
        
        reflection_coeff = np.where(np.abs(reflection_coeff_1) < 1,
                                    reflection_coeff_1, reflection_coeff_2)
        self.reflection_coeff = reflection_coeff
    
    def solve(self):
        self.K = (self.S11**2 - self.S21**2 +1) / (2*self.S11)
        self._get_reflection_coeff()
        self.transmition_coeff = (self.S11 + self.S21 - self.reflection_coeff) / (1 - self.reflection_coeff*(self.S11 + self.S21))
        self.lambda_0g = self.lambda_0 / np.sqrt(1-(self.lambda_0 / self.cutoff)**2)
        self.LAMBDA - (2*np.pi*self.d) / np.log(self.transmition_coeff)
        self.mu = self.lambda_0g / self.LAMBDA * (1+self.reflection_coeff) / (1-self.reflection_coeff)
        self.epsilon = self.lambda_0**2 * (1/self.LAMBDA**2 + 1/self.cutoff**2) / self.mu

df = read_excel('NRW~1\MeasData(19-12-23 17-10-06)_CNT.xls')

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

gT = graphTool((2,2), (10,10))
gT.plot((0,0), x=freq, y=real_mu, title='real_mu')
gT.plot((0,1), x=freq, y=imag_mu, title='imag_mu')
gT.plot((1,0), x=freq, y=real_epsilon, title = 'real_eps')
gT.plot((1,1), x=freq, y=imag_epsilon, title = 'imag_eps')
gT.show()

