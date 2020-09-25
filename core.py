import numpy as np
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
        self.LAMBDA = 1 / (1J / (2*np.pi*self.d) * np.log(self.transmition_coeff)) 
        self.mu = self.lambda_0g / self.LAMBDA * (1+self.reflection_coeff) / (1-self.reflection_coeff)
        self.epsilon = self.lambda_0**2 * (1/self.LAMBDA**2 + 1/self.cutoff**2) / self.mu



