import pandas as pd
import numpy as np
from utils import read_excel

df = read_excel('NRW~1\MeasData(19-12-23 17-10-06)_CNT.xls')

def decibel_to_linear(decibel):
    return 10**(decibel / 20)
    
def magnitude_phase_to_real_imaginary(magnitude, phase):
    real = magnitude * np.cos(phase * np.pi / 180)
    imag = magnitude * np.sin(phase * np.pi / 180)
    return real + 1J* imag

def get_S_parameter(DB, phase):
    mag = decibel_to_linear(DB)
    scattering_parameter = magnitude_phase_to_real_imaginary(mag, phase)
    return scattering_parameter