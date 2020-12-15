import numpy as np
import matplotlib.pyplot as plt 
from preprocessing import get_S_parameter

class NRWsolver():
    def __init__(self, freq, S11, S21, cutoff=0.02285*2, d=0.00969):
        c = 299792458 # m/s
        self.cutoff = cutoff
        self.d = d
        self.freq = np.asarray(freq)
        self.lambda_0 = c / (self.freq * 10**6)
        self.S11 = np.asarray(S11)
        self.S21 = np.asarray(S21)
        
    def solve(self):
        self.V1 = self.S21 + self.S11
        self.V2 = self.S21 - self.S11
        
        self.X = (1 - self.V1 * self.V2) / (self.V1 - self.V2)
        self._get_reflection_coeff()
        
        self.transmission_coeff = (self.S11 + self.S21 - self.reflection_coeff) / (1 - self.reflection_coeff*(self.S11 + self.S21))
        self.lambda_0g = self.lambda_0 / np.sqrt(1-(self.lambda_0 / self.cutoff)**2)
        log_term = np.log(self.transmission_coeff)
        log_term[np.imag(log_term) > 0 ] -= 2J*np.pi
        
        self.one_over_LAMBDA = 1J / (2*np.pi*self.d) * log_term
        self.LAMBDA = 1 / self.one_over_LAMBDA
        
        self.mu = self.lambda_0g / self.LAMBDA * (1+self.reflection_coeff) / (1-self.reflection_coeff)
        self.epsilon = self.lambda_0**2 * (1/self.LAMBDA**2 + 1/self.cutoff**2) / self.mu

    def _get_reflection_coeff(self):
        reflection_coeff_1 = self.X + np.sqrt(self.X**2 - 1)
        reflection_coeff_2 = self.X - np.sqrt(self.X**2 - 1)
        reflection_coeff = np.where(np.abs(reflection_coeff_1) < 1,
                                    reflection_coeff_1, reflection_coeff_2)
        self.reflection_coeff = reflection_coeff
 

class ErrorPropSolver:
    def __init__(self, freq, S11_list, S21_list, dx=0.01):
        self.freq = freq
        self.num_sample = len(freq)
        self.S11_list = S11_list
        self.S21_list = S21_list
        self.S11 = np.mean(S11_list, axis=0)
        self.S21 = np.mean(S21_list, axis=0)
        self.S11 
        self.dx = dx
    
    def solve_mean(self):
        self.solver_mean = NRWsolver(self.freq, self.S11, self.S21)
        self.solver_mean.solve()
        
    def solve_next(self, mode='Sr1'):
        if mode == 'Sr1':
            self.solver_Sr1_next = NRWsolver(self.freq, self.S11 + self.dx, self.S21)
            self.solver_Sr1_next.solve()
        if mode == 'Si1':
            self.solver_Si1_next = NRWsolver(self.freq, self.S11 + 1J * self.dx, self.S21)
            self.solver_Si1_next.solve()
        if mode == 'Sr2':
            self.solver_Sr2_next = NRWsolver(self.freq, self.S11, self.S21 + self.dx)
            self.solver_Sr2_next.solve()
        if mode == 'Si2':
            self.solver_Si2_next = NRWsolver(self.freq, self.S11, self.S21 + 1J * self.dx)
            self.solver_Si2_next.solve()
            
    def solve_prev(self, mode='Sr1'):
        if mode == 'Sr1':
            self.solver_Sr1_prev = NRWsolver(self.freq, self.S11 - self.dx, self.S21)
            self.solver_Sr1_prev.solve()
        if mode == 'Si1':
            self.solver_Si1_prev = NRWsolver(self.freq, self.S11 - 1J * self.dx, self.S21)
            self.solver_Si1_prev.solve()
        if mode == 'Sr2':
            self.solver_Sr2_prev = NRWsolver(self.freq, self.S11, self.S21 - self.dx)
            self.solver_Sr2_prev.solve()
        if mode == 'Si2':
            self.solver_Si2_prev = NRWsolver(self.freq, self.S11, self.S21 - 1J * self.dx)
            self.solver_Si2_prev.solve()    
    
    def cal_Jacobian_mu(self):
        modes = ['Sr1', 'Si1', 'Sr2', 'Si2']
        for mode in modes:
            self.solve_next(mode)
            self.solve_prev(mode)
        
        self.J11 = np.real(self.solver_Sr1_next.mu) - np.real(self.solver_Sr1_prev.mu) / (2*self.dx)
        self.J12 = np.real(self.solver_Si1_next.mu) - np.real(self.solver_Si1_prev.mu) / (2*self.dx)
        self.J13 = np.real(self.solver_Sr2_next.mu) - np.real(self.solver_Sr2_prev.mu) / (2*self.dx)
        self.J14 = np.real(self.solver_Si2_next.mu) - np.real(self.solver_Si2_prev.mu) / (2*self.dx)
        self.J21 = np.imag(self.solver_Sr1_next.mu) - np.imag(self.solver_Sr1_prev.mu) / (2*self.dx)
        self.J22 = np.imag(self.solver_Si1_next.mu) - np.imag(self.solver_Si1_prev.mu) / (2*self.dx)
        self.J23 = np.imag(self.solver_Sr2_next.mu) - np.imag(self.solver_Sr2_prev.mu) / (2*self.dx)
        self.J24 = np.imag(self.solver_Si2_next.mu) - np.imag(self.solver_Si2_prev.mu) / (2*self.dx)
        
        Jacobian = []
        for idx in range(self.num_sample):
            Jacobian.append(np.array([[self.J11[idx], self.J12[idx], 
                                       self.J13[idx], self.J14[idx]], 
                                      [self.J21[idx], self.J22[idx], 
                                       self.J23[idx], self.J24[idx]]]))
        Jacobian = np.asarray(Jacobian)
        self.Jacobian = Jacobian
    
    def cal_V_mu(self):
        V_list = []
        for idx in range(self.num_sample):
            S11_real = np.real(self.S11_list)[:, idx]
            S11_imag = np.imag(self.S11_list)[:, idx]
            S21_real = np.real(self.S21_list)[:, idx]
            S21_imag = np.imag(self.S21_list)[:, idx]
            freq_measure = np.c_[S11_real, S11_imag, S21_real, S21_imag]
            V = np.cov(freq_measure.T)
            V_list.append(V)
        self.V_list = V_list
    
    def cal_mu_error(self):
        self.cal_Jacobian_mu()
        self.cal_V_mu()
        error_list = []
        for idx in range(self.num_sample):
            error_matrix = np.matmul(self.Jacobian[idx], 
                                     np.matmul(self.V_list[idx], self.Jacobian[idx].T))
            error_list.append(error_matrix)
        error_list = np.asarray(error_list)
        return error_list 

    def cal_Jacobian_epsilon(self):
        modes = ['Sr1', 'Si1', 'Sr2', 'Si2']
        for mode in modes:
            self.solve_next(mode)
            self.solve_prev(mode)
        
        self.J11 = np.real(self.solver_Sr1_next.epsilon) - np.real(self.solver_Sr1_prev.epsilon) / (2*self.dx)
        self.J12 = np.real(self.solver_Si1_next.epsilon) - np.real(self.solver_Si1_prev.epsilon) / (2*self.dx)
        self.J13 = np.real(self.solver_Sr2_next.epsilon) - np.real(self.solver_Sr2_prev.epsilon) / (2*self.dx)
        self.J14 = np.real(self.solver_Si2_next.epsilon) - np.real(self.solver_Si2_prev.epsilon) / (2*self.dx)
        self.J21 = np.imag(self.solver_Sr1_next.epsilon) - np.imag(self.solver_Sr1_prev.epsilon) / (2*self.dx)
        self.J22 = np.imag(self.solver_Si1_next.epsilon) - np.imag(self.solver_Si1_prev.epsilon) / (2*self.dx)
        self.J23 = np.imag(self.solver_Sr2_next.epsilon) - np.imag(self.solver_Sr2_prev.epsilon) / (2*self.dx)
        self.J24 = np.imag(self.solver_Si2_next.epsilon) - np.imag(self.solver_Si2_prev.epsilon) / (2*self.dx)
        
        Jacobian = []
        for idx in range(self.num_sample):
            Jacobian.append(np.array([[self.J11[idx], self.J12[idx], 
                                       self.J13[idx], self.J14[idx]], 
                                      [self.J21[idx], self.J22[idx], 
                                       self.J23[idx], self.J24[idx]]]))
        Jacobian = np.asarray(Jacobian)
        self.Jacobian = Jacobian
    
    def cal_V_epsilon(self):
        V_list = []
        for idx in range(self.num_sample):
            S11_real = np.real(self.S11_list)[:, idx]
            S11_imag = np.imag(self.S11_list)[:, idx]
            S21_real = np.real(self.S21_list)[:, idx]
            S21_imag = np.imag(self.S21_list)[:, idx]
            freq_measure = np.c_[S11_real, S11_imag, S21_real, S21_imag]
            V = np.cov(freq_measure.T)
            V_list.append(V)
        self.V_list = V_list
    
    def cal_epsilon_error(self):
        self.cal_Jacobian_epsilon()
        self.cal_V_epsilon()
        error_list = []
        for idx in range(self.num_sample):
            error_matrix = np.matmul(self.Jacobian[idx], 
                                     np.matmul(self.V_list[idx], self.Jacobian[idx].T))
            error_list.append(error_matrix)
        error_list = np.asarray(error_list)
        return error_list 

