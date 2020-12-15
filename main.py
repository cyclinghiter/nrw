import pandas as pd
import numpy as np
from preprocessing import get_S_parameter
from utils import read_excel, graphTool
from core import NRWsolver, ErrorPropSolver
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

def func(xdata, constant):
    return constant 

def get_S(fname):
    df = read_excel(fname, header=[0])
    df.head()

    freq = np.asarray(df.index)

    dx = 0.01
    # gT = graphTool((2,1), (10,10))


    S11_DB = df[('log_s11')]
    S11_phase = df[(('phase_s11'))]
    S11 = get_S_parameter(S11_DB, S11_phase)
    S11 = S11.reshape(1,801)

    S21_DB = df[('log_s21')]
    S21_phase = df[(('phase_s21'))]
    S21 = get_S_parameter(S21_DB, S21_phase)
    S21 = S21.reshape(1,801)
    return freq, S11, S21

freq, S11, S21 = get_S('NRW~1/epoxy/epoxy.xlsx')
freq, S11_2, S21_2 = get_S('NRW~1 (2)/epoxy.xlsx')

S11_list = np.r_[S11, S11_2]
S21_list = np.r_[S21, S21_2]

std = np.std(np.real(S11_list), axis=0)
imag_std =np.std(np.imag(S11_list), axis=0)

errorsolver = ErrorPropSolver(freq, S11_list, S21_list, dx=0.01)
errorsolver.solve_mean()
epsilon = np.real(errorsolver.solver_mean.epsilon)
error = np.sqrt(errorsolver.cal_epsilon_error()[:,0,0])

epsilon_imag = np.imag(errorsolver.solver_mean.epsilon)
error_imag = np.sqrt(errorsolver.cal_epsilon_error()[:,1,1])

popt , pcov = curve_fit(func, xdata=freq, ydata=epsilon, sigma = error)
print(np.sum((epsilon - popt)**2 / (epsilon)))

popt , pcov = curve_fit(func, xdata=freq, ydata=epsilon_imag, sigma = error_imag)
print(np.sum((epsilon_imag - popt)**2 / (epsilon_imag)))

fig, axes = plt.subplots(2,2, figsize=(10,10))

axes[0,0].plot(freq/1000, np.real(errorsolver.S11))
axes[0,0].fill_between(freq/1000, np.real(errorsolver.S11) - std/2 , np.real(errorsolver.S11) + std/2, alpha=0.3, edgecolor='white')
axes[0,0].set_xlabel("frequency (GHz)")
axes[0,0].set_ylabel("real(S11)")
axes[0,0].set_title("S11 real")

axes[0,1].plot(freq/1000, np.imag(errorsolver.S11))
axes[0,1].fill_between(freq/1000, np.imag(errorsolver.S11) - imag_std/2 , np.imag(errorsolver.S11) + imag_std/2, alpha=0.3, edgecolor='white')
axes[0,1].set_xlabel("frequency (GHz)")
axes[0,1].set_ylabel("imag(S11)")
axes[0,1].set_title("S11 imag")

std = np.std(np.real(S21_list), axis=0)
imag_std =np.std(np.imag(S21_list), axis=0)

axes[1,0].plot(freq/1000, np.real(errorsolver.S21))
axes[1,0].fill_between(freq/1000, np.real(errorsolver.S21) - std/2 , np.real(errorsolver.S21) + std/2, alpha=0.3, edgecolor='white')
axes[1,0].set_xlabel("frequency (GHz)")
axes[1,0].set_ylabel("real(S21)")
axes[1,0].set_title("S21 real")

axes[1,1].plot(freq/1000, np.imag(errorsolver.S21))
axes[1,1].fill_between(freq/1000, np.imag(errorsolver.S21) - imag_std/2 , np.imag(errorsolver.S21) + imag_std/2, alpha=0.3, edgecolor='white')
axes[1,1].set_xlabel("frequency (GHz)")
axes[1,1].set_ylabel("imag(S21)")
axes[1,1].set_title("S21 imag")

plt.savefig("epoxy.")

plt.style.use("classic")
fig, axes = plt.subplots(1,2, figsize=(12,5))
axes[0].plot(freq/1000, epsilon, c='red', label='real')
axes[0].plot(freq/1000, epsilon_imag, c='green', label='imag')
axes[0].set_xlabel("frequency (GHz)")
axes[0].set_ylabel("epsilon_r")
axes[0].fill_between(freq/1000, epsilon-error/2, epsilon+error/2, 
                     fc='red', alpha=0.5, edgecolor = 'white')
axes[0].fill_between(freq/1000, epsilon_imag-error_imag/2, epsilon_imag+error_imag/2, 
                     fc='green', alpha=0.9, edgecolor = 'white')

axes[0].set_ylim(-2,10)
axes[0].set_title("Permittivity of CNT")
axes[0].legend()


mu = np.real(errorsolver.solver_mean.mu)
error = np.sqrt(errorsolver.cal_mu_error()[:,0,0])

mu_imag = np.imag(errorsolver.solver_mean.mu)
error_imag = np.sqrt(errorsolver.cal_mu_error()[:,1,1])

axes[1].plot(freq/1000, mu, c='red', label='real')
axes[1].plot(freq/1000, mu_imag, c='green', label='imag')
axes[1].set_xlabel("frequency (GHz)")
axes[1].set_ylabel("mu_r")
axes[1].fill_between(freq/1000, mu-error/2, mu+error/2, 
                     fc='red', alpha=0.5, edgecolor = 'white')
axes[1].fill_between(freq/1000, mu_imag-error_imag/2, mu_imag+error_imag/2, 
                     fc='green', alpha=0.9, edgecolor = 'white')

axes[1].set_ylim(-2,10)
axes[1].set_title("Permeability of CNT")
axes[1].legend()
plt.savefig("CNT.jpg")


popt , pcov = curve_fit(func, xdata=freq, ydata=mu, sigma = error)
print(np.sum((mu - popt)**2 / (mu)))


popt , pcov = curve_fit(func, xdata=freq, ydata=mu_imag, sigma = error_imag)
print(np.sum((mu_imag - popt)**2 / (mu_imag)))


# # df = read_excel('NRW~1\CNT\CNT.xlsx', header=[0])
# # df.head()

# # freq = np.asarray(df.index)

# # S11_DB = df[('log_s11')]
# # S11_phase = df[(('phase_s11'))]
# # S11 = get_S_parameter(S11_DB, S11_phase)


# # S11 = S11 + np.random.normal(0, scale=100*1e-4, size=S11.shape)

# # S21_DB = df[('log_s21')]
# # S21_phase = df[(('phase_s21'))]
# # S21 = get_S_parameter(S21_DB, S21_phase)

# # solver = NRWsolver(freq, S11, S21)
# # solver.solve()

# # real_mu = np.real(solver.mu)
# # imag_mu = np.imag(solver.mu)
# # real_epsilon = np.real(solver.epsilon)
# # imag_epsilon = np.imag(solver.epsilon)

# # gT = graphTool((2,1), (10,10))
# # gT.plot((0,0), x=freq, y=real_mu, c='r', label='real_mu')
# # gT.plot((0,0), x=freq, y=imag_mu, c='b', label='imag_mu')
# # gT.plot((1,0), x=freq, y=real_epsilon, c='r', label = 'real_eps')
# # gT.plot((1,0), x=freq, y=imag_epsilon, c='b', label = 'imag_eps')
# # gT.save("rog_lin")
