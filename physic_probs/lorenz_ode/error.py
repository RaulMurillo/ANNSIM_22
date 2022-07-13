from lib2to3.pgen2.pgen import DFAState
import numpy as np

LOGDIR = './build/'
rows = 4001
cols = 1+3

d_a   = np.genfromtxt(LOGDIR+'lorenz_ode_data_double.txt')#, delimiter=' ')
time = d_a[:,0]
d_a = d_a[:,1:]

f_a   = np.genfromtxt(LOGDIR+'lorenz_ode_data_float.txt')[:,1:]
hf_a  = np.genfromtxt(LOGDIR+'lorenz_ode_data_float16.txt')[:,1:]
bf_a  = np.genfromtxt(LOGDIR+'lorenz_ode_data_bfloat.txt')[:,1:]
p32_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_32_2.txt')[:,1:]
p28_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_28_2.txt')[:,1:]
p24_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_24_2.txt')[:,1:]
p20_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_20_2.txt')[:,1:]
p18_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_18_2.txt')[:,1:]
p16_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_16_2.txt')[:,1:]
p14_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_14_2.txt')[:,1:]
p12_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_12_2.txt')[:,1:]
p10_a = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_10_2.txt')[:,1:]
p8_a  = np.genfromtxt(LOGDIR+'lorenz_ode_data_posit_8_2.txt')[:,1:]


# print(d_a.shape, d_a[0], d_a[-1])

format = ['float', 'half float', 'Bfloat', 'posit_32', 'posit_28', 'posit_24', 'posit_20', 'posit_18', 'posit_16', 'posit_14', 'posit_12', 'posit_10', 'posit_8']
for i, a in enumerate([f_a, hf_a, bf_a, p32_a, p28_a, p24_a, p20_a, p18_a, p16_a, p14_a, p12_a, p10_a, p8_a]):
    # error = 0
    # for n in range(rows):
    #     dist = np.linalg.norm(d_a[n]-a[n])
    #     error += dist**2
    # error /= rows
    error = np.mean((d_a-a)**2)

    print(f'{format[i]} error: {error:.2f}')
