import numpy as np

LOGDIR = './results/'
Ls = 1024
Fs = 2000

d_a   = np.genfromtxt(LOGDIR+'double.log', delimiter=',')
f_a   = np.genfromtxt(LOGDIR+'float.log', delimiter=',')
hf_a  = np.genfromtxt(LOGDIR+'float_16_5.log', delimiter=',')
bf_a  = np.genfromtxt(LOGDIR+'float_16_8.log', delimiter=',')
p32_a = np.genfromtxt(LOGDIR+'posit_32.log', delimiter=',')
p28_a = np.genfromtxt(LOGDIR+'posit_28.log', delimiter=',')
p24_a = np.genfromtxt(LOGDIR+'posit_24.log', delimiter=',')
p20_a = np.genfromtxt(LOGDIR+'posit_20.log', delimiter=',')
p16_a = np.genfromtxt(LOGDIR+'posit_16.log', delimiter=',')
p14_a = np.genfromtxt(LOGDIR+'posit_14.log', delimiter=',')
p12_a = np.genfromtxt(LOGDIR+'posit_12.log', delimiter=',')
p10_a = np.genfromtxt(LOGDIR+'posit_10.log', delimiter=',')
p8_a  = np.genfromtxt(LOGDIR+'posit_8.log', delimiter=',')


dd = (d_a[2]**2+d_a[3]**2)/Ls
d_PS = dd[1:Ls//2+1]
d_PS[2:-1] = 2*d_PS[2:-1]

format = ['float', 'half float', 'Bfloat', 'posit_32', 'posit_28', 'posit_24', 'posit_20', 'posit_16', 'posit_14', 'posit_12', 'posit_10', 'posit_8']
for i, a in enumerate([f_a, hf_a, bf_a, p32_a, p28_a, p24_a, p20_a, p16_a, p14_a, p12_a, p10_a, p8_a]):
    aa = (a[2]**2+a[3]**2)/Ls
    a_PS = aa[1:Ls//2+1]
    a_PS[2:-1] = 2*a_PS[2:-1]

    error = np.mean((d_PS-a_PS)**2)

    print(f'{format[i]} error: {error:.2e}')