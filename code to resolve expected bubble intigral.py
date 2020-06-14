import scipy as sp
from scipy import special

import scipy.integrate as integrate

import matplotlib.pyplot as plt
import math
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

#define constants
a_inf =  1.091*10**(-29)
a_o = 1
H_o = 1.5*10**(-42)
m = 1.8*10**13
bar = 3.46e38

abs_lambda = 0.01

H = 1.1937e8 

eta_o_eta_inf = 3.21/(a_o * H_o)

def eta_inf_eta_n(n):
    coeff = 1/(a_inf* m)*(3 * sp.pi/2*sp.e)**0.5
    erfi_1 = (n + 1/2 )**0.5
    erfi_2 = sp.sqrt(0.5)
    return coeff*(special.erfi(erfi_1) - special.erfi(erfi_2))

def eta_o_eta_N(n):
    return eta_o_eta_inf + eta_inf_eta_n(n)

#this is a fucntion descibing the inside of the intigral
#the constants outside the initgral are not concidered

def integral(N):
    return (1+2*N)**(3/2) * math.exp(-3*N) * (eta_o_eta_N(N))**3

#%%
#Resolve the intigral - as in terms of e-foldings (N) we only concider from 0 to 60

result_0_60 = sp.integrate.quad(integral, 0, 60)

print(result_0_60[0])

#%%
#find the xi value, set the expected no. to 1 and take the log for 90pi^2xi^2/lambda

e_exponent = 4*3**0.5/27 * m**3 * a_inf**3 * result_0_60[0]

exponent = sp.log(e_exponent) 

xi_squared = exponent/(96*sp.pi**2/abs_lambda)

xi = xi_squared**0.5

print(xi)
