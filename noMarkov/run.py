from src.tools import *
import sys as s
from sympy import *
from scipy.linalg import sqrtm
import numpy as np
init_printing(use_unicode=True)
from matplotlib import pyplot as plt
#%matplotlib inline
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum import TensorProduct
import scipy.interpolate
import platform

def werner_state(c1, c2, c3):
    # c = [-0.8,-0.8,-0.8]
    c = [c1, c2, c3]
    index = 0
    rho = np.zeros((4,4),dtype=complex)
    for i in range(len(rho)):
        rho[i,i] = 1
    for i in c:
        index += 1
        rho += TensorProduct(i*Pauli(index),Pauli(index))
    #print(np.array(rho,dtype=complex))
    return rho
werner = werner_state(-0.8, -0.8, -0.8)
#stical modeling 
Mais=(cb(2,0)+cb(2,1))/sqrt(2)
Menos=(cb(2,0)-cb(2,1))/sqrt(2)
#werner

'''phase flip'''
def K_0(J):
    return sqrt(1-J/2)*Pauli(0)

def K_1(J):
    return sqrt(J/2)*Pauli(3)

def TP(a,b):
    return TensorProduct(a,b)

def proj(psi):
    z = Dagger(psi)
    return psi*z

# função pra obter o estado evoluído
def RHO_t_NM(state,J):
    tp1 = TP(K_0(J),K_1(J))
    tp2 = TP(K_1(J),K_0(J))
    return tp1*proj(state)*tp1.T + tp2*proj(state)*tp2.T

# Função para calcular o emaranhamento
def calculate_entanglement(rho):
    rho_sqrt = rho.applyfunc(sympify)  # Convert all matrix elements to sympy expressions
    eigenvalues = rho_sqrt.eigenvals()  # Calculate eigenvalues using SymPy's eigenvals method

    eigenvalues_real = [val for val in eigenvalues if val.is_real]
    eigenvalues_complex = [val for val in eigenvalues if not val.is_real]
    
    def custom_max(iterable):
        max_val = None
        for val in iterable:
            if max_val is None or val > max_val:
                max_val = val
        return max_val
    
    max_real = custom_max(eigenvalues_real)
    max_complex = custom_max(eigenvalues_complex)

    entanglement = max(0, max_real - sum([sqrt(val) for val in eigenvalues_complex]))
    return entanglement

def get_list_p_noMarkov(list_p, type):
    lamb = 0.01
    gamma_0 = 1
    list_p_noMarkov = []
    if type == 'Bellomo':
        def non_markov_list_p(lamb,gamma_0,t):
            d = sqrt(2*gamma_0*lamb-lamb**2)
            result = exp(-lamb*t)*(cos(d*t/2)+(lamb/d)*sin(d*t/2))**2
            return result
    if type == 'Ana':
        def non_markov_list_p(lamb,gamma_0,t):
            result = 1-exp(-lamb*t)*(cos(t/2)+(lamb)*sin(t/2))
            return result
    for p in list_p:
        list_p_noMarkov.append(non_markov_list_p(lamb,gamma_0,p))
    return list_p_noMarkov

#list_p = np.linspace(0.1,100,10)
T = np.linspace(0.01,100,20)
t_A = get_list_p_noMarkov(T, 'Ana')
# t_B = get_list_p_noMarkov(T, 'Bellomo')

#T = np.logspace(0.01, 100, 100) 
t_A = T
# t_B = T

state = werner_state(-0.8,-0.8,-0.8)
print(RHO_t_NM(state, 14))
print(type(RHO_t_NM(state, 14)))
print(calculate_entanglement(RHO_t_NM(state, 14)))
# y1 = [coh_l1(RHO_t_NM(state, i)) for i in t_A]

y1 = [calculate_entanglement(RHO_t_NM(state, i)) for i in t_A]
# y2 = [coh_l1(RHO_t_NM(state, i)) for i in t_B]

y3 = [concurrence(RHO_t_NM(state, i)) for i in t_A]
#y4 = [concurrence(RHO_t_NM(state, i)) for i in t_B]

#T = [ np.log(i) for i in t_A]
plt.plot(t_A,y1,label='coh_l1 - Ana')
# plt.plot(t_B,y2,label='coh_l1 - Bellomo')
plt.plot(t_A,y3,label='concurrence - Ana')
# plt.plot(t_B,y4,label='concurrence - Bellomo')
plt.ylabel('concurrence')
plt.xscale('log')
plt.xlabel('log(t)')

# plt.xlim(0.01, 200)
plt.grid(True)
plt.legend()
plt.show()