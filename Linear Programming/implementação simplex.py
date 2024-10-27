import numpy as np
import pandas as pd
from scipy.linalg import expm
import seaborn as sns
import matplotlib.pyplot as plt
import random

# Digitar as matrizes na forma padrão

cT = np.array ([-3, -2,0,0])
A = np.array ([[2,1,1,0],
              [1,2,0,1]])
bT= np.array ([10,12])
x= np.zeros([len(cT)])

# Selecionar uma base inicial

B = np.zeros((len(A),len(A)))


while np.linalg.matrix_rank(B) < len(A):
        v = np.random.choice(range(len(A[0])), size=len(A), replace=False)
        B = A[:, v]

# simplex

epsilon = np.zeros([len(bT)])
crT = np.zeros([len(cT)])

maxit = 1000
it = 0

while it<=maxit:
    crT = np.zeros([len(cT)])
    cBT = cT[v]
    B = A[:,v]
    # passo 1
    xBT = np.linalg.inv(B)@np.transpose(bT)
    #passo 2
    lambdaT = cBT@np.linalg.inv(B)
    crT = cT - lambdaT@A
    if np.all(crT>=0):
        if np.any(xBT<0): #se encontrar uma solução infactivel, encontrar uma nova base inicial
            B = np.zeros((len(A),len(A)))
            while np.linalg.matrix_rank(B) < len(A):
                v = np.random.choice(range(len(A[0])), size=len(A), replace=False)
                B = A[:, v]
        else:
            for m in range(len(v)):
                x[v[m]] = xBT[m]
            valot = xBT@np.transpose(cBT)
            print("Solução ótima: ", x)
            print("Valor mínimo: ", valot)
            it = maxit
    else:
        k = np.argmin(crT)
        y = np.linalg.inv(B)@A[:,k]
        if np.all(y<=0):
            print("O problema não tem solução finita")
            it = maxit
        else:
            for i in range(len(xBT)):
                if y[i]>0:
                    epsilon[i] = xBT[i]/y[i]
                else:
                    epsilon[i] = 1000
            l = v[np.argmin(epsilon)]
            for j in range(len(v)):
                if v[j] == l:
                    v[j] = k
    it+=1
