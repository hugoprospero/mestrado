#Nota: O método das potências inversas não está convergindo para autovalores complexos, não consegui descobrir o porquê disso, se for possível um feedback eu agradeço (se eu inicializo com um 
# vetor complexo aleatório, o método também deixa de convergir pra autovalores reais)

import numpy as np
np.set_printoptions(precision=4)

#Houselholder

def householder(A,m,n):
    Qt = np.eye(m)

    for k in range(m-2):
        x = A[k+1:m,k]
        e1 = np.zeros(n-(k+1))
        e1[0] = 1

        vk = np.linalg.norm(x)*np.sign(x[0])*e1 + x
        vk = vk/np.linalg.norm(vk)

        Qt = np.eye(m)
        Qt[k+1:m,k+1:m] = Qt[k+1:m,k+1:m] - 2*np.outer(vk,vk)
        
        A = Qt@A@np.transpose(Qt)
    return A,Qt

#Transformação de uma matriz para a forma de Hessenberg

def Hessenberg (A, m, n):
    if m != n:
        print("A matriz dada não é quadrada")
    else:
        H , Qt= householder(A,m,n)
        H[np.abs(H) < 1e-6] = 0
        print("Matriz na forma Hessenberg: \n", H)
    return H

#fatoração QR

def QR(H,m,n,maxit = 100):
    tol = 1e-6
    dec = 0
    it = 0

    while dec == 0 and it<maxit:
        Q, R = np.linalg.qr(H)
        H = R@Q
        dec = 1
        for i in range(m):
            for j in range(n):
                if i>j and np.abs(H[i,j]) > tol:
                    dec = 0
        it = it + 1
    H[np.abs(H) < 1e-4] = 0

    print("Número de iterações: ", it)
    print("Matriz após transformação QR: \n", H)
    return H
        

#fatoração QR double shift

def QRds(H,m,n,maxit = 101):
    tol = 1e-6
    dec = 0
    it = 0
    while dec == 0 and it<maxit:
        K=H[n-2:n,n-2:n]
        B = H[n-2:n,n-2:n]
        trB = B[0,0]+B[1,1]
        detB = B[0,0]*B[1,1] - B[0,1]*B[1,0]
        
        s1 = (trB + np.sqrt((trB)**2-4*detB+0j))/2
        s2 = (trB - np.sqrt((trB)**2-4*detB+0j))/2

        shift = H-s1*np.eye(n)
        Q, R = np.linalg.qr(shift)
        H = R@Q+s1*np.eye(n)

        it+=1
        if np.linalg.norm(K-H[n-2:n,n-2:n]) < tol:
            break

        #K=H[n-2:n,n-2:n]
        shift = H-s2*np.eye(n)
        Q, R = np.linalg.qr(shift)
        H = R@Q+s2*np.eye(n)      
       
        it+=1
        if np.linalg.norm(K-H[n-2:n,n-2:n]) < tol:
            break

    H.real[np.abs(H.real) < 1e-4] = 0
    H.imag[np.isclose(H.imag, 0)] = 0

    print("Número de iterações: ", it)
    print("Matriz após transformação QR com double shift: \n", H)
    return H

#Método das potências inversas

def invpower (A,m,n, lam, tol=1e-6,maxit=1000):
    v = np.ones(m, dtype = complex)
    I = np.eye(n, dtype = complex)
    con=1000
    it = 0

    while con >tol and it < maxit:
        vk = v.copy()
        v=np.linalg.inv(A-(lam+0.1)*I)@v
        v = v/np.linalg.norm(v)
        con = np.linalg.norm(np.sign(v)*v-np.sign(vk)*vk)
        it+=1
    v[np.abs(v) < 1e-6] = 0
    print("Número de iterações: ",it)
    return v

#1) Matriz com autovalores reais e distintos

print("1) Matriz com autovalores reais e distintos: ")

A = np.array([[6, -25, 2],
            [3, 86, -34],
            [-2, 12, -10]])

print(A)

m, n = np.shape(A)

#H = Hessenberg(A,m,n)
H = QR(A,m,n)
lam = np.diag(H)

print("Autovalores: ", lam)

x=np.zeros((m,m), dtype = complex)

for k in range(m):
    x[:,k] = invpower (A,m,n, lam[k])

print("Autovetores unitários associados a cada autovalor (colunas): \n", x)

#2) Matriz com autovalores todos de módulo 1

print("\n\n\n 2) Matriz com autovalores todos de módulo 1: ")

theta1 = np.pi / 3
theta2 = np.pi / 4

R1 = np.array([
    [np.cos(theta1), -np.sin(theta1)], 
    [np.sin(theta1),  np.cos(theta1)]
], dtype = complex)

R2 = np.array([
    [np.cos(theta2), -np.sin(theta2)], 
    [np.sin(theta2),  np.cos(theta2)]
], dtype = complex)

# Construir a matriz 4x4
R = np.block([
    [R1, np.zeros((2, 2), dtype = complex)],
    [np.zeros((2, 2), dtype = complex), R2]
])

print(R)

m, n = np.shape(R)

H = QRds(R,m,n) #A matriz resultante ja está na forma de Hessenberg


H1 = H[0:2,0:2]
H2 = H[n-2:n,n-2:n]

del1 = np.trace(H1)**2 - 4*np.linalg.det(H1)
del2 = np.trace(H2)**2 - 4*np.linalg.det(H2)

a1 = np.trace(H1)/2
a2 = np.trace(H2)/2
b1 = np.sqrt(np.sign(del1)*del1)/2
b2 = np.sqrt(np.sign(del2)*del2)/2

lam = np.array([a1+b1*1j,a1-b1*1j,a2+b2*1j,a2-b2*1j])

print("Autovalores: ", lam)

x=np.zeros((m,m), dtype = complex)

for k in range(m):
    x[:,k] = invpower (R,m,n, lam[k])

print("Autovetores unitários associados a cada autovalor (colunas): \n", x)

#3) Matriz real com dois pares de autovalores distintos

print("\n\n\n 3)Matriz real com dois pares de autovalores distintos")

A = np.array([[0, -1, 2, 1],
              [1, 0, 1, 2],
              [-2, -1, 0, 1],
              [-1, -2, -1, 0]])

print(A)

m, n = np.shape(A)

H = Hessenberg(A,m,n)
H = QRds(H,m,n)

H1 = H[0:2,0:2]
H2 = H[n-2:n,n-2:n]

del1 = np.trace(H1)**2 - 4*np.linalg.det(H1)
del2 = np.trace(H2)**2 - 4*np.linalg.det(H2)

a1 = np.trace(H1)/2
a2 = np.trace(H2)/2
b1 = np.sqrt(np.sign(del1)*del1)/2
b2 = np.sqrt(np.sign(del2)*del2)/2

lam = np.array([a1+b1*1j,a1-b1*1j,a2+b2*1j,a2-b2*1j])

print("Autovalores: ", lam)

x=np.zeros((m,m), dtype = complex)

for k in range(m):
    x[:,k] = invpower (A,m,n, lam[k])

print("Autovetores unitários associados a cada autovalor (colunas): \n", x)

#4) Matriz simétrica

print("\n\n\n 4) Matriz simétrica: ")

A = np.array([[2,1,0,1],
              [1,3,1,0],
              [0,1,4,1],
              [1,0,1,5]])

print(A)

m, n = np.shape(A)

H = Hessenberg(A,m,n)
H = QR(H,m,n)
lam = np.diag(H)

print("Autovalores: ", lam)

x=np.zeros((m,m), dtype = complex)

for k in range(m):
    x[:,k] = invpower (A,m,n, lam[k])


print("Autovetores unitários associados a cada autovalor (colunas): \n", x)