import numpy as np

# Ações de controle: 1 = inclui o k-ésimo item; 0 = não inclui
peso = np.array([12, 1, 6])
valor = np.array([24, 2, 6])
peso_max = 15
dx = np.gcd.reduce(peso)

Gx = np.arange(0, peso_max+1, dx)  # Grid dos pesos
T = len(peso) + 1  # Horizonte
N = len(Gx)  # Cardinalidade do espaço de estado

# Criação da matriz V e U
V = np.zeros((T, N), dtype = int)
U = np.zeros((T, N), dtype = int)

# Cálculo dos custos terminais
V[T - 1, -1] = 1000

# Cálculo dos custos c(k, x, u)
c = np.zeros((T - 1, N, 2))
for k in range(T - 1):
    for i in range(N):
        for u in range(2):
            c[k, i, u] = -valor[k] * (u == 1)

# Função f da dinâmica do sistema: f(k, x, u) = x + peso[k] * (u == 1)
f = np.zeros((T - 1, N, 2), dtype=int)
for k in range(T - 1):
    for i in range(N):
        for u in range(2):
            prox_peso = Gx[i] + peso[k] * (u == 1)
            prox_peso = min(prox_peso, Gx[-1])
            f[k, i, u] = np.where(prox_peso == Gx)[0][0]

# Loop principal
for k in range(T - 2, -1, -1):
    for i in range(N):
        Vaux = np.zeros(2)
        for u in range(2):
            Vaux[u] = c[k, i, u] + V[k + 1, f[k, i, u]]
            V[k, i] = min(Vaux)
            U[k, i] = np.argmin(Vaux)

# Simulação do sistema a partir da condição inicial x(0) = 0
i = np.zeros(T, dtype=int)
i[0] = 0
uot = np.zeros(T - 1, dtype=int)
for k in range(T - 1):
    i[k + 1] = f[k, i[k], U[k, i[k]]]
    uot[k] = U[k, i[k]]


i+=1

print(f)
print("Estado ótimo:", i)
print("Decisões ótimas:", uot)