import numpy as np
# acoes de controle: 1=cima, 2=dir, 3=baixo, 4=esq., 5=nada

# custos c(i,u):
c = np.full((16,5),np.inf)
c[0][0]=np.inf; c[0][1]=10;      c[0][2]=10;      c[0][3]=np.inf;   c[0][4]=0
c[1][0]=np.inf; c[1][1]=10;      c[1][2]=12;      c[1][3]=10;       c[1][4]=0
c[2][0]=np.inf; c[2][1]=8;       c[2][2]=10;      c[2][3]=10;       c[2][4]=0
c[3][0]=np.inf; c[3][1]=np.inf;  c[3][2]=22;      c[3][3]=8;        c[3][4]=0
c[4][0]=10;     c[4][1]=15;      c[4][2]=5;       c[4][3]=np.inf;   c[4][4]=0
c[5][0]=12;     c[5][1]=8;       c[5][2]=2;       c[5][3]=15;       c[5][4]=0
c[6][0]=10;     c[6][1]=15;      c[6][2]=8;       c[6][3]=8;        c[6][4]=0
c[7][0]=22;     c[7][1]=np.inf;  c[7][2]=10;      c[7][3]=15;       c[7][4]=0
c[8][0]=5;      c[8][1]=12;      c[8][2]=12;      c[8][3]=np.inf;   c[8][4]=0
c[9][0]=2;      c[9][1]=30;      c[9][2]=10;      c[9][3]=12;       c[9][4]=0
c[10][0]=8;     c[10][1]=22;     c[10][2]=15;     c[10][3]=30;      c[10][4]=0
c[11][0]=10;    c[11][1]=np.inf; c[11][2]=10;     c[11][3]=22;      c[11][4]=0
c[12][0]=12;    c[12][1]=5;      c[12][2]=np.inf; c[12][3]=np.inf;  c[12][4]=0
c[13][0]=10;    c[13][1]=8;      c[13][2]=np.inf; c[13][3]=5;       c[13][4]=0
c[14][0]=15;    c[14][1]=10;     c[14][2]=np.inf; c[14][3]=8;       c[14][4]=0
c[15][0]=10;    c[15][1]=np.inf; c[15][2]=np.inf; c[15][3]=10;      c[15][4]=0



# função f da dinamica do sistema:
f = np.zeros((16,5), dtype = int)

for i in range(16):
    f[i][4] = i
    f[i][1] = i+1
    f[i][3] = i-1
    f[i][0] = i-4
    f[i][2] = i+4

f[3][1]=3;   f[7][1]=7;   f[11][1]=11; f[15][1]=15
f[0][3]=0;   f[4][3]=4;   f[8][3]=8;   f[12][3]=12
f[0][0]=0;   f[1][0]=1;   f[2][0]=2;   f[3][0]=3
f[12][2]=12; f[13][2]=13; f[14][2]=14; f[15][2]=15

# custo terminal:

V=np.zeros((11,16))

for i in range(16):
    V[10][i] = 1000

V[10][11]=0


# loop principal

J=np.zeros((16,5))
U=np.zeros((10,16), dtype = int)

for k in range(9,-1,-1):
    for i in range(16):
        for u in range(5):
            J[i][u] = c[i][u] + V[k+1][f[i][u]]
        V[k][i] = np.min(J[i])
        U[k][i] = np.argmin(J[i])

# simulando o sistema a partir de cada condição inicial

x=np.zeros((16,11), dtype = int)

for i in range(16):
    x[i][0]=i

for i in range(16):
    for k in range (10):
        x[i][k+1] = f[x[i][k]][U[k][x[i][k]]]


x = x + 1

print(x)
print(V)
