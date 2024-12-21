import numpy as np
import math

u1 = np.linspace(-100, 100, num=((100 - (-100)) // 20) + 1, dtype=int)
u2 = np.linspace(-100, 100, num=((100 - (-100)) // 20) + 1, dtype=int)
V = np.zeros((101,101,3,3,7))
a = np.zeros((101,101,3,3,6))
aa = np.zeros((101,101,3,3,6))
c = np.zeros((101,101,3,3,len(u1),len(u2)))


Q = np.zeros((3,3))

Q[0,0] = 0.2
Q[0,1] = 0.3
Q[0,2] = 0.5

Q[1,0] = 0.15
Q[1,1] = 0.45
Q[1,2] = 0.4

Q[2,0] = 0.35
Q[2,1] = 0.5
Q[2,2] = 0.15

P = np.zeros((3,3))

P[0,0] = 0.05
P[0,1] = 0.75
P[0,2] = 0.2

P[1,0] = 0.15
P[1,1] = 0.7
P[1,2] = 0.15

P[2,0] = 0.15
P[2,1] = 0.8
P[2,2] = 0.05

dt = np.zeros ((2,3))

dt[0,0] = -0.05
dt[0,1] = 0
dt[0,2] = 0.06
dt[1,0] = -0.08
dt[1,1] = 0
dt[1,2] = 0.14

#print("teste")

for i in range(101):
    for j in range(101):
        for k in range(3):
            for l in range(3):
                if (i+k)>100:
                    V[i,j,k,l,6] = -1000

#print("teste")



for w in range(5,-1,-1):
    for i in range(101):
        for j in range(101):
            for k in range(3):
                for l in range(3):
                    for m in range(len(u1)):
                        for n in range(len(u2)):
                            if not (0 <= (i+u1[m]) < 101 and 0 <= (j+u2[n]) < 101 and (i+j+u1[m]+u2[n])<101):
                                        continue
                            else:
                                c[i,j,k,l,m,n] = (i+u1[m])*(Q[k,2]*dt[0,2]+Q[k,0]*dt[0,0])+(j+u2[n])*(P[k,2]*dt[1,2]+P[k,0]*dt[1,0])+(100-(i+u1[m])-(j+u2[n]))*0.01
                                J = c[i,j,k,l,m,n]
                                for o in range(3):
                                    for p in range(3):
                                        g1 = (i+u1[m])*(1+dt[0,o])
                                        g2 = (j+u2[n])*(1+dt[1,p])
                                        gs = (100-(i+u1[m]+j+u2[n]))*1.01
                                        gt = g1+g2+gs
                                        f1 = round(100*g1/gt)
                                        f2 = round(100*g2/gt)
                                        #print (i+u1[m])
                                        #print ("e", j+u2[n])
                                        #print(g1)
                                        #print(g2)
                                        #print(gs)
                                        #print(gt)
                                        J = J + V[f1,f2,o,p,w+1]*Q[k,o]*P[l,p]
                            if J > V[i,j,k,l,w]:
                                V[i,j,k,l,w] = J
                                a[i,j,k,l,w] = u1[m]
                                aa[i,j,k,l,w] = u2[n]    
        #print(i)
    #print(w)

c=50
d=50
theta1 = [1,1,2,0,0,2,0]
theta2 = [1,2,1,1,0,0,0]

for w in range(6):
    print(V[c,d,theta1[w],theta2[w],w])
    print ("Decisão ótima do mês ",w+1,": ", a[c,d,theta1[w],theta2[w],w], "de ações A e ", aa[c,d,theta1[w],theta2[w],w], "de ações B")
    g1 = (c+a[c,d,theta1[w],theta2[w],w])*(1+dt[0,theta1[w+1]])
    g2 = (d+aa[c,d,theta1[w],theta2[w],w])*(1+dt[1,theta2[w+1]])
    gt = 100+(c+a[c,d,theta1[w],theta2[w],w])*dt[0,theta1[w+1]]+(d+aa[c,d,theta1[w],theta2[w],w])*dt[1,theta2[w+1]]+(100-(c+a[c,d,theta1[w],theta2[w],w])-(d+aa[c,d,theta1[w],theta2[w],w]))*0.01
    c = round(100*g1/gt)
    d = round(100*g2/gt)