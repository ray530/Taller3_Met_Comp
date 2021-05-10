import numpy as np
import matplotlib.pyplot as plt



def solucion(grados):
	M=1.5
	r1=0.001 
	r2=0.002 
	e=0.8
	k=0.5
	u1=2.0
	theta=np.radians(grados)
	A=np.sin(theta)
	B=np.cos(theta)
	eq10=[A,B,-A,-B,1,1]
	eq20=[0,M,0,1,0,0]
	eq30=[M,0,1,0,0,0]
	eq40=[-B,A,B,-A,0,0]
	eq50=[A,B,0,0,-k,0]
	eq60=[0,0,A,B,0,k]

	P_t=np.array([eq10,eq20,eq30,eq40,eq50,eq60])
	y=np.array([0,0,M*u1,e*u1*B,u1*A,0])
	P=np.transpose(P_t)
	a=np.matmul(P_t,P)
	b=np.matmul(P_t,y)
	V=np.matmul(np.linalg.inv(a),b)
	
	return V

v=[]
theta=np.linspace(0.0,90,8)
for i in range(len(theta)):
	v.append(solucion(theta[i]))

sol=np.array(v[:][:])
a=np.transpose(sol)

plt.figure(figsize=(10,8))
plt.title("Variacion de la colision cuasi-elastica respecto al angulo")
plt.plot(theta,a[0][:],label=r'$V_{1x}$')
plt.plot(theta,a[1][:],label=r'$V_{1y}$')
plt.plot(theta,a[2][:],label=r'$V_{2x}$')
plt.plot(theta,a[3][:],label=r'$V_{2y}$')
plt.plot(theta,a[4][:],label=r'$r_1\omega_1$')
plt.plot(theta,a[5][:],label=r'$r_2\omega_2$')
plt.xlabel("Angulo °")
plt.legend()
plt.show()


