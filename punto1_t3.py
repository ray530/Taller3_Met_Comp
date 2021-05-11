import numpy as np
import matplotlib.pyplot as plt

r1=0.001 
r2=0.002 

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

vsol=np.vectorize(solucion)

theta=np.linspace(0.0,90,10)
v=[]

for i in range(len(theta)):
	v.append(vsol(theta[i]))
print(v)



sol=np.vstack((v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9]))
print(sol)

np.savetxt('dat_solucionp4.dat',sol)


v1x, v1y, v2x, v2y, w1, w2 = np.loadtxt('dat_solucionp4.dat', unpack=True)
#print(w1)

plt.figure(figsize=(10,8))
plt.plot(theta,v1x,label=r'$v_{1x}$')
plt.plot(theta,v1y,label=r'$v_{1y}$')
plt.plot(theta,v2x,label=r'$v_{2x}$')
plt.plot(theta,v2y,label=r'$v_{2y}$')
plt.plot(theta,w1*r1,label=r'$r_1 w_{1}$')
plt.plot(theta,w2*r2,label=r'$r_2 w_{2}$')
plt.xlabel("Grados")
plt.grid()
plt.legend()
plt.show()


