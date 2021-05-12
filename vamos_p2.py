import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import integrate 

#G = const.G.to('au^3/(M_sun year^2)')
G=4*np.pi**2
print(G)

M=1.
alfa=1.1*10**(-8)

def ODE(t,q0):

	x = q0[0]
	y = q0[1]
	u = q0[2] 
	v = q0[3]
	f = np.zeros(4)
	f[0] = u 
	f[1] = v 
	f[2] = -( G*M*x/(x**2 + y**2)**(3/2) * (1 + alfa/(x**2+y**2)**1))
	f[3] = - ( G*M*y/(x**2 + y**2)**(3/2) * (1 + alfa/(x**2+y**2)**1))
	return f

def cond(ODE,h,t,q0):
	x = q0[0]
	y = q0[1]
	u = q0[2] 
	v = q0[3]
	r = np.zeros(2)
	r[0] = x+ODE(t,q0)[0]*h+((h**2)/2)*ODE(t,q0)[2]
	r[1] = y+ODE(t,q0)[1]*h+((h**2)/2)*ODE(t,q0)[3] 
	return r

def RK4(ODE, h, t0, q0):

	k11 = h*ODE(t0, q0)
	
	k21 = h*ODE(t0 + h/2, q0 + k11/2)

	k31 = h*ODE(t0 + h/2, q0 + k21/2)

	k41 = h*ODE(t0 + h, q0 + k31)

	q1 = q0 + (k11 + 2*k21 + 2*k31 + k41)/6

	return q1

a=0.387
e=0.20563

r0 = [a*(1+e), 0.0]
v0 = [0.0, np.sqrt((G*1/a) * ((1-e)/(1+e)))]

# tiempo
t_0 = 0.
t_f = 2.5 #tiempo de 10 orbitas, 88 dias por orbita. 

# no. de pasos
n = 150000

# intervalo 
h = (t_f - t_0)/n #parecido al orden de alfa

# inicializar vectores 
t = np.linspace(t_0, t_f, n) 
QR = np.zeros([4,n]) 
cond=np.zeros([2,n])

# Condiciones iniciales
QR[0,0] = r0[0]
QR[1,0] = 0.
QR[2,0] = 0.
QR[3,0] = v0[-1]

for i in range(1,n):
    q0 = QR[:,i-1]
    qf = RK4(ODE, h, 0, q0)
    #cond=cond(qf,h,t,q0)
    QR[:,i] = qf[:]

wy=QR[3,:]/np.sqrt(QR[0,:]**2+QR[1,:]**2)
theta_y=integrate.cumulative_trapezoid(wy,t)


plt.figure(figsize=(10,6))
plt.plot(QR[0,:], QR[1,:], color='cornflowerblue', label=f'$h=$ {h:.2e}')
plt.title('RK4 Method')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend()

plt.figure(figsize=(10,6))
plt.plot(t[:-1], -theta_y/2062.6481 , color='cornflowerblue', label=f'$h=$ {h:.2e}')
plt.title('variacion angular')
plt.xlabel(r'$t$')
plt.ylabel(r'$\theta_y \ ar/siglo$')
plt.legend()

fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')
ax.plot(QR[0,:], QR[1,:], t, label='orbita de Mercurio')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$t$')
ax.legend()

'''
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')
ax.plot(QR[2,:], QR[3,:], t, label='velocidades')
ax.set_xlabel(r'$V_x$')
ax.set_ylabel(r'$V_y$')
ax.set_zlabel(r'$t$')
ax.legend()
'''
plt.show()
