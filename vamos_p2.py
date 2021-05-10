
import astropy.constants as const
import astropy.units as units
import matplotlib.pyplot as plt

import numpy as np

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

# Creation of the time grid (in years)
t_0 = 0.
t_f = 2.5

# Number of steps in the grid
n = 1500000

# Constant stepsize defined by the number of steps in the grid
h = (t_f - t_0)/n

# Arrays to store the solution
t = np.linspace(t_0, t_f, n) # Time information
QR = np.zeros([4,n]) # RK4's Method information
cond=np.zeros([2,n])

# Initial Conditions
QR[0,0] = r0[0]
QR[1,0] = 0.
QR[2,0] = 0.
QR[3,0] = v0[-1]

for i in range(1,n):
    q0 = QR[:,i-1]
    qf = RK4(ODE, h, 0, q0)
    #cond=cond(qf,h,t,q0)
    QR[:,i] = qf[:]

'''
def omega(a,e,o,t):
	d=a*(1-e)*(1+1/e)
	c=(e**2*d**2)/(np.sqrt(1-e**2)*(e**2-1))
	g=(e**2*d**2)/(e**2-1)
	return 2*np.atan(((e-1)/(np.sqrt(1-e**2)))*np.tan(o/2))*c + g*np.sin(o)/(1+e*np.cos(o))-t
from scipy import optimize

root = optimize.newton(omega(a,e,0,t), 0.)
print(root)
'''


plt.figure(figsize=(10,6))
plt.plot(QR[0,:], QR[1,:], color='cornflowerblue', label=f'$h=$ {h:.2e}')
plt.title('RK4 Method')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend()


from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')
ax.plot(QR[0,:], QR[1,:], t, label='parametric curve')
ax.legend()
plt.show()