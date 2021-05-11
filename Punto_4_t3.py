
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate
from tqdm import tqdm



def f(t):
    #return t
    
    return t**2

def GetFourier(t,T,N):
    
    y = t

    
    a0, error0 = integrate.quad(f, -np.pi, np.pi )
    a0 *= 2./T
    
    y = a0
    
    for n in tqdm(range(1, N)):
        an, erroran = integrate.quad(lambda t : f(t)*np.cos( (2*np.pi*n*t)/T ), -np.pi, np.pi )
        bn, errorbn = integrate.quad(lambda t : f(t)*np.sin( (2*np.pi*n*t)/T ), -np.pi, np.pi )
      
        y += 2*an*np.cos( (2*np.pi*n*t)/T )/T + 2*bn*np.sin( (2*np.pi*n*t)/T )/T 
	

    
    return a0, y-np.pi



t = np.arange(-np.pi,np.pi,0.001)
a0, y = GetFourier(t,2*np.pi,150)

valor=np.trapz(y,t)
print("Valor aprox ")
print(valor/(a0*np.pi)) #por normalización de integral y(t)dt
print("Valor real ")
print(np.pi**6 /945)



fig = plt.figure()
plt.plot(t,f(t),'k',label="F(t)")
plt.plot(t,y, label='Serie de Fourier')
plt.legend(loc=0)
plt.grid()
plt.show()

