import scipy as sp
import numpy as np
def TrapecioInt(a,b,N,f):
    "Integra por el método del trapecio, Es una función de 4 entradas (a,b,N,f) donde a y b es el intervalo de integración, N es el número de partición de mi integral y f es la función a integrar:"
    "Nota: definir tu función aparte"
    
    h = abs(b-a)/N
    I = (f(a)+f(b))/2
    for i in range (1,N):
        I = I + f(a + i*h)
    I = h*I   
    return I

def Simpson1(a,b,N,f):
    "Integra por el metodo de Simpson, Es una función de 4 entradas (a,b,N,f) donde a y b es el intervalo de integración, N es el número de partición de mi integral y f es la función a integrar: "
    "Nota: definir la función a integrar, usar un valor N par"
    if (N%2 == 0):
        
        h = abs(b-a)/N
        I = f(a)+f(b)+4*f(a+(N-1)*h)
        for i in np.arange(1,N/2,1):
            I = I + 4*f(a+(2*i-1)*h) + 2*f(a+(2*i)*h)
        
        I = (h/3)*I
    else:
        print("Tu número no es divisible entre 2")
    
    return I    

def CuadGauss(a,b,N,f):
    "Integra por el metodo de Cuadratura, Es una función de 4 entradas (a,b,N,f) donde a y b es el intervalo de integración, N es el número de partición de mi integral y f es la función a integra"
    S = 0
    
    x,w = gauss_xw_ab(N,a,b)
    
    for i in range(N):
        S = S + w[i]*f(x[i])
        
    return float(S)

def gauss_xw_ab(N,a,b):
    #import scipy as sp
    "Esta función es necesaria para obtener valores de raices en un intervalo a y b que no sean -1 y 1"
    x,w = sp.special.roots_legendre(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

def dcentral(x,h,f):
    "Calcula la derivada de la función f utilizando el método de derivación central, se recomienda el parametro h <= 0.01 y x puede ser un número o "
    " arreglo(s)"
    Derivada = (f(x+h/2)-f(x-h/2))/h
    return Derivada

def Parcialx(x,y,h,f): # definimos la función de derivada parcial respecto a x
    "Calcula la derivada parcial respecto de x, con el proceso de derivada central"
    px = (f(x+h/2,y)-f(x-h/2,y))/h
    return px

def Parcialy(x,y,h,f): # definimos la función de derivada parcial respecto a y
    "Calcula la derivada parcial respecto de y, con el proceso de derivada central"
    py = (f(x,y+h/2)-f(x,y-h/2))/h
    return py

def Integral2D(N,a1,b1,a2,b2,f):
    "Calcula la integral de una función de 2 dimensiones f(x1,x2) donde pedimos una N que calculara los raices del polinomio de legendre de grado N"
    " a1 y b1 son los limites inferior y superior de x1; a2 y b2 son los limites inferior y superior de x2"
    x1,w1 = gauss_xw_ab(N,a1,b1)
    x2,w2 = gauss_xw_ab(N,a2,b2)
    I = 0
    for i in range(N):
        for j in range(N):
            I = I + w1[i]*w2[j]*f(x1[i],x2[j])
    return I

def Parcialxarreglos(M,i,j,h):
    if i == 0:
        d = (M[j][i+1] - M[j][i])/h
    elif i == len(M[j])-1:
        d = (M[j][i]-M[j][i-1])/h
    else:
        d = (M[j][i+1]-M[j][i-1])/2*h
    return d

def Parcialyarreglos(M,i,j,h):
    if j == 0:
        d = (M[j+1][i] - M[j][i])/h
    elif j == np.shape(M)[0]-1:
        d = (M[j][i]-M[j-1][i])/h
    else:
        d = (M[j+1][i]-M[j-1][i])/2*h
    return d

### Metodos númericos ###

def RungeKutta4_2D(f,x0,t0,tf,h):
    
    T=[]
    X1=[]
    X2=[]
    
    X1.append(x0[0])
    X2.append(x0[1])
    T.append(t0)

    t=t0
    x=x0

    while (t<tf):
        
        k1 = f(x)*h
        k2 = f(x+(k1/2))*h
        k3 = f(x+(k2/2))*h
        k4 = f(x+k3)*h
        
        x = x + (k1+(2*k2)+(2*k3)+k4)/6
        t = t + h/2
        
        X1.append(x[0])
        X2.append(x[1])
        T.append(t)

    return X1,X2,T