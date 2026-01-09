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
    "Integra por el metodo de Simpson, Es una función de 4 entradas (a,b,N,f) donde a y b es el intervalo de integración, N es el número de partición de mi integral y f es la función a integra"
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

    