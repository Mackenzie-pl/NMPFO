from math import log
from numpy import arange, exp, asarray
from scipy.linalg import solve
from scipy.sparse import diags

def calcdx(S,E,n):
    return log(S/E)/n

def calcdtau(dx,vol,T,n):
    dtau=1000
    a=1
    while dtau>=(pow(dx,2)/2):
        dtau=(pow(vol,2)*T)/(2*n*a)
        a=a+1
    return dtau

def calca(dtau,dx):
    return dtau/pow(dx,2)

def calck(r,vol):
    return 2*r/pow(vol,2)

def ui0(k,dx,i,otype):
    if otype==0:
        return max([exp(0.5*(k+1)*i)-exp(0.5*(k-1)*i),0])
    elif otype==1:
        return max([exp(0.5*(k-1)*i)-exp(0.5*(k+1)*i),0])

def uXT(S,E,k,j,dtau,otype,dx,ostyle):
    if ostyle==0:
        if otype==0:
            return exp(((k+1)*(log(S/E)+100*dx)/2)+((pow((k+1)/2,2)*j*dtau)))
        elif otype==1:
            return 0
    elif ostyle==1:
        if otype==0:
            return exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)+100*dx))-exp(0.5*(k-1)*(log(S/E)+100*dx)),0])
        elif otype==1:
            return exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k-1)*(log(S/E)+100*dx))-exp(0.5*(k+1)*(log(S/E)+100*dx)),0])

def u0T(S,E,dx,k,j,dtau,otype,ostyle):
    if ostyle==0:
        if otype==0:
            return 0
        elif otype==1:
            return exp(-0.5*((k-1)*(log(S/E)+100*dx))+(pow((k+1)/2,2)-k)*j*dtau)
    elif ostyle==1:
        if otype==0:
            return exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-100*dx))-exp(0.5*(k-1)*(log(S/E)-100*dx)),0])
        elif otype==1:
            return exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k-1)*(log(S/E)-100*dx))-exp(0.5*(k+1)*(log(S/E)-100*dx)),0])

def uijnew(Sj,a,i,S,E,dx,otype,ostyle,k):
    if ostyle==0:
        return a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i]
    elif ostyle==1:
        if otype==0:
            return max([a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
        elif otype==1:
            return max([a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k-1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k+1)*(log(S/E)-(100-i)*dx)),0])])

def calc_S_0(S,E,k,dx):
    S0=[]
    for i in arange(log(S/E)-100*dx,log(S/E)+100*dx,dx):
        S0.append(ui0(k,dx,i,otype))
    return S0

def calc_S_jnew(Sj,a,S,E,dx,j,dtau,k,otype,ostyle):
    Sjnew=[u0T(S,E,dx,k,j,dtau,otype,ostyle)]
    for i in range(1,199):
        Sjnew.append(uijnew(Sj,a,i,S,E,dx,otype,ostyle,k))
    Sjnew.append(uXT(S,E,k,j,dtau,otype,dx,ostyle))
    return Sjnew

while True:
    try:
        ostyle=int(input("If the option is European, enter 0. If it is American, enter 1: "))
    except ValueError:
        print("Option style must be described as an integer.")
        continue
    else:
        if ostyle==0 or ostyle==1:
            break
        print("Option style must be 0 or 1.")

while True:
    try:
        otype=int(input("If the option is a call, enter 0. If it is a put, enter 1: "))
    except ValueError:
        print("Option type must be described as an integer.")
        continue
    else:
        if otype==0 or otype==1:
            break
        print("Option type must be 0 or 1.")

while True:
    try:
        n=int(input("Enter the number of steps as an integer: "))
    except ValueError:
        print("Number of steps must be an integer number.")
        continue
    else:
        if n>0:
            break
        print("Number of steps cannot be less than or equal to zero.")

while True:
    try:
        r=float(input("Please enter the risk-free rate as a decimal: "))
    except ValueError:
        print("Risk free-rate must be a decimal number.")
        continue
    else:
        break

while True:
    try:
        vol=float(input("Please enter the implied volatility as a decimal: "))
    except ValueError:
        print("Implied volatility must be a decimal.")
        continue
    else:
        if vol>=0:
            break
        print("Implied volatility must be greater than or equal to 0.")

while True:
    try:
        S=float(input("Please enter the starting stock value as a decimal: "))
    except ValueError:
        print("Starting stock value must be a decimal.")
        continue
    else:
        break

while True:
        try:
            E=float(input("Please enter the options strike price as a decimal: "))
        except ValueError:
            print("Strike price must be a decimal.")
            continue
        else:
            break

while True:
    try:
        T=float(input("Please enter the length of the option in years: "))
    except ValueError:
        print("Length of the option must be a decimal.")
        continue
    else:
        if T>0:
            break
        print("Length of the option must be positive.")

dx=calcdx(S,E,n)
dtau=calcdtau(dx,vol,T,n)
k=calck(r,vol)
a=calca(dtau,dx)
ldiag=[]
odiag=[]
u=[]

for i in arange(log(S/E)-100*dx,log(S/E)+100*dx,dx):
    u.append(ui0(k,dx,i,otype))

for j in range(0,200):
    ldiag.append(1+(2*a))
for j in range(0,199):
    odiag.append(-a)

diagonals=[odiag,ldiag,odiag]
A=diags(diagonals,[-1,0,1]).toarray()
A[0,1]=1
A[len(A)-2,len(A)-1]=1

for j in range(0,int((pow(vol,2)*T)/(2*dtau))):
    b=[u0T(S,E,dx,k,j,dtau,otype,ostyle),u[1]+a*u0T(S,E,dx,k,j,dtau,otype,ostyle)]
    for i in range(2,198):
        b.append(u[i])
    b.append(u[len(u)-2]+a*uXT(S,E,k,j,dtau,otype,dx,ostyle))
    b.append(uXT(S,E,k,j,dtau,otype,dx,ostyle))
    b=asarray(b)
    Iu=solve(A,b)
    if ostyle==1:
        for i in range(1,199):
            if otype==0:
                Iu[i]=max([Iu[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
            elif otype==1:
                Iu[i]=max([Iu[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
    Eu=calc_S_jnew(u,a,S,E,dx,j,dtau,k,otype,ostyle)
    for i in range(0,200):
        u[i]=0.5*(Eu[i]+Iu[i])
    

ux=u[int(len(u)/2)]
C=pow(E,0.5*(1+k))*pow(S,0.5*(1-k))*exp((-pow(k+1,2)*pow(vol,2)*T)/8)*ux
print("Option Priced at: ",C)
