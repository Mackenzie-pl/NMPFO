from numpy import arange, exp, asarray
from scipy.linalg import solve
from scipy.sparse import diags
from math import exp, log, sqrt
from tkinter import *
try:
    from statistics import NormalDist
except ImportError:
    print("Please update to python 3.8")
    exit()

def cncalcdx(S,E,n):
    return log(S/E)/n

def cncalcdtau(dx,vol,T,n):
    dtau=1000
    a=1
    while dtau>=(pow(dx,2)/2):
        dtau=(pow(vol,2)*T)/(2*n*a)
        a=a+1
    return dtau

def cncalca(dtau,dx):
    return dtau/pow(dx,2)

def cncalck(r,vol):
    return 2*r/pow(vol,2)

def cnui0(k,dx,i,otype):
    if otype==0:
        return max([exp(0.5*(k+1)*i)-exp(0.5*(k-1)*i),0])
    elif otype==1:
        return max([exp(0.5*(k-1)*i)-exp(0.5*(k+1)*i),0])

def cnuXT(S,E,k,j,dtau,otype,dx,ostyle):
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

def cnu0T(S,E,dx,k,j,dtau,otype,ostyle):
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

def cnuijnew(Sj,a,i,S,E,dx,otype,ostyle,k):
    if ostyle==0:
        return a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i]
    elif ostyle==1:
        if otype==0:
            return max([a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
        elif otype==1:
            return max([a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k-1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k+1)*(log(S/E)-(100-i)*dx)),0])])

def cncalc_S_0(S,E,k,dx):
    S0=[]
    for i in arange(log(S/E)-100*dx,log(S/E)+100*dx,dx):
        S0.append(cnui0(k,dx,i,otype))
    return S0

def cncalc_S_jnew(Sj,a,S,E,dx,j,dtau,k,otype,ostyle):
    Sjnew=[cnu0T(S,E,dx,k,j,dtau,otype,ostyle)]
    for i in range(1,199):
        Sjnew.append(cnuijnew(Sj,a,i,S,E,dx,otype,ostyle,k))
    Sjnew.append(cnuXT(S,E,k,j,dtau,otype,dx,ostyle))
    return Sjnew

def ecalcdx(S,E,n):
    return log(S/E)/n

def ecalcdtau(dx,vol,T,n):
    dtau=1000
    a=1
    while dtau>=(pow(dx,2)/2):
        dtau=(pow(vol,2)*T)/(2*n*a)
        a=a+1
    return dtau

def ecalca(dtau,dx):
    return dtau/pow(dx,2)

def ecalck(r,vol):
    return 2*r/pow(vol,2)

def eui0(k,dx,i,otype):
    if otype==0:
        return max([exp(0.5*(k+1)*i)-exp(0.5*(k-1)*i),0])
    elif otype==1:
        return max([exp(0.5*(k-1)*i)-exp(0.5*(k+1)*i),0])

def euXT(S,E,k,j,dtau,otype,dx,ostyle):
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

def eu0T(S,E,dx,k,j,dtau,otype,ostyle):
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

def euijnew(Sj,a,i,S,E,dx,otype,ostyle,k):
    if ostyle==0:
        return a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i]
    elif ostyle==1:
        if otype==0:
            return max([a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
        elif otype==1:
            return max([a*(Sj[i-1]+Sj[i+1])+(1-2*a)*Sj[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k-1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k+1)*(log(S/E)-(100-i)*dx)),0])])

def ecalc_S_0(S,E,k,dx):
    S0=[]
    for i in arange(log(S/E)-100*dx,log(S/E)+100*dx,dx):
        S0.append(eui0(k,dx,i,otype))
    return S0

def ecalc_S_jnew(Sj,a,S,E,dx,j,dtau,k,otype,ostyle):
    Sjnew=[eu0T(S,E,dx,k,j,dtau,otype,ostyle)]
    for i in range(1,199):
        Sjnew.append(euijnew(Sj,a,i,S,E,dx,otype,ostyle,k))
    Sjnew.append(euXT(S,E,k,j,dtau,otype,dx,ostyle))
    return Sjnew

def icalcdx(S,E,n):
    return log(S/E)/n

def icalcdtau(vol,T,n):
    return (pow(vol,2)*T)/(2*n)

def icalck(r,vol):
    return 2*r/pow(vol,2)

def icalca(dtau,dx):
    return dtau/pow(dx,2)

def iui0(k,dx,i,otype):
    if otype==0:
        return max([exp(0.5*(k+1)*i)-exp(0.5*(k-1)*i),0])
    elif otype==1:
        return max([exp(0.5*(k-1)*i)-exp(0.5*(k+1)*i),0])

def iuXT(S,E,k,j,dtau,otype,dx,ostyle):
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

def iu0T(S,E,dx,k,j,dtau,otype,ostyle):
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

def calcd1(S,E,r,vol,div,T):
    return (log(S/E)+(r-div+0.5*pow(vol,2))*T)/(vol*sqrt(T))

def calcd2(d1,vol,T):
    return d1-(vol*sqrt(T))

def calcV(otype,d1,d2,S,div,T,E,r):
    if otype==0:
        return (S*exp(-div*T)*NormalDist().cdf(d1))-(E*exp(-r*T)*NormalDist().cdf(d2))
    elif otype==1:
        return (E*exp(-r*T)*NormalDist().cdf(-d2))-(S*exp(-div*T)*NormalDist().cdf(-d1))

def calcu(r,div,vol,n,T):
    return exp((r-div)*T/n+vol*sqrt(T/n))

def calcd(r,div,vol,n,T):
    return exp((r-div)*T/n-vol*sqrt(T/n))

def calcpstar(u,d,T,n,r,div):
    return ((exp((r-div)*T/n)-d))/(u-d)

def calc_S_n(u,d,S,n):
    S_n=[]
    for i in range(0,n+1):
        S_n.append(S*pow(u,n-i)*pow(d,i))
    return S_n

def calc_C_n(E,S_n,otype,B,E1,E2):
    C_n=[]
    for S_i in S_n:
        if otype==0:
            C_n.append(max([S_i-E,0]))
        elif otype==1:
            C_n.append(max([E-S_i,0]))
        elif otype==2:
            if S_i>E:
                C_n.append(B)
            elif S_i<E:
                C_n.append(0)
        elif otype==3:
            C_n.append(max([S_i-E1,0])-max([S_i-E2,0]))
        elif otype==4:
            if S_i>E:
                C_n.append(S_i)
            elif S_i<E:
                C_n.append(0)
        elif otype==5:
            C_n.append(max([S_i-E,0])+max([E-S_i,0]))
    return C_n

def calc_next(C_n,pstar,r,T,n,ostyle,otype,u,d,S,E,j):
    C_nnew=[]
    for i in range(0,len(C_n)-1):
        if ostyle==0:
            C_nnew.append(exp(-r*T/n)*(pstar*C_n[i]+(1-pstar)*C_n[i+1]))
        elif ostyle==1:
            S_n=calc_S_n(u,d,S,n-j-1)
            if otype==0:
                C_nnew.append(max([exp(-r*T/n)*(pstar*C_n[i]+(1-pstar)*C_n[i+1]),S_n[i]-E]))
            elif otype==1:
                C_nnew.append(max([exp(-r*T/n)*(pstar*C_n[i]+(1-pstar)*C_n[i+1]),E-S_n[i]]))
            elif otype==2:
                if S_n[i]>E:
                    C_nnew.append(B)
                else:
                    C_nnew.append(exp(-r*T/n)*(pstar*C_n[i]+(1-pstar)*C_n[i+1]))
    return C_nnew

def bse_pressed():
    global method
    global ostyle
    global otype
    method=0
    bse_button.config(relief="sunken")
    binomial_button.config(relief="raised")
    implicit_button.config(relief="raised")
    explicit_button.config(relief="raised")
    cn_button.config(relief="raised")
    ostyle=2
    otype=6
    european_button.config(state="normal")
    american_button.config(state="disabled",relief="raised")
    call_button.config(state="normal")
    put_button.config(state="normal")
    cashonothing_button.config(state="disabled",relief="raised")
    bullishverticalspread_button.config(state="disabled",relief="raised")
    assetonothing_button.config(state="disabled",relief="raised")
    straddle_button.config(state="disabled",relief="raised")
    
    
def binomial_pressed():
    global method
    method=1
    bse_button.config(relief="raised")
    binomial_button.config(relief="sunken")
    implicit_button.config(relief="raised")
    explicit_button.config(relief="raised")
    cn_button.config(relief="raised")
    european_button.config(state="normal")
    american_button.config(state="normal")
    call_button.config(state="normal")
    put_button.config(state="normal")
    cashonothing_button.config(state="normal")
    bullishverticalspread_button.config(state="normal")
    assetonothing_button.config(state="normal")
    straddle_button.config(state="normal")
    
    
def implicit_pressed():
    global method
    global otype
    method=2
    bse_button.config(relief="raised")
    binomial_button.config(relief="raised")
    implicit_button.config(relief="sunken")
    explicit_button.config(relief="raised")
    cn_button.config(relief="raised")
    european_button.config(state="normal")
    american_button.config(state="normal")
    otype=6
    call_button.config(state="normal",relief="raised")
    put_button.config(state="normal",relief="raised")
    cashonothing_button.config(state="disabled",relief="raised")
    bullishverticalspread_button.config(state="disabled",relief="raised")
    assetonothing_button.config(state="disabled",relief="raised")
    straddle_button.config(state="disabled",relief="raised")

def explicit_pressed():
    global method
    global otype
    method=3
    bse_button.config(relief="raised")
    binomial_button.config(relief="raised")
    implicit_button.config(relief="raised")
    explicit_button.config(relief="sunken")
    cn_button.config(relief="raised")
    european_button.config(state="normal")
    american_button.config(state="normal")
    otype=6
    call_button.config(state="normal",relief="raised")
    put_button.config(state="normal",relief="raised")
    cashonothing_button.config(state="disabled",relief="raised")
    bullishverticalspread_button.config(state="disabled",relief="raised")
    assetonothing_button.config(state="disabled",relief="raised")
    straddle_button.config(state="disabled",relief="raised")

def cn_pressed():
    global method
    global otype
    method=4
    bse_button.config(relief="raised")
    binomial_button.config(relief="raised")
    implicit_button.config(relief="raised")
    explicit_button.config(relief="raised")
    cn_button.config(relief="sunken")
    european_button.config(state="normal")
    american_button.config(state="normal")
    otype=6
    call_button.config(state="normal",relief="raised")
    put_button.config(state="normal",relief="raised")
    cashonothing_button.config(state="disabled",relief="raised")
    bullishverticalspread_button.config(state="disabled",relief="raised")
    assetonothing_button.config(state="disabled",relief="raised")
    straddle_button.config(state="disabled",relief="raised")

def euro_pressed():
    global ostyle
    ostyle=0
    european_button.config(relief="sunken")
    american_button.config(relief="raised")

def amer_pressed():
    global ostyle
    ostyle=1
    american_button.config(relief="sunken")
    european_button.config(relief="raised")

def call_pressed():
    global otype
    otype=0
    call_button.config(relief="sunken")
    put_button.config(relief="raised")
    cashonothing_button.config(relief="raised")
    bullishverticalspread_button.config(relief="raised")
    assetonothing_button.config(relief="raised")
    straddle_button.config(relief="raised")
    E1_entry.config(state="disabled")
    E1_entry.delete(0, END)
    E2_entry.config(state="disabled")
    E2_entry.delete(0, END)
    B_entry.config(state="disabled")
    B_entry.delete(0, END)
    
def put_pressed():
    global otype
    otype=1
    call_button.config(relief="raised")
    put_button.config(relief="sunken")
    cashonothing_button.config(relief="raised")
    bullishverticalspread_button.config(relief="raised")
    assetonothing_button.config(relief="raised")
    straddle_button.config(relief="raised")
    E1_entry.config(state="disabled")
    E1_entry.delete(0, END)
    E2_entry.config(state="disabled")
    E2_entry.delete(0, END)
    B_entry.config(state="disabled")
    B_entry.delete(0, END)
    
def cashonothing_pressed():
    global otype
    otype=2
    call_button.config(relief="raised")
    put_button.config(relief="raised")
    cashonothing_button.config(relief="sunken")
    bullishverticalspread_button.config(relief="raised")
    assetonothing_button.config(relief="raised")
    straddle_button.config(relief="raised")
    E1_entry.config(state="disabled")
    E1_entry.delete(0, END)
    E2_entry.config(state="disabled")
    E2_entry.delete(0, END)
    B_entry.config(state="normal")

def bullishverticalspread_pressed():
    global otype
    otype=3
    call_button.config(relief="raised")
    put_button.config(relief="raised")
    cashonothing_button.config(relief="raised")
    bullishverticalspread_button.config(relief="sunken")
    assetonothing_button.config(relief="raised")
    straddle_button.config(relief="raised")
    E1_entry.config(state="normal")
    E2_entry.config(state="normal")
    B_entry.config(state="disabled")
    B_entry.delete(0, END)

def assetonothing_pressed():
    global otype
    otype=4
    call_button.config(relief="raised")
    put_button.config(relief="raised")
    cashonothing_button.config(relief="raised")
    bullishverticalspread_button.config(relief="raised")
    assetonothing_button.config(relief="sunken")
    straddle_button.config(relief="raised")
    E1_entry.config(state="disabled")
    E1_entry.delete(0, END)
    E2_entry.config(state="disabled")
    E2_entry.delete(0, END)
    B_entry.config(state="disabled")
    B_entry.delete(0, END)

def straddle_pressed():
    global otype
    otype=5
    call_button.config(relief="raised")
    put_button.config(relief="raised")
    cashonothing_button.config(relief="raised")
    bullishverticalspread_button.config(relief="raised")
    assetonothing_button.config(relief="raised")
    straddle_button.config(relief="sunken")
    E1_entry.config(state="disabled")
    E1_entry.delete(0, END)
    E2_entry.config(state="disabled")
    E2_entry.delete(0, END)
    B_entry.config(state="disabled")
    B_entry.delete(0, END)

def display(text):
    print(text)

def GO_pressed():
    global r
    global vol
    global div
    global T
    global S
    global E
    global n
    global B
    global E1
    global E2
    global prog_stringvar
    if method>4:
        print(method)
        display("Please select a method of calculation.")
        return 0
    if ostyle>1:
        display("Please select an option style.")
        return 0
    if otype>5:
        display("Please select an option type.")
        return 0
    if otype==0 or otype==1 or otype==4:
        r=r_stringvar.get()
        vol=vol_stringvar.get()
        div=div_stringvar.get()
        T=T_stringvar.get()
        S=S_stringvar.get()
        E=E_stringvar.get()
        n=n_stringvar.get()
        
    elif otype==2:
        r=r_stringvar.get()
        vol=vol_stringvar.get()
        div=div_stringvar.get()
        T=T_stringvar.get()
        S=S_stringvar.get()
        E=E_stringvar.get()
        n=n_stringvar.get()
        B=B_stringvar.get()
        try:
            B=float(B)
        except:
            display("Cash value must be decimal.")
            return 0
        if B<=0:
            display("Cash value must be positive.")
            return 0
    elif otype==3:
        r=r_stringvar.get()
        vol=vol_stringvar.get()
        div=div_stringvar.get()
        T=T_stringvar.get()
        S=S_stringvar.get()
        E=E_stringvar.get()
        n=n_stringvar.get()
        E1=E1_stringvar.get()
        E2=E2_stringvar.get()
        try:
            E1=float(E1)
        except:
            display("First strike price must be decimal.")
            return 0
        try:
            E2=float(E2)
        except:
            display("Second strike price must be decimal.")
            return 0            
    try:
        r=float(r)
    except:
        display("Risk-free rate not a decimal.")
        return 0
    try:
        vol=float(vol)
    except:
        display("Volatility not a decimal.")
        return 0
    try:
        div=float(div)
    except:
        display("Dividend not a decimal.")
        return 0
    try:
        T=float(T)
    except:
        display("Time period not a decimal.")
        return 0
    try:
        S=float(S)
    except:
        display("Asset value rate not a decimal.")
        return 0
    try:
        E=float(E)
    except:
        display("Strike price not a decimal.")
        return 0
    try:
        n=int(n)
    except:
        display("Number of iterations not an integer.")
        return 0
    if vol<0:
        display("Volatility must be non-negative.")
        return 0
    if div<0:
        display("Dividend yield must be non-negative.")
        return 0
    if n<=0:
        display("Number of iterations must be positive.")
        return 0
    if T<=0:
        display("Time period must be positive.")
        return 0
    display("All variable values validated. Beginning calculation")
    if method==0:
        d1 = calcd1(S,E,r,vol,div,T)
        d2 = calcd2(d1,vol,T)
        V=calcV(otype,d1,d2,S,div,T,E,r)
        display("Option priced at: "+str(V))
        
    elif method==1:
        d=calcd(r,div,vol,n,T)
        u=calcu(r,div,vol,n,T)
        pstar=calcpstar(u,d,T,n,r,div)
        S_n=calc_S_n(u,d,S,n)
        C_n=calc_C_n(E,S_n,otype,B,E1,E2)
        del S_n
        for j in range(0,n):
            prog_stringvar.set(str(round((j*(j+1))/(n*(n+1)),2))+"%")
            C_n=calc_next(C_n,pstar,r,T,n,ostyle,otype,u,d,S,E,j)
        C_n=C_n[0]
        display("Option priced at: "+str(C_n))

    elif method==2:
        dx=icalcdx(S,E,n)
        dtau=icalcdtau(vol,T,n)
        k=icalck(r,vol)
        a=icalca(dtau,dx)
        ldiag=[]
        odiag=[]
        u=[]
        for i in arange(log(S/E)-100*dx,log(S/E)+100*dx,dx):
            u.append(iui0(k,dx,i,otype))
        for j in range(0,200):
            ldiag.append(1+(2*a))
        for j in range(0,199):
            odiag.append(-a)
        diagonals=[odiag,ldiag,odiag]
        A=diags(diagonals,[-1,0,1]).toarray()
        A[0,1]=1
        A[len(A)-2,len(A)-1]=1
        for j in range(1,int((pow(vol,2)*T)/(2*dtau))):
            b=[iu0T(S,E,dx,k,j,dtau,otype,ostyle),u[1]+a*iu0T(S,E,dx,k,j,dtau,otype,ostyle)]
            for i in range(2,198):
                b.append(u[i])
            b.append(u[len(u)-2]+a*iuXT(S,E,k,j,dtau,otype,dx,ostyle))
            b.append(iuXT(S,E,k,j,dtau,otype,dx,ostyle))
            b=asarray(b)
            u=solve(A,b)
            if ostyle==1:
                for i in range(1,199):
                    if otype==0:
                        u[i]=max([u[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
                    elif otype==1:
                        u[i]=max([u[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
        ux=u[int(len(u)/2)]
        C=(1-div)*pow(E,0.5*(1+k))*pow(S,0.5*(1-k))*exp((-pow(k+1,2)*pow(vol,2)*T)/8)*ux
        display("Option Priced at: "+str(C))

    elif method==3:
        dx=ecalcdx(S,E,n)
        dtau=ecalcdtau(dx,vol,T,n)
        k=ecalck(r,vol)
        a=ecalca(dtau,dx)
        Sj=ecalc_S_0(S,E,k,dx)
        Sjnew=[]
        for j in range(0,int((pow(vol,2)*T)/(2*dtau))):
            Sjnew=ecalc_S_jnew(Sj,a,S,E,dx,j,dtau,k,otype,ostyle)
            Sj=Sjnew
        ux=Sj[int(len(Sj)/2)]
        C=(1-div)*pow(E,0.5*(1+k))*pow(S,0.5*(1-k))*exp((-pow(k+1,2)*pow(vol,2)*T)/8)*ux
        display("Option Priced at: "+str(C))

    elif method==4:
        dx=cncalcdx(S,E,n)
        dtau=cncalcdtau(dx,vol,T,n)
        k=cncalck(r,vol)
        a=cncalca(dtau,dx)
        ldiag=[]
        odiag=[]
        u=[]
        for i in arange(log(S/E)-100*dx,log(S/E)+100*dx,dx):
            u.append(cnui0(k,dx,i,otype))
        for j in range(0,200):
            ldiag.append(1+(2*a))
        for j in range(0,199):
            odiag.append(-a)
        diagonals=[odiag,ldiag,odiag]
        A=diags(diagonals,[-1,0,1]).toarray()
        A[0,1]=1
        A[len(A)-2,len(A)-1]=1
        for j in range(0,int((pow(vol,2)*T)/(2*dtau))):
            b=[cnu0T(S,E,dx,k,j,dtau,otype,ostyle),u[1]+a*cnu0T(S,E,dx,k,j,dtau,otype,ostyle)]
            for i in range(2,198):
                b.append(u[i])
            b.append(u[len(u)-2]+a*cnuXT(S,E,k,j,dtau,otype,dx,ostyle))
            b.append(cnuXT(S,E,k,j,dtau,otype,dx,ostyle))
            b=asarray(b)
            Iu=solve(A,b)
            if ostyle==1:
                for i in range(1,199):
                    if otype==0:
                        Iu[i]=max([Iu[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
                    elif otype==1:
                        Iu[i]=max([Iu[i],exp(pow((k+1)/2,2)*j*dtau)*max([exp(0.5*(k+1)*(log(S/E)-(100-i)*dx))-exp(0.5*(k-1)*(log(S/E)-(100-i)*dx)),0])])
            Eu=cncalc_S_jnew(u,a,S,E,dx,j,dtau,k,otype,ostyle)
            for i in range(0,200):
                u[i]=0.5*(Eu[i]+Iu[i])
        ux=u[int(len(u)/2)]
        C=pow(E,0.5*(1+k))*pow(S,0.5*(1-k))*exp((-pow(k+1,2)*pow(vol,2)*T)/8)*ux
        print("Option Priced at: ",C)

    

root=Tk()
left_widgets=Frame(root)

#title & notes
title_label=Label(left_widgets,text="Option Valuation Calculator",bg="white")
notes_label=Label(left_widgets,text="An option value calculator for multiple option styles, types and calculation methods.",bg="white")
#

#stringvars
r_stringvar=StringVar()
vol_stringvar=StringVar()
div_stringvar=StringVar()
T_stringvar=StringVar()
S_stringvar=StringVar()
E_stringvar=StringVar()
B_stringvar=StringVar()
E1_stringvar=StringVar()
E2_stringvar=StringVar()
n_stringvar=StringVar()
prog_stringvar=StringVar()
#

#algorithm method and frame
algo_button_frame=Frame(left_widgets)
bse_button=Button(algo_button_frame,text="BSE",command=bse_pressed)
binomial_button=Button(algo_button_frame,text="Binomial",command=binomial_pressed)
implicit_button=Button(algo_button_frame,text="Implicit",command=implicit_pressed)
explicit_button=Button(algo_button_frame,text="Explicit",command=explicit_pressed)
cn_button=Button(algo_button_frame,text="C-N",command=cn_pressed)
#

#option style buttons and frame
style_button_frame=Frame(left_widgets)
european_button=Button(style_button_frame,text="European",command=euro_pressed)
american_button=Button(style_button_frame,text="American",command=amer_pressed)
#

#option type buttons and frame
type_button_frame=Frame(left_widgets)
call_button=Button(type_button_frame,text="Call",command=call_pressed)
put_button=Button(type_button_frame,text="Put",command=put_pressed)
cashonothing_button=Button(type_button_frame,text="Cash-or-Nothing",command=cashonothing_pressed)
bullishverticalspread_button=Button(type_button_frame,text="Bullish Vertical Spread",command=bullishverticalspread_pressed)
assetonothing_button=Button(type_button_frame,text="Asset-or-Nothing",command=assetonothing_pressed)
straddle_button=Button(type_button_frame,text="Straddle",command=straddle_pressed)
#

#parameter labels and frame
parameter_frame=Frame(left_widgets)
r_label=Label(parameter_frame,text="Risk-free rate: ")
vol_label=Label(parameter_frame,text="Volatility: ")
div_label=Label(parameter_frame,text="Dividend yield: ")
T_label=Label(parameter_frame,text="Time period: ")
S_label=Label(parameter_frame,text="Underlying asset: ")
E_label=Label(parameter_frame,text="Strike price: ")
B_label=Label(parameter_frame,text="Cash : ")
E1_label=Label(parameter_frame,text="First strike price: ")
E2_label=Label(parameter_frame,text="Second strike price: ")
n_label=Label(parameter_frame,text="Number of iterations: ")
#

#parameter entries and frame
r_entry=Entry(parameter_frame,textvariable=r_stringvar)
vol_entry=Entry(parameter_frame,textvariable=vol_stringvar)
div_entry=Entry(parameter_frame,textvariable=div_stringvar)
T_entry=Entry(parameter_frame,textvariable=T_stringvar)
S_entry=Entry(parameter_frame,textvariable=S_stringvar)
E_entry=Entry(parameter_frame,textvariable=E_stringvar)
B_entry=Entry(parameter_frame,textvariable=B_stringvar)
E1_entry=Entry(parameter_frame,text=E1_stringvar)
E2_entry=Entry(parameter_frame,text=E2_stringvar)
n_entry=Entry(parameter_frame,text=n_stringvar)
#

#algorithm buttons
go_button=Button(left_widgets,text="GO",command=GO_pressed)
progress_label=Label(left_widgets,text="0%")
console_output_label=Label(root,text="output label content",bg="white")
#

#first offspring griding
left_widgets.grid(row=0,column=0)
title_label.grid(row=0,column=0)
notes_label.grid(row=1,column=0)
algo_button_frame.grid(row=2,column=0)
style_button_frame.grid(row=3,column=0)
type_button_frame.grid(row=4,column=0)
parameter_frame.grid(row=5,column=0)
#

#algorithm type griding
bse_button.grid(row=0,column=0)
binomial_button.grid(row=0,column=1)
implicit_button.grid(row=0,column=2)
explicit_button.grid(row=0,column=3)
cn_button.grid(row=0,column=4)
#

#style griding
european_button.grid(row=0,column=0)
american_button.grid(row=0,column=1)
#

#type griding
call_button.grid(row=0,column=0)
put_button.grid(row=0,column=1)
cashonothing_button.grid(row=0,column=2)
bullishverticalspread_button.grid(row=0,column=3)
assetonothing_button.grid(row=0,column=4)
straddle_button.grid(row=0,column=5)
#

#parameter griding
r_label.grid(row=0,column=0)
r_entry.grid(row=0,column=1)
vol_label.grid(row=0,column=2)
vol_entry.grid(row=0,column=3)
div_label.grid(row=0,column=4)
div_entry.grid(row=0,column=5)
T_label.grid(row=1,column=0)
T_entry.grid(row=1,column=1)
S_label.grid(row=1,column=2)
S_entry.grid(row=1,column=3)
E_label.grid(row=1,column=4)
E_entry.grid(row=1,column=5)
B_label.grid(row=2,column=0)
B_entry.grid(row=2,column=1)
E1_label.grid(row=2,column=2)
E1_entry.grid(row=2,column=3)
E2_label.grid(row=2,column=4)
E2_entry.grid(row=2,column=5)
n_label.grid(row=3,column=0)
n_entry.grid(row=3,column=1)
#

#algo griding
go_button.grid(row=6,column=0)
progress_label.grid(row=7,column=0)
console_output_label.grid(row=0,column=1)
#

method=5
ostyle=2
otype=6
r=vol=div=T=S=E=B=E1=E2=n="a"

root.mainloop()
