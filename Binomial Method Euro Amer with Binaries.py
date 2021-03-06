from math import exp,sqrt
from tkinter import *

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

top=Tk()
progSV=StringVar()
progress=Label(top,textvariable=progSV)
progSV.set("0%")
progress.pack()

print("----European and American style call/put option type approximator using binomial method. Created by Mackenzie Langdown 19/12/2019----")
print("To create approximations you must know the risk-free rate, implied volatility, length of the option, dividend payment, starting stock value, strike price, type of option and style of option.")

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
        otype=int(input("If the option is a call, enter 0. If put, enter 1. If Cash-Or-Nothing, enter 2. If Bullish-Vertical Spread, enter 3. If Asset-Or-Nothing, enter 4. If Straddle, enter 5: "))
    except ValueError:
        print("Option type must be described as an integer.")
        continue
    else:
        if otype>=0 and otype<=5:
            break
        print("Option type must be 0, 1, 2, 3, 4 or 5.")

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
        n=int(input("Please enter the number of steps: "))
    except ValueError:
        print("The number of steps must be an integer number.")
        continue
    else:
        if n>0:
            break
        print("The number of steps must be positive.")

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

while True:
    try:
        div=float(input("Please enter the dividend payment as a decimal: "))
    except ValueError:
        print("Dividend payment must be a decimal.")
        continue
    else:
        if div>=0:
            break
        print("Dividend payment must be greater than or equal to 0.")

while True:
    try:
        S=float(input("Please enter the starting stock value as a decimal: "))
    except ValueError:
        print("Starting stock value must be a decimal.")
        continue
    else:
        break

B=0
if otype!=3:
    while True:
        try:
            E=float(input("Please enter the options strike price as a decimal: "))
        except ValueError:
            print("Strike price must be a decimal.")
            continue
        else:
            break
if otype==2:
    while True:
        try:
            B=float(input("Please enter the Cash value for the option as a decimal: "))
        except ValueError:
            print("Cash value must be a decimal.")
            continue
        else:
            if B>0:
                break
            print("Cash value must be greater than zero.")
E1=0
E2=0
if otype==3:
    while True:
        try:
            E1=float(input("Please enter the options first strike price as a decimal: "))
        except ValueError:
            print("Strike price must be a decimal.")
            continue
        else:
            break
    while True:
        try:
            E2=float(input("Please enter the options second strike price as a decimal: "))
        except ValueError:
            print("Strike price must be a decimal.")
            continue
        else:
            break

print("Values accepted, starting backwards interations, a progress bar will be displayed in another window and the final result will be printed on the console.")
d=calcd(r,div,vol,n,T)
u=calcu(r,div,vol,n,T)
pstar=calcpstar(u,d,T,n,r,div)
S_n=calc_S_n(u,d,S,n)
C_n=calc_C_n(E,S_n,otype,B,E1,E2)
del S_n
for j in range(0,n):
    C_n=calc_next(C_n,pstar,r,T,n,ostyle,otype,u,d,S,E,j)
    progSV.set(str(round((j/n)*100,2))+"%")
    top.update()
progSV.set("100%")
top.update()
C_n=C_n[0]
print("Option Priced at: ",C_n)
