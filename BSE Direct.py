from math import exp, log, sqrt
try:
    from statistics import NormalDist
except ImportError:
    print("Please update to python 3.8")
    exit()

def calcd1(S,E,r,vol,div,T):
    return (log(S/E)+(r-div+0.5*pow(vol,2))*T)/(vol*sqrt(T))

def calcd2(d1,vol,T):
    return d1-(vol*sqrt(T))

def calcV(otype,d1,d2,S,div,T,E,r):
    if otype==0:
        return (S*exp(-div*T)*NormalDist().cdf(d1))-(E*exp(-r*T)*NormalDist().cdf(d2))
    elif otype==1:
        return (E*exp(-r*T)*NormalDist().cdf(-d2))-(S*exp(-div*T)*NormalDist().cdf(-d1))

while True:
    try:
        otype=int(input("If the option is a call, enter 0. If put, enter 1: "))
    except ValueError:
        print("Option type must be described as an integer.")
        continue
    else:
        if otype>=0 and otype<=5:
            break
        print("Option type must be 0 or 1.")

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

while True:
        try:
            E=float(input("Please enter the options strike price as a decimal: "))
        except ValueError:
            print("Strike price must be a decimal.")
            continue
        else:
            break

d1 = calcd1(S,E,r,vol,div,T)
d2 = calcd2(d1,vol,T)
print("BSE option cost is: ",str(calcV(otype,d1,d2,S,div,T,E,r)))
