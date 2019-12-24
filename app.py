import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from io import BytesIO
import streamlit as st

def fig1model(X, t, m):
    k1 = 0.04
    k2d = 0.04
    k2dd = 1
    k2ddd = 1
    k3d = 1
    k3dd = 10
    k4d = 2
    k4 = 35
    A = 0.0
    J3 = 0.04
    J4 = 0.04
    cdh1, cycb = X
    dcycb = k1- (k2d + k2dd * cdh1)*cycb
    dcdh1 = ((k3d + k3dd*A)*(1 - cdh1))/(J3 + 1 - cdh1) - (k4*m*cycb*cdh1)/(J4 + cdh1)
    return([dcdh1, dcycb])

def makeFig1():
    st.markdown("$\\frac{d[CycB]}{dt} = k_1 - (k_2' + k_2'' [Cdh1])[CycB]$")
    st.markdown("$\\frac{d[Cdh1]}{dt} = \\frac{(k_3' + k_3'' A)(1- [Cdh1])}{J_3 + 1 - [Cdh1]} - \\frac{k_4 m [CycB] [Cdh1]}{J_4 + [Cdh1]}$")
    cdh1 = np.logspace(-5,0.0,5000)
    k1 = 0.04
    k2d = 0.04
    k2dd = 1.
    k2ddd = 1.
    k3d = 1.
    k3dd = 10.
    k4d = 2.
    k4 = 35.
    A = 0.01
    J3 = 0.04
    J4 = 0.04
    def cycbnc(cdh1):
        beta = k1/k2dd
        J = k2d/k2dd
        cycb = [beta/(J + c) for c in cdh1]
        return cycb

    def cdh1nc(cdh1,m=0.03):
        p = (k3d + k3dd*A)/(k4*m)
        cycb = [p*((1-c)*(J4 + c))/(c*(J3 + 1 - c)) for c in cdh1]
        return cycb

    mval = st.slider(label='Mass', min_value=0.1, max_value=0.7, value=0.3,step=0.05)
    Cdh1_i = st.slider(label='Cdh1', min_value=0.0, max_value=1.0, value=0.1,step=0.1)
    CycB_i = st.slider(label='log(CycB)', min_value=-2., max_value=1., value=1e-1,step=1.0)

    t = np.linspace(0,100,500)
    y = odeint(fig1model, [Cdh1_i, 10**CycB_i],t,args=(mval,))

    cycb1 = cycbnc(cdh1)
    cycb2 = cdh1nc(cdh1,m=mval)
    plt.plot(cycb1, cdh1, 'b', label='CycB nullcline')
    plt.plot(cycb2, cdh1, 'r', label='Cdh1 nullcline')
    plt.plot(y[:,1], y[:,0], 'k--',alpha=0.5, lw=2.0)
    for i in range(len(cdh1)):
        if abs(cycb1[i] - cycb2[i]) < 10.**(round(np.log10(min(cycb1[i],cycb2[i])),0) -1.5) :
            plt.plot(cycb1[i], cdh1[i],'ko')
    plt.xscale('log')
    plt.xlabel('[CycB]')
    plt.ylabel('[Cdh1]')
    plt.ylim([-0.1,1.01]) 
    plt.xlim([0.01,10]) 
    plt.legend()

    st.pyplot()

def main():
    page = st.sidebar.selectbox('Choose Model',['Simplified Cdh1-CycB'])
    if page == 'Simplified Cdh1-CycB':
        st.header('This generate Figure 1 from Tyson and Novak')
        makeFig1()
if __name__ == '__main__':
    main()
