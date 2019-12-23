import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

cdh1 = np.logspace(-5,0.02,2000)
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
def cycbnc(cdh1):
    beta = k1/k2dd
    J = k2d/k2dd
    cycb = [beta/(J + c) for c in cdh1]
    return cycb

def cdh1nc(cdh1,m=0.03):
    rho = (k3d + k3dd*A)/(k4*m)
    cycb = [rho*((1-c)*(J4 + c))/(c*(J3 + 1 - c)) for c in cdh1]
    return cycb

mval = st.slider(label='Mass', min_value=0.1, max_value=0.7)
cycb1 = cycbnc(cdh1)
cycb2 = cdh1nc(cdh1,m=mval)
plt.plot(cycb1, cdh1, 'b', label='CycB nullcline')
plt.plot(cycb2, cdh1, 'r', label='Cdh1 nullcline')
for i in range(len(cdh1)):
    if abs(cycb1[i] - cycb2[i]) < 1e-3:
        plt.plot(cycb1[i], cdh1[i],'ko')
#plt.close()
# cycb_1 = cycbnc(cdh1)
# plt.plot(cycb_1,cdh1,'k')
# cycb_2 = cdh1nc(cdh1,m=0.3)
# plt.plot(cycb_2,cdh1,'g--',label='m=0.3')
# cycb_3 = cdh1nc(cdh1,m=0.6)
# plt.plot(cycb_3,cdh1,'r--',label='m=0.6')

plt.xscale('log')
plt.xlabel('[CycB]')
plt.ylabel('[Cdh1]')
plt.ylim([-0.01,1.01]) 
plt.xlim([0.01,10]) 
plt.legend()
st.pyplot()
#plt.savefig('fig2.png')
