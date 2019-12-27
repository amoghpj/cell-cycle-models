import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from io import BytesIO
import streamlit as st

def lotkavolterra(X, t, args):
    rabbit, fox = X
    alpha= 1.1
    beta = 0.4
    delta = 0.1
    gamma = 0.4
    drabbit = alpha*rabbit - beta*rabbit*fox
    dfox = delta*fox*rabbit  - gamma*fox
    return([drabbit, dfox])

def testLV():
    """Lotka-Volterra test"""
    f, ax = plt.subplots(1,1)
    x0 = [10, 5]
    rabbiti = st.slider('rabbit', min_value=0, max_value=20, step=1)
    foxi = st.slider('fox', min_value=0, max_value=20, step=1)
    X = odeint(lotkavolterra, [rabbiti, foxi], np.linspace(0,10,500), args = ({},))
    gridsize =20
    maxrabbit = 20
    maxfox = 10
    xlim = np.linspace(-1,maxrabbit,gridsize)
    ylim = np.linspace(-1,maxfox,gridsize)
    makePP(lotkavolterra, xlim, ylim,ax, gridsize, {})
    ax.plot(X[:,0],X[:,1])
    ax.plot(rabbiti, foxi, 'ko')
    ax.axhline(1.1/0.4)
    ax.axvline(0.1/0.4)
    ax.set_xlim(-1,maxrabbit)
    ax.set_ylim(-1,maxfox)
    ax.set_xlabel('rabbit')
    ax.set_ylabel('fox')
    st.pyplot()

def makePP(model, xlim, ylim, ax, gridsize,args):
    xstep = 0.5*max(xlim)/gridsize# length of vectors
    ystep = 0.5*max(ylim)/gridsize
    for x in xlim:
        for y in ylim:
            der = model([x,y],0.1, args)
            theta = np.arctan(abs(der[1]/der[0]))
            xsign = 1.
            ysign = 1.
            if der[0] < 0:
                xsign = -1
            if der[1] < 0:
                ysign = -1
            deltax = xstep*np.cos(theta)*xsign
            deltay = ystep*np.sin(theta)*ysign

            ax.arrow(x,np.log10(y),
                     deltax,
                     deltay,
                     head_width=0.5*min(xstep, ystep),
                     color='k',
                     alpha=0.3,
                     length_includes_head=False)

def getRoots(y1, y2):
    err = [abs(y1[i] - y2[i]) for i in range(len(y1))]
    winsize = 2
    differr = [(err[i+winsize] - err[i])/winsize for i in range(len(err)-winsize)]
    solutions = []

    for i in range(len(differr)-winsize):
        if err[i+winsize] < 0.05 and err[i] < 0.05:
            if differr[i + winsize] >0 and differr[i] < 0 :
                solutions.append(i+1)

    f, ax = plt.subplots(2,1)
    ax[0].plot(np.log10(err))
    signde = []
    for d in differr:
        if d > 0 :
            signde.append(1)
        elif d<0:
            signde.append(-1)
        else:
            signde.append(0)
    ax[1].plot(signde)
    extreme = [err[s] for s in solutions]
    #ax[0].axhline(np.log10(min(extreme)), color='b')
    for s in solutions:
        ax[1].plot(s,0,'r.')
        ax[0].axvline(s,color='r',alpha=0.1)
    ax[1].set_title('err')
    #st.pyplot()
    plt.close()
    return solutions

def integrate(func, x0, tspan, parameters, stepsize=0.01, method='rk45'):
    methoddict = {'rk45':rk45,
                  'euler':euler}
    xprev = x0
    t0 = min(tspan)
    tmax = max(tspan)
    size = int(tmax/stepsize)
    timecourse = np.zeros(shape=(size, len(x0)))
    t = t0
    counter = 0
    divide = False
    while counter < size:
        dX = func(xprev, t, parameters)
        x = []
        x = xprev + stepsize*(methoddict[method](func, xprev, t, stepsize, parameters))
        # cycb
        if x[1] < 0.1 and timecourse[counter-1,1] >0:# and division = False:
            divide = True
            #x = [v/2. for v in x]
        if divide is True:
            x[5] = 0.6 # mass
            divide = False
        xprev = x
        t += stepsize
        timecourse[counter,: ] = x
        counter += 1
    return(timecourse)

def euler(function, x, t, args):
    dx = function(x, t, args)
    return dx

def rk45(function, x, t, stepsize, args):
    k1 = function(x, t, args)
    k2 = function(x + k1*stepsize/2., t + stepsize/2, args)
    k3 = function(x + k2*stepsize/2., t + stepsize/2., args)
    k4 = function(x + k3*stepsize/2., t + stepsize, args)
    return(k1 + 2.*k2 + 2.*k3 + k4)/6.

def fig1model(X, t, args):
    m = args['m']
    k1 = args['k1']
    k2d = args['k2d']
    k2dd = args['k2dd']
    k2ddd = args['k2ddd']
    k3d = args['k3d']
    k3dd = args['k3dd']
    k4d = args['k4d']
    k4 = args['k4']
    A = args['A']
    J3 = args['J3']
    J4 = args['J4']
    cdh1, cycb = X
    dcycb = k1- (k2d + k2dd * cdh1)*cycb
    dcdh1 = ((k3d + k3dd*A)*(1 - cdh1))/(J3 + 1 - cdh1) - (k4*m*cycb*cdh1)/(J4 + cdh1)
    return([dcdh1, dcycb])

def makeFig3():
    return

def cycbnc_fig2(cdh1, parameters):
    beta = parameters['k1']/parameters['k2dd']
    J = parameters['k2d']/parameters['k2dd']
    cycb = [beta/(J + c) for c in cdh1]
    return cycb

def cdh1nc_fig2(cdh1,parameters):
    p = (parameters['k3d'] + parameters['k3dd']*parameters['A'])/(parameters['k4']*parameters['m'])
    cycb = [p*((1-c)*(parameters['J4'] + c))/(c*(parameters['J3'] + 1 - c)) for c in cdh1]
    return cycb

def makeFig2(parameters):
    st.markdown("$\\frac{d[CycB]}{dt} = k_1 - (k_2' + k_2'' [Cdh1])[CycB]$")
    st.markdown("$\\frac{d[Cdh1]}{dt} = \\frac{(k_3' + k_3'' A)(1- [Cdh1])}{J_3 + 1 - [Cdh1]} - \\frac{k_4 m [CycB] [Cdh1]}{J_4 + [Cdh1]}$")
    Cdh1_i = st.slider(label='Cdh1', min_value=0.0, max_value=1.0, value=0.1,step=0.1)
    CycB_i = st.slider(label='log(CycB)', min_value=-2., max_value=1., value=1e-1,step=0.1)
    mval = st.slider(label='Mass', min_value=0.3, max_value=0.7, value=0.3,step=0.01)    
    parameters['m'] = mval
    cdh1 = np.logspace(-6,0.1,8000)
    cycb1 = cycbnc_fig2(cdh1, parameters)
    cycb2 = cdh1nc_fig2(cdh1, parameters)
    parameters['m'] =  mval
    solutions = getRoots(np.log10(cycb1), np.log10(cycb2))
    plt.close()
    t = np.linspace(0,100,500)
    y = odeint(fig1model, [Cdh1_i, 10**CycB_i],t,args=(parameters,))

    # Root finding
    f, ax = plt.subplots(1,1)
    # # make phaseplane                                #
    # gridsize = 10                                    #
    # xlim= np.linspace(0,1., gridsize)                #
    # ylim = np.logspace(-2, 1., gridsize)             #
    # #makePP(fig1model,xlim, ylim, ax, gridsize,args) #
    # for i in range(len(cdh1)):
    #     if abs(cycb1[i] - cycb2[i]) < 10.**(np.ceil(np.log10(min(cycb1[i],cycb2[i])))-2):
    #         ax.plot(cdh1[i],np.log10(cycb1[i]), 'ko')

    ax.plot(cdh1,np.log10(cycb1), 'b', label='CycB nullcline')
    ax.plot(cdh1,np.log10(cycb2), 'r', label='Cdh1 nullcline')
    ax.plot(Cdh1_i,CycB_i,'ko')
    ax.plot(y[:,0], np.log10(y[:,1]), 'k--',alpha=0.5, lw=2.0)
    ax.plot(y[-1,0], np.log10(y[-1,1]), 'ro', lw=2.0)

    for s in solutions:
        ax.plot(cdh1[s], np.log10(cycb1[s]), 'go')

    ax.set_ylabel('log([CycB])')
    ax.set_xlabel('[Cdh1]')
    ax.set_xlim([-0.05,1.01]) 
    ax.set_ylim([-2,1]) 
    ax.legend()

    plt.tight_layout()
    st.pyplot()

    # f, ax = plt.subplots(1,1)
    # ax.plot(t, y[:,0], 'r', label='Cdh1')
    # ax.plot(t, np.log10(y[:,1]), 'b', label='CycB')
    # ax.set_xlabel('time')
    # ax.set_ylabel('[X]')       
    # ax.legend()
    # plt.tight_layout()
    # st.pyplot()

def makeFig3(parameters):       
    mvals = np.linspace(0.01,0.9,500)

    cdh1 = np.logspace(-4, 0.1, 7000)
    hyst = []
    pvals = []
    for m in mvals:
        p = (parameters['k3d'] + parameters['k3dd']*parameters['A'])/(parameters['k4']*m)
        parameters['m'] = m
        cycb1 = cycbnc_fig2(cdh1, parameters)
        cycb2 = cdh1nc_fig2(cdh1, parameters)
        solutions = getRoots(np.log10(cycb1), np.log10(cycb2))
        for s in solutions:
            hyst.append(cycb1[s])
            pvals.append(p)
    plt.plot(pvals, hyst,'k.')
    plt.xlabel('p')
    plt.ylabel('[CycB]')
    st.pyplot()

#def makeFig4(parameters):

def fullmodel(X, t, args):
    k1 = args['k1']
    k2d = args['k2d']
    k2dd = args['k2dd']
    k2ddd = args['k2ddd']
    k3d = args['k3d']
    k3dd = args['k3dd']
    k4d = args['k4d']
    k4 = args['k4']
    A = args['A']
    J3 = args['J3']
    J4 = args['J4']
    mu = args['mu']
    J5 = args['J5']
    Mad = args['Mad']
    k6 = args['k6']
    k7 = args['k7']
    k8 = args['k8']
    n = args['n']
    k5d = args['k5d']
    k5dd = args['k5dd']
    J7 = args['J7']
    J8 = args['J8']
    mstar = args['mstar']
    k9 = args['k9']
    k10 = args['k10']
    cdh1, cycb, cdc20t, cdc20a, iep, m = X
    # if cycb < 0.1:
    #     m = m/2.
    dcdh1 = ((k3d + k3dd*cdc20a)*(1 - cdh1))/(J3 + 1 - cdh1) - (k4*m*cycb*cdh1)/(J4 + cdh1)
    dcycb = k1- (k2d + k2dd * cdh1)*cycb    
    dcdc20t = k5d + k5dd*( (cycb*m/J5)**n /(1+ (cycb*(m/J5))**n )) - k6*cdc20t
    dcdc20a = (k7*iep*(cdc20t-cdc20a)/(J7 + cdc20t - cdc20a)) - (k8*Mad*cdc20a)/(J8+cdc20a) - k6*cdc20a
    diep = k9*m*cycb*(1-iep) - k10*iep
    dm = mu*m*(1-m/mstar)
    return np.array(([dcdh1, dcycb, dcdc20t, dcdc20a, diep, dm]))

def plottimecourses(parameters):
    x0 = [1.0, 0.5,1.5, 1.4, 0.7, 0.6]
    stepsize = 0.01
    tmax = 160
    t= np.linspace(0, tmax, int(tmax/stepsize))
    #y = odeint(fullmodel,x0, t, args=(parameters,))
    y = integrate(fullmodel, x0, t, parameters, stepsize=stepsize)
    f , ax = plt.subplots(3,1, figsize=(3,6))
    ax[0].plot(t,y[:,5], label='m')
    ax[0].legend()
    ax[1].plot(t,y[:,0], label='Cdh1')
    ax[1].plot(t,y[:,1], label='CycB')
    ax[1].legend()
    ax[2].plot(t,y[:,2], label='Cdc20T')
    ax[2].plot(t,y[:,3], label='Cdc20A')
    ax[2].plot(t,y[:,4], label='IEP')        
    ax[2].set_ylim([0,1.8])
    ax[2].legend()
    st.pyplot()

def makeIntroPage():
    with open('markdown/intro.md','r') as infile:
        introtext = ''.join(infile.readlines())
    st.markdown(introtext)
    # st.markdown("TLDR: This project seeks to make a series of abstract models of "
    # "the eukaryotic cell cycle accessble to the non-modelers. The content is organized "
    # "as per the ideas developed in [Tyson and Novak, 2001](https://www.ncbi.nlm.nih.gov/pubmed/11371178). "
    # "This interactive site is meant to be an educational tool, aimed at anyone who has been exposed to"
    # " the basic concepts of eukaroytic mitosis, and is curious about the utility of mathematical models "
    # "in making sense of complex biological processes.")

    # st.subheader("What are the cell cycle models all about?")

    # st.subheader("Why did you make this?")
    # st.markdown("The prototypical mathematical model of biological systems still seems to be the "
    # "Lotka-Volterra predator-prey model, from the 20th century. The curious student with an interest in molecular "
    # "biology *may* have come across the [reprissilator](https://en.wikipedia.org/wiki/Repressilator)."
    # " I believe that there is still a general lack of awareness of the success of mathematical models  of cellular processes, "
    # "ranging from the cell cycle, to circadian oscillations, to autophagy, and even dynamical models of cancer. "
    # "While there are general purpose tools [Cell Collective](https://cellcollective.org/#) that provide platforms to lower the barrier to entry"
    # " to these theoretical models, I have not come "
    # "across a curated, interactive resource exploring any of these models in depth. This is my attempt at "
    # "creating such a tool, focussed on the highly successful work by Tyson and Novak in the last couple of decades"
    # "on the yeast cell cycle."
    # "\n"
    # "\n"
    # "Please feel free to reach out with any feedback and comments at jamogh [at] vt [dot] edu, or on twitter [@amogh_jalihal](https://twitter.com/amogh_jalihal), "
    # "or by creating an issue on this project's [github repository](https://github.com/amoghpj/cell-cycle-models).")

def main():
    # parameterdict
    parameters = {
        'k1':0.04,
        'k2d':0.04,
        'k2dd':1.,
        'k2ddd':1.,
        'k3d':1.,
        'k3dd':10.,
        'k4d':2.,
        'k4':35.,
        'A':0.00,
        'J3':0.04,
        'J4':0.04,
        'k5d':0.005,
        'k5dd':0.2,
        'k6':0.1,
        'Mad':1.0,
        'k7':1.0,
        'k8':0.5,
        'k9':0.1,
        'k10':0.02,
        'k11':1.,
        'k12d':0.2,
        'k12dd':50,
        'k12ddd':100,
        'k13':1.,
        'k14':1.,
        'k15d':1.5,
        'k15dd':0.05,
        'k16d':1.0,
        'k16dd':3.0,
        'mu':0.01,
        'J5':0.3,
        'n':4,
        'J7':1e-3,
        'J8':1e-3,
        'Keq':1e3,
        'J15':0.01,
        'J16':0.01,
        'mstar':10,
    }

    page = st.sidebar.selectbox('Jump to...',['Introduction','Simplified Cdh1-CycB model', 'Hysteresis in transitions','Full Model'])
    if page == 'Introduction':
        st.header('Introduction')
        makeIntroPage()
    if page == 'Simplified Cdh1-CycB model':
        st.header('A simplified model of CycB/Cdk1-Cdh1/APC antagonism')
        makeFig2(parameters)
    if page == 'Hysteresis in transitions':
        st.header('Hystersis underlies cell state transitions')
        makeFig3(parameters)
    if page == 'Full Model':
        st.header('Full Model')
        plottimecourses(parameters)
if __name__ == '__main__':
    main()
