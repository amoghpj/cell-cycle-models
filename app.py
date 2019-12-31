import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from io import BytesIO
import streamlit as st
from PIL import Image
import time

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

from streamlit.ScriptRequestQueue import RerunData
from streamlit.ScriptRunner import RerunException
from streamlit.server.Server import Server
import streamlit.ReportThread as ReportThread


def rerun():
    """Rerun a Streamlit app from the top!"""
    widget_states = _get_widget_states()
    raise RerunException(RerunData(widget_states))


def _get_widget_states():
    # Hack to get the session object from Streamlit.

    ctx = ReportThread.get_report_ctx()

    session = None
    session_infos = Server.get_current()._session_infos.values()

    for session_info in session_infos:
        if session_info.session._main_dg == ctx.main_dg:
            session = session_info.session

    if session is None:
        raise RuntimeError(
            "Oh noes. Couldn't get your Streamlit Session object"
            'Are you doing something fancy with threads?')
    # Got the session object!

    return session._widget_states
# def test(counter = 1):
#     if counter == 8:
#         counter = 1
#     plt.close()
#     img = np.imread
#     with open('./data/budding_0.5x/0' + str(counter) + '.png', 'r') as infile:
#         img = np.array(infile.readlines())
#     # f, ax = plt.subplots(1,1)
#     # s = time.time()
#     # ax.set_xlim(0,10)
#     # ax.set_ylim(0,10)
#     # ax.annotate(s, (5,5))
#     st.image(img)
#     time.sleep(1)        
#     counter += 1
#     rerun()

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
    #ax[1].plot()
    extreme = [err[s] for s in solutions]
    #ax[0].axhline(np.log10(min(extreme)), color='b')
    for s in solutions:
        ax[1].plot(s,0,'r.')
        ax[0].axvline(s,color='r',alpha=0.1)
    ax[1].set_title('err')
    plt.tight_layout()
    #st.pyplot()
    plt.close()
    return solutions

def goldbeter(va, vi, ja, ji):
    alpha = vi - va
    beta = vi - va + va*ji + vi*ja
    gamma = va*ji
    g = (2*gamma)/(beta + np.sqrt(beta**2 - 4*alpha*gamma))
    return g

def integrate(func, x0, tspan, parameters, massindex=5,stepsize=0.01, method='rk45'):
    methoddict = {'rk45':rk45,
                  'euler':euler}
    xprev = x0
    t0 = min(tspan)
    tmax = max(tspan)
    size = int(tmax/stepsize)
    timecourse = np.zeros(shape=(size, len(x0)))
    t = t0
    counter = 0
    growing = False
    while counter < size:
        dX = func(xprev, t, parameters)
        x = []
        x = xprev + stepsize*(methoddict[method](func, xprev, t, stepsize, parameters))
        # cycb
        if x[massindex]> 0.8:
            growing = True
        if x[1] < 0.1 and growing == True:
            x[massindex] = x[massindex]/2. # mass
            #x[massindex] = 0.4 # mass
            growing = False
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
    with open('./markdown/two-variable-timecourse.md','r') as infile:
        sec1text = ''.join(infile.readlines())
    st.markdown(sec1text)
    ####################################
    ### Time courses
    Cdh1_i = st.slider(label='Cdh1',key='cdh1fig2tc', min_value=0.0, max_value=1.0, value=0.87,step=0.051)
    CycB_i = st.slider(label='log(CycB)',key='cdh1fig2tc', min_value=-2., max_value=1., value=-1.8,step=0.1)
    mval = st.slider(label='Mass', key='massfig2tc',min_value=0.3, max_value=0.7, value=0.3,step=0.01)    
    f, ax = plt.subplots(1,1)
    t = np.linspace(0,100,500)
    cdh1 = np.append(np.logspace(-5,-3,1500), np.logspace(-3,0.1,6000))
    parameters['m'] =  mval
    y = odeint(fig1model, [Cdh1_i, 10**CycB_i],t,args=(parameters,))
    ax.plot(t, y[:,0],'k', label = 'Cdh1')
    ax.set_ylim(0,1.0)
    ax1 = ax.twinx()
    ax.set_ylabel('[Cdh1]')
    ax.set_xlabel('time (min)')
    ax1.set_ylabel('[CycB]')
    ax1.plot(t, y[:,1],'k--', label = 'CycB')
    ax1.set_ylim(1e-2,10.)
    ax1.set_yscale('log')
    ax.legend()
    ax1.legend()
    ax.set_title('mass = ' + str(round(mval, 2)))
    st.pyplot()
    plt.close()
    ###################################
    ### Nullclines
    with open('./markdown/two-variable-nullcline.md','r') as infile:
        sec1text = ''.join(infile.readlines())
    st.markdown(sec1text)
    Cdh1_i = st.slider(label='Cdh1',key='cdh1fig2nc', min_value=0.0, max_value=1.0, value=0.9,step=0.1)
    CycB_i = st.slider(label='log(CycB)', key='cycbfig2nc',min_value=-2., max_value=1., value=-1.9,step=0.1)
    mval = st.slider(label='Mass', key='massfig2nc',min_value=0.1, max_value=0.7, value=0.3,step=0.01)    
    parameters['m'] = mval

    cycb1 = cycbnc_fig2(cdh1, parameters)
    cycb2 = cdh1nc_fig2(cdh1, parameters)
    parameters['m'] =  mval
    solutions = getRoots(np.log10(cycb1), np.log10(cycb2))
    plt.close()

    y = odeint(fig1model, [Cdh1_i, 10**CycB_i],t,args=(parameters,))

    f, ax = plt.subplots(1,1)
    ax.plot(cdh1,np.log10(cycb1), 'b', label='CycB nullcline')
    ax.plot(cdh1,np.log10(cycb2), 'r', label='Cdh1 nullcline')
    ax.plot(Cdh1_i,CycB_i,'ko')
    ax.plot(y[:,0], np.log10(y[:,1]), 'k--',alpha=0.5, lw=2.0)
    ax.plot(y[-1,0], np.log10(y[-1,1]), 'ro', lw=2.0)

    for s in solutions:
        ax.plot(cdh1[s], np.log10(cycb1[s]), 'go')
    ax.annotate("G1",(0.9,-1))
    ax.annotate("S/G2/M",(0.01,0.1))
    ax.set_ylabel('log([CycB])')
    ax.set_xlabel('[Cdh1]')
    ax.set_xlim([-0.05,1.01]) 
    ax.set_ylim([-2,1]) 
    ax.legend()

    plt.tight_layout()
    st.pyplot()
    #####################################3
    ### Conclusions
    with open('./markdown/two-variable-conclusion.md','r') as infile:
        sec1text = ''.join(infile.readlines())
    st.markdown(sec1text)

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

def makeFig3(parameters):       
    regenerate = False
    if regenerate == True:
        mvals = np.append(np.linspace(0.06,0.2,100), np.linspace(0.2,0.6,100))
        cdh1 = np.append(np.logspace(-5,-1.,1000), np.logspace(-1.0,0.1,6000))
        hyst = []
        pvals = []
        for m in mvals:
            p = (parameters['k3d'] + parameters['k3dd']*parameters['A'])/(parameters['k4']*m)
            parameters['m'] = m
            cycb1 = cycbnc_fig2(cdh1, parameters)
            cycb2 = cdh1nc_fig2(cdh1, parameters)
            solutions = getRoots(np.log10(cycb1), np.log10(cycb2))
            if len(solutions) >=1 and p>=0.265:
               solutions = [s for s in solutions if s < 0.4] 
            for s in solutions:
                hyst.append(cycb1[s])
                pvals.append(p)
        pval_s, hyst_s = zip(*sorted(zip(pvals, hyst)))
        vals = np.array([[p, h] for p,h in zip(pval_s, hyst_s)])
        df = pd.DataFrame(vals,columns=['p','cycb'])
        df.to_csv('data/hyst.dat')
    with open('./markdown/hysteresis-1.md','r') as infile:
        hyst1 = ''.join(infile.readlines())
    st.markdown(hyst1)
    df = pd.read_csv('data/hyst.dat')
    df = df.sort_values(by='cycb')
    f, ax = plt.subplots(1,1)
    ax.plot(df['p'], df['cycb'],'k',lw=4)
    ax.annotate('G1',(0.15,0.09))
    ax.annotate('S/G2/M',(0.2,0.8))
    t = np.linspace(0,100,1000)
    Cdh1_i = st.slider(label='Cdh1',key='cdh1fig2nc', min_value=0.0, max_value=1.0, value=0.9,step=0.1)
    CycB_i = st.slider(label='log(CycB)', key='cycbfig2nc',min_value=-2., max_value=1., value=-1.1,step=0.1)
    mval = st.slider(label='Mass', key='massfig2nc',min_value=0.1, max_value=0.7, value=0.3,step=0.01)    
    A = st.slider(label='A', key='afig2nc',min_value=0.0, max_value=0.6, value=0.0,step=0.01)    
    parameters['m'] = mval
    parameters['A'] = A
    p = (parameters['k3d'] + parameters['k3dd']*parameters['A'])/(parameters['k4']*mval)
    y = odeint(fig1model, [Cdh1_i, 10**CycB_i],t,args=(parameters,))
    ax.plot(p, y[0,1],'k.')
    ax.annotate('', xytext=(p,  10**CycB_i), xy=(p , y[-1,1]),arrowprops=dict(facecolor='black', arrowstyle='->'),)#, width=0.0025)
    #ax.plot([p for _ in range(len(t))], y[:,1])
    ax.set_xlim([0.0,0.3])
    ax.set_xlabel('p')
    ax.set_ylabel('[CycB]')

    st.pyplot()
    with open('./markdown/hysteresis-2.md','r') as infile:
        hyst2 = ''.join(infile.readlines())
    st.markdown(hyst2)

def cdc20ncfig4(cycb, m, parameters):
    cdc20 = (parameters['k5d'] + parameters['k5dd']*(cycb*m/parameters['J5'])**parameters['n']\
             /(1+(cycb*m/parameters['J5'])**parameters['n']))/parameters['k6']
    return cdc20

def cdh1ncfig4(cycb, p, parameters):
    cdh1 = goldbeter(p, cycb, parameters['J3'], parameters['J4'])
    return cdh1

def cycbncfig4(cdh1, parameters):
    beta = parameters['k1']/parameters['k2dd']
    J = parameters['k2d']/parameters['k2dd']
    cycb = beta/(J + cdh1) 
    return cycb

def threevariable(X, t, args):
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
    k5d = args['k5d']
    k5dd = args['k5dd']
    k6 = args['k6']
    J5 = args['J5']
    n = args['n']
    cdh1, cycb, cdc20t = X
    dcdh1 = ((k3d + k3dd*cdc20t)*(1 - cdh1))/(J3 + 1 - cdh1) - (k4*m*cycb*cdh1)/(J4 + cdh1)
    dcycb = k1- (k2d + k2dd * cdh1)*cycb
    dcdc20t = k5d + k5dd*((cycb*m/J5)**n/(1 + (cycb*m/J5)**n)) - k6*cdc20t
    return(np.array([dcdh1, dcycb, dcdc20t]))

def makeFig4(parameters):
    with open('./markdown/cdh1-activation-1.md','r') as infile:
        cdh1text = ''.join(infile.readlines())
    st.markdown(cdh1text)
    cycbvals = np.append(np.logspace(-6,-3,300),np.logspace(-3,1.1,600))
    # Sliders
    mval = st.selectbox(label='Mass', options=[0.4, 0.8, 1.0])
    Cdh1_i = st.slider(label='Cdh1',key='cdh1fig4tc', min_value=0.0, max_value=1.0, value=0.87,step=0.05)
    Cdc20_i = st.slider(label='Cdc20',key='cdc20fig4tc', min_value=0.0, max_value=1.0, value=0.01,step=0.01)
    CycB_i = st.slider(label='log(CycB)',key='cdh1fig4tc', min_value=0.0, max_value=1.0, value=0.01,step=0.01)
    ## 
    cdc20 = [cdc20ncfig4(c, mval, parameters) for c in cycbvals ]
    pvals = (parameters['k3d'] + parameters['k3dd']*np.array(cdc20))/(parameters['k4']*mval)
    cdh1 = [cdh1ncfig4(c, p, parameters) for c,p in zip(cycbvals, pvals)]
    cycb = cycbnc_fig2(cdh1, parameters)
    #mval = 0.8
    regenerate = False
    cdh1vals = np.linspace(0.0,1.0, 600)
    if regenerate == True:
        hyst = []
        pvals = []
        cdc20vals = []
        for c in cdh1vals:
            parameters['m'] = mval
            parameters['A'] = c
            cycb1 = cycbnc_fig2(cdh1, parameters)
            cycb2 = cdh1nc_fig2(cdh1, parameters)
            solutions = getRoots(np.log10(cycb1), np.log10(cycb2))
            for s in solutions:
                hyst.append(cycb1[s])
                cdc20vals.append(c)
        vals = np.array([[c, h] for c,h in zip(cdc20vals, hyst)])
        df = pd.DataFrame(vals,columns=['cdc20','cycb'])
        df.to_csv('data/hyst-cdc20-mid-m.dat')

    f, ax = plt.subplots(1,1)
    ax.plot(cdc20, cycbvals, 'k--',label='Cdh1 nullcline')
    ax.set_xlabel('Cdc20_T')
    ax.set_ylabel('CycB')
    ax.set_ylim(0,1.0)

    x0 = [Cdh1_i, CycB_i, Cdc20_i]
    tmax = 50
    stepsize=0.01
    t = np.linspace(0 ,tmax, int(tmax/stepsize))
    parameters['m'] = mval

    #y = integrate(threevariable, x0, t, parameters, stepsize=stepsize)
    y = odeint(threevariable, x0, t, args=(parameters,))

    fname = './data/hyst-cdc20-lo-m.dat'
    xmax = 0.5
    settings = {0.4:{'fname':'./data/hyst-cdc20-lo-m.dat','xmax':0.5},
                0.8:{'fname':'./data/hyst-cdc20-mid-m.dat','xmax':1.0},
                1.0:{'fname':'./data/hyst-cdc20-hi-m.dat','xmax':1.0}}

    df = pd.read_csv(settings[mval]['fname'])
    ax.plot(df['cdc20'], df['cycb'], 'k.', label='CycB nullcline')
    ax.plot(y[:,2], y[:,1],'r--')
    ax.set_xlim([0,settings[mval]['xmax']])
    ax.annotate('S/G2/M',(0.6,0.15))
    ax.annotate('G1',(0.051,0.05))
    ax.legend()
    st.pyplot()
    plt.close()

    with open('./markdown/cdh1-activation-2.md','r') as infile:
        cdh1conclude = ''.join(infile.readlines())
    st.markdown(cdh1conclude)

    # f, ax = plt.subplots(1,1)
    # ax.plot(t, y[:,0],label='cdh1')
    # ax.plot(t, y[:,1],label='cycb')
    # ax.plot(t, y[:,2],label='cdc20')
    # ax.legend()
    # st.pyplot()

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
    with open('./markdown/primitive.md', 'r') as infile:
        primitive = ''.join(infile.readlines())
    st.markdown(primitive)
    x0 = [1.0, 0.5,1.5, 1.4, 0.7, 0.6]
    stepsize = 0.01
    tmax = 160
    t= np.linspace(0, tmax, int(tmax/stepsize))
    #y = odeint(fullmodel,x0, t, args=(parameters,))
    y = integrate(fullmodel, x0, t, parameters, stepsize=stepsize)
    f , ax = plt.subplots(3,1)#, figsize=(1,3))
    ax[0].plot(t,y[:,5], label='m')
    ax[0].legend()
    ax[1].plot(t,y[:,0], 'k',label='Cdh1')
    axc = ax[1].twinx()
    axc.plot(t,y[:,1], 'r--',label='CycB')
    axc.set_ylim(0.,0.7)
    axc.legend()
    ax[1].legend()
    ax[2].plot(t,y[:,2], label='Cdc20T')
    ax[2].plot(t,y[:,3], label='Cdc20A')
    ax[2].plot(t,y[:,4], label='IEP')        
    ax[2].set_ylim([0,2.0])
    ax[2].legend()
    plt.tight_layout()
    st.pyplot()

def yeastmodel(X, t, args):
    parameters = args
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
    k11 = args['k11']
    k12d = args['k12d']
    k12dd = args['k12dd']
    k12ddd = args['k12ddd']
    k13 = args['k13']
    k14 = args['k14']
    k15d = args['k15d']
    k15dd = args['k15dd']
    k16d = args['k16d']
    k16dd = args['k16dd']
    Keq = args['Keq']
    J15 = args['J15']
    J16 = args['J16']
    cdh1, cycbt, ckit, sk, cdc20t, cdc20a, iep, m = X
    cycb = cycbt - trimer(cycbt, ckit, parameters)
    dcdh1 = ((k3d + k3dd*cdc20a)*(1 - cdh1))/(J3 + 1 - cdh1) - ((k4d*sk + k4*m*cycb)*cdh1)/(J4 + cdh1)
    dcycbt = k1- (k2d + k2dd*cdh1 + k2ddd*cdc20a)*cycbt    
    dckit = k11 - (k12d + k12dd*sk + k12ddd*m*cycb)*ckit
    dsk = k13*goldbeter(k15d*m + k15dd*sk, k16d + k16dd*m*cycb, J15, J16) - k14*sk
    dcdc20t = k5d + k5dd*( (cycb*m/J5)**n /(1+ (cycb*(m/J5))**n )) - k6*cdc20t
    dcdc20a = (k7*iep*(cdc20t-cdc20a)/(J7 + cdc20t - cdc20a)) - (k8*Mad*cdc20a)/(J8+cdc20a) - k6*cdc20a
    diep = k9*m*cycb*(1-iep) - k10*iep
    dm = mu*m*(1-m/mstar)
    return np.array(([dcdh1, dcycbt, dckit, dsk, dcdc20t, dcdc20a, diep, dm]))

def cdh1ncfig7(cycbt,ckit, cdc20a, m, sk, parameters):
    cycb = cycbt - trimer(cycbt, ckit, parameters)
    cdh1 = goldbeter(parameters['k3d'] + parameters['k3dd']*cdc20a,
                     parameters['k4d']*sk + parameters['k4']*m*cycb,
                     parameters['J3'],
                     parameters['J4'])
    return cdh1
def cycbncfig7(cdh1, cdc20a, parameters):
    cycbt = parameters['k1']/(parameters['k2d'] + parameters['k2dd']*cdh1 + parameters['k2ddd']*cdc20a)
    return cycbt

def ckitncfig7(cycbt, m, sk, parameters):
    cycb = cycbt - trimer(cycbt, )
    ckit = parameters['k11']/(parameters['k12d'] + parameters['k12dd']*sk + parameters['k12ddd']*m*cycb)
    return ckit

def trimer(cycbt, ckit, parameters):
    tri = (2*cycbt*ckit)/(cycbt + ckit + 1./parameters['Keq'] +\
                             np.sqrt((cycbt + ckit + 1./parameters['Keq'])**2 -4*cycbt*ckit))
    return tri

# def makeYeastPP(parameters):
#     # TODO: understand and implement nullclines in fig 7
    # mval = st.slider('mass', min_value=0.1, max_value=1.0,step=0.1,value=0.2)
    # mval = 1.0
    # parameters['m'] = mval
    # skval = st.slider('SK', min_value=-2., max_value=0.0,step=0.5,value=-2.)    
    # # cycbt values over which to compute cdh1 nullcline
    # cycbvals = np.append(np.logspace(-6,-3,300),np.logspace(-3,0.0,600))
    # ckitnc = [ckitncfig7(c, mval, 10**skval, parameters) for c in cycbvals]
    # ckitvals = np.linspace(0., 1., 100)
    # cdc20vals = [cdc20ncfig4(c, mval, parameters) for c in cycbvals]
    # pvals = (parameters['k3d'] + parameters['k3dd']*np.array(cdc20vals))/(parameters['k4']*mval)
    # cdh1 = [cdh1ncfig4(c, p, parameters) for c,p in zip(cycbvals, pvals)]
    # cycbnc = []
    # ckitx = []
    # for ckit in ckitvals:
    #     cdh1vals = [cdh1ncfig7(c, ckit, cdc20, mval, skval, parameters) for c,cdc20 in zip(cycbvals,cdc20vals)]
    #     cycb1 = cycbncfig7(np.array(cdh1vals), np.array(cdc20vals), parameters)
    #     cycb2 = cdh1nc_fig2(np.array(cdh1vals), parameters)
    #     roots = getRoots(np.log10(cycb1), np.log10(cycb2))
    #     for r in roots:
    #         cycbnc.append(r)
    #         ckitx.append(ckit)
    # #cycb = [cycbncfig7(c, 0, parameters) for c in cdh1vals]
    # f, ax = plt.subplots(1,1)
    # ax.plot(ckitnc, cycbvals)
    # ax.plot(ckitx, cycbnc,'k.')
    # ax.plot(ckitnc, [cycbt - trimer(cycbt,ckit,parameters ) for cycbt, ckit in zip(cycbvals, ckitnc)])
    # ax.set_title(str(len(cycbnc)))
    # # ax.set_xlim(0,1)
    # # ax.set_ylim(0,1)
    # ax.set_ylabel('[CycB$_T$]')    
    # ax.set_xlabel('[CKI$_T$]')

def makeFig7(parameters):
    ## time courses
    parameters['mu'] = 0.005
    x0 = [1.0, 0.01, 0.01, 0.04, 0.2,0.1, 0.2, 1.0]
    tmax = 300
    stepsize = 0.01
    t = np.linspace(0, tmax , int(tmax/stepsize))
    image = Image.open('data/yeast-model.png')
    st.image(image, caption="Wiring diagram of the yeast cell cycle, "\
             "taken from Tyson and Novak, 2001.", width=500)#use_column_width=True)
    with open('markdown/yeast-model-1.md', 'r') as infile:
        yeastdescription1 = ''.join(infile.readlines())

    st.markdown(yeastdescription1)
    strain  = st.selectbox(label='Select a yeast strain', options=['WT', 'SK-deletion','SK, CKI double deletion'])

    if strain == 'WT':
        pars = dict(parameters)
    elif strain == 'SK-deletion':
        pars = dict(parameters)
        pars['k13'] = 0
    elif strain == 'SK, CKI double deletion':
        pars = dict(parameters)
        pars['k13'] = 0
        pars['k11'] = 0
    y = integrate(yeastmodel,x0, t, pars, massindex=7, stepsize=stepsize)
    f, ax = plt.subplots(2,1)
    ax[0].plot(t, y[:,0],label='Cdh1')
    ax[0].plot(t, y[:,1],label='Cycb$_T$')
    ax[1].plot(t, y[:,2],label='Cki$_T$')
    ax[1].plot(t, y[:,3],label='SK')
    ax[1].plot(t, y[:,5],label='Cdc20$_A$')
    ax[0].legend()
    ax[1].legend()
    ax[0].set_title(strain)
    st.pyplot()
    with open('markdown/yeast-model-2.md', 'r') as infile:
        yeastdescription2 = ''.join(infile.readlines())
    st.markdown(yeastdescription2)

def makeIntroPage():
    # gifpath = r''' <img src="./data/budding.gif">''' # 
    # st.write(gifpath, unsafe_allow_html=True)
    with open('markdown/intro.md','r') as infile:
        introtext = ''.join(infile.readlines())
    #with open
    #st.image()
    st.markdown(introtext)

def makePages(parameters):
   page = st.sidebar.selectbox('Jump to...',['Introduction',
                                             'Cdh1-CycB Antagonism',
                                             'Hysteresis in transitions',
                                             'Regulation of Cdh1/APC',
                                             'A primitive model',
                                             'The yeast cell cycle', 'Summary'])
   if page == 'Introduction':
       st.header('Introduction')
       makeIntroPage()
   if page == 'Cdh1-CycB Antagonism':
       # st.header('A simplified model of CycB/Cdk1-Cdh1/APC antagonism')
       makeFig2(parameters)
   if page == 'Hysteresis in transitions':
       st.header('Hystersis underlies cell state transitions')
       makeFig3(parameters)
   if page == 'Regulation of Cdh1/APC':
       st.header('Activating the Cdh1/APC')
       makeFig4(parameters)
   if page == 'A primitive model':
       st.header('Primitive Model')
       plottimecourses(parameters)
   if page == 'The yeast cell cycle':
       st.header('The budding yeast cell cycle')
       makeFig7(parameters)
   if page == 'Summary':
       st.header("...But there's more!")
       summary()

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
    makePages(parameters)
        #plottimecourses(parameters)
    # if page == 'test':
    #     st.header('test')
        #test()
if __name__ == '__main__':
    main()
