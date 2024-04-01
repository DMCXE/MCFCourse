import numpy as np
import matplotlib.pyplot as plt
from CF88 import O16O16,DT,O16mbar,mbar2rate

kb = 1.38e-23
Na = 6.02e23

# K -> eV
def K2eV(K):
    return K/11604.51812

# 聚变功率P_fus
# def P_fus(n1,n2 ,T9,):
#     Q,f  = O16O16(T9)
#     return n1*n2*f*Q

def ntau(T9,delta,Zi,f,Yp):
    '''
    T9:    温度,单位GK
    delta: 相同反应物 delta=1,不同反应物 delta=0
    Zi:    离子电荷数
    f:     反应速率   ， 来着CF88
    Yp:    带电产物能量，MeV
    '''
    kbev = 1/11604.51812  #eV/K
    
    return 1.5*kbev*Zi*(1+Zi)*(1+delta)*T9*1e9*1e-3*Na/((1e-6)*f*(Yp*1e3))

def ntau_with_Q(T9e,T9i,delta,Zi,f,Q,fion,x1,x2,Qsci,ni,ne):
    '''
    T9i:   i,e表示为离子或电子，温度,单位GK
    delta: 相同反应物 delta=1,不同反应物 delta=0
    Zi:    离子电荷数
    f:     反应速率   ， 来着CF88
    Q:     反应Q值，MeV
    '''
    kbev  = 1/11604.51812    #eV/K
    kbTe  = kbev*T9e*1e9*1e-3 #keV
    kbTi  = kbev*T9i*1e9*1e-3
    C_B   = 5.34*1e-37
    me    = 9.10938215*1e-31
    c     = 3e8
    e     = 1.6*1e-19
    Zeff  = 2*ni*Zi**2 / ne  # 有效电荷数 
    Zeff  = 1.0
    mec2  = (me * c**2)/e * 1e-3 #kev
    geff  = (1+0.7936*kbTe/mec2 + 1.874*(kbTe/mec2)**2)\
            +(3/np.sqrt(2))*(kbTe/mec2)/Zeff

    up    = 1.5*(kbTe+kbTi/Zi)
    down1 = ((1/Qsci+fion)/(1+delta))*(x1*x2/(Zi**2))*(Q*1e3)*((1e-6)*f/Na)
    down2 = -C_B*np.sqrt(kbTe)*Zeff*geff
    return up/(down1+down2)

def p_fusion(ni,f,Q,delta):
    e = 1.6*1e-19
    return (1/(1+delta))*(ni*ni)* ((1e-6)*f/Na) *(Q*1e6*e) 

def p_bremsstrahlung(T9,Zi,ni,ne):
    '''
    W m-3
    '''
    C_B  = 5.34*1e-37
    me   = 9.10938215*1e-31
    c    = 3e8
    e    = 1.6*1e-19
    mec2 = (me * c**2)/e * 1e-3 #kev
    
    kbev = 1/11604.51812    #eV/K
    kbTe = kbev*T9*1e9*1e-3 #keV
    Zeff = 2*ni*Zi**2 / ne  # 有效电荷数 
    Zeff = 1
    return C_B* ne**2 * np.sqrt(kbTe)\
            *(Zeff*(1+0.7936*kbTe/mec2 + 1.874*(kbTe/mec2)**2)+(3/np.sqrt(2))*(kbTe/mec2))




if __name__ == '__main__':
    plt.figure()
    T9  = np.linspace(0.1, 10, 1000,endpoint=True)
    print('eV',K2eV(1e9)/1e3)


    '''DT'''
    Yp = 3.52
    delta = 0
    Zi = 1
    Q,f = DT(T9)
    nt = ntau(T9,delta,Zi,f,Yp)

    ni = 1e21
    ne = ni
    for qsci in np.array([1,10,100,1000,10000,1e5,1e6,1e100]):
        ntQ = ntau_with_Q(T9e=T9,T9i=T9,delta=delta,Zi = 2,
                        f=f,Q=Q,fion = Yp/Q, x1=0.5,x2=0.5,Qsci=qsci,ni=ni,ne=ne)
        #对数坐标
        #plt.loglog(K2eV(T9*1E9)/1E3,ntQ,label = 'DT in Q={}'.format(qsci))





    #plt.loglog(K2eV(T9*1E9)/1E3,nt*K2eV(T9*1E9)/1E3,label = 'DT')


    '''O16'''
    Q,f = O16O16(T9)
    Yp = Q-9.593
    Zi = 2
    delta = 1
    Qsci  = 10000
    nt = ntau(T9,delta,Zi,f,Yp)
    ni = 2e20
    ne = 2*ni
    # for qsci in np.array([1,10,100,1000,10000,1e5,1e6,1e100]):
    #     ntQ = ntau_with_Q(T9e=T9,T9i=T9,delta=delta,Zi = 2,
    #                     f=f,Q=Q,fion = Yp/Q, x1=1,x2=1,Qsci=qsci,ni=ni,ne=ne)
        
    #     #对数坐标
    #     plt.loglog(K2eV(T9*1E9)/1E3,ntQ,label = 'O16O16 in Q={}'.format(qsci))
    #     print(ntQ)

    ni = 2e27
    ne = 2*ni
    W = p_bremsstrahlung(T9,Zi,ni,ne)
    P = p_fusion(ni,f,Q,delta)
    # plt.loglog(K2eV(T9*1E9 )/1E3,W,label='bremsstrahlung')
    # plt.loglog(K2eV(T9*1E9 )/1E3,P,label='fusion')
    plt.plot(K2eV(T9*1E9 )/1E3,f/Na,label='fusion rate')

    MeV, MB = O16mbar()
    mev = np.linspace(MeV[0],MeV[-1],1000,endpoint=True)
    rate = np.array([])
    #plt.plot(K2eV(T9*1E9)/1E3,f)
    # for i in mev:
    #     print(i)
    #     rate = np.append(rate,mbar2rate(16,16,i,i,MeV,MB))
    # plt.loglog(K2eV(T9*1E9 )/1E3,rate)

 
    
    
    plt.xlabel('KeV')
    plt.ylabel('nt')
    plt.legend()
    plt.grid()
    plt.title('nt')
    plt.show()
