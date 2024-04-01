import numpy as np
import matplotlib.pyplot as plt


def power_bremsstrahlung(Te,Zi,ne):
    '''
    Te: 电子温度，单位keV
    Zi: 离子电荷数,同种离子反应密度自然相同
    ne: 电子密度，单位m^-3, ne = Zi ni
    '''
    C_B = 5.34*1e-37
    kbTe = Te #keV
    e    = 1.6*1e-19
    c    = 3*1e8
    me   = 9.10938215*1e-31
    mec2 = (me * c**2)/e * 1e-3 #kev
    Zeff = Zi
    return C_B* ne**2 * np.sqrt(kbTe)\
            *(Zeff*(1+0.7936*kbTe/mec2 + 1.874*(kbTe/mec2)**2)+(3/np.sqrt(2))*(kbTe/mec2))

def power_fusion(ni,rate,Y,delta):
    '''
    Y:MeV
    '''
    e    = 1.6*1e-19
    return ni**2 * rate * Y*1e6 *e /(1+delta)


def lawson_with_Q(Q,Te,Ti,Zi,rate,Y,fion,delta):
    '''
    只考虑聚变反应物是同种或同位素的情况
    '''
    kbev  = 1/11604.51812    #eV/K
    kbTe  = Te #keV
    kbTi  = Ti
    C_B   = 5.34*1e-37
    me    = 9.10938215*1e-31
    c     = 3e8
    e     = 1.6*1e-19
    Zeff  = Zi
    mec2  = (me * c**2)/e * 1e-3 #kev
    geff  = (1+0.7936*kbTe/mec2 + 1.874*(kbTe/mec2)**2)\
            +(3/np.sqrt(2))*(kbTe/mec2)/Zeff
    x1 = 0.5
    x2 = 0.5
    up    = 1.5*(kbTe+kbTi/Zi)
    down1 = ((1/Q+fion)/(1+delta))*(x1*x2/(Zi**2))*(Y*1e3)*(rate)
    down2 = -C_B*np.sqrt(kbTe)*Zeff*geff
    return up/(down1+down2)

def compare_power_88(ne):
    import CF88
    T9 = np.linspace(1e-1,1e1,1000)
    kbev  = 1/11604.51812    #eV/K
    T9kev = T9*1e9*kbev/1e3
    Y,f = CF88.O16O16(T9)
    Zi = 2
    Na = 6.02e23
    p88 = power_fusion(ne/2,f*1e-6/Na,Y,1)

    plt.loglog(T9kev,p88,label = 'CF88_O')

def compare_power_O16O16SI(ne):
    import EXFOR
    mr = 0.5*15.9994*1.6726*1e-27
    '''
    EXFOR
    '''
    MEVA,mBA = EXFOR.O16O16SI28_A0381006()
    MEVC,mBC = EXFOR.O16O16SI28_C1269002()
    MEVE,mBE = EXFOR.O16O16SI28_E1147029()

    KEVA = 1e3  * MEVA
    KEVC = 1e3  * MEVC
    KEVE = 1e3  * MEVE
    BA   = 1e-3 * mBA
    BC   = 1e-3 * mBC
    BE   = 1e-3 * mBE

    TrA = np.linspace(KEVA[0],KEVA[-1],100)
    TrC = np.linspace(KEVC[0],KEVC[-1],100)
    TrE = np.linspace(KEVE[0],KEVE[-1],100)
    T = np.linspace(1e1,20e4,1000)
    TrA = T
    TrC = T
    TrE = T

    rates_A  = EXFOR.B2R(KEVA,BA,TrA,mr)
    rates_C  = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_E  = EXFOR.B2R(KEVE,BE,TrE,mr)

    Zi = 2
    Y  = 9.593
    p_A = power_fusion(ne/Zi,rates_A,Y,1)
    p_C = power_fusion(ne/Zi,rates_C,Y,1)
    p_E = power_fusion(ne/Zi,rates_E,Y,1)

    plt.loglog(TrA,p_A,label = 'O16O16SI28_A0381006')
    plt.loglog(TrC,p_C,label = 'O16O16SI28_C1269002')
    plt.loglog(TrE,p_E,label = 'O16O16SI28_E1147029')
    
def compare_power_O16fu(ne):
    import EXFOR
    mr = 0.5*15.9994*1.6726*1e-27
    EXFOR.O16O16FUS_C1269009()
    EXFOR.O16O16FUS_D0769002()
    EXFOR.O16O16FUS_O1531004()
    MEVC,mBC = EXFOR.O16O16FUS_C1269009()
    MEVD,mBD = EXFOR.O16O16FUS_D0769002()
    MEVO,mBO = EXFOR.O16O16FUS_O1531004()

    KEVC = 1e3  * MEVC
    KEVD = 1e3  * MEVD
    KEVO = 1e3  * MEVO

    BC   = 1e-3 * mBC
    BD   = 1e-3 * mBD
    BO   = 1e-3 * mBO

    TrC = np.linspace(KEVC[0],KEVC[-1],100)
    TrD = np.linspace(KEVD[0],KEVD[-1],100)
    TrO = np.linspace(KEVO[0],KEVO[-1],100)

    T = np.linspace(1e1,20e4,1000)
    TrC = T
    TrD = T
    TrO  =T

    rates_C = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D = EXFOR.B2R(KEVD,BD,TrD,mr)
    rates_O = EXFOR.B2R(KEVO,BO,TrO,mr)

    Zi = 2
    Y  = 13.431 #!!!!!!!!!!!!!!!!
    p_C = power_fusion(ne/Zi,rates_C,Y,1)
    p_D = power_fusion(ne/Zi,rates_D,Y,1)
    p_O = power_fusion(ne/Zi,rates_O,Y,1)

    plt.loglog(TrC,p_C,label = 'O16O16FUS_C1269009')
    plt.loglog(TrD,p_D,label = 'O16O16FUS_D0769002')
    plt.loglog(TrO,p_O,label = 'O16O16FUS_O1531004')

def compare_power_SIfu(ne):
    import EXFOR
    mr = 0.5*28.0855*1.6726*1e-27
    '''
    EXFOR
    '''
    MEVC  ,mBC   = EXFOR.SI28SI28FUS_C0415002()
    MEVD75,mBD75 = EXFOR.SI28SI28FUS_D0750002()
    MEVD77,mBD77 = EXFOR.SI28SI28FUS_D0772002()
    MEVE  ,mBE   = EXFOR.SI28SI28FUS_E1435016()

    KEVC   = 1e3  * MEVC
    KEVD75 = 1e3  * MEVD75
    KEVD77 = 1e3  * MEVD77
    KEVE   = 1e3  * MEVE

    BC     = 1e-3 * mBC
    BD75   = 1e-3 * mBD75
    BD77   = 1e-3 * mBD77
    BE     = 1e-3 * mBE

    
    TrC   = np.linspace(KEVC[0],KEVC[-1],100)
    TrD75 = np.linspace(KEVD75[0],KEVD75[-1],100)
    TrD77 = np.linspace(KEVD77[0],KEVD77[-1],100)
    TrE   = np.linspace(KEVE[0],KEVE[-1],100)

    T = np.linspace(1e-3,20e4,1000)
    TrC = T
    TrD75 = T
    TrD77 = T
    TrE  = T

    rates_C    = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D75  = EXFOR.B2R(KEVD75,BD75,TrD75,mr)
    rates_D77  = EXFOR.B2R(KEVD77,BD77,TrD77,mr)
    rates_E    = EXFOR.B2R(KEVE,BE,TrE,mr)
    Zi = 2
    Y  = 49.431 #!!!!!!!!!!!!!!!!
    delta = 1
    p_C = power_fusion(ne/Zi,rates_C,Y,delta)
    p_D75 = power_fusion(ne/Zi,rates_D75,Y,delta)
    p_D77 = power_fusion(ne/Zi,rates_D77,Y,delta)
    p_E = power_fusion(ne/Zi,rates_E,Y,delta)

    
    plt.plot(TrC,p_C,label = 'SI28SI28FUS_C1269002')
    plt.plot(TrD75,p_D75,label = 'SI28SI28FUS_D0750002')
    plt.plot(TrD77,p_D77,label = 'SI28SI28FUS_D0772002')
    plt.plot(TrE,p_E,label = 'SI28SI28FUS_E1435016')




def compare_lawson_CF88_O(Q):
    import CF88
    T9 = np.linspace(1e-1,1e1,1000)
    kbev  = 1/11604.51812    #eV/K
    T9kev = T9*1e9*kbev/1e3
    Y,f = CF88.O16O16(T9)
    Zi = 2
    Na = 6.02e23
    nt = lawson_with_Q(Q,T9kev,T9kev,Zi,f*1e-6/Na,Y,1,1)
    #plt.figure()
    plt.loglog(T9kev,nt,label = 'CF88_O')
   # plt.show()
    
def comapre_lawson_O16O16SI(Q):
    import EXFOR
    mr = 0.5*15.9994*1.6726*1e-27
    '''
    EXFOR
    '''
    MEVA,mBA = EXFOR.O16O16SI28_A0381006()
    MEVC,mBC = EXFOR.O16O16SI28_C1269002()
    MEVE,mBE = EXFOR.O16O16SI28_E1147029()

    KEVA = 1e3  * MEVA
    KEVC = 1e3  * MEVC
    KEVE = 1e3  * MEVE
    BA   = 1e-3 * mBA
    BC   = 1e-3 * mBC
    BE   = 1e-3 * mBE

    TrA = np.linspace(KEVA[0],KEVA[-1],100)
    TrC = np.linspace(KEVC[0],KEVC[-1],100)
    TrE = np.linspace(KEVE[0],KEVE[-1],100)

    T = np.linspace(1e3,20e4,1000)
    TrA = T
    TrC = T
    TrE = T

    rates_A  = EXFOR.B2R(KEVA,BA,TrA,mr)
    rates_C  = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_E  = EXFOR.B2R(KEVE,BE,TrE,mr)

    Zi = 2
    Y  = 9.593
    #Q  = 1
    nt_A = lawson_with_Q(Q,TrA,TrA,Zi,rates_A,Y,1,1)
    nt_C = lawson_with_Q(Q,TrC,TrC,Zi,rates_C,Y,1,1)
    nt_E = lawson_with_Q(Q,TrE,TrE,Zi,rates_E,Y,1,1)

    plt.loglog(TrA,nt_A,label = 'O16O16SI28_A0381006')
    plt.loglog(TrC,nt_C,label = 'O16O16SI28_C1269002')
    plt.loglog(TrE,nt_E,label = 'O16O16SI28_E1147029')

def compare_lawson_O16fus(Q):
    import EXFOR
    mr = 0.5*15.9994*1.6726*1e-27
    EXFOR.O16O16FUS_C1269009()
    EXFOR.O16O16FUS_D0769002()
    EXFOR.O16O16FUS_O1531004()
    MEVC,mBC = EXFOR.O16O16FUS_C1269009()
    MEVD,mBD = EXFOR.O16O16FUS_D0769002()
    MEVO,mBO = EXFOR.O16O16FUS_O1531004()

    KEVC = 1e3  * MEVC
    KEVD = 1e3  * MEVD
    KEVO = 1e3  * MEVO

    BC   = 1e-3 * mBC
    BD   = 1e-3 * mBD
    BO   = 1e-3 * mBO

    TrC = np.linspace(KEVC[0],KEVC[-1],100)
    TrD = np.linspace(KEVD[0],KEVD[-1],100)
    TrO = np.linspace(KEVO[0],KEVO[-1],100)

    T = np.linspace(1e1,20e4,1000)
    TrC = T
    TrD = T
    TrO  =T

    rates_C = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D = EXFOR.B2R(KEVD,BD,TrD,mr)
    rates_O = EXFOR.B2R(KEVO,BO,TrO,mr)

    Zi = 2
    Y  = 13.431 #!!!!!!!!!!!!!!!!
    #Q  = 1
    nt_C = lawson_with_Q(Q,TrC,TrC,Zi,rates_C,Y,1,1)
    nt_D = lawson_with_Q(Q,TrD,TrD,Zi,rates_D,Y,1,1)
    nt_O = lawson_with_Q(Q,TrO,TrO,Zi,rates_O,Y,1,1)
    plt.loglog(TrC,nt_C,label = 'O16O16FUS_C1269009')
    plt.loglog(TrD,nt_D,label = 'O16O16FUS_D0769002')
    plt.loglog(TrO,nt_O,label = 'O16O16FUS_O1531004')

def compare_lawson_SIfus(Q):
    import EXFOR
    mr = 0.5*28.0855*1.6726*1e-27
    '''
    EXFOR
    '''
    MEVC  ,mBC   = EXFOR.SI28SI28FUS_C0415002()
    MEVD75,mBD75 = EXFOR.SI28SI28FUS_D0750002()
    MEVD77,mBD77 = EXFOR.SI28SI28FUS_D0772002()
    MEVE  ,mBE   = EXFOR.SI28SI28FUS_E1435016()

    KEVC   = 1e3  * MEVC
    KEVD75 = 1e3  * MEVD75
    KEVD77 = 1e3  * MEVD77
    KEVE   = 1e3  * MEVE

    BC     = 1e-3 * mBC
    BD75   = 1e-3 * mBD75
    BD77   = 1e-3 * mBD77
    BE     = 1e-3 * mBE

    
    TrC   = np.linspace(KEVC[0],KEVC[-1],100)
    TrD75 = np.linspace(KEVD75[0],KEVD75[-1],100)
    TrD77 = np.linspace(KEVD77[0],KEVD77[-1],100)
    TrE   = np.linspace(KEVE[0],KEVE[-1],100)

    T = np.linspace(1e3,20e4,1000)
    TrC = T
    TrD75 = T
    TrD77 = T
    TrE  = T

    rates_C    = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D75  = EXFOR.B2R(KEVD75,BD75,TrD75,mr)
    rates_D77  = EXFOR.B2R(KEVD77,BD77,TrD77,mr)
    rates_E    = EXFOR.B2R(KEVE,BE,TrE,mr)

    Zi = 4
    Y  = 13.525625
    #Y  = 49

    nt_C = lawson_with_Q(Q,TrC,TrC,Zi,rates_C,Y,1,1)
    nt_D75 = lawson_with_Q(Q,TrD75,TrD75,Zi,rates_D75,Y,1,1)
    nt_D77 = lawson_with_Q(Q,TrD77,TrD77,Zi,rates_D77,Y,1,1)
    nt_E = lawson_with_Q(Q,TrE,TrE,Zi,rates_E,Y,1,1)

    plt.loglog(TrC,nt_C,label = 'SI28SI28FUS_C0415002')
    plt.loglog(TrD75,nt_D75,label = 'SI28SI28FUS_D0750002')
    plt.loglog(TrD77,nt_D77,label = 'SI28SI28FUS_D0772002')
    plt.loglog(TrE,nt_E,label = 'SI28SI28FUS_E1435016')








if __name__ == '__main__':
    import EXFOR
    import CF88
    MEV,mB = EXFOR.O16O16SI28_E1147029()
    KEV    = 1e3  * MEV
    B      = 1e-3 * mB

    KEV,B = EXFOR.DT()
    
    Ti = np.linspace(KEV[0]*0.5,KEV[-1],1000)
    ne = 1e20
    Zi = 1
    mr = 1.2*1.66*1e-27 #DT体系的相对质量m,kg
    rates  = EXFOR.B2R(KEV,B,Ti,mr)
    P  = power_bremsstrahlung(Ti,Zi,ne)
    Pf = power_fusion(ne/Zi,rates,3.52,1)
    nt = lawson_with_Q(1000000000000,Ti,Ti,Zi,rates,17.571,3.52/17.571,1)
    plt.figure()
    #plt.loglog(Ti,P/Zi)
    #plt.loglog(Ti,Pf)
    plt.loglog(Ti,nt)
    #plt.ylim(bottom = 1e3)
    plt.show()


    Q = 1
    plt.figure()
    compare_lawson_CF88_O(Q)
    comapre_lawson_O16O16SI(Q)
    compare_lawson_O16fus(Q)
    compare_lawson_SIfus(Q)
    plt.legend()
    plt.xlim(left=3e2)
    plt.grid()
    plt.xlabel('Energy KeV')
    plt.ylabel('neτE m-3s')
    #plt.savefig('compareQ1.png')
    plt.show()

    ne =  1e21
    plt.figure(figsize=(16,9))
    #compare_power_88(ne)
    compare_power_O16O16SI(ne)
    compare_power_O16fu(ne)
    compare_power_SIfu(ne)
    pb = power_bremsstrahlung(np.linspace(1e1,10e4,1000),2,ne)
    plt.loglog(np.linspace(1e1,10e4,1000),pb,'--'
               ,c='black',label = 'Bremsstrahlung')
    plt.legend()
    plt.grid()
    plt.xlabel('Energy KeV')
    #plt.ylabel('<σv> W/m3')
    plt.ylabel('P W/m3')
    plt.ylim(bottom = 1e1)
    #plt.savefig('compare_power.png')
    plt.show()