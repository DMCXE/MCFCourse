import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio

import CF88
import EXFOR

'''
数据验证CF88/EXFOR/h.s.Xie的数据准确性
'''


'''
DT——CF88/EXFOR/H.S.Xie
'''

def verify_DT():
     
    '''
    EXFOR_1985
    '''
    KEV,B = EXFOR.DT()
    Tr = np.linspace(KEV[0],KEV[-1],100,endpoint=True)
    mr = 1.2*1.66*1e-27 #DT体系的相对质量m,kg
    rates = EXFOR.B2R(KEV,B,Tr,mr)

    '''
    EXFOR_1947_A1172.003
    '''
    MEV,mB = EXFOR.DT2()
    KEV2 = 1e3*MEV
    B2   = 1e-3*mB
    rates2 = EXFOR.B2R(KEV2,B2,Tr,mr)

    '''
    H.S.Xie from book
    '''
    data = scio.loadmat('xhsDTdata.mat')
    Teff = data['Teff']
    sigmadt = data['sgmv1']
    Teff = Teff[Teff<=KEV[-1]]
    Teff = Teff[Teff>=KEV[0]]
    sgmv = sigmadt[0][0:len(Teff)]

    '''
    CF88
    '''
    Na   = 6.022140857*1e23
    kbev = 1/11604.51812    #eV/K
    T9 = Tr*1e3/kbev/1e9
    Q,f = CF88.DT(T9)
    f /= Na*1e6

    plt.figure()
    plt.plot(Tr,rates,label='EXFOR_1985')
    plt.plot(Tr,rates2,label='EXFOR_1947')
    plt.plot(Teff,sgmv,label = 'H.S.Xie')
    plt.plot(Tr,f,label = 'CF88')
    plt.grid()
    plt.legend()
    plt.xlabel('KeV')
    plt.ylabel('<σv> m3/s')
    plt.savefig('Verify_DT.png')
    plt.show()


def verify_O16O16Si28(save=True,show=True):
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

    rates_A  = EXFOR.B2R(KEVA,BA,TrA,mr)
    rates_C  = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_E  = EXFOR.B2R(KEVE,BE,TrE,mr)

    '''
    CF88
    '''
    Na   = 6.022140857*1e23
    kbev = 1/11604.51812    #eV/K
    T9 = np.linspace(0.1,10,100)
    Q,f = CF88.O16O16(T9)
    f /= Na*1e6

    plt.figure()
    plt.plot(TrA/1e3,rates_A,label = 'EXFOR_A0381006')
    plt.plot(TrC/1e3,rates_C,label = 'EXFOR_C1269002')
    plt.plot(TrE/1e3,rates_E,label = 'EXFOR_E1147029')

    plt.plot(T9*1e9*kbev/1e6,f,label = 'CF88')

    plt.legend()
    plt.grid()
    plt.xlabel('MeV')
    plt.ylabel('<σv> m3/s')
    if save:
        plt.savefig('O16(O16,A)SI28.png')
    if show:
     plt.show()

def verify_O16O16FUS(save=True,show=True):
    mr = 0.5*15.9994*1.6726*1e-27
    '''
    EXFOR
    '''
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

    rates_C = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D = EXFOR.B2R(KEVD,BD,TrD,mr)
    rates_O = EXFOR.B2R(KEVO,BO,TrO,mr)

    plt.figure()
    plt.plot(TrC/1e3,rates_C,label = 'EXFOR_C1269009_FUS')
    plt.plot(TrD/1e3,rates_D,label = 'EXFOR_D0769002_FUS')
    plt.plot(TrO/1e3,rates_O,label = 'EXFOR_O1531004_FUS')

    plt.legend()
    plt.grid()
    plt.xlabel('MeV')
    plt.ylabel('<σv> m3/s')
    if save:
        plt.savefig('O16O16FUS.png')
    if show:
        plt.show()




def verify_SI28FUS():
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

    rates_C    = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D75  = EXFOR.B2R(KEVD75,BD75,TrD75,mr)
    rates_D77  = EXFOR.B2R(KEVD77,BD77,TrD77,mr)
    rates_E    = EXFOR.B2R(KEVE,BE,TrE,mr)

    plt.figure()
    
    plt.plot(TrC/1e3,rates_C,label = 'EXFOR_C1269002_FUS')
    plt.plot(TrD75/1e3,rates_D75,label = 'EXFOR_D0750002_FUS')
    plt.plot(TrD77/1e3,rates_D77,label = 'EXFOR_D0772002_FUS')
    plt.plot(TrE/1e3,rates_E,label = 'EXFOR_E1435016_FUS')


    plt.legend()
    plt.grid()
    plt.xlabel('MeV')
    plt.ylabel('<σv> m3/s')
    plt.savefig('SI28FUS.png')
    plt.show()

def verfiy_SIOFUS():
    mr = 0.5*15.9994*1.6726*1e-27
    '''
    EXFOR
    '''
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

    rates_C = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D = EXFOR.B2R(KEVD,BD,TrD,mr)
    rates_O = EXFOR.B2R(KEVO,BO,TrO,mr)

    plt.figure()
    plt.plot(TrC/1e3,rates_C,label = 'EXFOR_C1269009_FUS_O')
    plt.plot(TrD/1e3,rates_D,label = 'EXFOR_D0769002_FUS_O')
    plt.plot(TrO/1e3,rates_O,label = 'EXFOR_O1531004_FUS_O')

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

    rates_C    = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D75  = EXFOR.B2R(KEVD75,BD75,TrD75,mr)
    rates_D77  = EXFOR.B2R(KEVD77,BD77,TrD77,mr)
    rates_E    = EXFOR.B2R(KEVE,BE,TrE,mr)

    #plt.figure()
    
    plt.plot(TrC/1e3,rates_C,label = 'EXFOR_C1269002_FUS_Si')
    plt.plot(TrD75/1e3,rates_D75,label = 'EXFOR_D0750002_FUS_Si')
    plt.plot(TrD77/1e3,rates_D77,label = 'EXFOR_D0772002_FUS_Si')
    plt.plot(TrE/1e3,rates_E,label = 'EXFOR_E1435016_FUS_Si')

    plt.legend()
    plt.grid()
    plt.xlabel('MeV')
    plt.ylabel('<σv> m3/s')
    plt.savefig('compare.png')
    plt.show()

def power_fusion(ni1,ni2,rates,Q,delta):
    '''
    ni: ion density
    f: fusion rate
    Q: Q value
    delta: efficiency
    '''
    return ni1*ni2*rates*Q/(1+delta)

def power_fusion_Si(ni1,ni2,Q,delta):
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

    rates_C    = EXFOR.B2R(KEVC,BC,TrC,mr)
    rates_D75  = EXFOR.B2R(KEVD75,BD75,TrD75,mr)
    rates_D77  = EXFOR.B2R(KEVD77,BD77,TrD77,mr)
    rates_E    = EXFOR.B2R(KEVE,BE,TrE,mr)

    P_C = power_fusion(ni1,ni2,rates_C,Q,delta)
    P_D75 = power_fusion(ni1,ni2,rates_D75,Q,delta)
    P_D77 = power_fusion(ni1,ni2,rates_D77,Q,delta)
    P_E = power_fusion(ni1,ni2,rates_E,Q,delta)

    plt.figure()
    
    plt.plot(TrC/1e3,P_C,label = 'EXFOR_C1269002_FUS')
    plt.plot(TrD75/1e3,P_D75,label = 'EXFOR_D0750002_FUS')
    plt.plot(TrD77/1e3,P_D77,label = 'EXFOR_D0772002_FUS')
    plt.plot(TrE/1e3,P_E,label = 'EXFOR_E1435016_FUS')


    plt.legend()
    plt.grid()
    plt.xlabel('MeV')
    plt.ylabel('Power_fusion MW')
    plt.savefig('Power_SI28FUS.png')
    plt.show()



if __name__ == '__main__':
    #verify_DT()
    #verify_O16O16Si28(save=True,show=True)
    #verify_O16O16FUS(save=True,show=True)
    
    #verify_SI28FUS()
    #verfiy_SIOFUS()
    ni = 1e20
    power_fusion_Si(1e20,1e20,13.431,1)
    