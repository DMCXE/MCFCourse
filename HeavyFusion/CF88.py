import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

"""
CF88数据库
http://www.nuclear.csdb.cn/data/CF88/

Q    反应Q值，单位MeV
T9   反应温度，单位GK=10^9K
T9nm 反应温度，单位(T9)^(n/m)
rates experss as NA<σv>, NA for A
       units of (cm^3/mol)^(N-1)/s, N for number of body reaction

"""


def O16O16(T9):
    """    
    O16 + O16 -> Ne20 + He4
    O16+O16 (S32)
    Q159 = 16.542
    f159 = 7.10e+36/T923*exp(-135.93/T913-0.629*T923-0.445*T943
    |       +0.0103*T9**2)
    """
    Q = 16.542
    T9nm = lambda n,m: T9**(n/m)
    f = 7.10e+36/T9nm(2,3)*np.exp(-135.93/T9nm(1,3)-0.629*T9nm(2,3)-0.445*T9nm(4,3)
                                     +0.0103*T9**2)
    return Q, f


def DT(T9):
    """   
    Q11 = 17.589
    f11 = 8.09e+10/T923*exp(-4.524/T913-(T9/0.120)**2)*(1.0+0.092
    |      *T913+1.80*T923+1.16*T9+10.52*T943+17.24*T953)
    |      +8.73e+08/T923*exp(-0.523/T9)
    """
    Q = 17.589
    T9nm = lambda n,m: T9**(n/m)
    f = 8.09e+10/T9nm(2,3)*np.exp(-4.524/T9nm(1,3)-(T9/0.120)**2)*(1.0+0.092
                                     *T9nm(1,3)+1.80*T9nm(2,3)+1.16*T9+10.52*T9nm(4,3)+17.24*T9nm(5,3)
                                     )+8.73e+08/T9nm(2,3)*np.exp(-0.523/T9)
    return Q, f

def O16mbar():
    '''
    data from EXFOR 
    https://www-nds.iaea.org/exfor/servlet/X4sSearch5?reacc=8-O-16(8-O-16%2CA)14-SI-28%2C%2CSIG
                MEV        MB
                7.95       0.0919
                8.25       0.235
                8.45       0.448
                8.66       0.719
                8.86       1.30
                9.13       1.91
                9.34       3.05
                9.56       4.16
                9.78       5.77
                10.00       7.37
                10.22       9.46
                10.45      12.4
                10.67      15.4
                10.91      18.7
                11.14      22.9
                11.37      25.4
                11.61      27.6
                11.79      31.6
                12.10      34.1
                12.33      36.2
                12.63      37.2
                12.87      39.6
                13.10      42.8
                13.34      42.8
                13.58      46.6
                13.83      51.5`
    '''
    MeV = np.array([7.95,8.25,8.45,8.66,8.86,9.13,9.34,9.56,9.78,
                    10.00,10.22,10.45,10.67,10.91,11.14,11.37,11.61,
                    11.79,12.10,12.33,12.63,12.87,13.10,13.34,13.58,13.83])
    MB  = np.array([0.0919,0.235,0.448,0.719,1.30,1.91,3.05,4.16,
                    5.77,7.37,9.46,12.4,15.4,18.7,22.9,25.4,27.6,
                    31.6,34.1,36.2,37.2,39.6,42.8,42.8,46.6,51.5])
    # 多项式插值
    #MeV = np.linspace(7.95,13.83,10000)
    # MBi = interp1d(MeV,MB,kind='cubic')
    # # 通过多项式拟合
    # p = np.polyfit(MeV,MB,3)

    # M = np.linspace(7.95,13.83,10000)
    # plt.figure()
    # plt.plot(M,MBi(M),label='inter')
    # plt.plot(M,np.polyval(p,M),label='polyfig')
    # plt.scatter(MeV,MB,label='data')
    # plt.show()
    return MeV, MB

def mbar2rate(m1:float,m2:float,T1:float,T2:float,MeV,MB):
    '''
    T1, T2 -> MeV
    mBar -> reactionrate
    means <σv>
    [Nevins00]
    <σv> = sqrt(8/(pi*mr)) * (kbTr)**(-3/2) * integral(0,inf)(sigma(E)*E*exp(-E/kbTr)dE)
    mr = m1*m2/(m1+m2)
    Tr = (m1*T1 + m2*T2)/(m1+m2)

    '''
    kbev = 1/11604.51812    #eV/K
    Tr = (m1*T1 + m2*T2)/(m1+m2)
    kbTr = Tr #MeV
    mr = m1*m2/(m1+m2)
    # 处理离散数据 MeV, MB
    sigma = interp1d(MeV,MB,kind='cubic')
    # integral
    # for tr in Tr:
    #     integrand = lambda E: sigma(E)*E*np.exp(-E/kbTr)
    #     rate, _ = quad(integrand,0,np.inf)
    #     rate = np.sqrt(8/(np.pi*mr)) * (kbTr)**(-3/2) * rate
    integrand = lambda E: sigma(E)*(1e-3)*E*np.exp(-E/kbTr)
    rate, _ = quad(integrand,MeV[0],MeV[-1])
    rate = np.sqrt(8/(np.pi*mr)) * (kbTr)**(-3/2) * rate
    return rate

    
   

if __name__ == '__main__':
    # T9  = np.linspace(1e-3, 10, 100,endpoint=True)
    # #Q159, f159 = O16O16(T9)
    # Q159, f159 = DT(T9)
    # print(f159)

    # plt.figure()
    # plt.loglog(T9, f159)
    # plt.xlabel('T9')
    # plt.ylabel('rate')
    # plt.title('O16+O16')
    # plt.show()
    MeV, MB = O16mbar()
    MB *= 1e-3 * 1e-28
    mev = np.linspace(MeV[0],MeV[-1],100,endpoint=True)
    rate = np.array([])
    plt.figure()
    for i in mev:
        print(i)
        rate = np.append(rate,mbar2rate(16,16,i,i,MeV,MB))
    plt.loglog(mev*1e3,rate)
    plt.show()
        