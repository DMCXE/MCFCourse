import numpy as np
import scipy.io as scio


def DT(): 
    '''
    https://www-nds.iaea.org/exfor/servlet/X4sSearch5?reacc=1-H-3(D%2CN)2-HE-4%2C%2CSIG
    1985
    KEV        B
    '''
    KEV = np.array([13.8110,14.0760,18.0560,20.0880,23.220,
                    25.7230,25.8550,25.9780,26.8330,27.1110,
                    28.8220,28.9880,30.1210,33.1170,35.6250,
                    36.1240,39.0980,42.1020,45.1160,48.1330,
                    48.3020,51.0910,54.1350,57.1020,60.0860,
                    63.1180,66.1880,72.0920,73.6040,75.070,
                    78.2370,80.9910,84.0040,87.2040,89.0940,
                    91.0310,92.9050,96.6310,99.0120,101.090,
                    105.010,108.180,113.270,114.310])
    B = np.array([0.010, 0.012,0.037,0.059,0.109,
                    0.168,0.174,0.173,0.184,0.196,
                    0.246,0.247,0.268,0.382,0.468,
                    0.503,0.635,0.793,1.0,1.170,
                    1.250,1.390,1.660,1.960,2.170,
                    2.460,2.780,3.40,3.480,3.690,
                    3.950,4.050,4.310,4.490,4.540,
                    4.520,4.750,4.820,4.880,4.940,
                    4.950,5.070,4.980,4.890])
    return KEV,B

def DT2():
    '''
    https://www-nds.iaea.org/exfor/servlet/X4sSearch5?reacc=1-H-3(D%2CN)2-HE-4%2C%2CSIG
    1957
    : http://www-nds.iaea.org/EXFOR/A1172.003
    MEV        mB
    '''
    
    MEV = np.array([0.0481, 0.0537, 0.0671, 0.0733, 0.0867, 0.0931, 0.0954, 0.1011, 0.1062, 0.1130, 
                0.1141, 0.1211, 0.1260, 0.1321, 0.1332, 0.1412, 0.1458, 0.1532, 0.1659, 0.1773, 
                0.1859, 0.1936, 0.2062, 0.2151, 0.2265, 0.2350, 0.2548, 0.2600, 0.2675, 0.2880, 
                0.2957, 0.3098, 0.3348, 0.3854, 0.4602, 0.4822, 0.5350, 0.5570, 0.6108, 0.6300, 
                0.7045, 0.7336])
    

    mB = np.array([1400.0, 1650.0, 3080.0, 3440.0, 4510.0, 4770.0, 5100.0, 5200.0, 5100.0, 5050.0,
               5140.0, 4790.0, 4760.0, 4270.0, 4410.0, 3860.0, 3970.0, 3400.0, 3320.0, 3000.0,
               2630.0, 2490.0, 2290.0, 2190.0, 1980.0, 1910.0, 1690.0, 1680.0, 1510.0, 1370.0,
               1270.0, 1220.0, 1160.0, 850.0, 750.0, 670.0, 590.0, 570.0, 510.0, 470.0, 420.0,
               360.0])
    return MEV,mB


def O16O16SI28_A0381006():
    '''
    https://www-nds.iaea.org/exfor/servlet/X4sGetSubent?reqx=7151&subID=100381006&plus=1
    MEV  mB
    '''
    MEV = np.array([34.7, 39.7, 44.8, 47.7, 50.1, 54.3, 60.2, 65.2])
    mB  = np.array([47.0, 70.0, 85.0, 86.0, 65.0, 57.0, 27.0, 24.0])
    return MEV,mB

def O16O16SI28_E1147029():
    '''
    MEV mB
    '''
    MEV= np.array([3.527e+01, 4.013e+01, 4.523e+01, 4.767e+01, 
                   5.036e+01, 5.438e+01, 6.038e+01, 6.538e+01])
    mB = np.array([4.987e+01, 6.610e+01, 8.529e+01, 8.690e+01, 
                   6.391e+01, 5.687e+01, 2.736e+01, 2.265e+01])
    return MEV,mB

def O16O16SI28_C1269002():
    '''
    MEV mB
    '''
    MEV = np.array([7.00, 7.25, 7.50, 7.75, 8.00, 
                    8.25, 8.50, 8.75, 9.00, 9.25, 
                    9.50, 9.75, 10.00, 10.25, 10.50, 
                    10.75, 11.00, 11.25, 11.50, 11.75, 12.00])
    mB  = np.array([0.00309, 0.00816, 0.0297, 0.0712, 0.170, 
                    0.326, 0.575, 1.097, 1.70, 3.00, 4.40, 
                    6.29, 8.29, 11.4, 16.0, 18.1, 
                    22.0, 27.8, 34.2, 39.9, 46.3])
    return MEV,mB

def O16O16FUS_D0769002():
    '''
    MEV mB
    '''
    MEV = np.array([8.27,9.27,10.77,12.27])
    mB = np.array([1.64,21.1,152,407])
    return MEV,mB

def O16O16FUS_C1269009():
    '''
    MEV mB
    '''
    MEV = np.array([7.00,7.25,7.50,7.75,8.00,8.25,8.50,8.75,9.00,9.25,
                    9.50,9.75,10.00,10.25,10.50,10.75,11.00,11.25,11.50,11.75,12.00])
    mB  = np.array([0.0104,0.0332,0.109,0.289,0.747,1.54,2.08,5.24,8.78,16.3,
                    25.6,38.8,56.9,78.1,117,145,196,249,308,379,438])
    return MEV,mB

def O16O16FUS_O1531004():
    '''
    MEV mB
    '''
    MEV = np.array([13.48,13.98,14.47,14.92,15.4,16.0,
                    16.48,16.85,17.44,18.01,18.48,19.0,
                    19.35,19.97,20.39,20.82,22.44,22.86,
                    23.47,24.01,24.44,24.9,25.5,25.89,
                    26.5,26.82,27.31,27.99,28.47,29.23,
                    29.81,30.78,31.4,32.21,32.92])
    mB  = np.array([420,450,500,618,680,700,760,720,750,
                    820,900,910,900,870,900,900,940,970,
                    970,930,940,900,920,980,980,1020,940,
                    820,880,920,960,1070,1060,1010,990])
    return MEV,mB




def SI28SI28FUS_D0772002():
    '''
    MEV mB
    '''
    import numpy as np

    MEV = np.array([29.977, 30.228, 30.479, 30.730, 30.980, 31.231, 31.482, 31.733, 31.984, 32.235,
                        32.486, 32.737, 32.988, 33.239, 33.490, 33.741, 33.992, 34.242, 34.493, 34.744,
                        34.995, 35.246, 35.497, 35.744, 35.999, 36.124, 36.375, 36.626, 36.752, 36.877,
                        37.003, 37.128, 37.254, 37.379, 37.504, 37.630, 37.755, 37.881, 38.006, 38.132,
                        38.257, 38.508, 38.759])

    mB = np.array([88.967, 102.35, 111.88, 124.57, 138.74, 150.88, 162.75, 178.25, 191.62, 202.58,
                        222.40, 233.52, 248.18, 269.84, 283.83, 299.15, 317.80, 328.89, 342.26, 360.40,
                        377.73, 383.42, 402.30, 420.11, 438.51, 446.78, 464.45, 472.14, 496.36, 499.59,
                        506.69, 516.64, 524.54, 534.39, 536.21, 540.98, 538.71, 550.14, 557.26, 563.04,
                        562.86, 577.46, 581.55])
    return MEV,mB

def SI28SI28FUS_D0750002():
    '''
    MEV mB
    '''
    MEV = np.array([22.769, 23.071, 23.573, 24.076, 24.577, 25.080, 25.582, 26.084, 26.587, 27.089,
                        27.591, 28.093, 28.595, 29.097, 30.101, 31.106, 32.109, 33.113, 34.117, 35.120,
                        36.124, 37.003, 38.006, 39.010, 40.014])

    mB  = np.array([0.00062580, 0.00080790, 0.0045860, 0.010034, 0.027261, 0.062439, 0.16150, 0.40000,
                        1.0398, 2.9850, 6.8934, 15.957, 31.650, 49.240, 96.418, 145.26, 209.26, 252.44,
                        306.15, 367.21, 433.68, 506.69, 557.26, 592.27, 627.56])
    return MEV,mB

def SI28SI28FUS_C0415002():
    '''
    MEV mB
    '''
    MEV = np.array([30.0, 32.0, 34.0, 36.0, 38.0, 40.5, 42.5, 44.0, 45.0, 46.5, 48.0, 49.5, 52.5, 55.0, 57.5, 60.5, 66.3])
    mB = np.array([115., 260., 330., 446., 562., 715., 765., 828., 852., 922., 938., 970., 1072., 1068., 1064., 1058., 1050.])
    return MEV,mB

def SI28SI28FUS_E1435016():
    '''
    MEV mB
    '''
    MEV = np.array([30.0, 32.0, 34.0, 36.0, 38.0, 40.5, 42.5, 44.0, 45.0, 46.5, 48.0, 49.5, 52.5, 55.0, 57.5, 60.5, 66.3])
    mB  = np.array([115.0, 260.0, 330.0, 446.0, 562.0, 715.0, 765.0, 828.0, 852.0, 922.0, 938.0, 970.0, 1072.0, 1068.0, 1064.0, 1058.0, 1050.0])
    return MEV,mB



def from_xhs_DT_data():
    data = scio.loadmat('xhsDTdata.mat')

    return data

def B2R(keV,b,Tr,mr):
    '''
    keV: 截面数据的能量单位
    b:  截面数据的截面单位,m2
    Tr: 有效温度 keV ,ndarray
    mr: 有效质量 kg
    '''
    e  = 1.6022e-19 
    A1 = np.sqrt(8/(np.pi*mr)) # kg*(-1/2)
    A2 = (e*Tr*1e3)**(-3/2) # J ** -3/2
    #b *= *1e-28              # m2,
    Td = 1000 #keV
    # 不对数据进行拟合,通过原始数据进行离散积分
    dkev = np.diff(keV)
    Inter = np.zeros_like(Tr)
    for i,t in enumerate(Tr):
        for j,E in enumerate(keV):
            if j == 0 :
                Inter[i] +=0
            else:
                # Inter[i] += ((b[j]*1e-28)*(e*E*1e3)*np.exp(-E/t)+(b[j-1]*1e-28)*(e*keV[j-1]*1e3)*np.exp(-keV[j-1]/t))\
                #             *((E-keV[j-1])*e*1e3)/2
                #Inter[i] += (b[j]*1e-28)*(e*E*1e3)*np.exp(-E/t)*dkev[j-1]*e*1e3
                # 非maxwell模式
                Inter[i] += (b[j]*1e-28)*np.sqrt(e*E*1e3)*np.exp(-(E+Td)/t)*np.sinh(2*np.sqrt(Td*E)/t)*dkev[j-1]*e*1e3
    
    A1 = np.sqrt(2/(np.pi*mr*e*Tr*1e3*e*Td*1e3))
    A2 = 1
    return A1 * A2 *Inter


    


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    KEV,B = DT()
    
    Tr = np.linspace(KEV[0]*0.8,KEV[-1],100,endpoint=True)
    print(Tr)
    mr = 1.2*1.6726*1e-27 #DT体系的相对质量m,kg
    rates = B2R(KEV,B,Tr,mr)
    plt.figure()
    plt.scatter(Tr,rates)
    plt.show()

    plt.figure()
    plt.loglog(KEV,B*1e-28)
    plt.show()

    data = from_xhs_DT_data()
    Teff = data['Teff']
    sigmadt = data['sgmv1']
    Teff = Teff[Teff<KEV[-1]]
    sgmv = sigmadt[0][0:len(Teff)]

    print(len(Teff))
    print(len(sgmv))

    plt.figure()
    plt.scatter(Tr,rates)
    plt.scatter(Teff,sgmv)
    #plt.show()

    MEV,mB = DT2()
    KEV = 1e3*MEV
    B   = 1e-3*mB
    Tr = np.linspace(KEV[0],KEV[-1],100,endpoint=True)
    rates = B2R(KEV,B,Tr,mr)
    #plt.figure()
    plt.scatter(Tr,rates)
    #plt.show()

    #混合
    KEV1,B1 = DT()
    KEVm = np.append(KEV1,KEV)
    Bm   = np.append(B1,B)

    arg = np.argsort(KEVm)
    KEVm = KEVm[arg]
    Bm   = Bm[arg]
    plt.figure()
    plt.loglog(KEVm,Bm*1e-28)
    Tr = np.linspace(KEVm[0],KEVm[-1],100,endpoint=True)
    rates = B2R(KEVm,Bm,Tr,mr)
    #plt.scatter(Tr,rates,c='b')
    plt.show()


