import numpy as np
import scipy.optimize 
import matplotlib.pylab as plt
from matplotlib import rc
import pandas as pd
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'font.size': 25})
rc('text', usetex=True)
rc('legend',numpoints=1)
import matplotlib.patches as patches 
from scipy.optimize import curve_fit

class CO2:
    """This class contains some functions related to
    CO2 dissolved in the magma
    """
    
    def fO2_IW(self,Tk,Pbar,buf=0):
        """Calculates the Oxyge fugacity
        as a function of input temperature in K, pressure in bars
        and buffer units above or below the IW buffer
        Uses Eq A8 of Holloway et al. 1992
        returns the log10 value

        default value of buf is zero, equal to IW
        """
        logfO2=6.899-27714.0/Tk+0.05*(Pbar-1.0)/Tk + buf
        
        return(logfO2)
    
    def PCO2_calculate_new(self,T0,P0=3.0e4,buf=0):
        """
        This function calculates PCO2 and XCO3 in melt for a
        given fO2 buffer. The default is IW buffer. The expression for
        equilibrium coefficients only work for up to 3 GPa and 1700 C. 
        
        The default input pressure is 3 GPa (3e4 bars)
        buf is the unit above or below IW

        Input temperature in Celsius, pressure in bar
        Need to convert input temperature from celsius to K
        """
        #First calculate the  temperature
        P=P0 #Pressure in bar
        T=T0+273.0 #self.adiabat(r) #Temperature at the freezing front
        #Calculate log10K1
        logK1=40.07639-T*2.53932*1.0e-2+T*T*5.27096*1.0e-6 + 0.0267*(P-1.0)/T
        #calculate log10K2
        logK2=-6.24763-282.56/T-0.119242*(P-1.0e3)/T
        #Calculate the equilibrium constant Kp from P and T
        #Calculate the equilibrium constant K' eq. A1 of Holloway et al. 1992
        logKp=33.82876-(163.291+0.092542*P)/T-T*2.53932*1.0e-2+T*T*5.27096*1.0e-6
        logKP2=logK1+logK2
        
        #-(163.291+0.092542*P)/T-(2.53932*1.0e-2)*T + T*T*(5.27096*1.0e-6)
        #Calculate fO2 from O'neil 1987, need to provide pressure in bars
        # and temperature in Kelvin
        
        # equation for log10fO2 for IW buffer
        # Eq A8 of Holloway et al., 1992
        logfO2=self.fO2_IW(T,P0,buf)
        #Calculate log10Kf
        #p. 112 of Holloway et al. 1992
        logKf=logKp+logfO2
        Kf=np.exp(2.303*logKf)
        #mole fraction of CO3 in melt
        XCO3_melt=Kf/(1.0+Kf)
        #wt% of CO3 dissolved in melt
        wco2=(120.2656*XCO3_melt)/(1.0-0.202656*XCO3_melt)
        #concentration in fraction
        c=wco2/100.0
        #To calculate CO2 partial pressure calculate K1
        #Equation A1(a) of Holloway et al.
        
        
        #40.07639-2.53932e-2*T+(5.27096e-6)*T*T+0.0267*(P-1000.0)/T
        K1=np.exp(2.303*logK1)
        fO2=np.exp(2.303*logfO2) #Eq 3 from Holloway et al 1992
        PCO2=K1*fO2 #fCO2=PCO2 in bars
        #Calculate PCO2 from Pawley et al. 1992, page 223
        #XCO2_melt(ppm)= 0.492*fCO2
        # Assuming PCO2 is fCO2
        #PCO2=(XCO3_melt*1.0e6)/0.492
        return(PCO2,c,fO2)

def hirschmann_withers(mycO2):
    """This function plots figures 3
    and 4 from Hirschmann and Withers"""

    P1,c1,fO21= myCO2.PCO2_calculate_new(myCO2.T0,buf=1.0)
    P2,c2,fO22= myCO2.PCO2_calculate_new(mycO2.T0,buf=0.0)
    P3,c3,fO23= myCO2.PCO2_calculate_new(myCO2.T0,buf=-1.0)
    
    plt.figure(figsize=(12,16))
    plt.subplot(2,1,1)
    plt.semilogy(np.log10(fO21),c1*100.0,'-b',lw=4)
    plt.semilogy(np.log10(fO22),c2*100.0,'--r',lw=4)
    plt.semilogy(np.log10(fO23),c3*100.0,'-.k',lw=4)
    plt.xlabel(r'log$_{10} f_{O2}$ (bar)')
    plt.ylabel(r'CO$_2$ wt\% in melt ')
    plt.legend(['IW+1','IW','IW-1'],fancybox=True,framealpha=0.7,loc=4)
    plt.text(-11.6,0.4,'P = 3 GPa',fontsize=30)
    plt.text(-11.6,0.2,r'T = $1200-1600^{o}$C',fontsize=30)
    
    plt.subplot(2,1,2)
    
    PGPa=np.linspace(0.0,3.0)
    Pbar=PGPa*1.0e4
    T1=1300.0+273.0 #input temp in K
    #IW
    IW_bar=myCO2.fO2_IW(T1,Pbar,buf=0.0)
    #IW+1
    IWp1bar=myCO2.fO2_IW(T1,Pbar,buf=1.0)
    #Iw+2
    IWp2bar=myCO2.fO2_IW(T1,Pbar,buf=2.0)
    #Iw+3
    IWp3bar=myCO2.fO2_IW(T1,Pbar,buf=3.0)
    plt.plot(PGPa,IWp3bar,lw=3)
    plt.plot(PGPa,IWp2bar,'--',lw=3)
    plt.plot(PGPa,IWp1bar,'-.',lw=3)
    plt.plot(PGPa,IW_bar,':',lw=5)
   
    plt.legend(['IW+3','IW+2','IW+1','IW'],fancybox=True,loc=4,framealpha=0.6)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel(r'log$_{10}$($f_{O2}$) (bar)')
    plt.text(0.2,-6.5,r'T = 1573 K',fontsize=25)
    plt.ylim(-11.0,-6.0)

    ## Create a subset of Figure 4 from Hirschmann and Withers
    Tsol=np.array([1200.0,1320.0,1470.0,1520.0,1600.0])
    T1=np.linspace(Tsol[0], Tsol[0]+500.0)
    P1,c1,fO21= myCO2.PCO2_calculate_new(T1,P0=1.0e4,buf=0.0)
    T2=np.linspace(Tsol[1], Tsol[1]+500.0)
    P2,c2,fO22= myCO2.PCO2_calculate_new(T2,P0=2.0e4,buf=0.0)
    T3=np.linspace(Tsol[2], Tsol[2]+500.0)
    P3,c3,fO23= myCO2.PCO2_calculate_new(T3,P0=3.0e4,buf=0.0)
    T4=np.linspace(Tsol[3], Tsol[3]+500.0)
    P4,c4,fO24= myCO2.PCO2_calculate_new(T4,P0=4.0e4,buf=0.0)
    T5=np.linspace(Tsol[4], Tsol[4]+500.0)
    P5,c5,fO25= myCO2.PCO2_calculate_new(T5,P0=5.0e4,buf=0.0)

    cIW=np.array([c1[0],c2[0],c3[0],c4[0],c5[0]])
    
    plt.figure(figsize=(10,12))
    
    plt.semilogy(T1,c1*100.0,lw=3)
    plt.semilogy(T2,c2*100.0,lw=3)
    plt.semilogy(T3,c3*100.0,lw=3)
    plt.semilogy(T4,c4*100.0,lw=3)
    plt.semilogy(T5,c5*100.0,lw=3)
    
    plt.legend(['1 GPa', '2 GPa','3 GPa', '4 GPa', '5 GPa'], loc=4, fancybox=True, framealpha=0.5)
    plt.semilogy(Tsol,cIW*100,'--',lw=2)
    plt.xlabel(r'T ($^o$C)')
    plt.ylabel(r'CO$_3$ in melt wt$\%$')
    plt.xlim(1200,1800)

    
def temperature_plot(myCO2):
    """Plots various things as a function
    of temperature"""
    T0=myCO2.T0
    P1,c1,fO21= myCO2.PCO2_calculate_new(T0,buf=1.0)
    P2,c2,fO22= myCO2.PCO2_calculate_new(T0,buf=0.0)
    P3,c3,fO23= myCO2.PCO2_calculate_new(T0,buf=-1.0)
    plt.figure(figsize=(12,18))

    plt.subplot(3,1,1)
    plt.semilogy(T0,c1*100.0,'-b')
    plt.semilogy(T0,c2*100.0,'-r')
    plt.semilogy(T0,c3*100.0,'-k')
    plt.ylabel(r'CO$_3$ in melt wt$\%$')
    #plt.ylim(1.0e-4,1.0)
    
    plt.subplot(3,1,2)
    plt.semilogy(T0,P1,'-b')
    plt.semilogy(T0,P2,'-r')
    plt.semilogy(T0,P3,'-k')
    plt.ylabel('PCO2 (bar)')
    
    plt.subplot(3,1,3)
    plt.semilogy(T0,fO21,'-b')
    plt.semilogy(T0,fO22,'-r')
    plt.semilogy(T0,fO23,'-k')
    
    plt.legend(['IW+1','IW','IW-1'],fancybox=True,framealpha=0.7,loc=4)
    plt.ylabel('fO2 (bar)')
    plt.xlabel('Temperature ($^o$C)')


myCO2=CO2()
myCO2.T0=np.linspace(1200.0,1600.0)


hirschmann_withers(myCO2)
temperature_plot(myCO2)







#plt.ylim(1.0e-4,1)
plt.show()
