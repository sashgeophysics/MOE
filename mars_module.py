import numpy as np
import scipy.optimize 
import matplotlib.pylab as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'font.size': 25})
rc('text', usetex=True)
rc('legend',numpoints=1)
import matplotlib.patches as patches 
from scipy.optimize import curve_fit
def fit_func(x,m):
    """Fitting function"""
    return m*x

class CO2:
    """This class contains the functions
    necessary for calculating the concentration
    of CO2 in the magma ocean"""

    omega_CO2=44.0 #constant g/mole
    omegaM=36.0 #constant g/mol
    K0=np.exp(-14.83)*1.0e-5 #constant 1/Pa
    T0=1200 #reference temperature degree C
    dH = 5.2e3 #J enthalpy
    R= 8.314 #Universal gas constant
    DCO2=0.001 # Partitioning coefficient
    #def calc_kco2(self,T):
    #    """Calculates the equilibrium constant
    #    of CO2 partitioning from Pan et al. 1991"""
    #    self.KCO2=self.K0*np.exp(-self.dH*(1.0/T -1.0/self.T0)/self.R)
            
class H2O:
    """This class contains functions needed to 
    calculate the concentration of H2O"""
   
    KH2O=2.1e-7*1.0e-5 #Equilibrium constant 1/Pa
    omega_H2O = 18.0  #constant g/mole
    omega_Moore=62.5 #constant g/mole
class REE:
    """This class contains the partition coefficients
    of REE from Borg and Draper 2003, Table 5. 
    For each element, the array contains partition
    coefficients in the order 
    ol-opx-cpx-gt-ilm-titanite-phlogopite-amphibole
    """
    Rb=np.array([0.00018, 0.0006, 0.011, 0.001, 1E-05, 1E-05, 3.06, 0.31])
    Ba=np.array([0.00018, 0.0006, 0.0165, 0.0015, 1E-05, 1E-05, 1.09, 0.436])
    Th=np.array([0.0001, 0.001, 0.01, 0.0001, 1E-05, 1E-05, 0.31, 0.11])
    U=np.array([0.0001, 0.001, 0.02, 0.0001, 1E-05, 1E-05, 0.21, 0.15])
    Ta=np.array([0.01, 0.025, 0.05, 0.08, 1E-05, 19.0, 1E-05, 0.38])
    K=np.array([7.7E-05, 0.00026, 0.011, 0.018, 1E-05, 1E-05, 2.7, 0.4])
    La=np.array([0.0004, 0.002, 0.054, 0.01, 1E-05, 2.1, 0.036, 0.17])
    Ce=np.array([0.0005, 0.003, 0.0646, 0.021, 1E-05, 4.72, 0.034, 0.26])
    P=np.array([0.055, 0.03, 0.03, 0.261, 1E-05, 1E-05, 1E-05, 1E-05])
    Sr=np.array([0.00019, 0.007, 0.12, 0.0011, 1E-05, 1E-05, 0.081, 0.12])
    Nd=np.array([0.001, 0.0068, 0.107, 0.036, 1E-05, 9.96, 0.032, 0.44])
    Sm=np.array([0.0013, 0.01, 0.17, 0.103, 1E-05, 11.0, 0.031, 0.76])
    Zr=np.array([0.01, 0.01, 0.233, 0.268, 0.44, 0.01, 2.1, 0.7])
    Hf=np.array([0.01, 0.01, 0.233, 0.322, 0.44, 0.01, 2.1, 0.7])
    Eu=np.array([0.0016, 0.013, 0.28, 0.22, 0.001, 11.0, 0.03, 0.88])
    Gd=np.array([0.0015, 0.016, 0.32, 0.43, 0.001, 11.0, 0.03, 0.86])
    Tb=np.array([0.0015, 0.019, 0.31, 0.335, 0.001, 11.0, 0.032, 0.78])
    Dy=np.array([0.0017, 0.022, 0.33, 0.47, 0.001, 11.0, 0.034, 0.77])
    Y=np.array([0.0016, 0.026, 0.31, 1.01, 0.001, 11.0, 0.036, 0.73])
    Er=np.array([0.0015, 0.03, 0.3, 0.79, 0.001, 9.0, 0.038, 0.68])
    Tm=np.array([0.0015, 0.04, 0.29, 1.18, 0.01, 8.0, 0.04, 0.64])
    Yb=np.array([0.0015, 0.049, 0.29, 1.59, 0.01, 7.0, 0.042, 0.59])
    Lu=np.array([0.0015, 0.06, 0.28, 1.93, 0.01, 6.0, 0.044, 0.51])
    #########################################
    ## Element concentration in the CI chondrite
    ## From Lodder 2003, Table 3
    # Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu
    CI_Rb = 2.13e-6
    CI_Ba = 2.31e-6
    CI_Th = 0.0309e-6
    CI_U  = 0.0084e-6
    CI_Ta = 0.0144e-6
    CI_K  = 530.0e-6
    CI_La = 0.232e-6
    CI_Ce = 0.621e-6
    CI_P  = 920.0e-6
    CI_Sr = 7.74e-6
    CI_Nd = 0.457e-6
    CI_Sm = 0.145e-6
    CI_Zr = 3.96e-6
    CI_Hf = 0.115e-6
    CI_Eu = 0.0546e-6
    CI_Gd = 0.198e-6
    CI_Tb = 0.0356e-6
    CI_Dy = 0.238e-6
    CI_Y  = 1.53e-6
    CI_Er = 0.162e-6
    CI_Tm = 0.0237e-6
    CI_Yb = 0.163e-6
    CI_Lu = 0.0237e-6
    CI_REE=np.array([CI_Rb,CI_Ba,CI_Th,CI_U,CI_Ta,\
                     CI_K, CI_La, CI_Ce,CI_P, CI_Sr,\
                     CI_Nd, CI_Sm, CI_Zr, CI_Hf, CI_Eu,\
                     CI_Gd, CI_Tb,CI_Dy, CI_Y, CI_Er,\
                     CI_Tm, CI_Yb,  CI_Lu])
    ## CI REE from Anders and Grevesse, 1989
    CI_Rb_AG = 7.09e-6
    CI_Ba_AG = 4.49e-6
    CI_Th_AG = 0.0335e-6
    CI_U_AG  = 0.009e-6
    CI_Ta_AG = 0.0207e-6
    CI_K_AG  = 3770.0e-6
    CI_La_AG = 0.446e-6
    CI_Ce_AG = 1.136e-6
    CI_P_AG  = 1.04e-2
    CI_Sr_AG = 23.5e-6
    CI_Nd_AG = 0.8279e-6
    CI_Sm_AG = 0.2582e-6
    CI_Zr_AG = 11.4e-6
    CI_Hf_AG = 0.154e-6
    CI_Eu_AG = 0.0973e-6
    CI_Gd_AG = 0.33e-6
    CI_Tb_AG = 0.0603e-6
    CI_Dy_AG = 0.3942e-6
    CI_Y_AG  = 4.64e-6
    CI_Er_AG = 0.2508e-6
    CI_Tm_AG = 0.0378e-6
    CI_Yb_AG = 0.2479e-6
    CI_Lu_AG = 0.0367e-6
    CI_REE_AG=np.array([CI_Rb_AG,CI_Ba_AG,CI_Th_AG,CI_U_AG,\
                        CI_Ta_AG, CI_K_AG, CI_La_AG, CI_Ce_AG,\
                        CI_P_AG,CI_Sr_AG,CI_Nd_AG,CI_Sm_AG,CI_Zr_AG,\
                        CI_Hf_AG,CI_Eu_AG,CI_Gd_AG,CI_Tb_AG,CI_Dy_AG,
                        CI_Y_AG,CI_Er_AG,CI_Tm_AG,CI_Yb_AG,CI_Lu_AG])
    # NWA 1068 from Barrat 2002
    NWA1068=1.0e-6*np.array([5.75,127.0,0.409,0.1,0.2,1327.6,2.25,5.38,0.0,67.0,\
                      3.82,1.49,62.14,1.58,0.552,2.14,0.414,2.8,17.19,1.63,0.0,1.37,0.198])
    NWA2737 = 1.0e-6*np.array([1.28,37.5,0.13,0.056,0.07,414.9,1.17,2.87,436.4,27.2,\
                               1.43,0.266,5.01,0.14,0.0721,0.237,0.367,0.213,1.21,0.117,\
                               0.0,0.105,0.017])
    shergotty =1.0e-6*np.array([6.22,29.4,0.0,0.0,0.0,1392.0,1.5,3.51,0.0,51.0,2.6,1.01,0.0,0.0,\
                                0.43,1.64,0.0,2.16,0.0,1.33,0.0,1.19,0.176])
    zagami =1.0e-6*np.array([5.69,25.3,0.0,0.0,0.0,1243.0,1.6,3.75,0.0,45.9,2.89,1.17,0.0,0.0,\
                             0.476,0.0,0.0,2.66,0.0,1.6,0.0,1.38,0.201])
        ##########################################
class Earth(CO2,H2O):
    """This class contains functions for the Earth"""
    
    def solidus_liquidus(self,r):
        """Calculates silicate solidus enter pressure
        in GPa, these parameters are for the Earth
        The functions for solidus and liquidus are 
        from HMH2017"""
        
        PGPa=self.depth2PGPa(r)
        self.PGPa=PGPa
        if type(r) is float:
            if PGPa<20.0:
                dum1=1661.2*(PGPa/1.336 +1.0)**(1.0/7.437)
                dum2 = 1982.1*(PGPa/6.594+1.0)**(1.0/5.374)
            else:
                dum1=2081.8*(PGPa/101.69 + 1.0)**(1.0/1.226)
                dum2=78.74*(PGPa/4.054e-3 + 1.0)**(1.0/2.44)
        else:
            n=r.shape[0]        
            dum1=np.zeros(n)
            dum2=np.zeros(n)
            for ii in range(0,n):
                if PGPa[ii]<20.0:
                    dum1[ii]=1661.2*(PGPa[ii]/1.336 +1.0)**(1.0/7.437)
                    dum2[ii] = 1982.1*(PGPa[ii]/6.594+1.0)**(1.0/5.374)
                else:
                    dum1[ii]=2081.8*(PGPa[ii]/101.69 + 1.0)**(1.0/1.226)
                    dum2[ii]=78.74*(PGPa[ii]/4.054e-3 + 1.0)**(1.0/2.44)
        return(dum1,dum2)
    

class Mars(CO2,H2O,REE):
    """This class contains functions essential for Mars"""
    def __init__(self,tfinal=9.0e6*365.0*24.0*3600.0,noceans=1.0,\
                 HoverC=0.4,nsteps=1000,redox_factor=1.0):
        """Initial definitions, enter the number of
        oceans of H2O in the initial composition
        and the H:C ratio during initiation. Default values
        are 1 ocean and 0.4, respectively.
        nsteps is the number of steps in evolution"""
        self.area=144.4e12            # m^2 surface area of Mars
        self.radius=3.39e6            # m radius of Mars
        self.g=3.711                  # m/s^2, surface gravity of Mars
        self.mass = 0.624e24          # kg mass of Mars, Carr, encyclopedia of SS
        self.rho=4287.0               # kg/m^3 average density of Mars
        self.Cp=1.0e3                 # J/kg/K heat capacity
        self.alpha=5.0e-5             # 1/K coefficient of thermal expansion
        self.gam=1.3e-4               # K/m adiabatic gradient from Elkins-Tanton
        self.core=self.radius-1.995e6 # m Radius of the Martian core
        self.Fs = 0.0                 # W/m^2 solar heat flux
        self.tfinal = tfinal          # Final value of time
        self.dS     = 300.0           # Entropy of melting J/kg/K
        self.sigma  = 5.67e-8         # Steffan boltzman constant W/m^2/K^4
        self.PGPa=[]
        self.redox_factor = redox_factor # ratio fCO2/fO2
        ## Thermal evolution terms
        self.T_surf=self.freezing_front(self.radius)
        self.T_CMB=self.freezing_front(self.core)
        #Set initial temperature
        self.T_init=self.T_CMB-(self.radius-self.core)*self.gam+200.0
        #self.fit_coefficients=self.radius_temperature_analytical()
        #a-T polynomial fit coefficients for 1820<T<2070 deg C
        self.aT_fit_hi=[ -1.46402889e+00,   2.96827884e+03,   1.79228532e+05]
        #[ -1.19014770e+00,   2.03249128e+03,   1.09391926e+06]
        #a-T polynomial fit coefficient for 1350 < T <=1820 deg G
        self.aT_fit_lo= [ -1.52094133e-02,   6.59633338e+01,  \
                          -9.63920910e+04,   4.93403358e+07]
        #[ -1.52094133e-02,   6.82235562e+01,  -1.03039120e+05,  5.42788842e+07]
        self.dadT=0.0                 #da/dT in m/K
        
        #######Volatile masses in terms of oceans
        ocean=1.6e21 #kg one ocean mass
        
        self.CO2overH2O=1.0/(HoverC*2.45)
        self.H2OMASS=noceans*ocean
        self.CO2MASS=self.CO2overH2O*self.H2OMASS
        self.noceansH2O=noceans
        self.noceansCO2=self.CO2MASS/ocean
        self.HoverC=HoverC
        
        # Initial Martian REE concentrations
        # from Lodder and Fegley, 1997,  Table II
        # This includes the core, don't use for
        # post core formation
        ###########################################
        self.Rb_crust_mantle_core_conc = 2.7e-6
        self.Ba_crust_mantle_core_conc = 4.3e-6
        self.Th_crust_mantle_core_conc = 44.0e-9
        self.U_crust_mantle_core_conc  = 12.6e-9
        self.Ta_crust_mantle_core_conc = 0.023e-6
        self.K_crust_mantle_core_conc  = 730.0e-6
        self.La_crust_mantle_core_conc = 0.32e-6
        self.Ce_crust_mantle_core_conc = 0.89e-6
        self.P_crust_mantle_core_conc  = 1100.0e-6
        self.Sr_crust_mantle_core_conc = 10.7e-6
        self.Nd_crust_mantle_core_conc = 0.67e-6
        self.Sm_crust_mantle_core_conc = 0.2e-6
        self.Zr_crust_mantle_core_conc = 6.5e-6
        self.Hf_crust_mantle_core_conc = 0.18e-6
        self.Eu_crust_mantle_core_conc = 0.078e-6
        self.Gd_crust_mantle_core_conc = 0.31e-6
        self.Tb_crust_mantle_core_conc = 0.055e-6
        self.Dy_crust_mantle_core_conc = 0.36e-6
        self.Y_crust_mantle_core_conc  = 2.2e-6
        self.Er_crust_mantle_core_conc = 0.24e-6
        self.Tm_crust_mantle_core_conc = 0.04e-6
        self.Yb_crust_mantle_core_conc = 0.22e-6
        self.Lu_crust_mantle_core_conc = 0.035e-6
        ###############################################
        # Initial Martian crust and mantle abundance
        # From Lodder and Fegley, 1997, Table VI
        ##############################################
        self.Rb_init_conc = 3.5e-6
        self.Ba_init_conc = 5.4e-6
        self.Th_init_conc = 56.0e-9
        self.U_init_conc  = 16.0e-9
        self.Ta_init_conc = 29.0e-9
        self.K_init_conc  = 920.0e-6
        self.La_init_conc = 0.4e-6
        self.Ce_init_conc = 1.12e-6
        self.P_init_conc  = 740.0e-6
        self.Sr_init_conc = 13.5e-6
        self.Nd_init_conc = 0.85e-6
        self.Sm_init_conc = 0.25e-6
        self.Zr_init_conc = 8.3e-6
        self.Hf_init_conc = 0.229e-6
        self.Eu_init_conc = 0.099e-6
        self.Gd_init_conc = 0.395e-6
        self.Tb_init_conc = 0.069e-6
        self.Dy_init_conc = 0.454e-6
        self.Y_init_conc  = 2.8e-6
        self.Er_init_conc = 0.3e-6
        self.Tm_init_conc = 0.05e-6
        self.Yb_init_conc = 0.277e-6
        self.Lu_init_conc = 0.044e-6
        self.all_REE_init_conc=np.array([self.Rb_init_conc, self.Ba_init_conc, \
                    self.Th_init_conc, self.U_init_conc, self.Ta_init_conc,\
                    self.K_init_conc,self.La_init_conc,self.Ce_init_conc, \
                    self.P_init_conc, self.Sr_init_conc, self.Nd_init_conc, \
                    self.Sm_init_conc, self.Zr_init_conc, self.Hf_init_conc, \
                    self.Eu_init_conc, self.Gd_init_conc, self.Tb_init_conc, \
                    self.Dy_init_conc, self.Y_init_conc, self.Er_init_conc,  \
                    self.Tm_init_conc, self.Yb_init_conc,  self.Lu_init_conc ])
        ###########################
        # Trapped melt fraction related terms
        # compaction time for length 50 km and velocity 0.1m/yr
        #self.tau=1.577e13            #s
        # compaction time for length 10 km and velocity 0.01m/yr
        self.tau=3.15e13             #s
        self.phic=0.3                #dimensionless disaggregation melt fraction
        self.deltaT=100              #K difference between solidus and freezing front
        self.Ftl=0.0                 #Initial trapped melt fraction
        
        ##########################################################
        # Arrays for storing values of time marching
        # Start by declaring a large array size, these can later be resized
        # using the function discard zeros
        self.n=nsteps
        self.T=np.zeros(nsteps)      #oC temperature
        self.a=np.zeros(nsteps)      #m Radius of RM
        self.t=np.zeros(nsteps)      #s time
        self.MMO=np.zeros(nsteps)    #kg Mass of MO
        self.MRM=np.zeros(nsteps)    #kg Mass of residual mantle
        self.Ftl=np.zeros(nsteps)    #dimensionless trapped melt fraction
        self.dadT=np.zeros(nsteps)   #m/K rate of growth       
        self.MH2OMO=np.zeros(nsteps) #kg mass of H2O in MO
        self.MCO2MO=np.zeros(nsteps) #kg mass of CO2 in MO
        self.MH2ORM=np.zeros(nsteps) #kg mass of H2O in RM
        self.MCO2RM=np.zeros(nsteps) #kg mass of CO2 in RM
        self.CH2OMO=np.zeros(nsteps) #dimensionless conc. of H2O in MO
        self.CCO2MO=np.zeros(nsteps) #dimensionless conc. of CO2 in MO
        self.CH2ORM=np.zeros(nsteps) #dimensionless conc. of H2O in RM
        self.CCO2RM=np.zeros(nsteps) #dimensionless conc. of CO2 in RM
        self.PH2O=np.zeros(nsteps)   #PA Pressure of H2O in PA
        self.PCO2=np.zeros(nsteps)   #PA Pressure of CO2 in PA
        self.MH2OPA=np.zeros(nsteps) #kg Mass of H2O in PA
        self.MCO2PA=np.zeros(nsteps) #kg Mass of CO2 in PA
        self.REE_conc_RM = np.zeros((23,nsteps)) # Array for REE concentration
        # The REES are stored in the order:
        # Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu
        self.MREE   = np.zeros(23)          #Total Mass of REEs
        self.MREERM = np.zeros((23,nsteps)) # REE masses in the RM
        self.CREERM = np.zeros((23,nsteps)) # REE concentration in the RM
        self.MREEMO = np.zeros((23,nsteps)) # REE masses in the MO
        self.CREEMO = np.zeros((23,nsteps)) # REE concentration in the MO
        
       
        ############################################################
        #Set the first value of temperature to initial temperature
        self.T[0]=self.T_init
        self.MRM[0],self.MMO[0]=self.masses(self.a[0])
        #Initially all volatiles are in the MO
        self.MH2OMO[0]  = self.H2OMASS
        self.MCO2MO[0]  = self.CO2MASS
        self.CH2OMO[0]  = self.H2OMASS/self.MMO[0]
        self.CCO2MO[0]  = self.CO2MASS/self.MMO[0]
        self.CH2ORM[0]  = 0.0
        self.CCO2RM[0]  = 0.0
        # All REE in the MO
        self.MREE[:]    = self.all_REE_init_conc*self.MMO[0]
        self.MREEMO[:,0]= self.MREE[:]
        self.CREEMO[:,0]= self.all_REE_init_conc
        self.CREERM[:,0]= 0.0
        self.MREERM[:,0]= 0.0

    def kg2GELm(self,M,rho_fluid=1.0e3):
        """
        This function converts mass of water in kg to
        Global Equivalent Layer in m, height of a water
        column covering the surface of Mars in meters. The
        volume of a thin fluid shell of height h can be approximated
        by the formula h = mass/area/rho_fluid
        Input:
        M          : mass of fluid envelop
        rho_fluid  : density of fluid envelop
        Output:
        h          : Global Equivalent Layer (m)
        """
        h = M/self.area/rho_fluid
        return(h)
        
    def discard_zeros(self,ind):
        """This function resizes all arrays in the object to the value ind"""
        self.T=np.resize(self.T,ind)   
        self.a=np.resize(self.a,ind)
        self.heatflux=np.resize(self.heatflux,ind)
        self.t=np.resize(self.t,ind)
        #self.tma = np.resize(self.tma,ind)
        self.MMO=np.resize(self.MMO ,ind) 
        self.MRM=np.resize(self.MRM ,ind) 
        self.Ftl=np.resize(self.Ftl ,ind) 
        self.dadT=np.resize(self.dadT ,ind)
        self.MH2OMO=np.resize(self.MH2OMO ,ind) 
        self.MCO2MO=np.resize(self.MCO2MO ,ind) 
        self.MH2ORM=np.resize(self.MH2ORM ,ind) 
        self.MCO2RM=np.resize(self.MCO2RM ,ind) 
        self.CH2OMO=np.resize(self.CH2OMO ,ind) 
        self.CCO2MO=np.resize(self.CCO2MO ,ind) 
        self.CH2ORM=np.resize(self.CH2ORM ,ind) 
        self.CCO2RM=np.resize(self.CCO2RM ,ind) 
        self.PH2O=np.resize(self.PH2O ,ind)   
        self.PCO2=np.resize(self.PCO2 ,ind)   
        self.MH2OPA=np.resize(self.MH2OPA ,ind)
        self.MCO2PA=np.resize(self.MCO2PA ,ind)
        self.REE_conc_RM = np.resize(self.REE_conc_RM,(23,ind))
        #self.MREEMO= np.resize(self.MREEMO.T,(23,ind))
        #self.CREEMO= np.resize(self.CREEMO.T,(23,ind))
        #self.CREERM= np.resize(self.CREERM.T,(23,ind))
        #self.MREERM= np.resize(self.MREERM.T,(23,ind))
    def depth2PGPa(self,r):
        """Converts radius in m to pressure in GPa
        enter radius in meters
        Input  : 
        r      : radius in m
        Output : Pressure in GPa
        """
        z=self.radius-r
        P=self.rho*self.g*z
        PGPa=P*1.0e-9
        return(PGPa)
    def solidus(self, r1):
        """Returns solidus temperature
        enter radius of mantle in m
        Equation 2 of Elkins-Tanton, 2008"""
        Xliq=0.3
        r=r1*1.0e-3
        T = -1.963e-10*r**4+1.694e-6*r**3-0.00533*r*r+6.884*r\
            -830.0-(6/(0.2*Xliq+0.025))
        return(T)
    def liquidus(self,r):
        """This function outlines
        the liquidus of the Earth modified
        for Mars. We assume that the liquidus
        of Mars is depressed from that of the 
        Earth by a uniform amount due to the 
        presence of Fe in the innterior of
        Mars. The uniform amount is the difference
        between surface values of the HMH17 and
        ET08 solidi"""
        PGPa=self.depth2PGPa(r)
        self.PGPa=PGPa
        dT=1661.2-self.solidus(self.radius)
        if type(r) is float:
            if PGPa<20.0:
                dum2 = 1982.1*(PGPa/6.594+1.0)**(1.0/5.374)-dT
            else:
                dum2=78.74*(PGPa/4.054e-3 + 1.0)**(1.0/2.44)-dT
        else:
            n=r.shape[0]        
            dum2=np.zeros(n)
            for ii in range(0,n):
                if PGPa[ii]<20.0:
                    dum2[ii] = 1982.1*(PGPa[ii]/6.594+1.0)**(1.0/5.374)-dT
                else:
                    dum2[ii]=78.74*(PGPa[ii]/4.054e-3 + 1.0)**(1.0/2.44)-dT
        return(dum2)
    def freezing_front(self,r):
        """Calculates the freezing front between the
        solidus and the liquidus, enter radius in m"""
        Tsol=self.solidus(r)
        Tliq=self.liquidus(r)
        Tfront=Tsol+0.3*(Tliq-Tsol)
        return(Tfront)
        
    def adiabat(self,T0,r):
        """Returns the adiabatic temperature profile
        enter potential temperature T0 and radius r
        in m"""
        depth=self.radius-r #m depth from surface
        T=T0+self.gam*depth
        return(T)
    def dadT_analytical(self,T):
        """Calculates the value of da/dT from the fit"""
        #p=self.fit_coefficients
        p_hi=self.aT_fit_hi
        p_lo=self.aT_fit_lo
        if T<1350.0 or T>2070.0:
            dadT=0.0
        elif T>1350.0 and T<1820.0:
             #a_low=p_low[0]*T_low**3+p_low[1]*T_low**2+p_low[2]*T_low+p_low[3]
            dadT=3.0*p_lo[0]*T**2+2.0*p_lo[1]*T+p_lo[2]
        else:
            #a_high=p_high[0]*T_high**2+p_high[1]*T_high+p_high[2]
            dadT=2.0*p_hi[0]*T+p_hi[1]
            
        return(dadT)
    def dadT_analytical_old(self,T):
        """Calculates the value of da/dT from the fit"""
        p=self.fit_coefficients
        if T<1590.0 or T>2120.0:
            dadT=0.0
        else:
            dadT=5.0*p[0]*T**4+4.0*p[1]*T**3+3.0*p[2]*T**2\
                   +2.0*p[3]*T+p[4]
        return(dadT)
    
    def solid_find(self,r,T0):
        """This function defines the equation for intersection
        between a temperature profile
        enter radius r in m
        """       
        Tsol=self.freezing_front(r)   #Freezing front array with CMB as first index
        T=self.adiabat(T0,r)
        temp=Tsol-T
        return(temp)
        
    def solid_radius(self,T0,guess=1.8e6):
        """Calculates the solid radius for a potential temperature
        Returns the radius of residual mantle in m"""
        # WARNING!! The result for a is sensitive to the initial
        # guess in the fsolve call. 1.5e6-1.9e6 seems to work, a higher
        # or lower value can give negative a
        a=scipy.optimize.fsolve(self.solid_find,guess,T0)-self.core
        if a < 0.0:
            a=0.0
        return(a)
    def masses(self,a):
        """Calculates total mass of two reservoirs
        in kg. Enter radius of the RM in m"""        
        MRM=4.0*np.pi*self.rho*((a+self.core)**3-self.core**3)/3.0
        MMO=4.0*np.pi*self.rho*(self.radius**3-(a+self.core)**3)/3.0
        return(MRM,MMO)
    def partition_coefficients(self,a):
        """Calculates the partition coefficient
        of water between solid and melt. Enter
        RM radius in m"""
        #Depth of the FF in km
        depth=(self.radius-(a+self.core))*1.0e-3
        if depth>=991.0:
            self.DH2O=1.0e-3
            self.DCO2=1.0e-3
        else:
            self.DH2O=0.1
            self.DCO2=1.0e-3
    def REE_D(self, r):
        """Calculates the mineral modes in the crystallizing
        solid for a given radius of the resdiual mantle.
        The mineral modes are taken from Borg and Draper(2003), Table 3.
        The solid composition is an array consisting of wt% of the minerals
        in the following order
        ol-opx-cpx-gt-ilm-titanite-phlogopite-amphibole
        This order is the same as the partition coefficients of REE
        Then the partition coefficient of the REEs are calculated and returned
        """
        PGPa=self.depth2PGPa(r)
        if PGPa >= 15.0:
            self.bulk_composition=np.array([0.327,0.388,0.2,0.082,0.004,0.0,0.0,0.0])
        elif PGPa >10.0 and PGPa<15.0:
            self.bulk_composition=np.array([0.235,0.306,0.38,0.07,0.009,0.0,0.0,0.0])
        else:
            self.bulk_composition=np.array([0.699,0.117,0.173,0.006,0.005,0.0,0.0,0.0])
        D_Rb=np.dot(self.bulk_composition,self.Rb)
        D_Ba=np.dot(self.bulk_composition,self.Ba)
        D_Th=np.dot(self.bulk_composition,self.Th)
        D_U=np.dot(self.bulk_composition,self.U)
        D_Ta=np.dot(self.bulk_composition,self.Ta)
        D_K=np.dot(self.bulk_composition,self.K)
        D_La=np.dot(self.bulk_composition,self.La)
        D_Ce=np.dot(self.bulk_composition,self.Ce)
        D_P=np.dot(self.bulk_composition,self.P)
        D_Sr=np.dot(self.bulk_composition,self.Sr)
        D_Nd=np.dot(self.bulk_composition,self.Nd)
        D_Sm=np.dot(self.bulk_composition,self.Sm)
        D_Zr=np.dot(self.bulk_composition,self.Zr)
        D_Hf=np.dot(self.bulk_composition,self.Hf)
        D_Eu=np.dot(self.bulk_composition,self.Eu)
        D_Gd=np.dot(self.bulk_composition,self.Gd)
        D_Tb=np.dot(self.bulk_composition,self.Tb)
        D_Dy=np.dot(self.bulk_composition,self.Dy)
        D_Y=np.dot(self.bulk_composition,self.Y)
        D_Er=np.dot(self.bulk_composition,self.Er)
        D_Tm=np.dot(self.bulk_composition,self.Tm)
        D_Yb=np.dot(self.bulk_composition,self.Yb)
        D_Lu=np.dot(self.bulk_composition,self.Lu)
        D_REE=np.array([D_Rb,D_Ba,D_Th,D_U,D_Ta,D_K,D_La,D_Ce,D_P,\
                        D_Sr,D_Nd,D_Sm,D_Zr, D_Hf, D_Eu,D_Gd,D_Tb, D_Dy,\
                        D_Y,D_Er, D_Tm, D_Yb,D_Lu])
        return(D_REE)
                             
    def REE_masses(self,r,Ftl,dMRM,MREE_prev,CREE_prev):
        """Calculates the concentration
        of REE for a given radius of RM,r,
        and the trapped melt fraction, Ftl
        The formula is
        M_RM^Z = M^Z_prev+((1-Ftl)*DZ+Ftl)*c^Z_prev*dM_RM
        where 
        M_RM^Z    = mass of Z (REE) in the RM at the current time step
        M^Z_prev  = mass of Z (REE) in the RM at the last time step
        DZ = Partitioning coefficient for the given depth
        c0 = Concentration of Z in MO at the current time step
        Ftl = trapped melt fraction
        
        """
        D_REE=self.REE_D(r)
        dMRMREE=np.zeros(23)
        MRMREE_now=np.zeros(23)
        
        for kk in range(0,23):
            dMRMREE[kk]=dMRM*((1.0-Ftl)*D_REE[kk]+Ftl)*CREE_prev[kk]
            MRMREE_now[kk]=MREE_prev[kk]+dMRMREE[kk]

        return (MRMREE_now)
        
    def emmissivity_initial(self,H2O=0.4,CO2=0.5):
        """This function calculates the thermal
        emmissivity of the atmosphere from the current
        mass of H2O and CO2 in the atmosphere, from 
        Elkins-Tanton (2008)
        Input   :
        H2O     : conc of H2O in the atmosphere kg
        CO2     : conc of CO2 in the atmosphere kg
        
        Output  :
        epsilon : Thermal emmissivity SI units
        """
        
        p0     = 101325.0 #reference pressure, Pa
        a=np.array([0.3, 0.05]) #;
        n=np.array([0.52,0.45])#;
        k0=np.array([0.01,0.05])# m^2/kg absorption coefficient
        c= np.array([H2O,CO2]) #H2O and CO2 concentrations, wt%
        p=np.array([0.0,0.0])
        M=np.array([0.0,0.0])
        tau=np.array([0.0,0.0])
        for ii in range (0,2):        
            p[ii]=((c[ii]-a[ii])/2.08e-4)**(1.0/n[ii])
            M[ii]=4.0*np.pi*p[ii]*(self.radius**2)/self.g 
            tau[ii]=(3.0*M[ii]/8/np.pi/(self.radius**2))\
                *np.sqrt(k0[ii]*self.g/3.0/p0)
            
        epsilon = 2.0/(tau[0]+tau[1]+2.0)
        return(epsilon)
    def emmissivity_var(self,MCO2,MH2O):
        """
        This function calculates the variable emmissivity as a function
        of H2O and CO2 in the protoatmosphere, using the formula from 
        Elkins-Tanton (2008)
        Input:
        MH2O    : Mass of H2O in the PA kg
        MCO2    : Mass of CO2 in the PA kg
        Output:
        epsilon : Thermal emmissivity SI units
        """
        p0     = 101325.0 #reference pressure, Pa
        a=np.array([0.3, 0.05]) #;
        n=np.array([0.52,0.45])#;
        k0=np.array([0.01,0.05])# m^2/kg absorption coefficient
        M=np.array([MH2O,MCO2])
        tau=np.array([0.0,0.0])
        for ii in range (0,2):        
            tau[ii]=(3.0*M[ii]/8/np.pi/(self.radius**2))\
                *np.sqrt(k0[ii]*self.g/3.0/p0)
                        
        epsilon = 2.0/(tau[0]+tau[1]+2.0)
        return(epsilon)
        
    def PH2O_calculate(self,dMRM,MRM_prev,F,D,T,MMO):
        """Calculates pressure of H2O from 
        the masses on input
        dMRM : incremental increase in RM mass
        MRM_prev : Mass of H2O in RM from last time step
        F    : Trapped melt fraction at the time step
        D    : Water partition coefficient at that time step
        T    : Temperature at the time step
        """
        #First update the reaction coefficient for temperature
        
        term1=self.H2OMASS-MRM_prev
        term2=self.KH2O*T*self.omega_H2O/self.omega_Moore
        term3=(1.0-F)*D+F
        term4=self.area/self.g
        term5=term2*(MMO+dMRM*term3)
        P= term1/(term4+term5)
        c=P*term2
        
        return(P,c)
    def PCO2_calculate(self,dMRM,MRM_prev,F,D,T,MMO):
        """Calculates pressure of CO2 from 
        the masses on input
        dMRM : incremental increase in RM mass
        MRM_prev : Mass of CO2 in RM from last time step
        F    : Trapped melt fraction at the time step
        D    : CO2 partition coefficient at that time step
        T    : Temperature at the time step
        """
        #First update the reaction coefficient for temperature
        
        term1=self.CO2MASS-MRM_prev
        self.KCO2=self.redox_factor*self.K0\
                   *np.exp(-self.dH*(1.0/T -1.0/self.T0)/self.R)
        
        #self.calc_kco2(T)
        term2=self.KCO2*self.omega_CO2/self.omegaM
        term3=(1.0-F)*D+F
        term4=self.area/self.g
        term5=term2*(MMO+dMRM*term3)
        P= term1/(term4+term5)
        c=P*term2
        
        return(P,c)
    def heat_flux(self,Tm,e):
        """Calculates radiative heat flux"""
        Tinf = 0.0   # Temperature outside the planet
        F = self.sigma*e*(Tm**4 - (Tinf)**4)
        
        # IF the temperature is less than the solar flux
        # then set it to solar flux
        
        if F <  self.Fs:
            F=self.Fs
        return(F)
        
    def rhs(self,t,Tm,a,dadT,e):
        """This function evaluates the rhs of the thermal
        evolution equation:
        dT/dt = -R**2*(sigma*e*(T**4-Tinf**4)-Fs)\
        /(rho*cp*(R**3-a**3)/3 -rho*T*dS*dadT*a**2)

        Input   :
        t       : Time s
        Tm      : Mantle temperature C
        a       : Radius of residual mantle m
        da/dT   : Rate of change of radius with T, m/K
        e       : Thermal emmissivity SI units
        
        Output  :
        dTdt    : the value of dT/dt for the time step K/s
        """
        F     = self.heat_flux(Tm,e)
        term1 = -(self.radius**2)*(F-self.Fs)
        term2 = self.rho*Tm*self.dS*(a**2)*dadT
        term3 = self.rho*self.Cp*(self.radius**3 - a**3)/3.0
        dTdt  = term1/(term3-term2)
        return (dTdt)
    def rk4(self,t,Tm,a,dadT,e,dt):
        """Integrates the thermal evolution equation:
        dT/dt = -R**2*(sigma*e*(T**4-Tinf**4)-Fs)\
        /(rho*cp*(R**3-a**3)/3 -rho*T*dS*dadT*a**2)

        Input   : parameters from the i-th time step
        t       : Time s
        Tm      : Mantle temperature C
        a       : Radius of residual mantle m
        da/dT   : Rate of change of radius with T, m/K
        e       : Thermal emmissivity SI units
        dt      : Length of time step
        
        Output  :
        Tmip1   : Value of Tm at timestep i+1 
        """
        k1    = self.rhs(t,Tm,a,dadT,e)
        k2    = self.rhs(t+0.5*dt,Tm+0.5*k1*dt,a,dadT,e)
        k3    = self.rhs(t+0.5*dt,Tm+0.5*k2*dt,a,dadT,e)
        k4    = self.rhs(t+dt,Tm+k3*dt,a,dadT,e)
        Tmip1 = Tm + dt*(k1+2.0*(k2+k3)+k4)/6.0
        return(Tmip1)
    def output_filenames(self,const=False):
        """Generates file names for output
        The optional input const is for a constant trapped
        melt fraction"""
        if const == False:
            prefix = './data/Mars_noceansH2O_'+str(self.noceansH2O)+'_HC_'+str(self.HoverC)\
                 +'_redox_factor_'+str(self.redox_factor)
        else:
            prefix = './data/Mars_const_noceansH2O_'+str(self.noceansH2O)+'_HC_'+str(self.HoverC)\
                 +'_redox_factor_'+str(self.redox_factor)
            
        
        vars1  = '_t_T_rad_F_MMO_MRM_Ftl'
        vars2  = '_t_PCO2_PH2O'
        vars3  = '_t_MH2ORM_MH2OMO_MH2OPA_CH2ORM_CH2OMO'
        vars4  = '_t_MCO2RM_MCO2MO_MCO2PA_CCO2RM_CCO2MO'
        vars5  = 'CREERM'
        vars6  = 'CREEMO'
        vars7  = '_final'
        vars8  = 'CREE'
        vars_temp1='80pct'
        vars_temp2='90pct'
        #File extension
        ext    = '.csv'
        f1 = prefix+vars1+ext
        f2 = prefix+vars2+ext
        f3 = prefix+vars3+ext
        f4 = prefix+vars4+ext
        f5 = prefix+vars5+ext
        f6 = prefix+vars6+ext
        f7 = prefix+vars7+ext
        f8 = prefix+vars8+ext
        f9 = prefix+vars5+vars_temp1+ext
        f10= prefix+vars6+vars_temp2+ext
        return(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)
    def create_output(self,const=False):
        """This function writes the output
        of simulations into text files"""
        if const == False:
            fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,\
                fname9,fname10=self.output_filenames()
        else:
            fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8,\
            fname9,fname10=self.output_filenames(const=True)
            
        
        # Text for headers in the csv file
        txt1  = '# Time (Ma), Temperature (C), Radius of RM (km), Heat flux (W/m^2), Mass of MO (kg),Mass of RM (kg), Trapped melt fraction'
        txt2  = '# Time (Ma), CO2 pressure in PA (Pa), H2O pressure in PA (Pa)'
        txt3  = '# Time (Ma), H2O mass in RM (kg), H2O mass in MO (kg), H2O mass in PA (kg), H2O conc. in RM, H2O conc. in MO'
        txt4  = '# Time (Ma), CO2 mass in RM (kg), CO2 mass in MO (kg), CO2 mass in PA (kg), CO2 conc. in RM, CO2 conc. in MO'
        txt5  = '# Evolution of  REE concentration in RM in the order : Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu'
        txt6  = '#Evolution of REE concentration in MO in the order : Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu'
        txt7  = '# Time (Ma), CO2 pressure in PA (Pa), CO2 conc. in RM, CO2 conc. in MO, H2O pressure in PA (Pa), H2O conc. in RM, H2O conc. in MO'
        txt8  = '#Final REE concentration in MO(first col) and RM(second col) in the order : Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu'
        txt9  = '#80% crystallization REE concentration in MO(first col) and RM(second col) in the order : Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu'
        txt10= '#90% crystallization REE concentration in MO(first col) and RM(second col) in the order : Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu'
        
        temp1=np.array([self.tma,self.T,self.akm,self.heatflux,self.MMO,\
                        self.MRM,self.Ftl])
        np.savetxt(fname1,temp1.T,delimiter=",",header=txt1)

        temp2=np.array([self.tma,self.PCO2,self.PH2O])
        np.savetxt(fname2,temp2.T,delimiter=",",header=txt2)

        temp3=np.array([self.tma,self.MH2ORM,self.MH2OMO,self.MH2OPA,\
                        self.CH2ORM,self.CH2OMO])
        np.savetxt(fname3,temp3.T,delimiter=",",header=txt3)

        temp4=np.array([self.tma,self.MCO2RM,self.MCO2MO,self.MCO2PA,\
                        self.CCO2RM,self.CCO2MO])
        np.savetxt(fname4,temp4.T,delimiter=",",header=txt4)

        temp5=np.array([self.CREERM])
        #np.savetxt(fname5,temp5.T,delimiter=",",header=txt5)

        temp6=np.array([self.CREEMO])
        #np.savetxt(fname6,temp6.T,delimiter=",",header=txt6)

        n = np.shape(self.tma)[0]
        temp7= np.array([self.tma[n-1],self.PCO2[n-1],self.CCO2RM[n-1],\
                         self.CCO2MO[n-1],self.PH2O[n-1],self.CH2ORM[n-1]\
                         ,self.CH2OMO[n-1]])
        np.savetxt(fname7,temp7,delimiter=",",header=txt7)

        temp8=np.zeros((23,2))
        temp8[:,0]=self.CREEMO[:,n-1]
        temp8[:,1]=self.CREERM[:,n-1]
        np.savetxt(fname8,temp8,delimiter=",",header=txt8)

        temp9=np.zeros((23,2))
        temp9[:,0]=self.CREE_MO80
        temp9[:,1]=self.CREE_RM80
        np.savetxt(fname9,temp9,delimiter=",",header=txt9)

        temp10=np.zeros((23,2))
        temp10[:,0]=self.CREE_MO90
        temp10[:,1]=self.CREE_RM90
        np.savetxt(fname10,temp10,delimiter=",",header=txt10)

    def write_all_REE(self,const=False):
        """Writes REE concentrations for all time steps"""
        if const == False:
            fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8\
            =self.output_filenames()
        else:
            fname1,fname2,fname3,fname4,fname5,fname6,fname7,fname8\
            =self.output_filenames(const=True)
    
        txt5  = '# Evolution of  REE concentration in RM in the order :\
        Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu'
        txt6  = '#Evolution of REE concentration in MO in the order : \
        Rb,Ba,Th,U,Ta,K,La,Ce,P,Sr,Nd,Sm,Zr,Hf,Eu,Gd,Tb,Dy,Y,Er,Tm,Yb,Lu'
        temp5=np.array([self.CREERM])
        np.savetxt(fname5,temp5.T,delimiter=",",header=txt5)
        temp6=np.array([self.CREEMO])
        np.savetxt(fname6,temp6.T,delimiter=",",header=txt6) 
        
    def load_object_from_file(self):
       """This function loads data from simulated values.
       Call this only for visualization
       """
       fname1,fname2,fname3,fname4,fname5,fname6=self.output_filenames()
       data1=np.loadtxt(fname1,delimiter=',')
       self.tma      = data1[:,0]
       n=np.shape(self.tma)[0]
       self.T        = data1[:,1]
       self.akm      = data1[:,2]
       #self.heatflux = data1[:,3]
       self.MMO      = data1[:,4]
       self.MRM      = data1[:,5]
       #self.Ftl      = data1[:,6]
       data2=np.loadtxt(fname2,delimiter=',')
       self.PCO2     = data2[:,1]
       self.PH2O     = data2[:,2]
       #data3=np.loadtxt(fname3,delimiter=',')
       #self.MH2ORM   = data3[:,1]
       #self.MH2OMO   = data3[:,2]
       #self.MH2OPA   = data3[:,3]
       #self.CH2ORM   = data3[:,4]
       #self.CH2OMO   = data3[:,5]
       data4=np.loadtxt(fname4,delimiter=',')
       #self.MCO2RM   = data4[:,1]
       #self.MCO2MO   = data4[:,2]
       #self.MCO2PA   = data4[:,3]
       self.CCO2RM   = data4[:,4]
       self.CCO2MO   = data4[:,5]
       #data5=np.loadtxt(fname5,delimiter=',')
       #self.CREERM   = data5.T
       #data6=np.loadtxt(fname6,delimiter=',')
       #self.CREEMO   = data6.T
       return(n)
    def time_marching(self,tma,const_Ftl=False):
        """Time marching and calculation of masses
        enter final crystallization time in Ma"""
        
        nsteps=self.n
        #Convert crystallization time from Ma to s
        tfinal=tma*1.0e6*365*24*3600 #time in s
        dt=tfinal/nsteps
        dT=self.T_init-self.T_surf
        dTdt=dT/tfinal
        t=0
        e=self.emmissivity_initial()
        self.dTdt=0.0*self.T
        self.heatflux=0.0*self.T
        
        
        for ii in range (1,nsteps):
            self.t[ii]=t+ii*dt
            self.T[ii] = self.rk4(self.t[ii],self.T[ii-1],self.a[ii],self.dadT[ii],e,dt)
            # Solve for the RM radius
            # Use the value at the last step as initial
            # guess
            self.a[ii]=self.solid_radius(self.T[ii],guess=\
                                         self.a[ii-1]+self.core)
            self.dTdt[ii]=self.rhs(self.t[ii],self.T[ii-1],self.a[ii],self.dadT[ii],e)
            # Calculate masses of the reservoirs
            self.MRM[ii],self.MMO[ii]=self.masses(self.a[ii])
            self.dadT[ii]=self.dadT_analytical(self.T[ii])
            # Calculate the heat flux
            self.heatflux[ii]=self.heat_flux(self.T[ii],e)
            #calculate trapped melt fraction
            if const_Ftl == False:
                self.Ftl[ii]=self.phic*self.tau*dTdt/self.deltaT
            else:
                self.Ftl[ii]=0.01
                
            
            if self.Ftl[ii]>=0.3:
                self.Ftl[ii]=0.3
            #Get the partition coefficents
            #This function won't work if an array is passed as an argument
            self.partition_coefficients(self.a[ii])
            #change in residual mantle mass
            dMRM=self.MRM[ii]-self.MRM[ii-1]
            ################################################
            ## Calculate water concentration and masses
            ################################################
            #This following step calculates the pressure
            #of H2O from masses, but works for less than
            #3 oceans of water in the starting composition
            self.PH2O[ii],self.CH2OMO[ii]=self.PH2O_calculate(dMRM,\
             self.MH2ORM[ii-1],self.Ftl[ii],self.DH2O,self.T[ii],self.MMO[ii])
            #Now update the mass of H2O in self.CH2OMO[ii]each reservoir
            self.MH2OMO[ii]=self.MMO[ii]*self.CH2OMO[ii]
            self.MH2OPA[ii]=self.area*self.PH2O[ii]/self.g
            dMRMH2O=dMRM*((1.0-self.Ftl[ii])*self.DH2O+self.Ftl[ii])\
                     *self.CH2OMO[ii]
            self.MH2ORM[ii]=self.MH2ORM[ii-1]+dMRMH2O
            #################################################
            ## Calculate CO2 mass and pressures
            #################################################
            self.PCO2[ii],self.CCO2MO[ii]=self.PCO2_calculate(dMRM,\
             self.MCO2RM[ii-1],self.Ftl[ii],self.DCO2,self.T[ii],self.MMO[ii])
            #Now update the mass of CO2 in self.CH2OMO[ii]each reservoir
            self.MCO2MO[ii]=self.MMO[ii]*self.CCO2MO[ii]
            self.MCO2PA[ii]=self.area*self.PCO2[ii]/self.g
            dMRMCO2=dMRM*((1.0-self.Ftl[ii])*self.DCO2+self.Ftl[ii])\
                     *self.CCO2MO[ii]
            self.MCO2RM[ii]=self.MCO2RM[ii-1]+dMRMCO2
            #####################################################
            ## Calculate REE Mass in the RM
            #####################################################
            self.MREERM [:,ii] =self.REE_masses(self.a[ii], self.Ftl[ii],dMRM,\
                                self.MREERM[:,ii-1],self.CREEMO[:,ii-1])
            self.MREEMO[:,ii] = self.MREE[:]-self.MREERM [:,ii]
            #####################################################
            ## Calculate the CO2, H2O, and REE concentration in the RM
            ######################################################
            if self.MRM[ii]==0:
                self.CH2ORM[ii]=0.0
                self.CCO2RM[ii]=0.0
                self.CREERM[:,ii]=0.0
            else:
                self.CH2ORM[ii]=self.MH2ORM[ii]/self.MRM[ii]
                self.CCO2RM[ii]=self.MCO2RM[ii]/self.MRM[ii]
                self.CREERM[:,ii]=self.MREERM[:,ii]/self.MRM[ii]

            ##########################################################
            ## Calculate REE concentration in the MO
            #########################################################
            if self.MMO[ii]==0:
                self.CREEMO[:,ii]=0.0
            else:
                self.CREEMO[:,ii]=self.MREEMO[:,ii]/self.MMO[ii]

            #######################################################
            ## Update the emmissivity for next time step
            ######################################################
            
            e = self.emmissivity_var(self.MCO2PA[ii],self.MH2OPA[ii])
            ######################################################
            # check if solidified, then exit
            #####################################################
            percent_solid=100.0*self.a[ii]/(self.radius-self.core)
            #print '% MO solidified: ', percent_solid            
            # Create a special object for creating output
            # of REEs at 80% and 90% crystallization
            #if percent_solid >80.0 and percent_solid<81.0:
            #    self.CREE_RM80=self.CREERM[:,ii]
            #    self.CREE_MO80=self.CREEMO[:,ii]
               
            #Only uncomment the following for special cases    
            #if percent_solid>90.0 and percent_solid<91.0:
            #    self.CREE_RM90 =self.CREERM[:,ii]
            #    self.CREE_MO90=self.CREEMO[:,ii]

            if percent_solid > 99.5:
                # variables for print out
                temp1=self.CH2ORM[ii]*1.0e6
                temp2=self.CCO2RM[ii]*1.0e6
                temp3= self.PCO2[ii]*1.0e-5
                temp4= self.PH2O[ii]*1.0e-5
                temp5=self.t[ii]*1.0e-6/365/24/3600
                break
        print '80% shape and values'
        print np.shape(self.CREE_RM80)
        print self.CREE_RM80
        print '90% shape and values'
        print np.shape(self.CREE_RM90)
        print self.CREE_RM90
        self.discard_zeros(ii-1)    
        self.tma=self.t*1.0e-6/365/24/3600 #Time in Ma
        self.akm=self.a*1.0e-3
        self.create_output(const=const_Ftl)
        # Print out some useful numbers
        print '###########################################'
        print 'Simulation finished after iterations',ii
        print 'Total oceans of H2O:',self.noceansH2O
        print 'Total oceans of CO2:',self.noceansCO2
        print 'H:C ratio',self.HoverC
        print 'Redox factor',self.redox_factor
        print 'H2O conc. (ppm) in RM:',temp1
        print 'CO2 conc.(ppm) in RM:', temp2
        print 'Pressure of CO2 in atmosphere (bars):', temp3
        print 'Pressure of H2O in atmosphere (bars):', temp4
        print 'Time to crystallization (Ma):', temp5
        print '###########################################'
        return(ii-1)
 ###################################################
 #######Functions for testing the code
 ###################################################
        
    def radius_test(self,T0,fignum=1):
        """Tests for proper finding
        of the intersection between the adiabat and the solidus"""
        r=np.linspace(self.core,self.radius,100) #Mantle radius in km
        sol=self.solidus(r)
        Tad=self.adiabat(T0,r)
        akm=(r-self.core)*1.0e-3
        plt.figure()
        plt.subplot(1,2,1)
        plt.plot(sol,akm,'-b')
        plt.plot(Tad,akm,'-r')
        plt.legend(['Solidus','Adiabat'],loc=3)
        plt.xlabel('Solidus Temperature C')
        plt.ylabel('Planetary radius (km)')
        plt.subplot(1,2,2)
        rts=self.solid_find(r,T0)
        plt.plot(rts.T,akm,'-b')
        #plt.plot(r*0.0,akm,'-r')
        plt.xlabel(r'Residual from $T_{ad}-T_{sol}$')
        plt.ylabel('Planetary radius (km)')
    def radius_temperature_analytical(self,plot=False):
        """Calculates the radius as a function of temperature
        and gives a polynomial fit"""
        
        T=np.linspace(self.T_surf,self.T_init,100)
        a=0.0*T
        
        index=np.zeros(100)
        for ii in range (0,99):
            a[ii]=self.solid_radius(T[ii])
            
        for ii in range (0,99):
            # Fit the radius only for nonzero values
            if T[ii]>1350.0 and T[ii] < 1820.0:
                index[ii]=True
            else:
                index[ii]=False
    
        TT=np.extract(index,T)
        aa=np.extract(index,a)
        p_low=np.polyfit(TT,aa,3)
    
        index=np.zeros(100)
        for ii in range (0,99):
            # Fit the radius only for nonzero values
            if T[ii]>1820.0 and T[ii] < 2080.0:
                index[ii]=True
            else:
                index[ii]=False
    
        TT1=np.extract(index,T)
        aa1=np.extract(index,a)
        p_high=np.polyfit(TT1,aa1,2)
        
        
        if plot==True:
            temp_label=[1400,1800,2200,2600]
            plt.plot(T,a*1.0e-3,'s',color='skyblue',markersize=20,alpha=0.7)
            T_high=np.linspace(1821.0,2080.0)
            T_low=np.linspace(1350.0,1819.0)

            a_high=p_high[0]*T_high**2+p_high[1]*T_high+p_high[2]
            a_low=p_low[0]*T_low**3+p_low[1]*T_low**2+p_low[2]*T_low+p_low[3]
            plt.plot(T_high,a_high/1.0e3,'-',color='steelblue',lw=4)
            plt.legend(['Numerical','Fit'],fancybox=True,framealpha=0.7)
            plt.plot(2050,146.0,'o',color='tomato',markersize=20,alpha=0.7)
            plt.plot(1550,1776.0,'o',color='tomato',markersize=20,alpha=0.7)
            plt.plot(T_low,a_low/1.0e3,'-',color='steelblue',lw=4)
            plt.xticks(temp_label)
            plt.ylim(0.0,2000.0)
            plt.xlabel(r'Potential Temperature ($^\mathrm{o}$C)',fontsize=30)
            plt.ylabel(r'Residual mantle radius (km)',fontsize=30)
            print 'Coefficients between 1820 and 2080',p_high
            print 'Coefficients between 1350 and 1820',p_low
            
        
class Mars_read(REE):
    """
    A class for loading some data for a Mars
    object
    """
    def __init__(self,noceans=1.0,HoverC=0.4,redox_factor=1.0,Ftl_const=False,Fcl=False):
        self.HoverC=HoverC
        self.redox_factor = redox_factor # ratio fCO2/fO2
        self.noceansH2O=noceans
        self.H2OMASS=noceans*1.6e21
        self.CO2overH2O=1.0/(HoverC*2.45)
        self.CO2MASS=self.CO2overH2O*self.H2OMASS
        #########################
        f=self.output_filenames(const=Ftl_const)
        data=np.loadtxt(f[6],delimiter=',')
        self.area=144.4e12            # m^2 surface area of Mars
        self.mass = 0.624e24          # kg mass of Mars, Carr, encyclopedia of SS
        self.tma_final   = data[0]
        self.PCO2_final  = data[1]
        self.CCO2RM_final= data[2]
        self.CCO2MO_final= data[3]
        self.PH2O_final  = data[4]
        self.CH2ORM_final= data[5]                         ,
        self.CH2OMO_final= data[6]
        ###############################
        data2=np.loadtxt(f[7],delimiter=',')
        self.CREEMO_final=data2[:,0]
        self.CREERM_final=data2[:,1]
        ####################################
        # This is only for special cases
        if Fcl == True:            
            data3=np.loadtxt(f[8],delimiter=',')
            self.CREEMO_80=data3[:,0]
            self.CREERM_80=data3[:,1]        
            data4=np.loadtxt(f[9],delimiter=',')
            self.CREEMO_90=data4[:,0]
            self.CREERM_90=data4[:,1]
        ###############################################
        # Initial Martian crust and mantle abundance
        # From Lodder and Fegley, 1997, Table VI
        ##############################################
        self.Rb_init_conc = 3.5e-6
        self.Ba_init_conc = 5.4e-6
        self.Th_init_conc = 56.0e-9
        self.U_init_conc  = 16.0e-9
        self.Ta_init_conc = 29.0e-9
        self.K_init_conc  = 920.0e-6
        self.La_init_conc = 0.4e-6
        self.Ce_init_conc = 1.12e-6
        self.P_init_conc  = 740.0e-6
        self.Sr_init_conc = 13.5e-6
        self.Nd_init_conc = 0.85e-6
        self.Sm_init_conc = 0.25e-6
        self.Zr_init_conc = 8.3e-6
        self.Hf_init_conc = 0.229e-6
        self.Eu_init_conc = 0.099e-6
        self.Gd_init_conc = 0.395e-6
        self.Tb_init_conc = 0.069e-6
        self.Dy_init_conc = 0.454e-6
        self.Y_init_conc  = 2.8e-6
        self.Er_init_conc = 0.3e-6
        self.Tm_init_conc = 0.05e-6
        self.Yb_init_conc = 0.277e-6
        self.Lu_init_conc = 0.044e-6
        self.all_REE_init_conc=np.array([self.Rb_init_conc, self.Ba_init_conc, \
                    self.Th_init_conc, self.U_init_conc, self.Ta_init_conc,\
                    self.K_init_conc,self.La_init_conc,self.Ce_init_conc, \
                    self.P_init_conc, self.Sr_init_conc, self.Nd_init_conc, \
                    self.Sm_init_conc, self.Zr_init_conc, self.Hf_init_conc, \
                    self.Eu_init_conc, self.Gd_init_conc, self.Tb_init_conc, \
                    self.Dy_init_conc, self.Y_init_conc, self.Er_init_conc,  \
                    self.Tm_init_conc, self.Yb_init_conc,  self.Lu_init_conc ])
    def output_filenames(self,const=False):
        """Generates file names for output
        The optional input const is for a constant trapped
        melt fraction"""
        if const == False:
            prefix = './data/Mars_noceansH2O_'+str(self.noceansH2O)+'_HC_'+str(self.HoverC)\
                 +'_redox_factor_'+str(self.redox_factor)
        else:
            prefix = './data/Mars_const_noceansH2O_'+str(self.noceansH2O)+'_HC_'+str(self.HoverC)\
                 +'_redox_factor_'+str(self.redox_factor)
            
        
        vars1  = '_t_T_rad_F_MMO_MRM_Ftl'
        vars2  = '_t_PCO2_PH2O'
        vars3  = '_t_MH2ORM_MH2OMO_MH2OPA_CH2ORM_CH2OMO'
        vars4  = '_t_MCO2RM_MCO2MO_MCO2PA_CCO2RM_CCO2MO'
        vars5  = 'CREERM'
        vars6  = 'CREEMO'
        vars7  = '_final'
        vars8  = 'CREE'
        vars_temp1='80pct'
        vars_temp2='90pct'
        #File extension
        ext    = '.csv'
        f1 = prefix+vars1+ext
        f2 = prefix+vars2+ext
        f3 = prefix+vars3+ext
        f4 = prefix+vars4+ext
        f5 = prefix+vars5+ext
        f6 = prefix+vars6+ext
        f7 = prefix+vars7+ext
        f8 = prefix+vars8+ext
        f9 = prefix+vars5+vars_temp1+ext
        f10= prefix+vars6+vars_temp2+ext
        
        return(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)
    
    def load_evolution(self,Ftl_const=False):
        """Loads data from files"""
        f=self.output_filenames(const=Ftl_const)
        data1=np.loadtxt(f[0],delimiter=',')
        self.tma      = data1[:,0]        
        n1=np.shape(self.tma)[0]
        #self.niter=n1
        self.T        = data1[:,1]
        self.akm      = data1[:,2]
        self.heatflux = data1[:,3]
        self.MMO      = data1[:,4]
        self.MRM      = data1[:,5]
        self.Ftl      = data1[:,6]
        data1_p = np.loadtxt(f[1],delimiter=',')
        self.PCO2     = data1_p[:,1]
        self.PH2O     = data1_p[:,2]
        data1_CO2 = np.loadtxt(f[3],delimiter=',')
        #self.MCO2RM   = data1_CO2[:,1]
        #self.MCO2MO   = data1_CO2[:,2]
        #self.MCO2PA   = data1_CO2[:,3]
        self.CCO2RM   = data1_CO2[:,4]
        self.CCO2MO   = data1_CO2[:,5]
        #data5=np.loadtxt(f[4],delimiter=',')
        #self.CREERM   = data5.T
        #data6=np.loadtxt(f[5],delimiter=',')
        #self.CREEMO   = data6.T
        # return(n1)
    def evolution_plot(self):
        """plots the thermal evolution of an object
        """
        self.load_evolution()
        plt.figure(figsize=(12,10))
        ax1=plt.subplot(3,1,1)
        plt.plot(self.tma,self.T,'firebrick',linewidth=4)
        plt.ylabel(r'$T (^\mathrm{o}$C)',fontsize=30)
        
        GEL=self.kg2GELm(self.H2OMASS)/1000.0
        plt.text(0.8,2200,'{:3.2f} GEL (km)'.format(GEL),fontsize=30)
        plt.text(0.8,2000,'H:C ={:3.2f}'.format(self.HoverC),fontsize=30)
        plt.text(0.8,1800,r'$K^\ast$ ={:3.2f}'.format(self.redox_factor),fontsize=30)
        plt.text(0.03,1500,'(a)',fontsize=40,fontweight='bold')
        ax1.set_yticks([1600,2000,2400])
        ax1=plt.subplot(3,1,2)
        plt.plot(self.tma,self.MMO/1.0e23,'salmon',linewidth=4)
        plt.plot(self.tma,self.MRM/1.0e23,'steelblue',linewidth=4)
        plt.legend(['MO','RM'],fancybox=True,loc=4,framealpha=0.5)
        plt.text(0.03,1,'(b)',fontsize=40,fontweight='bold')
        ax1.yaxis.set_ticks_position('right') 
        ax1.yaxis.set_label_position('right') 
        ax1.spines['right'].set_position(('outward', 0))
        plt.ylabel(r'Mass ($10^{23}$ kg)',fontsize=30,labelpad=20)
        ax1.set_yticks([2.0,4.0,6.0])
        ax1=plt.subplot(3,1,3)
        plt.plot(self.tma,self.akm,'forestgreen',linewidth=4)
        plt.ylabel('RM radius (km)',fontsize=30)
        plt.xlabel('Time (Ma)',fontsize=30)
        plt.text(0.03,300,'(c)',fontsize=40,fontweight='bold')
        ax1.set_yticks([500,1000,1500])
    def kg2GELm(self,M,rho_fluid=1.0e3):
        """
        This function converts mass of water in kg to
        Global Equivalent Layer in m, height of a water
        column covering the surface of Mars in meters. The
        volume of a thin fluid shell of height h can be approximated
        by the formula h = mass/area/rho_fluid
        Input:
        M          : mass of fluid envelop
        rho_fluid  : density of fluid envelop
        Output:
        h          : Global Equivalent Layer (m)
        """
        h = M/self.area/rho_fluid
        return(h)
    def GELm2kg(self,h,rho_fluid=1.0e3):
        """
        This function converts mass of water in kg to
        Global Equivalent Layer in m, height of a water
        column covering the surface of Mars in meters. The
        volume of a thin fluid shell of height h can be approximated
        by the formula h = mass/area/rho_fluid
        Input:
        h          : Global Equivalent Layer (m)
        rho_fluid  : density of fluid envelop
        Output:
        M          : mass of fluid envelop (kg)
        """
        M = h*self.area* rho_fluid
        return(M)

