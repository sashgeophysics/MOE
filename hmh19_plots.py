from mars_module import*
def plot_profile(obj):
    """This function plots a profile of RM concentrations
    the input object should be of class Mars_read"""
    f2=obj.output_filenames(const=False)
    data2=np.loadtxt(f2[3],delimiter=',')
    obj.CCO2RM   = data2[:,4]
    del data2
    data2=np.loadtxt(f2[0],delimiter=',')
    obj.akm      = data2[:,2]
    del data2
    data2=np.loadtxt(f2[2],delimiter=',')
    obj.CH2ORM   = data2[:,4]
    del data2,f2



    plt.figure(figsize=(10,12))
    plt.plot(obj.CCO2RM*1.0e6,obj.akm,'dimgray',lw=4)
    plt.plot(obj.CH2ORM*1.0e6,obj.akm,'teal',lw=4,ls='dashed')
    plt.legend(['Total C',r'$\mathrm{H_2O}$'],fancybox=True\
               ,framealpha=0.7,loc=4)
    plt.xlabel(r'Concentration in the RM (ppm)')
    plt.ylabel('Height above CMB (km)')
def compare_Ftl_models():
    """Compares the thermal evolutions of two objects
    The second object should have a constant trapped
    melt fraction
    """
    obj1 = Mars_read(noceans=1.0,HoverC=0.55,redox_factor=1.0)
    obj2 = Mars_read(noceans=1.0,HoverC=0.55,redox_factor=1.0,Ftl_const=True)
    f=obj1.output_filenames()
    data3=np.loadtxt(f[2],delimiter=',')
    obj1.tma      = data3[:,0]
    obj1.MH2ORM   = data3[:,1]
    obj1.MH2OMO   = data3[:,2]
    #self.MH2OPA   = data3[:,3]
    obj1.CH2ORM   = data3[:,4]
    #self.CH2OMO   = data3[:,5]
    data2=np.loadtxt(f[1],delimiter=',')
    obj1.PCO2     = data2[:,1]
    obj1.PH2O     = data2[:,2]
    del data2,data3,f
    f=obj2.output_filenames(const=True)
    data3=np.loadtxt(f[2],delimiter=',')
    obj2.tma      = data3[:,0]
    obj2.MH2ORM   = data3[:,1]
    obj2.MH2OMO   = data3[:,2]
    #self.MH2OPA   = data3[:,3]
    obj2.CH2ORM   = data3[:,4]
    #self.CH2OMO   = data3[:,5]
    data2=np.loadtxt(f[1],delimiter=',')
    obj2.PCO2     = data2[:,1]
    obj2.PH2O     = data2[:,2]
    del data2,data3,f
    
    plt.figure(figsize=(18,12))
    plt.subplot(2,2,1)
    plt.plot(obj1.tma,obj1.PCO2/1.0e5,color='forestgreen',linewidth=4)
    plt.plot(obj2.tma,obj2.PCO2/1.0e5,ls='--',color='darkorange',linewidth=4)
    plt.ylim(275,325)
    plt.legend(['Dynamic $F_{tl}$',r'Constant $F_{tl}$'],loc=2,\
               fancybox=2,framealpha=0.7)
    plt.ylabel(r"${P_{CO2}}$ (bar)",fontsize=30)

    plt.subplot(2,2,2)
    plt.plot(obj1.tma,obj1.MH2ORM,color='forestgreen',linewidth=4)
    plt.plot(obj2.tma,obj2.MH2ORM,ls='--',color='darkorange',linewidth=4)
    plt.ylabel(r"$\mathrm{H_2O}$ in RM (kg)",fontsize=30)
    
    
    plt.subplot(2,2,3)
    plt.plot(obj1.tma,obj1.PH2O/1.0e5,color='forestgreen',linewidth=4)
    plt.plot(obj2.tma,obj2.PH2O/1.0e5,ls='--',color='darkorange',linewidth=4)
    plt.ylabel(r"${P_{H2O}}$ (bar)",fontsize=30)
    plt.xlabel('Time (Ma)',fontsize=30)
    
    plt.subplot(2,2,4)
    plt.plot(obj1.tma,obj1.MH2OMO,color='forestgreen',linewidth=4)
    plt.plot(obj2.tma,obj2.MH2OMO,ls='--',color='darkorange',linewidth=4)
    plt.ylabel(r"$\mathrm{H_2O}$ in MO (kg)",fontsize=30)
    
    
    plt.xlabel('Time (Ma)',fontsize=30)
def compare_evolution(obj1,obj2):
    """Compares the thermal evolutions of two objects"""
    obj1.load_evolution()
    obj2.load_evolution()
    plt.figure(figsize=(10,16))

    plt.subplot(3,1,1)
    plt.plot(obj1.tma,obj1.T,color='cornflowerblue',linewidth=4)
    plt.plot(obj2.tma,obj2.T,ls='--',color='dimgray',linewidth=4)
    plt.legend([r'$K^\ast = ${:.2f}'.format(obj1.redox_factor), r'$K^\ast = ${:.2f}'.format(obj2.redox_factor)],fancybox=True,framealpha=0.7)
    plt.ylabel(r'$T^\mathrm{o}$C',fontsize=30)
    
    plt.subplot(3,1,2)
    plt.plot(obj1.tma,obj1.PCO2/1.0e5,color='cornflowerblue',linewidth=4)
    plt.plot(obj2.tma,obj2.PCO2/1.0e5,ls='--',color='dimgray',linewidth=4)
    plt.ylabel(r"${P_{CO2}}$ (bar)",fontsize=30)
    
    plt.subplot(3,1,3)
    plt.plot(obj1.tma,obj1.PH2O/1.0e5,color='cornflowerblue',linewidth=4)
    plt.plot(obj2.tma,obj2.PH2O/1.0e5,ls='--',color='dimgray',linewidth=4)
    plt.ylabel(r"${P_{H2O}}$ (bar)",fontsize=30)
    plt.xlabel('Time (Ma)',fontsize=30)

# First, create instances of the object for visualization
k1=np.array([0.01,0.1,1.0,10.0,1.0e2])
print k1
# Create an array for labels for REE concentration
REE=np.zeros(23)
for ii in range (0,23):
    REE[ii]=ii
    
REE_label=['Rb', 'Ba', 'Th', 'U','Ta', 'K', 'La', 'Ce' ,'P' , 'Sr' , 'Nd' , 'Sm' , 'Zr' ,  'Hf' , 'Eu' ,'Gd' ,  'Tb' ,  'Dy' ,  'Y' ,  'Er' ,   'Tm' ,  'Yb' ,   'Lu']

mars1=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=k1[0])
mars2=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=k1[1])
mars3=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=k1[2])
mars4=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=k1[3])
mars5=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=k1[4])

#########################################################
# Plot: Show solidus liquidus and adiabat
# 
###########################################################
mars_sol=Mars()
print 'MArs core radius (km)',mars_sol.core/1.0e3
print 'Mars planetary radius (km)',mars_sol.radius/1.0e3
print 'Mars mantle thickness',(mars_sol.radius-mars_sol.core)/1.0e3
r=np.linspace(mars_sol.core,mars_sol.radius) #radius in m
r_mantle_km=(r-mars_sol.core)*1.0e-3
PGPa=mars_sol.depth2PGPa(r)
solidus=mars_sol.solidus(r)
liquidus=mars_sol.liquidus(r)
front=mars_sol.freezing_front(r)
rad1=mars_sol.solid_radius(2100.0)

T1=mars_sol.adiabat(2100.0,r)
rad2=mars_sol.solid_radius(1600.0)


T2=mars_sol.adiabat(1600.0,r)
temp_label=[1400,1800,2200,2600]

plt.figure(figsize=(19,12))
#
plt.subplot(1,2,1)
ax = plt.gca()

plt.plot(liquidus,r_mantle_km,'-',color='firebrick',linewidth=4)
plt.plot(front,r_mantle_km,'-',color='teal',linewidth=4)
plt.plot(solidus,r_mantle_km,'-',color='darkslategray',linewidth=4)
plt.xticks(temp_label)
ax.fill_betweenx(r_mantle_km,front, solidus, where=front>= solidus, alpha=0.7,  facecolor='teal',)

plt.legend(['Liquidus','Front','Solidus'],loc=1,fancybox=True,framealpha=0.7)
plt.plot(T1,r_mantle_km,'-k',linewidth=4)
plt.plot(T2,r_mantle_km,'-k',linewidth=4)
plt.plot(2340,146.0,'o',color='tomato',markersize=20,alpha=0.7)
plt.plot(1621,1776.0,'o',color='tomato',markersize=20,alpha=0.7)
#plt.ylim(1395.0,3390.0)
plt.ylim(0.0,2000.0)

plt.text(1550,1350,r'1400$^\mathrm{o}$ adiabat',rotation=-85,fontsize=30)
plt.text(2200,1350,r'2100$^\mathrm{o}$ adiabat',rotation=-85,fontsize=30)
plt.xlabel(r'Temperature ($^\mathrm{o}$C)',fontsize=30)
plt.ylabel('Height above CMB (km)',fontsize=30)
#
plt.subplot(1,2,2)
mars_sol.radius_temperature_analytical(plot=True)
plt.xticks(temp_label)
plt.ylim(0.0,2000.0)
#
#Plot dadT if neede
#plt.subplot(1,3,1)
#T=np.linspace(1300.0,2150.0,100)
#dadT=0.0*T
#for ii in range (0,99):
#            dadT[ii]=mars_sol.dadT_analytical(T[ii])
#plt.plot(T,dadT/1.0e3,'or')

del mars_sol
#########################################################
# Plot: Show thermal evolution related variables
# Uncomment the line below for plot, may take some time
###########################################################

#mars3.evolution_plot()

##############################################
# Plot: compare thermal evolutions
# for two extreme fO2 scenarios
# Uncomment the line below for plot, may take some time
##############################################
#compare_evolution(mars1,mars5)

##############################################
# Plot: compare end products for different 
# values of K*
##############################################
plt.figure(figsize=(10,16))
freezing_time=np.array([mars1.tma_final,mars2.tma_final,mars3.tma_final,mars4.tma_final,mars5.tma_final])
print freezing_time
CO2_pressure=np.array([mars1.PCO2_final,mars2.PCO2_final,mars3.PCO2_final,mars4.PCO2_final,mars5.PCO2_final])
RM_CO2 = np.array([mars1.CCO2RM_final,mars2.CCO2RM_final,mars3.CCO2RM_final,mars4.CCO2RM_final,mars5.CCO2RM_final])

plt.subplot(3,1,1)
plt.semilogx(k1,freezing_time,'s',markersize=20,color = 'indianred')
plt.ylabel('Freezing Time (Ma)')
plt.xlim(0.001,1000.0)
plt.ylim(0.1,1.5)
plt.subplot(3,1,2)
plt.semilogx(k1,CO2_pressure/1.0e5,'s',markersize=20,color = 'indianred')
plt.yticks([150,200,250,300,350])
plt.ylabel(r"${P_{CO2}}$ (bar)",fontsize=30)
plt.xlim(0.001,1000.0)
plt.subplot(3,1,3)
plt.semilogx(k1,RM_CO2*1.0e6,'s',markersize=20,color = 'indianred')
plt.yticks([100,300,500,700,900])
plt.xlim(0.001,1000.0)
plt.ylim(0.0,1000.0)
plt.ylabel('Mantle C (ppm)')
plt.xlabel(r'$K^\ast$')
##############################################
# Plot: compare REE concentrations
# 
##############################################

plt.figure(figsize=(16,12))

plt.semilogy(REE,(mars1.CREEMO_final/mars1.CI_REE),'s',color='darkgreen',markersize=20)
plt.semilogy(REE,(mars1.NWA1068/mars1.CI_REE),'o',color='moccasin',markersize=20)
plt.semilogy(REE,(mars1.shergotty/mars1.CI_REE),'o',color='tomato',markersize=20)
plt.semilogy(REE,(mars1.zagami/mars1.CI_REE),'o',color='#D4B1B1',markersize=20)
plt.semilogy(REE,(mars1.CREERM_final/mars1.CI_REE),'o',color='darkgreen',markersize=20)

plt.legend([ r'MO after 99.5$\%$ crystallization', 'NWA1068','Shergotty','Zagami','RM after 99.5$\%$ crystallization'],loc=4,fancybox=True,framealpha=0.7)
plt.ylabel('Concentration/Chondrite',fontsize=30)
plt.semilogy(REE,(mars1.CREEMO_final/mars1.CI_REE),color='darkgreen')
plt.semilogy(REE,(mars1.CREERM_final/mars1.CI_REE),color='darkgreen')
plt.semilogy(REE,(mars1.NWA1068/mars1.CI_REE),'-',color='#005668')
plt.semilogy(REE,(mars1.shergotty/mars1.CI_REE),'-',color='darkgreen')
plt.semilogy(REE,(mars1.zagami/mars1.CI_REE),'-',color='darkgreen')

ax=plt.gca()
plt.xticks(REE,REE_label,rotation='vertical',fontsize=25)


############################
### delete the used objects
##############################
del mars1,mars2,mars3,mars4,mars5
#################################
## Read in new objects
#################################
n_oceans=np.array([0.5,1.0,1.5,2.0,2.5,3.0])
H_over_C = 0.55
CO2overH2O=1.0/(H_over_C*2.45)
CO2_mass=n_oceans*1.6e21*CO2overH2O
H2O_mass=n_oceans*1.6e21
redox_fac=1.0        
mars1=Mars_read(noceans=n_oceans[0],HoverC=H_over_C,redox_factor=redox_fac)
mars2=Mars_read(noceans=n_oceans[1],HoverC=H_over_C,redox_factor=redox_fac)
mars3=Mars_read(noceans=n_oceans[2],HoverC=H_over_C,redox_factor=redox_fac)
mars4=Mars_read(noceans=n_oceans[3],HoverC=H_over_C,redox_factor=redox_fac)
mars5=Mars_read(noceans=n_oceans[4],HoverC=H_over_C,redox_factor=redox_fac)
mars6=Mars_read(noceans=n_oceans[5],HoverC=H_over_C,redox_factor=redox_fac)

mars1_const=Mars_read(noceans=n_oceans[0],HoverC=H_over_C,\
                      redox_factor=redox_fac,Ftl_const=True)
mars2_const=Mars_read(noceans=n_oceans[1],HoverC=H_over_C,\
                      redox_factor=redox_fac,Ftl_const=True)
mars3_const=Mars_read(noceans=n_oceans[2],HoverC=H_over_C,\
                      redox_factor=redox_fac,Ftl_const=True)
mars4_const=Mars_read(noceans=n_oceans[3],HoverC=H_over_C,\
                      redox_factor=redox_fac,Ftl_const=True)
mars5_const=Mars_read(noceans=n_oceans[4],HoverC=H_over_C,\
                      redox_factor=redox_fac,Ftl_const=True)
mars6_const=Mars_read(noceans=n_oceans[5],HoverC=H_over_C,\
                      redox_factor=redox_fac,Ftl_const=True)

redox_fac2=100.0
mars1_a=Mars_read(noceans=n_oceans[0],HoverC=H_over_C,redox_factor=redox_fac2)
mars2_a=Mars_read(noceans=n_oceans[1],HoverC=H_over_C,redox_factor=redox_fac2)
mars3_a=Mars_read(noceans=n_oceans[2],HoverC=H_over_C,redox_factor=redox_fac2)
mars4_a=Mars_read(noceans=n_oceans[3],HoverC=H_over_C,redox_factor=redox_fac2)
mars5_a=Mars_read(noceans=n_oceans[4],HoverC=H_over_C,redox_factor=redox_fac2)
mars6_a=Mars_read(noceans=n_oceans[5],HoverC=H_over_C,redox_factor=redox_fac2)
##############################################
# Plot: compare H2O and CO2 concentrations
# as a function of initial CO2 content
# For two different values of K*
##############################################
final_water=np.array([mars1.CH2ORM_final,mars2.CH2ORM_final,\
                      mars3.CH2ORM_final,mars4.CH2ORM_final,\
                      mars5.CH2ORM_final,mars6.CH2ORM_final])
final_CO2=np.array([mars1.CCO2RM_final,mars2.CCO2RM_final,\
                      mars3.CCO2RM_final,mars4.CCO2RM_final,\
                      mars5.CCO2RM_final,mars6.CCO2RM_final])
final_time=np.array([mars1.tma_final,mars2.tma_final,\
                      mars3.tma_final,mars4.tma_final,\
                      mars5.tma_final,mars6.tma_final])
final_water_a=np.array([mars1_a.CH2ORM_final,mars2_a.CH2ORM_final,\
                      mars3_a.CH2ORM_final,mars4_a.CH2ORM_final,\
                      mars5_a.CH2ORM_final,mars6_a.CH2ORM_final])
final_CO2_a=np.array([mars1_a.CCO2RM_final,mars2_a.CCO2RM_final,\
                      mars3_a.CCO2RM_final,mars4_a.CCO2RM_final,\
                      mars5_a.CCO2RM_final,mars6_a.CCO2RM_final])
final_time_a=np.array([mars1_a.tma_final,mars2_a.tma_final,\
                      mars3_a.tma_final,mars4_a.tma_final,\
                      mars5_a.tma_final,mars6_a.tma_final])
final_water_const=np.array([mars1_const.CH2ORM_final,mars2_const.CH2ORM_final,\
                      mars3_const.CH2ORM_final,mars4_const.CH2ORM_final,\
                      mars5_const.CH2ORM_final,mars6_const.CH2ORM_final])
final_CO2_const=np.array([mars1_const.CCO2RM_final,mars2_const.CCO2RM_final,\
                      mars3_const.CCO2RM_final,mars4_const.CCO2RM_final,\
                      mars5_const.CCO2RM_final,mars6_const.CCO2RM_final])


plt.figure(figsize=(10,16))

plt.subplot(3,1,1)
plt.plot(CO2_mass,final_water*1.0e6,'o',markersize=15,color='powderblue')
plt.plot(CO2_mass,final_water_a*1.0e6,'s',markersize=15,color='dimgray',alpha=0.5)
plt.ylabel(r'Mantle $\mathrm{H_2O}$ (ppm)')
plt.yticks([1000,2000,3000,4000,5000])
#plt.xlim(0.0,3.5)
plt.subplot(3,1,2)
plt.semilogy(CO2_mass,final_CO2*1.0e6,'o',markersize=15,color='powderblue')
plt.semilogy(CO2_mass,final_CO2_a*1.0e6,'s',markersize=15,color='dimgray')
plt.ylabel('Mantle C (ppm)')
#plt.xlim(0.0,3.5)
plt.subplot(3,1,3)
plt.plot(CO2_mass,final_time,'o',markersize=15,color='powderblue')
plt.plot(CO2_mass,final_time_a,'s',markersize=15,color='dimgray')
plt.legend([r'$K^\ast = ${:.2f}'.format(redox_fac), r'$K^\ast = ${:.2f}'.format(redox_fac2)],fancybox=True,loc=2,framealpha=0.7)
plt.yticks([0.5,1.5,2.5,3.5])
plt.ylabel('Time (Ma)')
plt.xlabel(r'Initial $\mathrm{CO_2}$ (kg)')
#plt.xlim(0.0,3.5)
###########################################################

##########################################################################
### Plot comparing the final CO2 and H2O contents for the two different
### types of trapped melt fraction. K* is constant
##########################################################################
plt.figure(figsize=(10,16))
plt.subplot(2,1,1)
plt.plot(H2O_mass,final_water*1.0e6,'o',markersize=20,color='darkgreen')
plt.plot(H2O_mass,final_water_const*1.0e6,'o',markersize=20,markeredgecolor='darkorange',markerfacecolor='none',markeredgewidth=3)
plt.legend(['Dynamic $F_{tl}$',r'Constant $F_{tl}$'],loc=2,fancybox=2,framealpha=0.7)
plt.ylabel(r'Mantle $\mathrm{H_2O}$ (ppm)')
plt.yticks([1000,2000,3000,4000,5000])
plt.xlabel(r'Initial $\mathrm{H_2O}$ (kg)')
#plt.xlim(0.0,3.5)
plt.subplot(2,1,2)
plt.semilogy(CO2_mass,final_CO2*1.0e6,'o',markersize=20,color='darkgreen')
plt.semilogy(CO2_mass,final_CO2_const*1.0e6,'o',markersize=20,markeredgecolor='darkorange', markerfacecolor='none',markeredgewidth=3)
plt.ylabel('Mantle C (ppm)')
#plt.xlim(0.0,3.5)
plt.xlabel(r'Initial $\mathrm{CO_2}$ (kg)')
##############################################
# Plot: compare REE concentrations
# for dyamic trapping and constant Ftl
##############################################
plt.figure(figsize=(16,10))


plt.semilogy(REE,(mars2.all_REE_init_conc/mars2.CI_REE),'s',color='lightsalmon',markersize=20)
plt.semilogy(REE,(mars2.CREERM_final/mars2.CI_REE),'o',color='darkgreen',markersize=20)

plt.semilogy(REE,(mars2_const.CREERM_final/mars2.CI_REE),'o',markeredgecolor='darkorange',markersize=20, markerfacecolor='none',markeredgewidth=3)


plt.legend([ r'Initial abundance','RM, dynamic trapping','RM, constant $F_{tl}$ '],loc=4,fancybox=True,framealpha=0.7)

plt.ylabel('Concentration/Chondrite',fontsize=30)

plt.semilogy(REE,(mars2.CREERM_final/mars2.CI_REE),color='darkgreen')
plt.semilogy(REE,(mars2.all_REE_init_conc/mars2.CI_REE),'-',color='lightsalmon')
plt.semilogy(REE,(mars2_const.CREERM_final/mars2_const.CI_REE),'--',color='darkorange')

ax=plt.gca()
plt.xticks(REE,REE_label,rotation='vertical',fontsize=25)
######################################################################

##########################################################################
### Compare the thermal evoloution of two cases
### one with dynamic melt trapping, the other with constant Ftl
### Can take some time
##########################################################################


compare_Ftl_models()
##########################################################################
#######################################################
## Create a profile of H2O and CO2 in the RM after
## crystallization ends
## May take some time to run
#####################################################
k1=np.array([0.01,0.1,1.0,10.0,1.0e2])

mars2=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=k1[4])
plot_profile(mars2)

############################################################################
plt.show()
