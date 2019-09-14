from mars_module import*
from plot_module import*
from scipy.optimize import curve_fit
def fit_func(x,m):
    """Fitting function"""
    return m*x


n_oceans=np.array([0.01,0.04,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0])
redox_fac= np.array([0.01,0.1,1.0,10.0,100.0,1000.0,1.0e4])
H_over_C = 0.55
CO2overH2O=1.0/(H_over_C*2.45)
CO2_mass=n_oceans*1.6e21*CO2overH2O
H2O_mass=n_oceans*1.6e21
########################################################################
## Plot show solidus and adiabat
########################################################################
#solidus_plot()
########################################################################
# Plot: Show thermal evolution related variables
# Uncomment the line below for plot, may take some time
########################################################################
#mars=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=redox_fac[2])
#print 'water mass, GEL (km)', mars.H2OMASS, mars.kg2GELm(mars.H2OMASS)/1.0e3
#mars.evolution_plot()
########################################################################
# Plot: compare thermal evolutions
# for two extreme fO2 scenarios
# Uncomment the line below for plot, may take some time
########################################################################
#mars1=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=redox_fac[0])
#mars2=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=redox_fac[4])
#compare_evolution(mars1,mars2)
########################################################################
# Plot 1: Show Final concentrations as a function of initial mass
# Plot 2: Mantle C and freezing time as a function of K*
########################################################################

#Final_concentrations(n_oceans)
########################################################################
# Plot: Show freezing time, PCO2, and mantle C as a function of K*
########################################################################
#redox_plots(redox_fac,n1=0.8,n2=0.4)
########################################################################
### Plot: Compare the thermal evoloution of two cases
### one with dynamic melt trapping, the other with constant Ftl
### Can take some time
########################################################################
#compare_Ftl_models()
########################################################################
########################################################################
### Plot: Create a plot for incompatible element concentrations
########################################################################
#mars1=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=redox_fac[2])
#REE_plots(mars1)
########################################################################
## Plot: Create a profile of H2O and CO2 in the RM after
## crystallization ends
## May take some time to run
########################################################################
#mars2=Mars_read(noceans=1.0,HoverC=0.55,redox_factor=redox_fac[4])
#print 'water mass, GEL (km)', mars2.H2OMASS, mars2.kg2GELm(mars2.H2OMASS)/1.0e3
#plot_profile(mars2)
############################################################################

########################################################
# Plot several time steps in REE plots
################################################
mars=Mars_read(noceans=0.8,HoverC=0.55,redox_factor=0.01, Fcl=True)
print 'Water in GEL m in the bulk',mars.kg2GELm(mars.H2OMASS)
REE_plots(mars)
plt.show()
