from mars_module import*

###############################################
# Large parameter space
# Uncomment the block below to explore a large
# parameter space
###############################################
#n_oceans=np.array([0.5,1.0,1.5,2.0,2.5,3.0])
#redox_fac=np.array([0.01,0.1,1.0,10.0,100.0])
#H_over_C=np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
#for ii in range(0,5):
#    for jj in range(0,4):
#        for kk in range (0,7):
            
#            mars=Mars(noceans=n_oceans[ii],HoverC=H_over_C[kk],\
#                       nsteps=1000000,redox_factor=redox_fac[jj])
#            mars.time_marching(4.0)
#            del mars
###############################################

n_oceans=0.8
#n_oceans=np.array([0.01,0.04,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0])

H_over_C = 0.55
#redox_fac= np.array([0.01,0.1,1.0,10.0,100.0,1000.0,1.0e4])
redox_fac=0.01      
for ii in range (0,1):
    #mars=Mars(noceans=n_oceans[ii],HoverC=H_over_C, nsteps=5000000,redox_factor=redox_fac)
    mars=Mars(noceans=n_oceans,HoverC=H_over_C, nsteps=5000000,redox_factor=redox_fac)
    mars.time_marching(4.0,const_Ftl=False)
    #mars.write_all_REE()
    
    del mars
