from mars_module import*

mars1=Mars()

r=np.linspace(mars1.core,mars1.radius) #radius in m
PGPa=mars1.depth2PGPa(r)
solidus=mars1.solidus(r)
liquidus=mars1.liquidus(r)
front=mars1.freezing_front(r)

plt.figure(1)

plt.subplot(1,2,1)
plt.plot(PGPa,r/1.0e3)
plt.xlabel('Pressure (GPa)')
plt.ylabel('Radius (km)')
#
plt.subplot(1,2,2)
plt.plot(solidus,PGPa,'-b',linewidth=4)
plt.plot(liquidus,PGPa,'-r',linewidth=4)
plt.plot(front,PGPa,'--g',linewidth=4)
plt.ylim=(np.max(PGPa),np.min(PGPa))
plt.gca().invert_yaxis()
plt.legend(['Solidus','Liquidus','Front'],loc=3)


plt.xlabel(r"Temperature $^\mathsf{o}$C")
plt.ylabel('Pressure (GPa)')


plt.figure(3)
mars1.radius_temperature_analytical(plot=True)
plt.plot(mars1.T,0.0*mars1.T)
plt.show()
