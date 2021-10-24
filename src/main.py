"""

Author: Quentin Michalski

"""
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
# import xlsxwriter as xl

def volumeConservation(pressure,volumes,gases,volumeChambre,fullOutput=False):
        compressedVolumes = np.zeros(len(volumes))
        # equilibrateSet = [0,0,0]
        equilibrateSet = [0,1,1]
        # print(volumes)
        for jj,gas in enumerate(gases):
            T,P,Y = gas.TPY
            uncompressedDensity = gas.density
            gas.SP = gas.entropy_mass, pressure # isentropic compression of every slice
            if equilibrateSet[jj]:
                gas.equilibrate('SP')
            compressedDensity = gas.density
            compressionRatio = uncompressedDensity / compressedDensity
            compressedVolumes[jj] = volumes[jj] * compressionRatio
            gas.TPY = T,P,Y
        if fullOutput:
            for jj,gas in enumerate(gases):
                gas.SP = gas.entropy_mass, pressure
                if equilibrateSet[jj]:
                    gas.equilibrate('SP')
            return(compressedVolumes)
        else:
            return (volumeChambre-sum(compressedVolumes))

# =========================================
# Initialization
# =========================================

rayonChamber = 0.078
volumeChambre = 4/3*np.pi*rayonChamber**3  # volume de la chambre

T0 = 300
p0 = 1e5
ER = 1

mech = 'gri30.yaml'

volumes = np.array([volumeChambre, 0, 0]) # fresh gas / burning slice / burnt gases
gases = [ct.Solution(mech) for m in volumes]  # gases[0]: gaz frais / gases[1]: front de flamme / gases[2]: gaz brulés
gases[0].TP = T0, p0
gases[0].set_equivalence_ratio(phi=ER, fuel={'CH4': 1}, oxidizer={'O2': 1.0, 'N2': 0.79/0.21})
# =============================================================================
# Get adiabatic values
# =============================================================================
gases[0].equilibrate('UV')
T_adb = gases[0].T
p_adb = gases[0].P

gases[0].TP = T0, p0
gases[0].set_equivalence_ratio(phi=ER, fuel={'CH4': 1}, oxidizer={'O2': 1.0, 'N2': 0.79/0.21})
masses = np.array([volumes[0]*gases[0].density, 0, 0])
m0 = sum(masses)

dt = 1e-5  # (s)
t = 0
vf = 1 # m/s
rf = 0 # tracking flame radius

p = []
U = []
x = []
time = []
VSet = []
m_fSet = []
m_bSet = []
V_fSet = []
V_bSet = []
T_fSet = []
T_bSet = []

p.append(p0)
U.append(gases[0].int_energy_mass)
x.append(rf)
time.append(t)
VSet.append(volumes)
m_fSet.append(masses[0])
m_bSet.append(masses[2])
V_fSet.append(volumes[0])
V_bSet.append(volumes[2])
T_fSet.append(T0)
T_bSet.append(T0)

ii = 0

while rf < rayonChamber:
    print('Step {} at p = {} bar '.format(ii,p[ii]/1e5))
    
# =============================================================================
#     Flame burn step
# =============================================================================
    
    # calculating volume then mass of fresh gases that burns
    # note that for a random flame geometry
    # you need to specify an evolution of the flame surface that sweeps
    # the entire chamber so that along the propagation vector ds
    # the integral of the surface is equal to the volume
    # here i just make it spherical
    if volumes[2] > 0:
        # adjusting flame position based on new burnt gases sphere (accounts for compression)
        rf = ( 3 * volumes[2] / ( 4 * np.pi ) )**(1/3)
    # print(rf)
    ds = vf * dt
    rf = rf + ds
    rf = min(rf,rayonChamber) #tidying the ends
    print(rf)
    Sf = 4 * np.pi * (rf**2) # flame surface assuming a spherical combustion
    dv = Sf * ds
    rho_f = gases[0].density
    dm = dv * rho_f
    if dm > masses[0]: #tidying the ends
        dm = masses[0]
        dv = volumes[0]
    masses[0] = masses[0] - dm
    masses[1] = dm
    volumes[0] = volumes[0] - dv
    volumes[1] = dv
    # getting the burning slice properties
    gases[1].TPY = gases[0].TPY
    gases[1].equilibrate('HP')
    rho_rb = gases[1].density
    # isobaric dilatation of the burnt slice by ratio of burnt / fresh density 
    volumes[1] = volumes[1] * rho_f / rho_rb
    # note that now sum(volumes) > volumeChambre
    
# =============================================================================
#     Calculation of the compression to make it so that sum(volumes) = volumeChambre
# ============================================================================
    
    P = fsolve(volumeConservation, p[ii], args=(volumes,gases,volumeChambre)) # could replace
    P = P[0]
    
    compressedVolumes = volumeConservation(P,volumes,gases,volumeChambre,fullOutput=True)
    volumes = compressedVolumes
    # for jj,volume in enumerate(volumes):
    #     volumes[jj]=compressedVolumes[ii]
    
# =============================================================================
#     Transfering burnt slice into mix of burnt gases
# =============================================================================

    hb = (gases[1].enthalpy_mass*masses[1]+gases[2].enthalpy_mass*masses[2])/(masses[2] + masses[1])
    Yb = (gases[1].Y*masses[1]+gases[2].Y*masses[2])/(masses[2] + masses[1])
    gases[2].HPY = hb,P,Yb
    masses[2] = masses[2] + masses[1]
    masses[1] = 0
    volumes[2] = volumes[2] + volumes[1]
    volumes[1] = 0

    # total internal energy
    Ut = (gases[0].int_energy_mass*masses[0]+gases[1].int_energy_mass*masses[1]+gases[2].int_energy_mass*masses[2])/(sum(masses))
    
    # exporting information
    p.append(P)
    U.append(Ut)
    x.append(rf)
    time.append(t)
    m_fSet.append(masses[0])
    m_bSet.append(masses[2])
    V_fSet.append(volumes[0])
    V_bSet.append(volumes[2])
    T_fSet.append(gases[0].T)
    T_bSet.append(gases[2].T)
    
    ii+=1 # increment counter
    t+=dt # increment time

print('Mass balance {}'.format(sum(masses)/m0))
print('Internal energy balance {}'.format(U[-1]/U[0]))
print('Adiabatic pressure {}'.format(p_adb/p[-1]))
print('Adiabatic temperature {}'.format(T_adb/T_bSet[-1]))

x = np.array(x)
U = np.array(U)
p = np.array(p)
m_fSet = np.array(m_fSet)
m_bSet = np.array(m_bSet)
VSet = np.array(VSet)

plt.plot(x/rayonChamber,p)
plt.plot([x[0]/rayonChamber,x[-1]/rayonChamber],[p_adb,p_adb])
plt.show()

plt.plot(x/rayonChamber,m_fSet)
plt.plot(x/rayonChamber,m_bSet)
plt.show()

plt.plot(x/rayonChamber,V_fSet)
plt.plot(x/rayonChamber,V_bSet)
plt.show()