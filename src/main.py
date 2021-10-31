"""

Author: Quentin Michalski

"""
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

def _volumeConservation(pressure,volumes,gases,volumeChambre,fullOutput=False):
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

def _radiusToFlameSurface(rf):
    # definition of a flame surface as a function of a radius from origin
    # this could be a tabulated or analytical function
    Sf = 4 * np.pi * (rf**2)
    return(Sf)

def _volumesToContactSurface(volumes,geometry):
    # definition of a flame surface as a function of a radius from origin
    # this could be a tabulated or analytical function
    S = np.zeros(len(volumes))
    if geometry['name'] == 'spherical':
        S_burnt = 0
        S_slice = 0
        S_unburnt =  geometry['totalSurface'] 
    S[0] = S_unburnt
    S[1] = S_slice
    S[2] = S_burnt
    return(S)

def _burntVolumeToRadius(Vb,initialKernelRadius=1e-3):
    # the reciprocal which is required after compression to compute the new
    # radius based on the recompressed volume of burnt gases
    # the initialKernelRadius avoid the ignition problem
    if Vb > 0:
        rf = ( 3 * Vb / ( 4 * np.pi ) )**(1/3)
    else: 
        rf = initialKernelRadius
    return(rf)

def _flameSpeed():
    # definition of the flameSpeed
    # however its done and whatever it requires
    vf = 1 # m/s
    return(vf)

def constantVolumeCombustionStep(volumes,masses,gases):
    # =============================================================================
    #     Flame burn step
    # =============================================================================
    
    rf = _burntVolumeToRadius(volumes[2])
        
    vf = _flameSpeed()
    
    ds = vf * dt # differential element of flame movement
    
    Sf = _radiusToFlameSurface(rf) # flame surface assuming a spherical combustion
    
    rf = rf + ds # flame moves of its differential element of movement
    
    dv = Sf * ds # differential element of volume that burns
    
    rho_f = gases[0].density # fetch density of fresh gases
    
    dm = dv * rho_f # differential element of mass that burns
    
    if dm > masses[0]: #tidying the ends so mass is conserved
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
    
    P = fsolve(_volumeConservation, p[ii], args=(volumes,gases,volumeChambre)) # could replace
    P = P[0]
    
    compressedVolumes = _volumeConservation(P,volumes,gases,volumeChambre,fullOutput=True)
    volumes = compressedVolumes
    
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
    
    return(volumes,masses)

def heatExchangeStep(volumes,masses,gases,geometry,correlation=0,fullOutput=False):
    # =============================================================================
    #   Heat exchange step
    # =============================================================================
    Twall = 300
    dqSet=np.zeros(len(volumes))
    if correlation == 0:
        # backing up gas for calculation
        Hf,Pf,Yf = gases[0].HPY
        Ymix = (gases[0].Y*masses[0]+gases[1].Y*masses[1]+gases[2].Y*masses[2])/(masses[0] + masses[1] + masses[2])
        Hmix = (gases[0].enthalpy_mass*masses[0]+gases[1].enthalpy_mass*masses[1]+gases[2].enthalpy_mass*masses[2])/(masses[0]+masses[1]+masses[2])
        gases[0].HPY = Hmix,Pf,Ymix 
        Tmix = gases[0].T
        gases[0].HPY = Hf,Pf,Yf
        h = 1.15 * ( ((Pf/1e6)**2) * Tmix )**(1/3)
        Sexchange = _volumesToContactSurface(volumes,geometry)
        
    for jj,gas in enumerate(gases):
        dQ = Sexchange[jj] * h * (gases[jj].T - Twall)
        dq = dQ / masses[jj]
        dq = np.nanmax([dq,0])
        intialDensity = gas.density
        gases[jj].HP = gases[jj].enthalpy_mass - dq, Pf
        finalDensity = gas.density
        expansionRatio = intialDensity / finalDensity
        volumes[jj] = volumes[jj] * expansionRatio
        dqSet[jj]=dq
    
    P = fsolve(_volumeConservation, p[ii], args=(volumes,gases,volumeChambre)) # could replace
    P = P[0]
    
    compressedVolumes = _volumeConservation(P,volumes,gases,volumeChambre,fullOutput=True)
    volumes = compressedVolumes
    if fullOutput:
        return(volumes,dqSet)
    else:    
        return(volumes)

# =========================================
# Initialization
# =========================================

rayonChamber = 0.078
volumeChambre = 4/3*np.pi*rayonChamber**3  # volume de la chambre

geometry = {}
geometry['name'] = 'spherical'
geometry['totalSurface'] = 4*np.pi*rayonChamber**2

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

dt = 1e-4  # (s)
t = 0
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
dqSet = []
    
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
dqSet.append([0,0,0])
ii = 0

while volumes[0] > 0:
    print('Step {} at p = {} bar '.format(ii,p[ii]/1e5))
    
    volumes,dq = heatExchangeStep(volumes,masses,gases,geometry,correlation=0,fullOutput=True)                            
    volumes,masses = constantVolumeCombustionStep(volumes,masses,gases)
    
    # =============================================================================
    # Exports
    # =============================================================================
    rf = _burntVolumeToRadius(volumes[2])
    
    # total internal energy
    Ut = (gases[0].int_energy_mass*masses[0]+gases[1].int_energy_mass*masses[1]+gases[2].int_energy_mass*masses[2])/(sum(masses))
    
    # exporting information
    p.append(gases[0].P)
    U.append(Ut)
    x.append(rf)
    time.append(t)
    m_fSet.append(masses[0])
    m_bSet.append(masses[2])
    V_fSet.append(volumes[0])
    V_bSet.append(volumes[2])
    T_fSet.append(gases[0].T)
    T_bSet.append(gases[2].T)
    dqSet.append(dq)
    
    ii+=1 # increment counter
    t+=dt # increment time

print('Mass balance {} %'.format((m0-sum(masses))*100/m0))
print('Internal energy balance {} %'.format((U[0]-U[-1])/U[0]*100))
print('Adiabatic pressure {} %'.format((p_adb-p[-1])*100/p_adb))
print('Adiabatic temperature {} %'.format((T_adb-T_bSet[-1])*100/T_adb))

x = np.array(x)
U = np.array(U)
p = np.array(p)
m_fSet = np.array(m_fSet)
m_bSet = np.array(m_bSet)
VSet = np.array(VSet)
dqSet = np.array(dqSet)

plt.plot(x/rayonChamber,p)
plt.plot([x[0]/rayonChamber,x[-1]/rayonChamber],[p_adb,p_adb])
plt.show()

plt.plot(x/rayonChamber,m_fSet)
plt.plot(x/rayonChamber,m_bSet)
plt.show()

plt.plot(x/rayonChamber,V_fSet)
plt.plot(x/rayonChamber,V_bSet)
plt.show()

plt.plot(x/rayonChamber,U/U[0])
plt.show()

plt.plot(x/rayonChamber,dqSet/geometry['totalSurface']/1e6)
plt.show()