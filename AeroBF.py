# -*- coding: utf-8 -*-
"""
AERONAUTICS BASIC FUNCTIONS
@author: ELC
"""

#Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Parameters
g = 9.80665 #[m/s**2]
airConstantR = 287 #[J/(kg*K)]
nauticMileToM = 1852 #[m]
gammaAir = 1.4
radiusEarth = 6384*10**3 #[m]
massMolarAir = 28.97 # [g/mol]
rhoAirStandard = 1.225 #[kg/m**3]
temperatureStandard = 15 + 273.15 #[K]
pressureStandard = 101325 #[pa]
kBoltzmann = 1.38 * 10**(-23) #[J/K]
viscocityDynAir = 1.789*10**(-5) #[pa*s]

toDegConv = 180/np.pi

# Functions


def calcTrueAirspeed(eqAirspeed,rho):
    """
    Description:
    Return true airspeed calculation for a given equivalent airspeed and
    local density

    Parameters: Local Cruise Density(kg/m**3)
                Equivalent Airspeed (m/s)

    Output: True airspeed(m/s)
    """      
    
    result = np.sqrt(rhoAirStandard/rho)*eqAirspeed
    return result
    
def calcEqAirspeed(trueAirspeed,rho):
    """
    Description:
    Return equivalent airspeed calculation for a given true airspeed and
    local density

    Parameters: Local Cruise Density(kg/m**3)
                True Airspeed (m/s)

    Output: Equivalent airspeed(m/s)
    """   
    result = np.sqrt(rho/rhoAirStandard)*trueAirspeed
    return result


def calcSonicSpeed(temperature,R=airConstantR,gamma=gammaAir):
    """
    Description:
    Return speed of sound in air for a given temperature in K

    Parameters: Temperature(K)
                R constant (J/kgK) (air by default)
                gamma (air by default)

    Output: sound speed (m/s) (for air by default)
    """        
    result = np.sqrt(gamma*R*temperature)
    return result
 

def calcNoMach(temperature,airspeed,regimeOut=False):
    """
    Description:
    Return mach number for a given Temperature and airspeed

    Parameters: Temperature(K)
                airspeed (m/s)

    Output: noMach and regime (if True)
    """      
    result = airspeed/calcSonicSpeed(temperature)
    if regimeOut == False:
    
        if result <= 0.8:
            regime = "Subsonic"
        elif result > 0.8 and result <= 1.2:
            regime = "Transonic"
        elif result > 1.2 and result <= 5:
            regime = "Supersonic"
        else:
            regime = "Hipersonic"
        return (result,regime)
    
    else:
        return result
    
def calcGeometricAltitude(altitude):
    """
    Description:
    Return geometric altitude for a given gepotential altitude (normal altitude)

    Parameters: Altitude (m)
             

    Output: Geometric Altitude (m)
        
    """
    result = radiusEarth*altitude*1/(radiusEarth - altitude)
    return result


def calcGeopotentialAltitude(altitude):
    """
    Description:
    Return geopotential altitude (normal altitude) for a given geometric altitude

    Parameters:Geometric Altitude (m)
             

    Output: Geopotential Altitude (m)
        
    """
    result = radiusEarth*altitude*1/(radiusEarth + altitude)
    return result


     
# AEROSTATICS --------------------------------------------

def calcGasEquivalentAirTemp(massMolar):
    """
    Description:
    Return required increase of temperature (deltaT) in for achieving
    other gas lift for a given volume

    Parameters: Molar mass of the desired gas (g/mol)
             
    Output: Temperature increase deltaT (K)
        
    """   

    rhoGas = rhoAirStandard*(massMolar/massMolarAir) #kg/m**3 
    deltaT =temperatureStandard*(rhoAirStandard-rhoGas)/rhoGas 
    return deltaT

     
def calcAirTempForLift(desVolume,desPayload,tempEnv = 288.15, pressure = 101325 ):
    """
    Description:
    Return required increase of temperature (deltaT) for generating enough lift
    for lifting a desired payload with a given volume

    Parameters: desired Volume (m**3)
                desired payload (kg)
                Temperature (optional) [K]
                pressure (optional) [Pa]
                
    Output: Temperatue increase deltaT (K)
        
    """  
    rhoAir = pressure/(airConstantR*tempEnv)

    deltaT = airConstantR*(airConstantR/(rhoAir*desVolume-desPayload)) #Relación de deltaT = f(m,rho_aire,V,Tenv)
    return deltaT

# FUSELAGE STRESS ----------------------------------------------------------

def calcFusHoopStress(pressureCabin,pressureOuter,radius,wallThk):
    """
    Description:
    Return the hoop stress in the fuselage (given diameter and wall thk) for a given pressure at cabin
    and for a given external pressure

    Parameters: Pressure at cabin (Pa)
                Outer Pressure (Pa)
                Fuselage Radius (m)
                Wall Thk (m)
                
    Output: Fuselage Stress (Pa)
        
    """    
    deltaPressure = pressureCabin - pressureOuter
    sigmaHoop = deltaPressure*radius*1/wallThk
    return sigmaHoop


def clacFusLongStress(pressureCabin,pressureOuter,radius,wallThk):
    """
    Description:
    Return the longitudinal stress in the fuselage (given diameter and wall thk) for a given pressure at cabin
    and for a given external pressure

    Parameters: Pressure at cabin (Pa)
                Outer Pressure (Pa)
                Fuselage Radius (m)
                Wall Thk (m)
                
    Output: Fuselage Stress (Pa)
        
    """    
    deltaPressure = pressureCabin - pressureOuter
    sigmaLong = deltaPressure*radius*1/(2*wallThk)
    return sigmaLong


def calcFusDesignThk(pressureCabin,pressureOuter,radius,stressUlt,safetyFactor = 1):
    """
    Description:
    Return the design thk of the fuselage (given diameter) for a given pressure at cabin,
    for a given external pressure condition and for a specified safety factor and Ultimate Stress

    Parameters: Pressure at cabin (Pa)
                Outer Pressure (Pa)
                Fuselage Radius (m)
                Safety Factor (optional)
                Ultimate Stress (Pa)
                
    Output: Design Wall Thk (m)
        
    """     
    deltaPressure = pressureCabin - pressureOuter
    stressAllow = stressUlt/safetyFactor
    
    thkDesign = deltaPressure*radius/stressAllow #[m]
    return thkDesign

# Aerodynamicas Fundamentals ---------------

def calcAirspeedForPitot(pressureTotal,pressureDynamic,rhoAir):
    """
    Description:
    Returns the calculated airspeed for the given pressure measurements

    Parameters: - Static pressure (Pa)
                - Dynamic pressure (pa) 
                - Air density at given conditions (kg/m**3)
                
    Output: Airspeed (m/s)
        
    """        
    airspeed = np.sqrt(2*(pressureTotal - pressureDynamic)/rhoAir)
    return airspeed

def calcDynPressure(airspeed,rhoAir=rhoAirStandard):
    """
    Description:
    Returns the dynamic pressure for a given airspeed and a given
    air density

    Parameters: - Airspeed (m/s)
                - Air density (kg/m**3) 
                
    Output: Dynamic Pressure (pa)
    """
    q = 1/2*rhoAir*airspeed**2
    return q
    
    

def calcIsentropicFormDeux(mode,givenValue,givenPoint,mach,gamma=1.4):
    # Calcula cualquier presión, temperatura o densidad para cualquier modo (givenReservoirCond, givenPointCond)
    
    """
    Description:
    Returns the temperature/pressure/density using the second form of the
    isentropic energy form involving a reservoir and a certain point in the 
    flowstream

    Parameters: - Temperature(K)/Pressure(pa)/Air Density(kg/m**3) Given (givenValue)
                - Mach Number
                - Mode (T,P,rho)
                - Gamma (optional)
                - given Point (reservoir or flow)
                
    Output: Temperature(K)/Pressure(pa)/Air Density(kg/m**3) Desired
    """    
    if mode == "T":
        if givenPoint == "flow":
            desiredPointValue = givenValue*(1 + (gamma-1)/2*mach**2) #To
        elif givenPoint == "reservoir":
            desiredPointValue = givenValue*1/(1 + (gamma-1)/2*mach**2) #T1
        else:
            raise ValueError("Not a valid input for point")
    elif mode == "P":
        if givenPoint == "flow":
            desiredPointValue = givenValue*(1 + (gamma-1)/2*mach**2)**(gamma/(gamma-1)) #P0
        elif givenPoint == "reservoir":
            desiredPointValue = givenValue*1/((1 + (gamma-1)/2*mach**2)**(gamma/(gamma-1))) #P1
        else:
            raise ValueError("Not a valid input for point")            
    elif mode == "rho":
        if givenPoint == "flow":
            desiredPointValue = givenValue*(1 + (gamma-1)/2*mach**2)**(1/(gamma-1)) #rho0
        elif givenPoint == "reservoir":
            desiredPointValue = givenValue*1/((1 + (gamma-1)/2*mach**2)**(1/(gamma-1))) #rho1
        else:
            raise ValueError("Not a valid input for point")
    else:
        raise ValueError("Not a valid input for mode")        
    return desiredPointValue
    
    
# Pressure distribution y flujo viscoso ------------------------------------

def calcReynoldsNo(rhoFluid,x,flowspeed,viscocityDyn=viscocityDynAir):
    """
    Description:
    Returns the Reynolds number of a flow for a given dynamic pressure

    Parameters: - Density (kg/m**3)
                - Viscocity (pa*s)
                - flowspeed (m/s)
                - length (m)
                
    Output: Reynolds Number
    """    
    Re = rhoFluid*flowspeed*x/viscocityDyn
    return Re
    
 
def calcHeightBLPrandtl(noReynolds,x):   
       """
       Description:
       Returns the Boundary layer height by Prandtl's Boundary layer
       expression (Flat thin Plate).
       Parameters: - Reynolds number 
       Output: Boundary layer height in the given point
       """
       delta = 5.2*x*1/(np.sqrt(noReynolds))
       return delta    
    
def calcHeightAndFricCoeffBL(noReynolds,x,mode="laminar"):
       """
       Description:
       Returns the calculated integrated expressions for friction coefficient and BL height for
       laminar flows and experimental calculated values for turbulent regime.
       Parameters: - Reynolds number 
                   - length (m)
                   - Mode ("laminar" or "turbulent")
       Output: Boundary layer height in the given point
               Friction Coefficient at given point
       """    
       if mode == "laminar":
           delta = 5.2*x*1/(np.sqrt(noReynolds))
           coeffFriction = 1.328/(np.sqrt(noReynolds))
       elif mode == "Turbulent":
           delta = 0.37*x/((noReynolds**0.2))
           coeffFriction = 0.074/((noReynolds**0.2))
       return (delta,coeffFriction)

def plotCriticalMachForCp(cpZero):
    """
    Description:
    Plot the graphs for critical Mach Number and corrected Cp
    Parameters: - Reynolds number 
                - Pressure Coefficient without Correction (m)
    Output: - None (Plot graph)
            - Intersection Values (matchVals)
    """  

    machLim = 1 #[m] # 50000
    step = 0.001
        
    machVec = np.arange(0.4,machLim,step) #[m]
    ecn = np.zeros(len(machVec))
    cpCrit = np.zeros(len(machVec))
    matchVals = []
    
    for i in range(len(machVec)):
        if machVec[i] >= 1:    
            beta = 1
        else:
            beta = (1-machVec[i]**2)**(1/2)
        
        cpCrit[i] = -cpZero/beta
    
        ecn[i] = 2/(gammaAir*-machVec[i]**2)*(((2+(gammaAir-1)*machVec[i]**2)/(gammaAir+1))**(gammaAir/(gammaAir-1))-1) 
        
        if abs(cpCrit[i] - ecn[i]) <= 0.01:
            matchVals.append(machVec[i])
    
    plt.plot(machVec,ecn,machVec,cpCrit)
    plt.show()
    return matchVals


def calcBetaCorrectionFactor(noMachStream):
    """
    Description:
    Calculates the beta correction factor for M < 0.7
    Parameters: - Mach number for stream flight conditions           
    Output: - Correction Factor Beta
    """       
    beta = np.sqrt(1-noMachStream**2)
    return beta

def calcLiftGradientWing(liftGradientAirfoil,aspectRatio,effSpan,beta=1):
    """
    Description:
    Calculates lift gradient for wing for a given set of aircraft parameters.
    Parameters: - Airfoil Lift Gradient (1/deg) 
                - Aspect Ratio 
                - Span Efficiency factor 
                - Beta Correction Factor (optional)
                
    Output: - Wing Lift Gradient (1/deg)
    """
     
    a = liftGradientAirfoil/beta*1/(1 + toDegConv*liftGradientAirfoil/(np.pi*aspectRatio*effSpan))
    return a    
    

def calcAirfoilLiftDragDynFactor(rho,airspeed,chordLength):
    """
    Description:
    Calculates proportional factor for L = factor*Cl or D = factor*Cd for
    an airfoil.
    Parameters: - Air density at given conditions [kg/m**3]
                - Air speed [m/s]
                - chord Length [m]
    Output: - Airfoil Dynamic Factor (N/m *1m)  
    """    
    dynFac = calcDynPressure(airspeed,rho)*chordLength
    return dynFac
    
def calcAirfoilMomentDynFactor(rho,airspeed,chordLength):
    """
    Description:
    Calculates proportional factor for M = factor*Cm for
    an airfoil.
    Parameters: - Air density at given conditions [kg/m**3]
                - Air speed [m/s]
                - chord Length [m]
    Output: - Airfoil Dynamic Factor (N *1m)  
    """    
    dynFac = calcDynPressure(airspeed,rho)*chordLength**2
    return dynFac

def calcWingLiftDragDynFactor(rho,airspeed,S):
    """
    Description:
    Calculates proportional factor for L = factor*Cl or D = factor*Cd for
    an airfoil.
    Parameters: - Air density at given conditions [kg/m**3]
                - Air speed [m/s]
                - Wetted Area S [m**2]
    Output: - Airfoil Dynamic Factor (N)  
    """    
    dynFac = calcDynPressure(airspeed,rho)*S
    return dynFac
    
def calcWingMomentDynFactor(rho,airspeed,chordLength,S):
    """
    Description:
    Calculates proportional factor for M = factor*Cm for
    an airfoil.
    Parameters: - Air density at given conditions [kg/m**3]
                - Air speed [m/s]
                - chord Length [m]
                - Wetted Area S [m**2]
    Output: - Airfoil Dynamic Factor (N*m)  
    """    
    dynFac = calcDynPressure(airspeed,rho)*chordLength*S
    return dynFac

def calcWingInducedDragCoeff(liftCoefficient,aspectRatio,effSpan):
    """
    Description: Calculates Induced Drag Coefficient Term
    Parameters: - lift Coefficient
                - Aspect Ratio 
                - Span Efficiency factor 
    Output: - Wing Induced Drag Coefficient Cdi (1/deg)
    """
    cdi = liftCoefficient**2*1/(np.pi*aspectRatio*effSpan)
    return cdi

def calcAirspeedForCL(lift,rho,wettedArea,coeffLiftDesired):
    """
    Description: Calculates necessary airpseed for a calculated CL depending on the circunstances.
    Parameters: - Lift [N]
                - Air Density rho [kg/m**3] 
                - Wetted Area [m**2]
                - Desired Coefficient of Lift (CL)
                
    Output: - Airspeed for desired CL [m/s]
    """  
    vel = abs((2*lift/(rho*wettedArea*coeffLiftDesired))**(1/2))
    return vel
    


