# -*- coding: utf-8 -*-
"""
Created on Fri May 16 22:36:06 2025

@author: ELC
"""

import numpy as np
import matplotlib.pyplot as plt


class flightConditions:
    
    
# Módulo que genera plot y carga vectores datos T vs h
# Vector altitud

    def __init__(self,step=0.1):
    
        
        heightLimit = 50000 #[m] # 50000
        stepVec = step #[m]
    
        hVec = np.arange(0,heightLimit,stepVec) #[m]
        TVec = np.zeros(len(hVec))
        PVec = np.zeros(len(hVec))
        rhoVec = np.zeros(len(hVec))
    
    
        # Set aLapse and altitude limits
        aLapseTrop = -6.5 #[K/km]
        aLapseStrato1 = 1 #[K/km]
        aLapseStrato2 = 2.8 #[K/km]
        
        tropSupLimit = 11000 #[m]
        stratoBotLimit = 20000 #[m]
        stratoTopLimit1 = 32000 #[m]
        stratoTopLimit2 = 47000#[m]
        
        # Constants and ref values
        To = 288.15 #[K]
        
        Po = 101325 #[Pa]
        airRSpecific = 287.05 #[J/kg * K]
        g = 9.81 #[m/s**2]
        rho_o = Po/(airRSpecific*To)   # debe dar 1.225 #[kg/m**3]
        
        
        # Units Conversion
        aLapseTrop = aLapseTrop *1/1000
        aLapseStrato1 =aLapseStrato1 *1/1000
        aLapseStrato2 = aLapseStrato2 *1/1000
        
        # Tempearature profile data
        # Temperature,pressure and rho for aLapseShiftPoints------
        # Temperature
        To2 = To + aLapseTrop*tropSupLimit
        To3 = To2 + aLapseStrato1*(stratoTopLimit1-stratoBotLimit)
        To4 = To3 + aLapseStrato2*(stratoTopLimit2-stratoTopLimit1)
        
        # Pressure
        Po2 = Po*(To2/To)**(-g/(aLapseTrop*airRSpecific))
        Po3 = Po2*np.exp(-g/(airRSpecific*To2)*(stratoBotLimit-tropSupLimit))
        Po4 = Po3*(To3/To2)**(-g/(aLapseStrato1*airRSpecific))
        Po5 = Po4*(To4/To3)**(-g/(aLapseStrato2*airRSpecific))
        
        # Rho
        
        rho_o2 = rho_o*(To2/To)**(-g/(aLapseTrop*airRSpecific)-1)
        rho_o3 = rho_o2*np.exp(-g/(airRSpecific*To2)*(stratoBotLimit-tropSupLimit))
        rho_o4 = rho_o3*(To3/To2)**(-g/(aLapseStrato1*airRSpecific)-1)
        rho_o5 = rho_o4*(To4/To3)**(-g/(aLapseStrato2*airRSpecific)-1)
        
        
        
        # Profile generation
        for i in range(len(hVec)):
            if hVec[i] <= tropSupLimit:
                TVec[i] = To + aLapseTrop*(hVec[i])
                PVec[i] = Po*(TVec[i]/To)**(-g/(aLapseTrop*airRSpecific))
                rhoVec[i] = rho_o*(TVec[i]/To)**(-g/(aLapseTrop*airRSpecific)-1)
        
            elif hVec[i] > tropSupLimit and hVec[i] <= stratoBotLimit:
                TVec[i] = To2
                PrhoRatio = np.exp(-g/(airRSpecific*TVec[i])*(hVec[i]-tropSupLimit))
                PVec[i] = Po2*PrhoRatio
                rhoVec[i] = rho_o2*PrhoRatio
        
            elif hVec[i] > stratoBotLimit and hVec[i] <= stratoTopLimit1:
                TVec[i] = To2 + aLapseStrato1*(hVec[i]-stratoBotLimit)
                PVec[i] = Po3*(TVec[i]/To2)**(-g/(aLapseStrato1*airRSpecific))
                rhoVec[i] = rho_o3*(TVec[i]/To2)**(-g/(aLapseStrato1*airRSpecific)-1)
        
            elif hVec[i] > stratoTopLimit1 and hVec[i] <= stratoTopLimit2:
                TVec[i] = To3 + aLapseStrato2*(hVec[i]-stratoTopLimit1)
                PVec[i] = Po4*(TVec[i]/To3)**(-g/(aLapseStrato2*airRSpecific))
                rhoVec[i] = rho_o4*(TVec[i]/To3)**(-g/(aLapseStrato2*airRSpecific)-1)
        
            else:
                TVec[i] = To4
                PrhoRatio2 = np.exp(-g/(airRSpecific*TVec[i])*(hVec[i]-stratoTopLimit2))
                PVec[i] = Po5*PrhoRatio2
                rhoVec[i] = rho_o5*PrhoRatio2
        
        self.altitudeVec = hVec
        self.pressureVec = PVec
        self.temperatureVec = TVec
        self.densityVec = rhoVec
            
    def profileDictioGen(self):
        
    # Generate data Dictionary

        tempDictio = {}
        pressureDictio = {}
        rhoDictio = {}
        
        
        for i in range(len(self.altitudeVec)):
            tempDictio[self.altitudeVec[i]] = self.temperatureVec[i]
            pressureDictio[self.altitudeVec[i]] = self.pressureVec[i]
            rhoDictio[self.altitudeVec[i]] = self.densityVec[i]

        self.tempDictio = tempDictio
        self.pressureDictio = pressureDictio
        self.rhoDictio = rhoDictio
        
    def tempForAltitude(self,altitude,units = "m"):
        
        if units == "m":
            altitude = round(altitude,1)
            return self.tempDictio[altitude]
        elif units == "ft":
            altitude = altitude*0.3048 # ft to m Convert
            altitude = round(altitude,1)
            return self.tempDictio[altitude]
        else:
            print("Not valid unit. Input must be either m or ft")
  
    
    def pressureForAltitude(self,altitude,units = "m"):
        
        if units == "m":
            altitude = round(altitude,1)
            return self.pressureDictio[altitude]
        elif units == "ft":
            altitude = altitude*0.3048 # ft to m Convert
            altitude = round(altitude,1)
            return self.pressureDictio[altitude]
        else:
            print("Not valid unit. Input must be either m or ft")
            
    def densityForAltitude(self,altitude,units = "m"):
         
        if units == "m":
            altitude = round(altitude,1)
            return self.rhoDictio[altitude]
        elif units == "ft":
            altitude = altitude*0.3048 # ft to m Convert
            altitude = round(altitude,1)
            return self.rhoDictio[altitude]
        else:
            print("Not valid unit. Input must be either m or ft")   
            
    def plotTempProfile(self):          
            
        # Plot
        plt.title('Temperature Profile Vs Altitude')
        plt.plot(self.temperatureVec,self.altitudeVec)
        plt.xlabel('Temperature [K]')
        plt.ylabel('Altitude [m]')
        plt.show()
        
    def plotDensityProfile(self):
       # Plot
       
       plt.title('Density Profile Vs Altitude')
       plt.plot(self.densityVec,self.altitudeVec)
       plt.xlabel('Density [kg/m**3]')
       plt.ylabel('Altitude [m]')
       plt.show()
               
    def altitudeForlocalRho(self,rhoInput):
        for i in range(len(self.densityVec)):
            if abs(self.densityVec[i] - rhoInput) <= 0.001:
                hOutput = self.altitudeVec[i]
                break
        return hOutput
        
      
    def giveHeightData(self):
        return self.altitudeVec
        
        
    # def checkForDataGen:
        
    # try:
    #     self.altitudeVec
    #     print("")