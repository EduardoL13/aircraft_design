# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 23:39:54 2025

@author: ELC
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
from aeronautics.unitsAE import *
from mathE.generalFunctions import *

class aircraft:
    
    def setWeight(self,mass):
        """
        Description:
        Sets Aircraft weight
        Parameters: - Mass (kg)
        Output: None (loads parameters for class object)  
        """          
        self.weight = mass*9.81
    
    def setDragPolarParams(self,cd0,k1,k2): # En el futuro, en caso de que haya más términos para la función de drag polar, se puede crear un vector que tenga todos los coeficientes
       """
       Description:
       Sets parameters for aircraft/wing drag polar characteristics (Assumes quadratic drag polar)
       Parameters: - Parasitic drag cd0
                   - Constant k1
                   - constant k2
       Output: None (loads parameters for class object)  
       """            
       self.cd0 = cd0 
       self.k1 = k1
       self.k2 = k2
       
       
    def setFlightConditionsParams(self,rhoAir): # Más adelante se puede mirar cómo implementar el módulo FC para esto, si vale la pena
       """
       Description:
       Sets parameters for aircraft flight conditions (air density for the moment being)
       Parameters: - Air Density (rhoAir) [kg/m**3]
       Output: None (loads parameters for class object)  
       """     
       self.rhoAir = rhoAir
       

    
    def setWingGeometryParams(self,wettedAreaS,bSpan=0,aspectRatio=0):   
  
       """
       Description:
       Sets parameters for aircraft wing geometry
       Parameters: - Wing Span (bSpan) Optional [m]
                   - Wetted Area (wetted Area S) Optional [m**2]
       Output: None (loads parameters for class object)  
       """    
       if bSpan != 0:
           self.aspectRatio = bSpan**2/wettedAreaS
           self.wettedAreaS = wettedAreaS
           self.bSpan = bSpan
       else:
           self.aspectRatio = aspectRatio
           self.wettedAreaS = wettedAreaS
           
    def velMin(self,coeffLiftMax):
        """
        Description:
        Returns propeller engine min velocity for existing parameters within the propeller class.
        Parameters: - Max Lift Coefficient 
 
        Output: Aircraft Min. Velocity (velResult) [m/s] 
        """             
        velResult = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, coeffLiftMax)
        return velResult 
    
    
    def minGlideAngle(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns aircraft optimal angle for optimal range in descending operation without thrust (glide)
        Parameters: - None 
 
        Output: Aircfaft Optimal glide angle (gamma) [deg]
                
        """      
        
        # (velOpt,clOpt) = self.maxRange()
        # coeffDrag = self.cd0 + self.k1*clOpt + self.k2*clOpt**2
        # gammaGlide = arcSinDeg((coeffDrag)/(clOpt))
        # return (gammaGlide,clOpt) 
        
        optimumGlideAngleFactor = 2*self.cd0/(np.sqrt(self.cd0/self.k2)) + self.k1 # SALE DE LA CONDICIÓN DE DRAG MÍNIMO, APLICADA DIRECTAMENTE EN LA ECUACIÓN (EOM) DE GLIDE OPTIMO
        optimumGlideAngle1 = arcSinDeg(optimumGlideAngleFactor)
        return optimumGlideAngle1
    
    def minROD(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns propeller optimal minimal (rate of descend)
        Parameters: - None 
 
        Output: Aircraft ROD [m/s]
                Aircraft Optimal Velocity for max endurance (velResult) [m/s] 
                Aircraft Lift Coefficient  (clOpt)

        """      
        
        quadFactor = np.sqrt(self.k1**2 + 12*self.k2*self.cd0) 
        clOpt = (self.k1 + quadFactor)/(2*self.k2)
        velOpt = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, clOpt)
        cd = self.cd0 + self.k1*clOpt + self.k2*clOpt**2
        dragDynFactor = calcWingLiftDragDynFactor(self.rhoAir, velOpt, self.wettedAreaS)
        drag = cd*dragDynFactor
        

        ROC = (-drag*velOpt)/self.weight #
        
        
        
        return (ROC,velOpt,clOpt)  
  

class propellerAircraft(aircraft):
    
    def setEngineParams(self,powerMax,effProp=1):
       """
       Description:
       Sets parameters for aircraft engine curve relevant for cruise conditions
       Parameters: - Engine/Propeller shaft Output (powerMax) [W]
                   - Propeller efficiency (effProp) (depending on the given info)
       Output: None (loads parameters for class object)  
       """             
       self.powerMax = effProp*powerMax
    
    def velMax(self,solType="graphical"):
        """
        Description:
        Returns propeler engine max velocity for existing parameters within the propeller class.
        Parameters: - Solution Type (optional: graphical by default)
 
        Output: Aircraft Max. Velocity (velResult) [m/s]  
        """  
        if solType == "graphical":
            limCL = 20 #[m] # 50000
            step = 0.01
            
            CLVec = np.arange(0.1,limCL,step) #[m]
            ecnGraph = np.zeros(len(CLVec))
            powerMaxVec = np.zeros(len(CLVec))
            velVec = np.zeros(len(CLVec))
            matchVals = []
            
            for i in range(len(CLVec)):

                powerMaxVec[i] = self.powerAvailable

                ecnGraph[i] = np.sqrt(self.weight**3/self.wettedAreaS*2/self.rhoAir*(self.cd0 + self.k2*CLVec[i]**2)**2*1/(CLVec[i]**3)) 
            
                velVec[i] = np.sqrt(self.weight/self.wettedAreaS*2/self.rhoAir*1/CLVec[i])
                
                if abs(ecnGraph[i] - velVec[i]) <= 0.01:
                    matchVals.append(machVec[i])
    
            plt.plot(velVec,ecnGraph,velVec,powerMaxVec)
            plt.show()
            velResult = min(matchVals)
            return velResult
        
        elif solType == "numerical":
            
            CL = smp.symbols('CL', real=True, positive=True)

            ecn = self.powerAvailable**2*self.wettedAreaS*self.rhoAir*1/(self.weight**3*2) - (self.cd0 + self.k2*CL**2)**2*1/CL**3

            res = smp.solve(ecn,CL)
            optCL = min(res) 

            velResult = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, optCL)
            return velResult
        else:
            print("No proper input was given for Propeller Max Velocity calculation")
            return 0
        
   
    
    def maxRange(self):
        """
        Description:
        Returns propeller engine optimal velocity fot max specific range for existing parameters within the propeller class.
        Parameters: - None 
 
        Output: Aircraft Optimal Velocity for max range (velResult) [m/s] 
        """        
        coeffLiftOpt = np.sqrt(self.cd0/self.k2)
        velResult = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, coeffLiftOpt)
        return velResult   

    def maxEndurance(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns propeller engine optimal velocity for max Endurance time for existing parameters within the propeller class.
        Parameters: - None 
 
        Output: Aircraft Optimal Velocity for max endurance (velResult) [m/s] 
        """        
        quadFactor = np.sqrt(self.k1**2 + 12*self.k2*self.cd0) 
        coeffLiftOpt = (self.k1 + quadFactor)/(2*self.k2)
        velResult = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, coeffLiftOpt)
        return (velResult,coeffLiftOpt)   
    
    
    def maxClimbAngle(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns propeller engine optimal for climb
        Parameters: - None 
 
        Output: Aircraft Max climb angle [deg]
                Aircraft Optimal Velocity for max endurance (velResult) [m/s] 
                Aircraft Optimal Lift Coeficient (CL)
        """      
        
        (velOpt,clOpt) = self.maxRange()
        coeffDrag = self.cd0 + self.k1*clOpt + self.k2*clOpt**2
        dragDynFactor = calcWingLiftDragDynFactor(self.rhoAir,velOpt,self.wettedAreaS)
        drag = coeffDrag*dragDynFactor

        gammaMax = arcSinDeg((self.powerMax-drag*velOpt)/(self.weight*velOpt))
        
        return (gammaMax,velOpt,clOpt)    
    
    
    def maxROC(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns propeller optimal ROC (rate of climb)
        Parameters: - None 
 
        Output:   [1] Aircraft Max ROC [m/s]
                  [2] Aircraft Optimal Velocity for max endurance (velResult) [m/s] 
                  [3] Aircraft Lift Coefficient  (clOpt)        
        """      
        
        (velOpt,clOpt) = self.maxEndurance()
        coeffDrag = self.cd0 + self.k1*clOpt + self.k2*clOpt**2
        dragDynFactor = calcWingLiftDragDynFactor(self.rhoAir,velOpt,self.wettedAreaS)
        drag = coeffDrag*dragDynFactor

        ROC = (self.powerMax - drag*velOpt)/self.weight
        
        return (ROC,velOpt,clOpt)       
    
    
    
    
    
    
    #REPETIR LO MISMO PARA JET (considerar bien qué casos se invierten)
    
class JetAircraft(aircraft):
        
    def setEngineParams(self,thrustMax,effProp=1):
        """
        Description:
        Sets parameters for aircraft jet engine curve relevant for cruise conditions
        Parameters: - Jet Output (powerMax) [N]
        Output: None (loads parameters for class object)  
        """             
        self.thrustMax = thrustMax
        
    def velMax(self):
        """
        Description:
        Returns aircraft engine max velocity for existing parameters within the jet class.
        Parameters: - Solution Type (optional: graphical by default)
     
        Output: Aircraft Max. Velocity (velResult) [m/s]  
        """ 
        quadFactor = np.sqrt((self.k1 - self.thrustMax/self.weight)**2 - 4*self.k2*self.cd0)
        optCL = (-(self.k1 - self.thrustMax/self.weight) + quadFactor)/(2*self.k2)
        velResult = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, optCL)
        return (velResult,optCL)
   
    def maxRange(self):
        """
        Description:
        Returns Jet engine optimal velocity fot max specific range for existing parameters within the Jet class.
        Parameters: - None 
     
        Output: Aircraft Optimal Velocity for max range (velResult) [m/s] 
        """        
        quadFactor = np.sqrt(self.k1**2 + 12*self.k2*self.cd0) 
        coeffLiftOpt = (self.k1 + quadFactor)/(6*self.k2)
        velResult = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, coeffLiftOpt)
        return (velResult,coeffLiftOpt)   
 
    def maxEndurance(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns jet engine optimal velocity for max Endurance time for existing parameters within the jet class.
        Parameters: - None 
 
        Output: Aircraft Optimal Velocity for max endurance (velResult) [m/s] 
        """               
        coeffLiftOpt = np.sqrt(self.cd0/self.k2)
        velResult = calcAirspeedForCL(self.weight, self.rhoAir, self.wettedAreaS, coeffLiftOpt)
        return (velResult,coeffLiftOpt)           
    
    def maxClimbAngle(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns jet engine optimal lift coefficient, airspeed and climb angle for jet performance.
        Parameters: - None 
 
        Output: Maximum climb angle (gammaMax) [deg] 
                Optimal Velocity (velOpt) [m/s]
                Optimal Lift Coefficinet (clOpt)
        """               
        (velOpt,clOpt) = self.maxEndurance()
        coeffDrag = self.cd0 + self.k1*clOpt + self.k2*clOpt**2
        dragDynFactor = calcWingLiftDragDynFactor(self.rhoAir,velOpt,self.wettedAreaS)
        drag = coeffDrag*dragDynFactor

        gammaMax = arcSinDeg((self.thrustMax-drag)/self.weight)
        
        return (gammaMax,velOpt,clOpt)           
        
    def maxROC(self):
        # Agregar aquí positive por defecto al discriminant. la idea es que si el valor que uno recibe de velocidad es raro,
        # Se pueda cambiar a "negative" para calcular discriminante negativ y comparar
        """
        Description:
        Returns jet optimal ROC (rate of climb)
        Parameters: - None 
 
        Output: Aircraft Optimal Velocity for max endurance (velResult) [m/s] 
                Aircraft Lift Coefficient  (clOpt)
                Aircraft ROC
        """      
        
        (velOpt,clOpt) = self.maxEndurance()
        coeffDrag = self.cd0 + self.k1*clOpt + self.k2*clOpt**2
        dragDynFactor = calcWingLiftDragDynFactor(self.rhoAir,velOpt,self.wettedAreaS)
        drag = coeffDrag*dragDynFactor

        ROC = (self.thrustMax*velOpt - drag*velOpt)/self.weight
        
        return (ROC,velOpt,clOpt)         
        
        
        