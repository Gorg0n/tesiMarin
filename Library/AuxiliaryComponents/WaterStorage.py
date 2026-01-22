from Library.AuxiliaryComponents.AuxiliaryComponent import AuxiliaryComponent
from Library.AuxiliaryComponents.Constants import *
import numpy as np


class WaterStorage(AuxiliaryComponent):
    def __init__(self, id, replacement_year, replacement_cost, sp_capex, sp_opex, capacity, soc_min=0.1, soc_max=1, tech='waterstorage'):
        """
        :param id: identification code
        :param tech:string --> 'waterstorage', 'hydrogenstorage', 'desalinator', 'compressor'
        :param replacement_year:float
        :param replacement_cost:float
        :param capacity:float --> tank volume [m3]
        :param soc_min, soc_max:float --> minimum and maximum state of charge of the storage [-]
        """
        super().__init__(id=id, tech=tech, replacement_year=replacement_year, replacement_cost=replacement_cost, sp_capex=sp_capex, sp_opex=sp_opex, capacity=capacity)

        self.capacity = capacity
        self.soc_min = soc_min
        self.soc_max = soc_max

    def compute_output(self, water_input, water_output, soc_in): 
        """
        :param water_input:float --> Input water flow rate [kg]
        :param water_output:float --> Output water flow rate [kg]
        OUTPUT:
            soc:float --> State of charge of the storage [-]
            water_surplus:float --> Surplus water [kg]
            water_deficit:float --> Deficit water [kg]
            delta_w:float --> variation of water (stored if positive, given if negative) [kg]
        """

        soc=soc_in

        water_balance = water_input - water_output # kg
        rho_water = 1000 # kg/m3
        water_storable = self.capacity * (self.soc_max - soc) * rho_water # kg
        water_available = self.capacity * (soc-self.soc_min) * rho_water # kg
        if water_balance >= 0: # water input >= water output  (CARICA)
            if water_balance >= water_storable: 
                soc = self.soc_max
                water_surplus = water_balance - water_storable # quello che entra effettivamente meno quello che viene storato 
                water_deficit = 0
                delta_w = self.capacity*rho_water * (soc-soc_in)
            else: 
                soc += water_balance / (self.capacity * rho_water)
                water_surplus = 0
                water_deficit = 0
                delta_w = self.capacity*rho_water * (soc-soc_in)
        
        else: # SCARICA 
            if abs(water_balance) >= water_available: # l'acqua effettivamente richiesta è maggiore di quella disponibile
                soc = self.soc_min
                water_surplus = 0
                water_deficit = abs(water_balance) - water_available
                delta_w = self.capacity*rho_water *(soc_in-soc)
            else: 
                soc -= abs(water_balance) / (self.capacity * rho_water)
                water_surplus = 0
                water_deficit = 0
                delta_w = self.capacity*rho_water *(soc_in-soc)
        
        delta_w = self.capacity*rho_water *(soc-soc_in) 

        return soc, water_surplus, water_deficit, delta_w

