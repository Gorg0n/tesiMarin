from Library.AuxiliaryComponents.AuxiliaryComponent import AuxiliaryComponent
from Library.AuxiliaryComponents.Constants import *
import numpy as np


class GasStorage(AuxiliaryComponent):
    def __init__(self, id, replacement_year, replacement_cost, sp_capex, sp_opex, fuel, capacity, pressure_max, tech='GasStorage'):
        """
        :param id: identification code
        :param tech:string --> 'waterstorage', 'hydrogenstorage', 'desalinator', 'compressor', 'GasStorage'
        :param replacement_year:float
        :param replacement_cost:float
        :param fuel:string --> air, LPG (liquefied gas) , propane (C3H8), methane (CH4), natural gas, hydrogen (H2), syngas, biogas, wood pellets, wood chips
        :param capacity:float --> tank volume [m3]
        :param pressure_max:float --> maximum admissible pressure [bar]
        """

        super().__init__(id=id, tech=tech, replacement_year=replacement_year, replacement_cost=replacement_cost, sp_capex=sp_capex, sp_opex=sp_opex, capacity=capacity)

        self.fuel = fuel
        self.capacity = capacity
        self.pressure_max = pressure_max

    def compute_output(self, gas_flow_rate_input, gas_flow_rate_output, initial_pressure_level, temperature_storage):
        """
        :param gas_flow_rate_input:float --> Input gas flow rate [Nm3/time] 
        :param gas_flow_rate_output:float --> Required output gas flow rate [Nm3/time] 
        :param initial_pressure_level:float --> initial pressure level in the tank=initial pressure/maximum allowable pressure ∈[0,1] (0:completely empty, 1:completely full)
        :param temperature_storage:float --> temperature of gas in the storage [°C] 
        OUTPUT:
            pressure:float --> pressure of gas in the storage tank [bar]
            pressure_level:float --> pressure level in the tank=pressure/maximum allowable pressure ∈[0,1]
            v:float --> volume of gas stored in tank [Nm3] 
            v_surplus:float --> surplus volume of gas [Nm3/time] 
            v_unmet:float --> deficit volume of gas [Nm3/time]
        """
        rho_ref = p_ref * conv_bar_pascal / (r_gas * (t_ref + conv_celsius_kelvin)) # mol/Nm3 
        # rho_ref = p_ref * conv_bar_pascal / (r_gas * 293.15) # mol/Nm3 in condizioni normali 
        p_initial = initial_pressure_level * self.pressure_max * conv_bar_pascal  # Pa
        mol = p_initial * self.capacity / (r_gas * (temperature_storage + conv_celsius_kelvin)) # mol 

        temperature_storage = temperature_storage + conv_celsius_kelvin # K
        mol_max = self.pressure_max * conv_bar_pascal * self.capacity / (r_gas * temperature_storage) # mol
        mol_storable = mol_max - mol
        mol_available = mol
        mol_input = gas_flow_rate_input * rho_ref # mol
        mol_output = gas_flow_rate_output * rho_ref # mol

        mol_balance = mol_input - mol_output
        if mol_balance >= 0: # CARICA 
            if mol_balance >= mol_storable:
                mol = mol_max
                mol_surplus = mol_balance - mol_storable
                mol_unmet = 0
            else:
                mol += mol_balance
                mol_surplus = 0
                mol_unmet = 0
        else: # SCARICA
            if abs(mol_balance) >= mol_available:
                mol = 0
                mol_surplus = 0
                mol_unmet = abs(mol_balance) - mol_available
            else:
                mol -= abs(mol_balance)
                mol_surplus = 0
                mol_unmet = 0
        pressure = (mol * r_gas * temperature_storage / self.capacity) / conv_bar_pascal # bar
        pressure_level = pressure / self.pressure_max 
        v = mol / rho_ref # Nm3 
        v_surplus = mol_surplus / rho_ref # Nm3
        v_unmet = mol_unmet / rho_ref  # Nm3

        return pressure, pressure_level, v, v_surplus, v_unmet, gas_flow_rate_input, gas_flow_rate_output
