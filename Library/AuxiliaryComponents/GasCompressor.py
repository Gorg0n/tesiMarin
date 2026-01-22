from Library.AuxiliaryComponents.AuxiliaryComponent import AuxiliaryComponent
from Library.AuxiliaryComponents.Constants import *
import numpy as np


class GasCompressor(AuxiliaryComponent):
    def __init__(self, id, replacement_year, replacement_cost, sp_capex, sp_opex, rated_power, fuel, tech='compressor'):
        """
        :param id: identification code
        :param tech:string --> 'waterstorage', 'hydrogenstorage', 'desalinator', 'compressor'
        :param replacement_year:float
        :param replacement_cost:float
        :param fuel:string --> air, LPG (liquefied gas) , propane (C3H8), methane (CH4), natural gas, hydrogen (H2), syngas, biogas, wood pellets, wood chips
        :param sp_capex:float --> specific capital cost [€/kW]
        :param sp_opex:float --> specific operational cost [€/kW/year]
        :param rated_power:float --> Rated power of the compressor [kW]
        """
        super().__init__(id=id, tech=tech, replacement_year=replacement_year, replacement_cost=replacement_cost, sp_capex=sp_capex, sp_opex=sp_opex, capacity=rated_power)

        self.fuel = fuel

    def compute_output(self, time, temperature_gas_input, pressure_gas_input, gas_flow_rate_input, pressure_gas_output):
        """
        :param time: 1 if hourly analysis, 0.25 if quarterly analysis
        :param temperature_gas_input:float --> Input gas temperature [°C]
        :param pressure_gas_input:float --> Input gas pressure [bar]
        :param gas_flow_rate_input:float --> Input gas flow rate [Nm3/time] (0°C, 1 atm)
        :param pressure_gas_output:float --> Output gas pressure required [bar]
        OUTPUTS: 
            :output power_el_abs:float --> Absorbed power by compressore [W]
            :output q_power_developed:float --> Heat power developed by compressor [W]
            :output temperature_gas_output:float --> Ouput gas temperature [°C]
        """
        r = 2 # compression ratio per stadio 
        n_stages=int(np.log(pressure_gas_output/pressure_gas_input)/np.log(r)) # numero di stadi
        n_parallel=1

        temperature_gas_input = temperature_gas_input + conv_celsius_kelvin
        rho_ref = p_ref * conv_bar_pascal / (r_gas * (t_ref + conv_celsius_kelvin)) # mol/Nm3 
        c_p =self.fuel.a+self.fuel.b*temperature_gas_input+self.fuel.c*temperature_gas_input**2+self.fuel.d*temperature_gas_input**3 # kJ/(K*mol K) 
        c_v = c_p - r_gas
        polytropic_coeff = c_p / c_v
        
        pressure_ratio = (pressure_gas_output / pressure_gas_input)**(1/n_stages)
        p_stages=np.zeros([n_stages+1,1])
        p_stages[0]=pressure_gas_input
        p_stages[n_stages]=pressure_gas_output
        for stage in range(1,n_stages):
            p_stages[stage]=pressure_ratio**(stage)*p_stages[0] 
        w_stage=(polytropic_coeff*r_gas*temperature_gas_input)/(polytropic_coeff-1)*(1-(p_stages[1]/p_stages[0])**((polytropic_coeff-1)/polytropic_coeff))        
        # print(polytropic_coeff,'deve essere positivo compreso tra 1 e 2')
        w_polytropic=w_stage*n_stages
        mol_input = gas_flow_rate_input * rho_ref # mol/h 
        w_iso = -r_gas*temperature_gas_input*np.log(pressure_gas_output/pressure_gas_input)
        p_iso = w_iso * mol_input / (conv_hour_sec * time) # W
        power_el_abs = n_parallel*abs(w_polytropic) * mol_input / (conv_hour_sec * time)
        if power_el_abs == 0:
            eta_iso = 0
        else:
            eta_iso = p_iso / power_el_abs
        # print('efficienza compressore:', eta_iso)
        temperature_gas_output = (temperature_gas_input * (pressure_ratio) ** ((polytropic_coeff - 1) / polytropic_coeff))
        q_power_developed = mol_input / (conv_hour_sec * time) * c_p * (temperature_gas_output - temperature_gas_input)*n_stages*n_parallel
        temperature_gas_output = temperature_gas_output - conv_celsius_kelvin

        return power_el_abs, q_power_developed, temperature_gas_output, eta_iso
    
