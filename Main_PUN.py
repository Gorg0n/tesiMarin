import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle

path_dir = os.getcwd()  # os.getcwd() prende la directory corrente, che è Library
sys.path.append(path_dir)

from Library.Utility.Funzioni import compute_h2_demand, compute_hydrogen_density 
from Library.Electrolyzer.PEMclass import PEM
from Library.ControlSystems.ControlSystem_PUN import Controller
from Library.AuxiliaryComponents.GasCompressor import GasCompressor
from Library.AuxiliaryComponents.GasStorage import GasStorage
from Library.AuxiliaryComponents.WaterStorage import WaterStorage
from Library.AuxiliaryComponents.Desalinator import Desalinator
from Library.AuxiliaryComponents.Fuel import Fuel
from Library.AuxiliaryComponents.Constants import *


# ------ OPERATING CONDITIONS ------
T = 80
p = 30
time = 1 
lifetime = 20

# ------ INPUT WIND POWER ------
# Import wind power data from excel 
xlsx_path = os.path.join(path_dir, 'Data_input', 'Wind_power.xlsx')
wind_power = pd.read_excel(xlsx_path) # Dataframe  [kW]
wind_power = wind_power.iloc[:, 0]  # Array [kW]

# ------ HYDROGEN DEMAND ------
hydrogen = Fuel(fuel='hydrogen')
steel_demand = 500000 #ton/year
yearly_h2_demand = steel_demand*60*0.78/1000 # ton/year
year_length = len(wind_power) 
h2_demand = compute_h2_demand(yearly_h2_demand=yearly_h2_demand, year_length=year_length, time='half-day')   
for i in range(len(h2_demand)):
    h2_demand[i] = h2_demand[i] / hydrogen.rho # [Nm3/h]

# ------ ELECTROLYZER PARAMETERS ------
h2_prod_nom = 1062 # kg/day
pem_n_parallel = 8
pem_n_series = 1
pem_n_groups = 8 
p_nom = 2500     
phi_mem = 0.0055 # cm  
sp_capex = 608 * pem_n_groups # [€/kW]
sp_opex = 608*0.015 * pem_n_groups # [€/kW/year]
repl_cost = 608*0.12 * pem_n_groups # [€/kW]
repl_year = 10 # [years]
pem = PEM(p_nom=p_nom, nominal_prod_rate=h2_prod_nom/24, V_module=30000, n_parallel=pem_n_parallel, n_series=pem_n_series, phi_mem=phi_mem, sp_capex=sp_capex, sp_opex=sp_opex, replacement_cost=repl_cost, replacement_year=repl_year) 

# ------ DESALINATOR PARAMETERS ------
nom_water_diss_m3_h = 30 # m3/h
nom_water_diss_kg = nom_water_diss_m3_h * rho_water # kg/h
power_diss = 37 # kW
sp_capex = 488 # [€/m3]
sp_opex = 0.18 # [€/m3/year]
diss_capacity = nom_water_diss_m3_h *24 #m3/day
desalinator = Desalinator(id=1, replacement_year=0, replacement_cost=0, sp_capex=sp_capex, sp_opex=sp_opex, capacity=diss_capacity, consumption=power_diss, nom_water_prod=nom_water_diss_kg)

# ------ WATER STORAGE PARAMETERS ------
margin = 1 # days
nom_wat_cons_pem = 15.9 * h2_prod_nom *pem_n_parallel*pem_n_series*pem_n_groups # l/day
nom_wat_prod_diss = nom_water_diss_m3_h*24 *1000 # l/day
wat_stor_capacity = (nom_wat_prod_diss - nom_wat_cons_pem + margin * nom_wat_cons_pem)/rho_water #m3
sp_capex = 236 # [€/m3]
watstorage = WaterStorage(id=2, replacement_year=0, replacement_cost=0, sp_capex=sp_capex, sp_opex=0, capacity=wat_stor_capacity)

# ------ COMPRESSOR PARAMETERS ------
repl_year = 10 # [years]
repl_cost = 2353 # [€/kW]
sp_capex = 2353 # [€/kW]
sp_opex = sp_capex*0.03 # [€/kW/year]
compressor = GasCompressor(id=3, replacement_year=repl_year, replacement_cost=repl_cost, sp_capex=sp_capex, sp_opex=sp_opex, rated_power=4032, fuel=hydrogen)

# ------ HYDROGEN STORAGE PARAMETERS ------
pressure = 500 # bar
module_hs_capacity = 6.5 #kg da scheda tecnica Mahytech
h2_stor_days = 0.5
n_modules_hs = round(h2_prod_nom*pem_n_parallel*pem_n_series*pem_n_groups * h2_stor_days / module_hs_capacity)
h2_capacity_kg = n_modules_hs * module_hs_capacity #kg
rho_h2_stor = massa_molare_h2 * pressure*1e5 / (r_gas*(T+273.15)) # [kg/m3]
h2_capacity_m3 = h2_capacity_kg / rho_h2_stor #m3/day
sp_capex = 455 * rho_h2_stor # [€/m3]
sp_opex = 0.02*sp_capex # [€/m3/year]
h2_storage = GasStorage(id=4, replacement_year=0, replacement_cost=0, sp_capex=sp_capex, sp_opex=sp_opex, capacity=h2_capacity_m3, pressure_max=pressure, fuel=hydrogen)

# ------ CONTROLLER ------
# Storage state of charge
soc_water_in = 0 
soc_water_limit = 0.8

# Price data
threshold = 98 # €/MWh
PUN_table = pd.read_excel(os.path.join(path_dir, 'Data_input', 'PUN_GME_hour_20192024.xlsx'))
PUN_table_2024 = PUN_table[pd.to_datetime(PUN_table.iloc[:,0]).dt.year == 2024] # estract 2024 data
PUN_2024 = PUN_table_2024.iloc[:,1] # €/MWh 
PUN_2024 = np.array(PUN_2024) # €/MWh ARRAY 

# Controller class
controller = Controller(PEM=pem, Desalinator=desalinator, WaterStorage=watstorage, Compressor=compressor, H2Storage=h2_storage, fuel=hydrogen, n_group=pem_n_groups, soc_water_in=soc_water_in, soc_water_limit=soc_water_limit)
output_controller = controller.energy_perfomance(production=wind_power, T=T, time=time, p=p, hourly_h2_demand=h2_demand, threshold=threshold, el_price=PUN_2024)

# Save output power flows
df_out_power=pd.DataFrame({'p_wind[kW]':output_controller['production_out'],'p_pem_cons[kW]':output_controller['p_pem_cons'], 
                     'p_diss[kW]':output_controller['consumption_diss'], 'p_compr[kW]':output_controller['power_compr'],
                     'surplus[kW]':output_controller['surplus_tot']})
df_out_power.to_excel(f'Data_output/Complete_system_PUN/Power_flows.xlsx', index=False)

# Save output hydrogen flows
df_out_hydrogen = pd.DataFrame({'p_wind[kW]':output_controller['production_out'],'H2_prod[kg/h]':output_controller['H2_production'], 'n_elec_on':output_controller['n_elec_on'], 
                                'h2_stor_input[Nm3/h]':output_controller['gas_flow_rate_input_storage'], 'h2_demand[Nm3/h]':output_controller['h2_stor_demand'], 
                                'P_stor[bar]':output_controller['pressure_gas_in_storage'], 'pressure_level':output_controller['pressure_level'],
                                'V_stor[Nm3/h]':output_controller['volume_gas_stored'], 'V_surplus[Nm3/h]':output_controller['v_surplus'], 'V_unmet[Nm3/h]':output_controller['v_unmet'], 
                                })
df_out_hydrogen.to_excel(f'Data_output/Complete_system_PUN/Hydrogen_flows.xlsx', index=False)

# Save output water flows
wat_av=output_controller['soc_wat']*wat_stor_capacity*rho_water-watstorage.soc_min*wat_stor_capacity*rho_water
df_out_water = pd.DataFrame({'p_wind[kW]':output_controller['production_out'],'water_cons[kg/h]':output_controller['water_consumption'], 
                             'wat_prod_diss':output_controller['wat_prod_diss'], 'soc_wat':output_controller['soc_wat'], 'water_av[kg]':wat_av,
                             'delta_wat':output_controller['delta_w'], 'water_surplus':output_controller['water_surplus'], 
                             'water_deficit':output_controller['water_deficit'], 'status_diss':output_controller['status_diss']})
df_out_water.to_excel(f'Data_output/Complete_system_PUN/Water_flows.xlsx', index=False)