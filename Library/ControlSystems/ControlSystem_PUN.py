import numpy as np
import copy
from Library.AuxiliaryComponents.Constants import *

#Control logic for multiple PEM and multiple electric storage

class Controller:
    def __init__(self, PEM, Desalinator, WaterStorage, Compressor, H2Storage, fuel, n_group, soc_water_in, soc_water_limit):
        """
        :param PEM: list obj by PEM
        :param Desalinator: obj by Desalinator
        :param WaterStorage: obj by WaterStorage
        :param Compressor: obj by GasCompressor
        :param H2Storage: obj by GasStorage
        :param fuel: obj by Fuel
        :param n_group:float number of pem which do not operate simultaneously, but instead are activated sequentially, one after another (series)
        :param soc_water_in:float initial state of charge of the water storage
        :param soc_water_limit:float limit state of charge of the water storage
        """
        self.PEM = PEM
        self.Desalinator = Desalinator
        self.WaterStorage = WaterStorage
        self.Compressor = Compressor
        self.H2Storage = H2Storage
        self.fuel = fuel
        self.n_group = n_group
        self.soc_water_in = soc_water_in
        self.soc_water_limit = soc_water_limit

        self.PEM_water_need = PEM.nominal_prod_rate*15 # [kg/h] NOMINALE 


    def energy_perfomance(self, production, T, time, p, hourly_h2_demand, threshold, el_price):
        '''
        :input production:array --> power from renewable source [kW]  
        :input T:float--> operative temperature of electrolyzer [°C]
        :input p:float--> operative pressure of electrolyzer [bar]
        :input time:float--> 1 if hourly analysis, 0.25 if quarterly analysis
        :input hourly_h2_demand:array --> hourly hydrogen demand [Nm3]
        :input threshold:float --> threshold for electricity price [€/MWh]
        :input el_price:array --> electricity price [€/MWh]
        OUTPUTS (arrays): 
            production_out: power from renewable source [kW]
            H2_production: hydrogen production [kg/h]
            p_pem_cons: power consumed by the PEM [kW]
            surplus_tot: surplus power [kW]
            water_consumption: water consumption [kg/h]
            n_elec_on: number of electrolyzer on
            case: case of operation
            soc: state of charge of the water storage
            status_diss: status of the desalination plant
            consumption_diss: energy consumption of the desalination plant [kWh]
            soc_surplus: surplus of the state of charge of the water storage
            case_diss: case of operation of the desalination plant
            power_compr: power absorbed by the compressor [kW]
            q_compr: heat absorbed by the compressor [kW]
            temp_h2_out: temperature of the hydrogen at the outlet of the compressor [°C]  
        '''
        PEM = self.PEM
        PEM_group = [copy.deepcopy(PEM) for i in range(self.n_group)]
        rho_water = 1000 # kg/m3
        len_ref = len(production)
        desalinator = self.Desalinator
        water_storage = self.WaterStorage
        compressor = self.Compressor
        h2_storage = self.H2Storage

        for electrolyzer in PEM_group:
            electrolyzer.en_perf_evolution = {}
            electrolyzer.en_perf_evolution['H2_production'] = np.zeros(len_ref)
            electrolyzer.en_perf_evolution['elec_consumption'] = np.zeros(len_ref)
            electrolyzer.en_perf_evolution['power_surplus'] = np.zeros(len_ref)
            electrolyzer.en_perf_evolution['system_efficiency'] = np.zeros(len_ref)
            electrolyzer.en_perf_evolution['status'] = np.zeros(len_ref)
            electrolyzer.en_perf_evolution['water_consumption'] = np.zeros(len_ref)

        production_out_ev = np.array(production)
        mode_ev  = np.zeros(len_ref)
        # electrolyzer 
        H2_production_ev = np.zeros(len_ref)
        p_pem_cons_ev = np.zeros(len_ref)
        surplus_tot_ev = np.zeros(len_ref)
        water_consumption_ev = np.zeros(len_ref)
        n_elec_on_ev = np.zeros(len_ref)
        # dissalatore 
        consumption_diss_ev = np.zeros(len_ref)
        status_diss_ev = np.zeros(len_ref)
        water_produced_ev = np.zeros(len_ref)
        # water storage 
        soc_wat_ev = np.zeros(len_ref)
        water_surplus_ev = np.zeros(len_ref)
        water_deficit_ev = np.zeros(len_ref)
        delta_w_ev = np.zeros(len_ref)
        soc_wat_tot = self.soc_water_in
        case_stor_wat_ev = np.zeros(len_ref)
        case_stor_h2_ev = np.zeros(len_ref)
        case_el_ev = np.zeros(len_ref)
        # compressor
        power_compr_ev = np.zeros(len_ref)
        q_compr_ev = np.zeros(len_ref) 
        temp_out_compr_ev = np.zeros(len_ref)
        # hydrogen storage
        p_h2_storage_ev = np.zeros(len_ref)
        p_level_h2_storage_ev = np.zeros(len_ref) 
        v_h2_storage_ev= np.zeros(len_ref) 
        v_h2_storage_surplus_ev=np.zeros(len_ref) 
        v_h2_storage_unmet_ev= np.zeros(len_ref) 
        flow_in_storage_ev = np.zeros(len_ref)
        flow_out_storage_ev = np.zeros(len_ref)
        initial_pressure_level = 0
        available_space_h2stor = h2_storage.capacity #m3
        rho_h2_stor = 38 # [kg/m3] densità dell'idrogeno nello storage

        for i in range(len_ref):
            production_out = production[i]
            h2_demand = hourly_h2_demand[i]
            pun = el_price[i] # €/MWh
            surplus_tot = production_out
            # electrolyzer 
            H2_production_tot = 0
            p_pem_cons_tot = 0
            water_consumption_tot = 0
            n_elec_on_tot = 0
            case_el = 0 
            # dissalatore
            status_diss = 0
            consumption_diss_tot = 0
            # water storage 
            water_surplus_tot = 0
            water_deficit_tot = 0 
            delta_w_tot = 0     
            case_stor_wat = 0    

            max_n_pem_series = np.ceil(production_out / PEM.p_nom)
            max_h2_prod = max_n_pem_series * PEM.nominal_prod_rate / self.fuel.rho # [Nm3/h]
            power_el_abs, q_power_developed, temperature_gas_output, eta_iso = compressor.compute_output(time=time, temperature_gas_input=T, pressure_gas_input=p, gas_flow_rate_input=max_h2_prod, pressure_gas_output=h2_storage.pressure_max)
            max_power_el_abs = float(power_el_abs*1e-3) # [kW] 
            
            if pun < threshold: # se il prezzo è basso, allora si può usare l'energia per produrre idrogeno
                if production_out > max_power_el_abs:  # non maggiore uguale perché se Pw=0 e Pcompr=0 devo andare al caso 5 
                    if soc_wat_tot == water_storage.soc_max:    
                        power_to_pem = production_out - max_power_el_abs
                        case = 1
                    elif self.soc_water_limit < soc_wat_tot < water_storage.soc_max:
                        power_to_pem = production_out - max_power_el_abs
                        case = 2
                    else: # 0.2 < soc_wat_tot < self.soc_water_limit
                        if production_out - desalinator.consumption >= max_power_el_abs:
                            power_to_pem = production_out - desalinator.consumption - max_power_el_abs
                            case = 3
                        else: 
                            power_to_pem = 0
                            case = 4
                else: 
                    power_to_pem = 0
                    case = 5

                PEM_sorted = sorted(PEM_group, key=lambda PEM_group: PEM_group.status, reverse=True)
                for electrolyzer in PEM_sorted:
                    if power_to_pem >= electrolyzer.min_power_input:
                        H2_production, p_pem_cons, power_surplus, system_efficiency, status, water_consumption = electrolyzer.energy_performance(
                        production=power_to_pem, time=time, p=p, T=T)
                        electrolyzer.status=status
                        memories = {'H2_production': H2_production, 'elec_consumption': p_pem_cons, 'power_surplus':power_surplus, 'system_efficiency': system_efficiency,
                                    'status': status, 'water_consumption': water_consumption}
                        for key, value in memories.items():
                            electrolyzer.en_perf_evolution[key][i] = value
                        case_el = 1

                    else: 
                        H2_production, p_pem_cons, power_surplus, system_efficiency, status, water_consumption = electrolyzer.energy_performance(
                        production=0, time=time, p=p, T=T)
                        electrolyzer.status=status
                        memories = {'H2_production': H2_production, 'elec_consumption': p_pem_cons, 'power_surplus':power_surplus, 'system_efficiency': system_efficiency,
                                    'status': status, 'water_consumption': water_consumption}
                        for key, value in memories.items():
                            electrolyzer.en_perf_evolution[key][i] = value
                        case_el = 2

                    soc_wat_pem = soc_wat_tot
                    available_water = water_storage.capacity*rho_water * (soc_wat_pem-water_storage.soc_min)  # [kg]
                    if available_water >= water_consumption:  
                        volume_h2_prod = H2_production / rho_h2_stor # [m3] volume calcolato con la densità che avrà nello storage di idrogeno
                        if available_space_h2stor >= volume_h2_prod:
                            power_to_pem -= p_pem_cons
                            H2_production_tot += H2_production
                            p_pem_cons_tot += p_pem_cons
                            water_consumption_tot += water_consumption
                            surplus_tot -= p_pem_cons                   
                            n_elec_on_tot += status

                            soc_wat_pem -= water_consumption / (water_storage.capacity*rho_water) 
                            available_space_h2stor -= volume_h2_prod

                            case_stor_h2 = 1 
                            case_stor_wat = 1

                        else: 
                            H2_production_tot += 0
                            p_pem_cons_tot += 0
                            water_consumption_tot += 0
                            n_elec_on_tot += 0

                            soc_wat_pem -= 0
                            available_space_h2stor -= 0

                            case_stor_h2 = 2
                            case_stor_wat = 2

                    else:
                        H2_production_tot += 0
                        p_pem_cons_tot += 0
                        water_consumption_tot += 0
                        n_elec_on_tot += 0
                        soc_wat_pem -= 0
                        available_space_h2stor -= 0 

                        case_stor_h2 = 3
                        case_stor_wat = 3     

                # Definizione della lookup table
                mode_lookup = {
                    (1, 1, 1):  1, (1, 1, 2):  2, (1, 1, 3):  3,
                    (1, 2, 1):  4, (1, 2, 2):  5, (1, 2, 3):  6,
                    (2, 1, 1):  7, (2, 1, 2):  8, (2, 1, 3):  9,
                    (2, 2, 1): 10, (2, 2, 2): 11, (2, 2, 3): 12,
                    (3, 1, 1): 13, (3, 1, 2): 14, (3, 1, 3): 15,
                    (3, 2, 1): 16, (3, 2, 2): 17, (3, 2, 3): 18,
                    (4, 1, 1): 19, (4, 1, 2): 20, (4, 1, 3): 21,
                    (4, 2, 1): 22, (4, 2, 2): 23, (4, 2, 3): 24,
                    (5, 1, 1): 25, (5, 1, 2): 26, (5, 1, 3): 27,
                    (5, 2, 1): 28, (5, 2, 2): 29, (5, 2, 3): 30
                }
                # Recupero del valore dal dizionario
                mode = mode_lookup.get((case, case_el, case_stor_wat), None)
                

                # DESALINATION AND WATER STORAGE 
                if case == 1: 
                    power_to_diss = 0
                elif case == 2: 
                    power_to_diss = production_out - p_pem_cons_tot - max_power_el_abs
                elif case == 3 or case == 4: 
                    power_to_diss = desalinator.consumption
                else: # case =5 
                    if production_out >= desalinator.consumption:
                        power_to_diss = production_out
                    else:
                        power_to_diss = 0

                water_produced, consumption, status = desalinator.compute_output(production=power_to_diss)
                soc_water, water_surplus, water_deficit, delta_w = water_storage.compute_output(water_input=water_produced, water_output=water_consumption_tot, soc_in=soc_wat_tot)

                consumption_diss_tot += consumption 
                status_diss = status
                soc_wat_tot = soc_water
                water_surplus_tot += water_surplus 
                water_deficit_tot += water_deficit 
                delta_w_tot += delta_w
                surplus_tot -= consumption

                # COMPRESSOR
                h2_flow_rate = H2_production_tot / self.fuel.rho # [Nm3/h]
                power_el_abs, q_power_developed, temperature_gas_output, eta_iso = compressor.compute_output(time=time, temperature_gas_input=T, pressure_gas_input=p, gas_flow_rate_input=h2_flow_rate, pressure_gas_output=h2_storage.pressure_max)
                surplus_tot -= power_el_abs*1e-3

                # HYDROGEN STORAGE
                pressure, pressure_level, v, v_surplus, v_unmet, gas_flow_rate_input, gas_flow_rate_output = h2_storage.compute_output(gas_flow_rate_input=h2_flow_rate, gas_flow_rate_output=h2_demand, initial_pressure_level=initial_pressure_level, temperature_storage=temperature_gas_output)  
                initial_pressure_level = pressure_level
                available_space_h2stor = h2_storage.capacity - v * 1.01325/h2_storage.pressure_max*(temperature_gas_output+273.15)/273.15   # così tiene conto del fatto che riaumenta perché si è liberato spazio grazie all'idrogeno uscito 
                rho_h2_stor  = massa_molare_h2 * h2_storage.pressure_max*1e5 / (r_gas*(temperature_gas_output+273.15)) # [kg/m3]

            else: # PUN > threshold
                H2_production, p_pem_cons, power_surplus, system_efficiency, status, water_consumption = electrolyzer.energy_performance(
                production=0, time=time, p=p, T=T)
                electrolyzer.status=status
                memories = {'H2_production': H2_production, 'elec_consumption': p_pem_cons, 'power_surplus':power_surplus, 'system_efficiency': system_efficiency,
                            'status': status, 'water_consumption': water_consumption}
                for key, value in memories.items():
                    electrolyzer.en_perf_evolution[key][i] = value 

                H2_production_tot += H2_production
                p_pem_cons_tot += p_pem_cons
                water_consumption_tot += water_consumption
                n_elec_on_tot += status         

                power_to_diss = 0
                water_produced, consumption, status = desalinator.compute_output(production=power_to_diss)
                consumption_diss_tot = consumption
                status_diss = status
                soc_water, water_surplus, water_deficit, delta_w = water_storage.compute_output(water_input=water_produced, water_output=water_consumption_tot, soc_in=soc_wat_tot)
                
                h2_flow_rate = H2_production_tot / self.fuel.rho # [Nm3/h]
                power_el_abs, q_power_developed, temperature_gas_output, eta_iso = compressor.compute_output(time=time, temperature_gas_input=T, pressure_gas_input=p, gas_flow_rate_input=h2_flow_rate, pressure_gas_output=h2_storage.pressure_max)

                pressure, pressure_level, v, v_surplus, v_unmet, gas_flow_rate_input, gas_flow_rate_output = h2_storage.compute_output(gas_flow_rate_input=h2_flow_rate, gas_flow_rate_output=h2_demand, initial_pressure_level=initial_pressure_level, temperature_storage=temperature_gas_output)  
                initial_pressure_level = pressure_level
                
                mode = 0
                case_stor_wat = 0
                case_stor_h2 = 0
                case_el = 0

               
            H2_production_ev[i] = H2_production_tot
            p_pem_cons_ev[i] = p_pem_cons_tot
            surplus_tot_ev[i] = surplus_tot
            water_consumption_ev[i] = water_consumption_tot
            n_elec_on_ev[i] = n_elec_on_tot
            status_diss_ev[i] = status_diss
            consumption_diss_ev[i] = consumption_diss_tot
            water_produced_ev[i] = water_produced
            soc_wat_ev[i] = soc_wat_tot
            water_surplus_ev[i] = water_surplus_tot
            water_deficit_ev[i] = water_deficit_tot 
            delta_w_ev[i] = delta_w_tot
            power_compr_ev[i] = power_el_abs
            q_compr_ev[i] = q_power_developed
            temp_out_compr_ev[i] = temperature_gas_output
            p_h2_storage_ev[i] = pressure
            p_level_h2_storage_ev[i] = pressure_level
            v_h2_storage_ev[i]= v
            v_h2_storage_surplus_ev[i]=v_surplus
            v_h2_storage_unmet_ev[i]= v_unmet
            flow_in_storage_ev[i] = gas_flow_rate_input
            flow_out_storage_ev[i] = gas_flow_rate_output
            mode_ev[i] = mode
            case_stor_wat_ev[i] = case_stor_wat
            case_stor_h2_ev[i] = case_stor_h2
            case_el_ev[i] = case_el


        return {
            'production_out':production_out_ev,
            'PUN': el_price,
            'H2_production':H2_production_ev, 
            'p_pem_cons':p_pem_cons_ev,
            'surplus_tot':surplus_tot_ev,
            'water_consumption':water_consumption_ev,
            'n_elec_on':n_elec_on_ev,
            'status_diss': status_diss_ev,
            'consumption_diss': consumption_diss_ev,
            'wat_prod_diss': water_produced_ev,
            'soc_wat': soc_wat_ev,
            'water_surplus': water_surplus_ev,
            'water_deficit': water_deficit_ev,
            'delta_w': delta_w_ev,
            'power_compr': power_compr_ev*1e-3,
            'q_compr': q_compr_ev*1e-3,
            'temp_out_compr': temp_out_compr_ev,
            'pressure_gas_in_storage': p_h2_storage_ev, 
            'pressure_level':p_level_h2_storage_ev, 
            'volume_gas_stored': v_h2_storage_ev, 
            'v_surplus': v_h2_storage_surplus_ev, 
            'v_unmet':v_h2_storage_unmet_ev, 
            'gas_flow_rate_input_storage':flow_in_storage_ev, 
            'h2_stor_demand':flow_out_storage_ev,
            'mode':mode_ev,
            'case_stor_wat': case_stor_wat_ev,
            'case_stor_h2': case_stor_h2_ev,
            'case_el': case_el_ev,
        }
