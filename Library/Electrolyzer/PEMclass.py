import numpy as np
import pandas as pd
import math 

# Class of a PEM group (composto da più moduli in parallelo e in serie)
class PEM:
    def __init__(self, p_nom, nominal_prod_rate, V_module, n_parallel, n_series, phi_mem, sp_capex, sp_opex, replacement_year, replacement_cost, A=100, alpha_a=2, alpha_c=0.5, status_in=0):

        '''
        :param p_nom:float -> nominal power of the PEM module [kW]
        :param nominal_prod_rate:float -> nominal production rate of the electrolyzer module [kg/h]
        :param n_parallel --> number of modules simultaneously on or off (parallel)
        :param n_series --> number of modules in series
        :param phi_mem --> membrane thickness [cm]
        :param A:float -> cell active area [cm^2]
        :param V_module:float -> operating voltage of the module [V]
        :param alpha_a:float -> anodic transfer coefficient (usually 2)
        :param alpha_c:float -> cathodic transfer coefficient (usually 0.5)
        :param status_in -> 1:on; 0:off
        :param sp_capex: specific capital cost [€/kW]
        :param sp_opex: specific operational cost [€/kW/year]
        ''' 

        self.n_parallel = n_parallel
        self.n_series = n_series
        self.p_nom = p_nom*self.n_parallel*self.n_series 
        self.nominal_prod_rate = nominal_prod_rate*self.n_parallel*self.n_series
        self.A = A
        self.phi_mem = phi_mem
        self.V_module = V_module
        self.alpha_a = alpha_a
        self.alpha_c = alpha_c
        self.status = status_in 
        self.sp_capex = sp_capex
        self.sp_opex = sp_opex
        self.replacement_year = replacement_year
        self.replacement_cost = replacement_cost

        self.F = 96485                # Costante di Faraday in C/mol
        self.R = 8.314                # Costante dei gas in J/(mol*K)
        self.n = 2                    # Numero di elettroni scambiati
        self.max_power_input = self.p_nom            # [kW]                                  
        self.min_power_input = self.p_nom  * 0.1     # [kW] 

    def power_input_controller(self, production):  
        '''
        :input production:float -> power output from renewable source [kW]
        OUTPUT:
            power_input_PEMEL:float -> power input to the PEM stack [kW]
        '''
        if production > self.max_power_input:
            power_input_PEMEL = self.max_power_input
        elif production < self.min_power_input:
            power_input_PEMEL = 0
        else:
            power_input_PEMEL =  production
        return power_input_PEMEL
    
    def polarization_curve(self, T, i, p):
        '''
        :input T:float -> operative temperature [°C]
        :input i:float -> current density [A/cm2]
        :input p:float -> operative pressure [bar]
        OUTPUT:
            V_cell:float -> cell voltage at the specific current density [V]

        Non è una vera e propria curva di polarizzazione, ma è la tensione della cella in funzione della corrente che metto io.
        Per farla essere una curva di polarizzazione avrei dovuto mettere un array di correnti standard.
        '''
        R = 8.314                # Costante dei gas in J/(mol*K)
        n = 2                    # Numero di elettroni scambiati
        F = 96485                # Costante di Faraday in C/mol
        T = T + 273.15           # Temperatura in K
        P =  p * 100             # Pressione in kPa
        pH2 = 2/3 * P            # Pressione parziale H2 kPa LEGGE DI DALTON
        pO2 = 1/3 * P            # Pressione parziale O2 kPa
        A = 8.07131
        B = 1730.63
        C = 233.426
        pH2O = 5065/38 * 10**(A-(B/(C-273+T))) / 1000  # kPa (Brezak)

        if i == 0:
            V_cell = 0
        else:
            i_0a = 6.6*1e-4  # A/cm2
            i_0c = 5.5*1e-4 # A/cm2
            if i <= 0 or i_0a <= 0 or i_0c <= 0:
                V_act = 0  # O assegna un valore predefinito
            else:
                V_act = R*T/(n*F) * (1/self.alpha_a * np.log(i/i_0a) + 1/self.alpha_c * np.log(i/i_0c))
            l = 0.08533*T - 6.77632
            sigma_mem = (0.005139*l - 0.00326) * np.exp(1268*(1/303-1/T)) # S/cm
            # print('sigma_mem:', sigma_mem)
            # print('resistvity:',1/sigma_mem)
            V_ohm_mem = self.phi_mem*10 / sigma_mem * i
            # resistivity = 101.8 # cm/S
            resistivity = 1/sigma_mem  # cm/S= ohm*S
            R_ele = resistivity*self.phi_mem/self.A
            V_ohm_ele = R_ele * i
            V_ohm = V_ohm_mem + V_ohm_ele

            Vrev0 = 1.229 - 0.9*1e-3*(T-298)
            V_rev = Vrev0 + R*T/(n*F) * np.log(pH2*pO2**(1/2)/pH2O)

            V_cell = V_rev + V_act + V_ohm

        return V_cell

    def faraday_efficiency(self, i, p):
        '''
        :input i:float -> current density [A/cm2]
        :input p:float -> operating pressure [bar]
        OUTPUT:
            etaf:float -> Faraday efficiency at the specific current density
        '''
        P_anodo = p # pressione in ingresso [bar]
        P_H2 = p # Electrolyzer H2 outlet pressure [bar]
        P_O2 = 2 # Electrolyzer O2 outlet pressure [bar]
        epsilon_h2_fick = 4.7 * 1e-11 # mol/(cm*s*bar)
        epsilon_h2_darcy = 2*1e-11 # mol/(cm*s*bar)
        S_h2 = 0.7 * 1e-7 # mol/(cm3*bar)
        n_h2o_eo = 3.5 # molh2o/molH+
        C_h2o = 40 * 1e3 # mol/m3
        P_h2_anodo = P_anodo # bar 
        P_h2_catodo = P_H2 # bar
        delta_P_h2 = P_h2_catodo-P_h2_anodo
        j_h2_fick = epsilon_h2_fick * (delta_P_h2/self.phi_mem) 
        P_o2_anodo = P_O2 
        j_h2_darcy = epsilon_h2_darcy * ((P_h2_catodo-P_o2_anodo)/self.phi_mem)
        j_h2_drag = n_h2o_eo * (P_h2_catodo*S_h2)/C_h2o
        j_h2_loss = j_h2_fick + j_h2_darcy - j_h2_drag

        epsilon_o2_fick = 2 * 1e-11 # mol/(cm*s*bar)
        S_o2 = 0.8 *1e-7 # mol/(cm3*bar) 
        P_o2_catodo = 0 
        delta_P_o2 =  P_o2_catodo -P_o2_anodo
        j_o2_fick = epsilon_o2_fick * (delta_P_o2/self.phi_mem)
        j_o2_drag = n_h2o_eo * (P_o2_anodo*S_o2)/C_h2o
        j_O2_loss = j_o2_fick + j_o2_drag 

        if i == 0:
            etaF = 0
        else:
            etaF = 1-(2*self.F*(j_h2_loss+2*j_O2_loss)/i)

        return etaF
    
    def cells_efficiency(self, i, T, p):
        '''
        :input i:float -> current density [A/cm2]
        :input T:float -> operating temperature [°C]
        :input p:float -> operating pressure [bar]
        OUTPUT:
            cell_efficiency:floar -> efficiency of the cell for the specific power input
        '''
        eta_f = self.faraday_efficiency(i=i, p=p)
        V_cell = self.polarization_curve(T=T, i=i, p=p)
        T = T + 273.15
        Vth = 1.4045 + (1.5784*1e-4*T) + 3.8037*1e-7*T**2

        if V_cell != 0: 
            eta_v = Vth / V_cell
        else: 
            eta_v = 0
        cells_efficiency = eta_f * eta_v

        return cells_efficiency
        
    def energy_performance(self, production, time, p, T): 
        '''
        :input production:float -> power output from renewable source [kW]  
        :input time:float -> 1 if hourly analysis, 0.25 if quarterly analysis
        :input T:float/array -> operating temperature [°C]
        :input p:float -> operating pressure [bar]
        OUTPUT: 
            H2_production:float -> H2 rate of production for the specific input power [kg/h]
            elec_consumption:float -> electricity consumption in one year [kWh]
            energy_surplus:float -> energy surplus in one year (the energy that doesn't go to the electrolysers) [kWh]
            system_efficiency:float -> efficiency of the stack for the specific power input 
            water_consumption:float -> water consumption [kg/h]
        '''
        power_input = self.power_input_controller(production)                                      # kW 
        I_module = power_input*1e3 / (self.V_module*self.n_parallel*self.n_series)                 # A
        I_cell = I_module                                                                          # A HP: celle del singolo modulo disposte in serie 
        i_cell = I_cell / self.A                                                                   # A/cm2
        H2_prod_nom = self.nominal_prod_rate                                                       # kg/h 
        I_cell_nom = self.max_power_input*1e3 / (self.V_module*self.n_parallel*self.n_series)      # A
        i_cell_nom = I_cell_nom / self.A                                                           # A/cm2
        if power_input == 0: 
            n_cells = 0
        else:
            n_cells = H2_prod_nom / (self.n_parallel*self.n_series*I_cell_nom * self.faraday_efficiency(i=i_cell_nom, p=p) / (2*self.F) * 2.016/1000 * 3600)  # [cells]

        H2_production = self.n_parallel*self.n_series * n_cells * I_cell * self.faraday_efficiency(i=i_cell, p=p) / (2*self.F) * 2.016/1000 * 3600  # [kg/h]
        if power_input==0: 
            elec_consumption = 0
            system_efficiency = 0
            energy_surplus = production * time
            self.status = 0
        else: 
            elec_consumption = power_input * time
            # system_efficiency = H2_production * 141.86 / power_input * 1e3/3600
            system_efficiency = self.cells_efficiency(i=i_cell, T=T, p=p)
            energy_surplus = (production - power_input) * time
            self.status = 1

        status = self.status         

        water_consumption = H2_production * 15  # [kg/h] (Giampieri)

        return H2_production, elec_consumption, energy_surplus, system_efficiency, status, water_consumption
    
    def compute_economics(self):

        capex=self.sp_capex*self.p_nom              # €
        opex=self.sp_opex*self.p_nom                # €/year 
        repl_cost=self.replacement_cost*self.p_nom  # €

        return capex, opex, repl_cost