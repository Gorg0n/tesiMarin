from Library.AuxiliaryComponents.AuxiliaryComponent import AuxiliaryComponent

class Desalinator(AuxiliaryComponent):
    def __init__(self, id, replacement_year, replacement_cost, sp_capex, sp_opex, capacity, consumption, nom_water_prod, tech='desalinator'):
        """
        :param id: identification code
        :param tech:string --> 'waterstorage', 'hydrogenstorage', 'desalinator', 'compressor'
        :param replacement_year:float
        :param replacement_cost:float
        :param consumption:float --> nominal power consumption [kW]
        :param nom_water_prod:float --> nominal water production [kg]
        """
        super().__init__(id=id, tech=tech, replacement_year=replacement_year, replacement_cost=replacement_cost, sp_capex=sp_capex, sp_opex=sp_opex, capacity=capacity)

        self.consumption = consumption
        self.nom_water_prod = nom_water_prod
    
    def compute_output(self, production):
        """
        :param production:float --> power output from renewable source [kW]
        OUTPUT:
            water_produced:float --> water produced [kg]
            status:float --> 1 if the desalinator is on, 0 if it is off
        """
        if production >= self.consumption:
            water_produced = self.nom_water_prod 
            consumption = self.consumption
            status = 1
        else: 
            water_produced = 0
            consumption = 0
            status = 0

        return water_produced, consumption, status
    

