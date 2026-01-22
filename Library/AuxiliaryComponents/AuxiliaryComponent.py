class AuxiliaryComponent:
    def __init__(self, id, tech, replacement_year, replacement_cost, sp_capex, sp_opex, capacity):
        """
        :param id: identification code
        :param tech:'waterstorage', 'hydrogenstorage', 'desalinator', 'compressor'
        :param replacement_year
        :param replacement_cost: replacement cost [€/capacity]
        :param sp_capex: specific capital cost [€/capacity]
        :param sp_opex: specific operational cost [€/capacity/year]
        :param capacity: 
        - water_storage m3
        - hydrogen_storage kgH2
        - desalinator m3 
        - compressor kW
        """
        self.id = id
        self.tech = tech
        self.replacement_year = replacement_year
        self.replacement_cost = replacement_cost
        self.sp_capex = sp_capex
        self.sp_opex = sp_opex
        self.capacity = capacity

    def compute_economics(self):
        '''
        OUTPUT: 
        - capex: capital expenditure [€]
        - opex: operational expenditure [€/year]
        - repl_cost: replacement cost [€]
        '''

        capex=self.sp_capex*self.capacity # €
        opex=self.sp_opex*self.capacity # €/year 
        repl_cost=self.replacement_cost*self.capacity # €
        
        return capex, opex, repl_cost