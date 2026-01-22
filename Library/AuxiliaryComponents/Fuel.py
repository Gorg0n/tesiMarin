
class Fuel:
    def __init__(self, fuel, state=None, lhv=None, rho=None, molecular_weight=None, c_p=None, tC02_emission=None,
                 cost=None):
        """
        :param fuel: air, LPG (liquefied gas) , propane (C3H8), methane (CH4), natural gas, hydrogen (H2), syngas, biogas, wood pellets, wood chips
        :param state:'solid' or 'liquid' or 'gas'
        :param lhv: lower heating value [kJ/kg]
        :param rho: density (rho) for liquid [kg/l] at 25°C and 1 atm, for gas [kg/Nm3] 
        :param molecular_weight [g/mol]
        :param c_p= pressure specific heat [J/(K*mol)] at 25°C 
        :param tC02_emission: tCO2_emission [tonC02/ton]
        :param cost: [€/ton]] 
        """ 

        self.library = {'LGP': {'lhv': 46100, 'state': 'gas', 'rho': 2.250, 'tCO2_emission': 100,
                                'molecular weight': 100, 'c_p': 100, 'cost': 100},
                        'diesel': {'lhv': 42500, 'state': 'liquid', 'rho': 0.835, 'tCO2_emission': 100,
                                   'molecular weight': 100, 'c_p': 100, 'cost': 100},
                        'propane': {'lhv': 46300, 'state': 'gas', 'rho': 2, 'tCO2_emission': 100,
                                    'molecular weight': 100, 'c_p': 100, 'cost': 100},
                        'methane': {'lhv': 50000, 'state': 'gas', 'rho': 0.720, 'tCO2_emission': 100,
                                    'molecular weight': 16.04, 'c_p': 35.8, 'cost': 0},
                        'natural gas': {'lhv': 100, 'state': 'gas', 'rho': 100, 'tCO2_emission': 100,
                                        'molecular weight': 100, 'c_p': 100, 'cost': 100},
                        'hydrogen': {'lhv': 120000, 'state': 'gas', 'rho': 0.08988, 'tCO2_emission': 100,
                                     'molecular weight': 2.016, 'c_p': 28.85, 'cost': 100,
                                     'a': 29.11, 'b': 0.1916*1e-2, 'c': 0.4003*1e-5, 'd': -0.8704*1e-9}, 
                        'syngas': {'lhv': 5580, 'state': 'gas', 'rho': 1.1, 'tCO2_emission': 100,
                                   'molecular weight': 20.38, 'c_p': 31, 'cost': 100},
                        'wood chips': {'lhv':17640 , 'state': 'solid', 'rho': 100, 'tCO2_emission': 100,
                                       'molecular weight': 100, 'c_p': 100, 'cost': 80},
                        'wood pellets': {'lhv': 17640, 'state': 'solid', 'rho': 100, 'tCO2_emission': 100,
                                         'molecular weight': 40, 'c_p': 100, 'cost':80}
                        }
        self.fuel=fuel
        if fuel not in self.library.keys():
            self.state = state
            self.lhv = lhv
            self.rho = rho
            self.tCO2_emission = tC02_emission
            self.molecular_weight = molecular_weight
            self.c_p = c_p
            self.cost = cost
        else:
            self.state = self.library[fuel]['state']
            self.lhv = self.library[fuel]['lhv']
            self.rho = self.library[fuel]['rho']
            self.tCO2_emission = self.library[fuel]['tCO2_emission']
            self.molecular_weight = self.library[fuel]['molecular weight']
            self.c_p = self.library[fuel]['c_p']
            self.cost = self.library[fuel]['cost']
            self.a = self.library[fuel]['a']
            self.b = self.library[fuel]['b']
            self.c = self.library[fuel]['c']
            self.d = self.library[fuel]['d']


