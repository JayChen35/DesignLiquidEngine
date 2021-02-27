# Classes for PropSimPython
# Project Caelus, Aphlex 1C Engine
# Jason Chen, 14 February, 2021


class Struct():
    """
    Generic data class to store named variables. Can be initialized recursively using an input dictionary
    to work with nested dictionary structures. Used to easily call variables with the "." operator.
    """
    def __init__(self, *args, **kwargs):
        super().__init__()
        if len(args) == 1 and type(args[0]) == dict: # Create from a dictionary
            for key, value in args[0].items():
                if isinstance(value, dict): # If the value is a dictionary (more nests to deal with)
                    setattr(self, key, Struct(value)) # Recursively create another Struct instance
                else:
                    setattr(self, key, value) # If the value is regular data, set the key, value pair
    
    def __getitem__(self, key): # Defines the "." operator to access a value in the internal dictionary.
        return self.__dict__[key]

    def __repr__(self): # __repr__ allows printing of this class to be more readible.
        return "{%s}" % str(", ".join("%s : %s" % (k, repr(v)) for (k, v) in self.__dict__.items()))
    
    def fieldnames(self):
        return [key for key in self.__dict__]  


class Gas():
    """
    Contains properties of a gas given a string input (currently only nitrogen is supported).
    """
    def __init__(self, type: str):
        super().__init__()
        if type == "Nitrogen":
            self.c_v = 0.743e3 # Specific volume (J/kg*K)
            self.mol_mass = 2*14.0067e-3 # Molecular mass (kg/mol)
        else:
            raise ValueError("Invalid gas type. Currently, only nitrogen is supported. Use \"Gas(Nitrogen)\".")
        self.R_spec = get_r_spec()
        self.c_p = get_c_p()
        self.gamma = get_gamma()

    def get_r_spec(self):
        R_u = 8.3144598 # J/K*mol
        return R_u/self.mol_mass

    def get_c_p(self):
        return self.c_v + self.R_spec

    def get_gamma(self):
        return self.c_p/self.c_v

    def to_dict(self):
        output = {
            "c_v": self.c_v,
            "mol_mass": self.mol_mass,
            "R_spec" : self.R_spec,
            "c_p": self.c_p,
            "gamma": self.gamma
        }
        return output

