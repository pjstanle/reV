import numpy as np
from PySAM import Windpower

class WindPlant():
    def __init__(self):

        # definitely need to define
        self.wind_turbine_powercurve_windspeeds = np.array([])
        self.wind_turbine_powercurve_powerout = np.array([])
        self.rotor_diameter = 0.0
        self.hub_height = 0.0
        self.turbine_x = np.array([0.0])
        self.turbine_y = np.array([0.0])
        self.wind_directions = np.array([])
        self.wind_frequencies = np.array([])
        self.wind_speeds = np.array([])

        # fine to remain with default
        # 0:Simple, 1:Park (WAsp) 2:Eddy Viscosity 3:Constant %
        self.wake_model = 1
        self.system_model = Windpower.default("WindPowerSingleOwner")

        # local variables
        self.n_turbines = 1
        self.plant_capacity = 0.0
        self.turbine_capacity = 0.0

        # outputs
        self.aep = 0.0
        self.capacity_factor = 0.0

    def initialize_plant(self):

        # turbine parameters
        self.system_model.Turbine.wind_turbine_powercurve_windspeeds = \
            self.wind_turbine_powercurve_windspeeds
        self.system_model.Turbine.wind_turbine_powercurve_powerout = \
            self.wind_turbine_powercurve_powerout
        self.system_model.Turbine.wind_turbine_rotor_diameter = \
            self.rotor_diameter
        self.system_model.Turbine.wind_turbine_hub_ht = \
            self.hub_height

        # plant parameters
        self.system_model.Farm.wind_farm_xCoordinates = self.turbine_x
        self.system_model.Farm.wind_farm_yCoordinates = self.turbine_y

        self.n_turbines = len(self.turbine_x)
        self.turbine_capacity = np.max(self.wind_turbine_powercurve_powerout)
        self.plant_capacity = self.n_turbines*self.turbine_capacity
        self.system_model.Farm.system_capacity = self.plant_capacity

        # wind resource
        self.system_model.Resource.wind_resource_model_choice = 2
        n_directions = len(self.wind_frequencies)
        wind_resource_distribution = [] #speed (m/s), dir (deg), freq
        for k in range(n_directions):
            wind_resource_distribution.append([self.wind_speeds[k], \
                self.wind_directions[k], self.wind_frequencies[k]])
        self.system_model.Resource.wind_resource_distribution = \
            wind_resource_distribution
        self.system_model.Farm.wind_farm_wake_model = self.wake_model

    def define_turbine_locations(self):
        # plant parameters
        self.system_model.Farm.wind_farm_xCoordinates = self.turbine_x
        self.system_model.Farm.wind_farm_yCoordinates = self.turbine_y

        self.n_turbines = len(self.turbine_x)
        self.turbine_capacity = np.max(self.wind_turbine_powercurve_powerout)
        self.plant_capacity = self.n_turbines * self.turbine_capacity
        self.system_model.Farm.system_capacity = self.plant_capacity

    def evaluate_aep(self):

        self.system_model.execute()
        self.aep = self.system_model.Outputs.annual_energy
        self.capacity_factor = self.system_model.Outputs.capacity_factor
