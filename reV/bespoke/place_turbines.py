# TODO need to fix the load exlusions function
# TODO need to figure out how to pass in wind data

import numpy as np
from evaluate_plant import WindPlant
from pack_turbs import PackTurbines
from layout_functions import load_exclusions
from shapely.geometry import Polygon
from gradient_free import GeneticAlgorithm

class PlaceTurbines():


    def __init__(self):

        # need to be assigned
        self.capex_function = None
        self.om_function = None
        self.wind_turbine_powercurve_windspeeds = None
        self.wind_turbine_powercurve_powerout = None
        self.rotor_diameter = 100.0
        self.hub_height = 100.0
        self.wind_directions = np.array([])
        self.wind_speeds = np.array([])
        self.wind_frequencies = np.array([])
        self.exclusions_filename = None
        self.min_spacing = 0.0
        self.ga_time = 20.0
        self.fcr = 0.063
        self.max_coe_multiplier = 1.05

        # internal variables
        self.plant = None
        self.packing = None
        self.x_locations = np.array([])
        self.y_locations = np.array([])
        self.safe_polygons = None
        self.max_coe = 1E6

        # outputs
        self.turbine_x = np.array([])
        self.turbine_y = np.array([])
        self.nturbs = 0
        self.capacity = 0.0
        self.area = 0.0
        self.capacity_density = 0.0
        self.coe = 0.0
        self.aep = 0.0

        self.turbine_x_coe = np.array([])
        self.turbine_y_coe = np.array([])
        self.nturbs_coe = 0
        self.capacity_coe = 0.0
        self.area_coe = 0.0
        self.capacity_density_coe = 0.0
        self.coe_coe = 0.0
        self.aep_coe = 0.0


    def initialize_plant(self):

        self.plant = WindPlant()
        self.plant.wind_turbine_powercurve_windspeeds = self.wind_turbine_powercurve_windspeeds
        self.plant.wind_turbine_powercurve_powerout = self.wind_turbine_powercurve_powerout
        self.plant.rotor_diameter = self.rotor_diameter
        self.plant.hub_height = self.hub_height
        self.plant.turbine_x = np.array([0])
        self.plant.turbine_y = np.array([0])
        self.plant.wind_directions = self.wind_directions
        self.plant.wind_speeds = self.wind_speeds
        self.plant.wind_frequencies = self.wind_frequencies
        self.plant.initialize_plant()


    def define_exclusions(self):

        safe_polygons,exclusions_x,exclusions_y,_ = load_exclusions(self.exclusions_filename,polygons=True)
        # add extra setback to cell boundary
        minx = np.min(exclusions_x)-45.0
        maxx = np.max(exclusions_x)+45.0
        miny = np.min(exclusions_y)-45.0
        maxy = np.max(exclusions_y)+45.0
        maxx = maxx-minx - self.min_spacing/2.0
        minx = minx-minx + self.min_spacing/2.0
        miny = miny-maxy + self.min_spacing/2.0
        maxy = maxy-maxy - self.min_spacing/2.0

        boundary_poly = Polygon(((minx,miny),(minx,maxy),(maxx,maxy),(maxx,miny)))
        self.safe_polygons = boundary_poly.intersection(safe_polygons)


    def initialize_packing(self):

        self.packing = PackTurbines()
        nturbs = 1E6
        self.packing.safe_polygons = self.safe_polygons
        mult = 1.0
        while nturbs > 300:
            self.packing.clear()
            self.packing.min_spacing = self.min_spacing*mult
            self.packing.pack_turbines_poly()
            nturbs = len(self.packing.turbine_x)
            mult *= 1.1
        self.x_locations = self.packing.turbine_x
        self.y_locations = self.packing.turbine_y


    def optimize_min_coe(self):
        nlocs = len(self.x_locations)
        ga_coe = GeneticAlgorithm()
        ga_coe.bits = np.ones(nlocs,dtype=int)
        bounds = np.zeros((nlocs,2),dtype=int)
        bounds[:,1] = 2
        variable_type = np.array([])
        for i in range(nlocs):
            variable_type = np.append(variable_type,"int")

        ga_coe.bounds = bounds
        ga_coe.variable_type = variable_type

        ga_coe.objective_function = self.objective_coe
        ga_coe.max_generation = 10000
        ga_coe.population_size = 25
        ga_coe.crossover_rate = 0.2
        ga_coe.mutation_rate = 0.01
        ga_coe.tol = 1E-6
        ga_coe.convergence_iters = 10000
        ga_coe.max_time = self.ga_time
        ga_coe.optimize_ga(print_progress=False)

        self.max_coe = ga_coe.optimized_function_value*self.max_coe_multiplier


        optimized_design_variables = ga_coe.optimized_design_variables
        optimized_design_variables = [bool(y) for y in optimized_design_variables]

        self.turbine_x_coe = self.x_locations[optimized_design_variables]
        self.turbine_y_coe = self.y_locations[optimized_design_variables]
        self.nturbs_coe = np.sum(optimized_design_variables)

        self.plant.turbine_x = self.x_locations[optimized_design_variables]
        self.plant.turbine_y = self.y_locations[optimized_design_variables]
        self.plant.define_turbine_locations()
        self.plant.evaluate_aep()
        self.capacity_coe = self.plant.plant_capacity
        self.aep_coe = self.plant.aep
        annual_cost = self.fcr*self.capex_function(self.capacity_coe) + self.om_function(self.capacity_coe)

        self.coe_coe = annual_cost/(self.aep_coe/1000) #$/MWh
        self.area_coe = self.safe_polygons.area
        self.capacity_density_coe = self.capacity_coe/self.area_coe*1E6


    def optimize_max_aep(self):
        nlocs = len(self.x_locations)
        ga_aep = GeneticAlgorithm()
        ga_aep.bits = np.ones(nlocs,dtype=int)
        bounds = np.zeros((nlocs,2),dtype=int)
        bounds[:,1] = 2
        variable_type = np.array([])
        for i in range(nlocs):
            variable_type = np.append(variable_type,"int")

        ga_aep.bounds = bounds
        ga_aep.variable_type = variable_type

        ga_aep.objective_function = self.objective_aep
        ga_aep.max_generation = 10000
        ga_aep.population_size = 25
        ga_aep.crossover_rate = 0.2
        ga_aep.mutation_rate = 0.01
        ga_aep.tol = 1E-6
        ga_aep.convergence_iters = 10000
        ga_aep.max_time = self.ga_time
        ga_aep.optimize_ga(print_progress=False)

        optimized_design_variables = ga_aep.optimized_design_variables
        optimized_design_variables = [bool(y) for y in optimized_design_variables]

        self.turbine_x = self.x_locations[optimized_design_variables]
        self.turbine_y = self.y_locations[optimized_design_variables]
        self.nturbs = np.sum(optimized_design_variables)

        self.plant.turbine_x = self.x_locations[optimized_design_variables]
        self.plant.turbine_y = self.y_locations[optimized_design_variables]
        self.plant.define_turbine_locations()
        self.plant.evaluate_aep()
        self.aep = self.plant.aep
        self.capacity = self.plant.plant_capacity
        annual_cost = self.fcr*self.capex_function(self.capacity) + self.om_function(self.capacity)

        self.coe = annual_cost/(self.aep/1000) #$/MWh
        self.area = self.safe_polygons.area
        self.capacity_density = self.capacity/self.area*1E6


    def place_turbines(self):

        self.initialize_plant()
        self.define_exclusions()
        self.initialize_packing()
        self.optimize_min_coe()
        self.optimize_max_aep()


    def objective_coe(self,x):

        if sum(x) == 0:
            return 1E6
        else:
            x = [bool(y) for y in x]
            self.plant.turbine_x = self.x_locations[x]
            self.plant.turbine_y = self.y_locations[x]
            self.plant.define_turbine_locations()
            self.plant.evaluate_aep()
            aep = self.plant.aep

            capacity = self.plant.plant_capacity
            annual_cost = self.fcr*self.capex_function(capacity) + self.om_function(capacity)

            coe = annual_cost/(aep/1000) #$/MWh
            return coe


    def objective_aep(self,x):

        x = [bool(y) for y in x]
        self.plant.turbine_x = self.x_locations[x]
        self.plant.turbine_y = self.y_locations[x]
        self.plant.define_turbine_locations()
        self.plant.evaluate_aep()
        aep = self.plant.aep

        capacity = self.plant.plant_capacity
        annual_cost = self.fcr*self.capex_function(capacity) + self.om_function(capacity)

        coe = annual_cost/(aep/1000) #$/MWh
        if coe > self.max_coe:
            adder = (coe-self.max_coe)*1E8
        else:
            adder = 0.0

        return -aep/1000.0 + adder
