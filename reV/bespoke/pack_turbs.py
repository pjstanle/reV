# import libraries
from json import load
from plotting_functions import plot_poly, plot_turbines, get_xy
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
# SCIPY functions
from scipy.optimize import minimize
from layout_functions import calc_spacing, check_exclusions, load_exclusions

from shapely.geometry import Polygon, MultiPolygon, Point, LineString
import time

class PackTurbines():

    def __init__(self):

        self.min_spacing = 0.0 # minimum spacing
        self.max_x = 0.0
        self.min_x = 0.0
        self.max_y = 0.0
        self.min_y = 0.0
        self.obj_scale = 1.0
        self.weight_x = 1.0

        # turbine locations
        self.turbine_x = np.array([])
        self.turbine_y = np.array([])

        # exclusions data
        # self.exclusions_x = np.array([])
        # self.exclusions_y = np.array([])
        # self.exclusions_data = np.array([])
        # self.dx = 0.0
        # self.dy = 0.0
        self.safe_polygons = None

    def plot_turbines(self,label=None,c='ro',fill=True):

        nTurbs = len(self.turbine_x)
        for i in range(nTurbs):
            x = np.linspace(self.turbine_x[i] - self.min_spacing, \
                self.turbine_x[i] + self.min_spacing, 1000)
            y1 = np.sqrt(self.min_spacing ** 2 - (x - self.turbine_x[i]) \
                ** 2) + self.turbine_y[i]
            y2 = -np.sqrt(self.min_spacing ** 2 - (x - self.turbine_x[i]) \
                ** 2) + self.turbine_y[i]
            plt.plot(x, y1, 'k--')
            plt.plot(x, y2, 'k--')
            if fill:
                plt.fill_between(x, y1, y2, color='k')
            if i == 0:
                plt.plot(self.turbine_x[i], self.turbine_y[i], c, markersize=5,label=label)
            else:
                plt.plot(self.turbine_x[i], self.turbine_y[i], c, markersize=5)

    def obj_func_spacing(self, x):

        new_x, new_y = x

        # TODO change minimize x+y to include some angle

        if len(self.turbine_x) < 1:
            return (new_x+new_y)/self.obj_scale
        else:
            temp_x = np.append(self.turbine_x,new_x)
            temp_y = np.append(self.turbine_y,new_y)
            spacing = np.min(calc_spacing(temp_x,temp_y))
            if np.min(spacing) < self.min_spacing:
                return 1E12
            else:
                return (new_x+new_y)/self.obj_scale

    def add_random_turbine(self):
        # add random feasible turbine
        count = 0
        dist = 0
        boundary = -1E12
        while (dist < self.min_spacing or (boundary < 0)) and (count < 1000) :
            x = (self.max_x-self.min_x)*np.random.rand(1) + self.min_x
            y = (self.max_y-self.min_y)*np.random.rand(1) + self.min_y

            if len(self.turbine_x) > 0:
                dist = np.min(np.sqrt( (x-self.turbine_x)**2 + (y-self.turbine_y)**2 ))
            else:
                dist = 2*self.min_spacing
            # excl = check_exclusions(x,y,self.exclusions_x,self.exclusions_y,
            #         self.exclusions_data,self.dx,self.dy)
            boundary = self.boundary_constraint((x,y))
            count = count + 1

        if count == 1000:
            return False
        else:
            return (x, y)

    def pack_turbines(self):
        # add bounds
        bnds = [(self.min_x,self.max_x),(self.min_y,self.max_y)]

        can_add_more = True

        while can_add_more == True:
            x0 = self.add_random_turbine()
            if x0 == False:
                can_add_more = False
            else:
                if len(self.turbine_x) > 0:
                    cons = ({'type': 'ineq', 'fun': lambda x: np.min(np.sqrt( (x[0] - self.turbine_x)**2 + (x[1] - self.turbine_y)**2)) - self.min_spacing},
                            {'type': 'ineq', 'fun': self.boundary_constraint})
                    res = minimize(self.obj_func_spacing,x0,method='SLSQP',bounds=bnds,constraints=cons)
                else:
                    cons = ({'type': 'ineq', 'fun': self.boundary_constraint})
                    res = minimize(self.obj_func_spacing,x0,method='SLSQP',bounds=bnds,constraints=cons)


                new_x, new_y = res.x
                if self.boundary_constraint((new_x,new_y)) > -1:
                    self.turbine_x = np.append(self.turbine_x,new_x)
                    self.turbine_y = np.append(self.turbine_y,new_y)

                # print(len(self.turbine_x))


        # TODO not necessarily want to add to min y point
        print("before: ", len(self.turbine_x))
        s = time.time()
        can_add_more = True
        while can_add_more:
            leftover = MultiPolygon(self.safe_polygons)
            for k in range(len(self.turbine_x)):
                turbine = Point(self.turbine_x[k],self.turbine_y[k]).buffer(min_spacing)
                leftover = leftover.difference(turbine)
            if type(leftover) == Polygon:
                leftover = MultiPolygon([leftover])
            nareas = len(leftover)
            if nareas > 0:
                areas = np.zeros(len(leftover))
                for i in range(nareas):
                    areas[i] = leftover[i].area
                smallest_area = leftover[np.argmin(areas)]
                exterior_coords = smallest_area.exterior.coords[:]
                x,y = get_xy(exterior_coords)
                self.turbine_x = np.append(self.turbine_x,x[np.argmin(x)])
                self.turbine_y = np.append(self.turbine_y,y[np.argmin(x)])
                print("added turbine")
            else:
                can_add_more = False
        print("after: ", len(self.turbine_x))


    def pack_turbines_poly(self):

        # TODO not necessarily want to add to min y point
        can_add_more = True
        leftover = MultiPolygon(self.safe_polygons)
        while can_add_more:
            if len(self.turbine_x) > 0:
                leftover = leftover.difference(new_turbine)
            if type(leftover) == Polygon:
                leftover = MultiPolygon([leftover])
            nareas = len(leftover)
            if nareas > 0:
                areas = np.zeros(len(leftover))
                for i in range(nareas):
                    areas[i] = leftover[i].area
                smallest_area = leftover[np.argmin(areas)]
                exterior_coords = smallest_area.exterior.coords[:]
                x,y = get_xy(exterior_coords)
                metric = self.weight_x*x+y
                self.turbine_x = np.append(self.turbine_x,x[np.argmin(metric)])
                self.turbine_y = np.append(self.turbine_y,y[np.argmin(metric)])
                new_turbine = Point(x[np.argmin(metric)],y[np.argmin(metric)]).buffer(self.min_spacing)
            else:
                can_add_more = False


    def boundary_constraint(self,x):
        point = Point(x)
        min_dist = -1E12
        for i in range(len(self.safe_polygons)):
            coords = self.safe_polygons[i].exterior.coords
            line = LineString(coords)
            dist = line.distance(point)
            if self.safe_polygons[i].contains(point) == False:
                dist = dist*-1
            if dist > min_dist:
                min_dist = dist

        return min_dist


    def clear(self):
        self.turbine_x = np.array([])
        self.turbine_y = np.array([])


if __name__=="__main__":

    import warnings
    warnings.filterwarnings("ignore")

    safe_polygons = load_exclusions(polygons=True)
    packing = PackTurbines()


    rotor_diameter = 150
    min_spacing = 4*rotor_diameter
    minx, miny, maxx, maxy = safe_polygons.bounds

    packing.min_spacing = min_spacing # minimum spacing
    packing.max_x = maxx
    packing.min_x = minx
    packing.max_y = maxy
    packing.min_y = miny

    packing.obj_scale = 0.1

    # packing.exclusions_x = exclusions_x
    # packing.exclusions_y = exclusions_y
    # packing.exclusions_data = exclusions_data
    # packing.dx = dx
    # packing.dy = dy

    packing.safe_polygons = safe_polygons

    import time
    start = time.time()
    # packing.pack_turbines()
    packing.pack_turbines_poly()
    print("time to run: ", time.time()-start)

    x = packing.turbine_x
    y = packing.turbine_y
    print("nturbs: ", len(x))
    # print(x)
    # print(y)

    # plt.figure(1)
    # I,J = np.shape(exclusions_x)
    # for i in range(I):
    #     for j in range(J):
    #         if exclusions_data[i,j] == True:
    #             color="blue"
    #         else:
    #             color="red"
    #         plt.plot(exclusions_x[i,j],exclusions_y[i,j],"o",markersize=5,color=color)
    # plt.axis("equal")

    # plt.plot(packing.turbine_x,packing.turbine_y,"o")

    plt.figure(1)
    packing.plot_turbines()
    plot_poly(safe_polygons,ax=plt.gca(),linewidth=2,color="blue")
    plt.axis("equal")

    plt.figure(2)
    plot_poly(safe_polygons,ax=plt.gca())
    plt.axis("equal")

    plot_turbines(x,y,rotor_diameter/2.0,"C0")

    plt.show()

    # plt.show()
