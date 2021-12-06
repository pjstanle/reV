import numpy as np
import matplotlib.pyplot as plt
import h5py
from numpy.core.fromnumeric import shape
from reV.supply_curve.exclusions import ExclusionMaskFromDict
from scipy.spatial import ConvexHull

def load_exclusions(filename="ri_exclusions.h5",polygons=False):

    excl_dict = {'ri_smod': {'include_values': [1, ], 'weight': 1,
                             'exclude_nodata': True}}

    with ExclusionMaskFromDict(filename, layers_dict=excl_dict) as f:
        truth = f.mask

    with h5py.File(filename, "r") as f:
        latitude = list(f["latitude"])
        longitude = list(f["longitude"])

    N = np.shape(latitude)[0]
    latitude = [[float(y) for y in x] for x in latitude]
    longitude = [[float(y) for y in x] for x in longitude]

    lat = np.mean(latitude)
    tx = [[y * 40075000 * np.cos(np.deg2rad(lat)) / 360 for y in x] \
        for x in longitude]
    ty = [[y * 111.32 * 1000 for y in x] for x in latitude]

    tx = np.array(tx)
    ty = np.array(ty)
    truth = np.array(truth)

    ncells = 128
    s1 = 500
    s2 = 500

    x = tx[s1:s1 + ncells, s2:s2 + ncells]
    y = ty[s1:s1 + ncells, s2:s2 + ncells]
    safe = truth[s1:s1 + ncells, s2:s2 + ncells]

    if polygons == True:
        from shapely.geometry import Polygon, MultiPolygon, Point
        import shapely.affinity
        import rasterio.features

        shapes = rasterio.features.shapes(safe)
        polygons = [Polygon(shape[0]["coordinates"][0]) for shape in shapes \
            if shape[1] == 1]

        for i in range(len(polygons)):
            polygons[i] = shapely.affinity.scale(polygons[i], xfact=90.0, \
                yfact=-90.0, origin=(0,0))

        safe_polygons = MultiPolygon(polygons)

        return safe_polygons, x, y, safe

    else:
        return x, y, safe

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def check_exclusions(x, y, exclusions_x, exclusions_y, \
    exclusions_data, dx, dy):

    if x > np.max(exclusions_x)+dx or x < np.min(exclusions_x)-dx \
        or y > np.max(exclusions_y)+dy or y < np.min(exclusions_y)-dy:
        return False
    else:
        xidx = find_nearest(exclusions_x[0, :], x)
        yidx = find_nearest(exclusions_y[:, 0], y)
        return exclusions_data[yidx, xidx]

def get_grid_locs(x_spacing, y_spacing, shear, rotation, center_x, \
    center_y, exclusions_x, exclusions_y, exclusions_data):
    # create grid
    width = np.max(exclusions_x)-np.min(exclusions_x)
    height = np.max(exclusions_y)-np.min(exclusions_y)
    nrows = int(np.max([width, height])/np.min([x_spacing, y_spacing])) * 2 + 1
    ncols = nrows
    xlocs = np.arange(0, ncols) * x_spacing
    ylocs = np.arange(0, nrows) * y_spacing
    row_number = np.arange(0, nrows)

    d = np.array([i for x in xlocs for i in row_number])
    layout_x = np.array([x for x in xlocs for y in ylocs]) + \
        d * y_spacing*np.tan(shear)
    layout_y = np.array([y for x in xlocs for y in ylocs])

    # rotate
    rotate_x = np.cos(rotation) * layout_x - np.sin(rotation) \
        * layout_y
    rotate_y = np.sin(rotation) * layout_x + np.cos(rotation) \
        * layout_y
    # move center of grid
    rotate_x = (rotate_x - np.mean(rotate_x)) + center_x
    rotate_y = (rotate_y - np.mean(rotate_y)) + center_y

    # get rid of points outside of boundary and violate setback constraints
    meets_constraints = np.zeros(len(rotate_x),dtype=bool)
    dx = (exclusions_x[0, 1] - exclusions_x[0, 0]) / 2.0
    dy = (exclusions_y[0, 0] - exclusions_y[1, 0]) / 2.0
    for i in range(len(rotate_x)):
        meets_constraints[i] = check_exclusions(rotate_x[i], \
            rotate_y[i], exclusions_x, exclusions_y, exclusions_data, \
            dx, dy)
    # print("constraints time: ", time.time()-start_constraints)

    # arrange final x,y points
    return_x = rotate_x[meets_constraints]
    return_y = rotate_y[meets_constraints]

    return return_x, return_y

def calc_spacing(turbine_x, turbine_y):
    #calculate the spacing between each turbine and every other turbine (without repeating)
    nturbs = len(turbine_x)
    npairs = int((nturbs * (nturbs - 1)) / 2)
    spacing = np.zeros(npairs)

    ind = 0
    for i in range(nturbs):
        for j in range(i, nturbs):
            if i != j:
                spacing[ind] = np.sqrt((turbine_x[i] - turbine_x[j])**2 \
                    + (turbine_y[i] - turbine_y[j])**2)
                ind += 1

    return spacing

def arbitrary_boundary(turbineX, turbineY, boundaryVertices, boundaryNormals):

    if type(turbineX) == np.float64:
        nTurbines = 1
    else:
        nTurbines = len(turbineX)
    locations = np.zeros([nTurbines, 2])
    if nTurbines == 1:
        locations[0] = np.array([turbineX, turbineY])
    else:
        for i in range(0, nTurbines):
            locations[i] = np.array([turbineX[i], turbineY[i]])

    # calculate distance from each point to each face
    boundaryDistances = calculate_distance(locations, boundaryVertices, \
        boundaryNormals)

    return boundaryDistances

def calculate_boundary(vertices):

    # find the points that actually comprise a convex hull
    hull = ConvexHull(list(vertices))

    # keep only vertices that actually comprise a convex hull and arrange in CCW order
    vertices = vertices[hull.vertices]

    # get the real number of vertices
    nVertices = vertices.shape[0]

    # initialize normals array
    unit_normals = np.zeros([nVertices, 2])

    # determine if point is inside or outside of each face, and distance from each face
    for j in range(0, nVertices):

        # calculate the unit normal vector of the current face (taking points CCW)
        if j < nVertices - 1:  # all but the set of point that close the shape
            normal = np.array([vertices[j + 1, 1]-vertices[j, 1],
                               -(vertices[j + 1, 0]-vertices[j, 0])])
            unit_normals[j] = normal/np.linalg.norm(normal)
        else:   # the set of points that close the shape
            normal = np.array([vertices[0, 1]-vertices[j, 1],
                               -(vertices[0, 0]-vertices[j, 0])])
            unit_normals[j] = normal/np.linalg.norm(normal)

    return vertices, unit_normals

def calculate_distance(points, vertices, unit_normals, return_bool=False):

    """
    :param points: points that you want to calculate the distance from to the faces of the convex hull
    :param vertices: vertices of the convex hull CCW in order s.t. vertices[i] -> first point of face for
           unit_normals[i]
    :param unit_normals: unit normal vector for each face CCW where vertices[i] is first point of face
    :param return_bool: set to True to return an array of bools where True means the corresponding point
           is inside the hull
    :return face_distace: signed perpendicular distance from each point to each face; + is inside
    :return [inside]: (optional) an array of zeros and ones where 1.0 means the corresponding point is inside the hull
    """

    # print points.shape, vertices.shape, unit_normals.shape

    nPoints = points.shape[0]
    nVertices = vertices.shape[0]

    # initialize array to hold distances from each point to each face
    face_distance = np.zeros([nPoints, nVertices])

    if not return_bool:
        # loop through points and find distance to each face
        for i in range(0, nPoints):

            # determine if point is inside or outside of each face,
            # and distance from each face
            for j in range(0, nVertices):

                # define the vector from the point of interest to
                # the first point of the face
                pa = np.array([vertices[j, 0] - points[i, 0], \
                    vertices[j, 1] - points[i, 1]])

                # find perpendicular distance from point to current
                # surface (vector projection)
                d_vec = np.vdot(pa, unit_normals[j]) * unit_normals[j]

                # calculate the sign of perpendicular distance from
                # point to current face (+ is inside, - is outside)
                face_distance[i, j] = np.vdot(d_vec, unit_normals[j])

        return face_distance

    else:
        # initialize array to hold boolean indicating whether a
        #  point is inside the hull or not
        inside = np.zeros(nPoints)

        # loop through points and find distance to each face
        for i in range(0, nPoints):

            # determine if point is inside or outside of each
            # face, and distance from each face
            for j in range(0, nVertices):

                # define the vector from the point of interest
                # to the first point of the face
                pa = np.array([vertices[j, 0] - points[i, 0], vertices[j, 1] \
                    - points[i, 1]])

                # find perpendicular distance from point to current
                # surface (vector projection)
                d_vec = np.vdot(pa, unit_normals[j]) * unit_normals[j]

                # calculate the sign of perpendicular distance from
                # point to current face (+ is inside, - is outside)
                face_distance[i, j] = np.vdot(d_vec, unit_normals[j])

            # check if the point is inside the convex hull by checking
            #  the sign of the distance
            if np.all(face_distance[i] >= 0):
                inside[i] = 1.0

        return face_distance, inside

if __name__=="__main__":
    exclusions_x,exclusions_y,exclusions_data = load_exclusions()
    exclusions_x = exclusions_x - np.mean(exclusions_x)
    exclusions_y = exclusions_y - np.mean(exclusions_y)

    plt.figure(1)
    I,J = np.shape(exclusions_x)
    for i in range(I):
        for j in range(J):
            if exclusions_data[i,j] == True:
                color="blue"
            else:
                color="red"
            plt.plot(exclusions_x[i,j],exclusions_y[i,j],"o",markersize=5,color=color)
    plt.axis("equal")

    x_spacing = 1000
    y_spacing = 2000
    shear = np.pi/4
    rotation = np.pi/4
    center_x = 0.0
    center_y = 0.0
    x,y = get_grid_locs(x_spacing,y_spacing,shear,rotation,center_x,center_y,exclusions_x,exclusions_y,exclusions_data)

    plt.figure(2)
    plt.plot(x,y,"o",markersize=5)
    plt.axis("equal")
    plt.show()

    # N = np.arange(10)*100+10
    # time_array1 = np.zeros(len(N))
    # time_array2 = np.zeros(len(N))
    # for i in range(len(N)):
    #     x = np.random.rand(N[i])
    #     y = np.random.rand(N[i])
    #     start_time = time.time()
    #     s1 = calc_spacing(x,y)
    #     time_array1[i] = time.time()-start_time

    #     start_time = time.time()
    #     s2 = calc_spacing2(x,y)
    #     time_array2[i] = time.time()-start_time

    #     # print(np.sum(s1-s2))

    # print(time_array1)
    # print(time_array2)
    # plt.plot(N,time_array1,label="1")
    # plt.plot(N,time_array2,label="2")
    # plt.legend()
    # plt.show()
