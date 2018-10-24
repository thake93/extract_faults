from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import os

def ba(num_arr, fname, aname, od, mn, Nx, Ny, res, absolute_min_strain, crop_x, crop_y):
    def remove_topo(i_x,topo,arr,arr2, min_value):
        #removes data above topography
        #division between brittle and ductile parts
        #by setting parts to min_value
        #returns the reduced array
        for i_y in range(np.size(arr,0)):
            if arr2[i_y][i_x] < 0.8 or np.isnan(arr[i_y][i_x]):
                arr[i_y][i_x] = min_value
            for pkt in topo:
                x, ymax_temp = pkt
                if x >= xi[i_x] and x < xi[i_x+1] and \
                yi[i_y] > ymax_temp:
                    arr[i_y][i_x] = min_value
        return arr
    if not os.path.isdir(od+mn):
        os.mkdir(od+mn)
    if not os.path.isdir(aname[:len(aname)-28]):
        os.mkdir(aname[:len(aname)-28])
    for i in range (num_arr):
        FILE_NAME = fname%i
        ARR_NAME = aname%i
        ######################
        ### Load pvtu data ###
        ######################
        # pvtu loading and gridding largely based on https://stackoverflow.com/questions/23138112/vtk-to-maplotlib-using-numpy
        # load a pvtu file as input
        print(FILE_NAME)
        reader = vtk.vtkXMLPUnstructuredGridReader()
        reader.SetFileName(FILE_NAME)
        reader.Update()

        # Get the coordinates of nodes in the mesh
        nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

        number_of_fields = reader.GetOutput().GetPointData().GetNumberOfArrays()
        for i in range(number_of_fields):
            if reader.GetOutput().GetPointData().GetArrayName(i) == "strain_rate":
                id_strain_rate = i
            elif reader.GetOutput().GetPointData().GetArrayName(i) == "plastic_yielding":
                id_plast_yield = i

        #The "strain_rate" field is scalar number 17 in this pvtu file
        #The "plastic_yielding" field is scalar number 21 in this pvtu file
        field_vtk_array = reader.GetOutput().GetPointData().GetArray(id_strain_rate)
        field_vtk_array_plast_yield = reader.GetOutput().GetPointData().GetArray(id_plast_yield)

        #Get the coordinates of the nodes and their temperatures
        nodes_numpy_array = vtk_to_numpy(nodes_vtk_array)
        field_numpy_array = vtk_to_numpy(field_vtk_array)
        field_numpy_array_plast_yield = vtk_to_numpy(field_vtk_array_plast_yield)

        x_y_crop = []
        z_strain_rate_crop = []
        z_plastic_yield_crop = []

        if crop_x == None:
            crop_x = (min(nodes_numpy_array[:,0]),max(nodes_numpy_array[:,0]))
            crop_y = (min(nodes_numpy_array[:,1]),max(nodes_numpy_array[:,1]))

        for j,elem in enumerate(nodes_numpy_array):
            if elem[0] >= crop_x[0] and elem[0] <= crop_x[1]\
            and elem[1] >= crop_y[0] and elem[1] <= crop_y[1]:
                x_y_crop.append(elem)
                z_strain_rate_crop.append(field_numpy_array[j])
                z_plastic_yield_crop.append(field_numpy_array_plast_yield[j])


        x_y_crop_arr = np.array(x_y_crop)
        z_strain_rate_crop_arr = np.array(z_strain_rate_crop)
        z_plastic_yield_crop_arr = np.array(z_plastic_yield_crop)

        x,y,z= x_y_crop_arr[:,0] , x_y_crop_arr[:,1] , x_y_crop_arr[:,2]
        Z = z_strain_rate_crop_arr
        Z_plast_yield = z_plastic_yield_crop_arr

        #############################################################
        ### Grid the data (images are arrays, so we need to grid) ###
        #############################################################

        #Draw contours
        xmin, xmax = min(x), max(x)
        ymin, ymax = min(y), max(y)

        # define grid
        xi = np.linspace(xmin, xmax, Nx)
        yi = np.linspace(ymin, ymax, Ny)
        #
        # remove the data above topography
        # grid the data into array arr
        arr = griddata((x, y), Z, (xi[None,:], yi[:,None]), method='cubic')
        arr_plast_yield = griddata((x, y), Z_plast_yield, (xi[None,:], yi[:,None]), method='cubic')

        #extract the topography
        grid = {}
        i = 0
        for col in x:
            if col in grid:
              grid[col].append(y[i])
            else:
                grid[col] = [y[i]]
            i += 1

        topo = []#max topography for each x-value
        for k, v in grid.items():
            if k % res == 0:
                topo.append((k,max(v)))

        for i in range(np.size(arr,1)-1):
            remove_topo(i,topo,arr, arr_plast_yield, absolute_min_strain)

        np.save(ARR_NAME, arr)
