#
#Exract Faults v. 2.0 by Tim Hake & Sascha Brune
#GFZ Potsdam Sec. 2.5
#31.08.18
#
'''
If you like to import another .pvtu-file change the 'FIlE_NAME'
and set import_new_data_set to T (True).
Since it takes a while to loop through all grid notes
to divide between brittle and ductile parts and remove the extrapolated data
above the topography the new data will be saved under ARR_NAME as an array.
Afterwards import_new_data_set can be reset to F (False) again and the new array
will be imported from ARR_NAME.
Note: Probably the strain rate and plastic yield fields vary in your pvtu file,
so adapt the field_strain_rate and field_plastic_yield value. As well as the
grid notes and resolution.


'''
T = True
F = False

#### EDIT HERE ################################

#Choose T if you import new data sets in order to save interpolated files
#Choose F if you want to reload previously saved files (saves a lot of time)
import_new_data_set = F

#Input/Output folders and names
M_NAME = 'FaultExtraction_asym'#a folder with M_NAME will be created in OUTPUT_DIRECTORY to save the output files
INPUT_DIRECTORY = 'D:/Uni/MA/praktikum_brune/Fortim/Model_asym/m_asym/output/solution/'
OUTPUT_DIRECTORY = 'D:/Uni/MA/praktikum_brune/Fortim/Model_asym/'
#INPUT_DIRECTORY = '/Users/brune/_prog/aspect/11_2DMenderes/omer/43a-500k/output/solution/'
#OUTPUT_DIRECTORY = '/Users/brune/_prog/aspect/11_2DMenderes/omer/43a-500k/'
FILE_NAME = INPUT_DIRECTORY + 'solution-%05d.pvtu'
ARR_NAME = OUTPUT_DIRECTORY + M_NAME + '/reduced_arrays/strain_rate_reduced_%05d.npy'

#v2.0 has time stepping. Currently it assumes 1 My time steps
num_files = 41#number of timesteps(e.g. from 0Ma to 40Ma = 41 timesteps)

#Pick your region of interest (model coordinates)
cropping = T# set to T (True) to crop the model during the creation of the Numpy-Arrays
crop_x = (150000, 350000)#outcrop of data in x-direction
crop_y = (100000, 170000)#outcrop of data in y-direction

#Insert here the ID number of strain_rate & plastic yielding field from the pvtu file
# field_strain_rate = 17
# field_plastic_yielding = 21

# set grid size
Nx = 800 # number of grid nodes in x-direction
Ny = 280 # number of grid nodes in y-direction
res = 250 # max resolution

#strain rates setting
absolute_min_strain = 1e-22#value for divison between brittle/ductile parts and removal of data above topography
cbar_range = (-21,-12)#range of colorbar for plotting the strain rate 10^-12 to 10^-21

#cut-off values for the strain rate
min_strain = 9e-16
max_strain = 2e-15

#parameter Probabilistic Hough Transform
#for explanation see skimage.transform.probabilistic_hough_line documentation
threshold = 0
min_line_length = 2
line_gap = 0

# parameter DBSCAN Cluster Analysis
#for explanation see sklearn.cluster.DBSCAN documentation
max_distance = 32
min_samples = 2

#weighting each axis by multiplying
stretch_x = 1
stretch_y = 1
stretch_z = 2.2

####
#this needs some adjustment or will be removed...
merge = F #deactivate (F) or activate the function 'merge_custer'
####
grad_tol = 15 #tolerance between the gradinet of two points and mean angle of the first cluster
alpha_tol = 10 #tolerance between the mean angle of two cluster
#see merge_cluster for further information

width = 25#the width +/- from a center linear function
#to create the envelopes to sort clusters from 'cluster_temp' to cluster

alpha_tol_2 = 25# onther tolerance between the mean angle of two cluster in line 528

#final polishing of the cluster
finalizing_cluster = T #set to T (True) to enable pop, join and squares
pop = [-1,8,11,21] #[-1,9,23,22,21,20,19,13,17] #None or [key-values, key-values, ...]
join = [[1],[4,5,15]] #None or [[key-values, key-values, ...], [key-values, key-values, ...]]
squares = [[(2,1),(225000,125000),(260000,105000),(225000,130000),(260000,110000)],
[(2,1), (200000,160000),(222000,135000),(215000,160000),(225000,135000)],
[(1,2), (230000,132000),(247000,118000),(230000,137000),(249000,122000)]]# [[(14,2),(220000,130000),(240000,120000),(240000,130000),(220000,140000)]]
#None or [[(clust_in, clust_out), p1, p2, p3, p4], [(clust_in, clust_out), p1, p2, p3, p4],...]
#see function describtion of 'disjoin_cluster' & 'create_final_cluster' below for further information

#the list pop contains key-values of the dictionary 'cluster' and those key-values will be removed from the cluster
#the list join contains list of key-values, each list will be joined together

#Plots
timestep = None #if None all timesteps will be saved
#a number allows to only saves the output of a certain timestep

#set to T (True)/ F (False) to save the corresponding figure
plot_data = F
plot_selected_strain = F
plot_skeleton = F
plot_hough_line_segments = F
plot_cluster_3D_points = F#not working jet
plot_cluster_2D_points = F
plot_all_cluster_line_segments = T
plot_new_cluster_line_segments = T
final_results = F

###############################################

def mean_tuple_lst(lst):
    #takes a list containing a start and end point and the gradient between the two points (p0,p1,z)
    #reduces the line with start (p0 = (x0,y0)) and end point (p1 = (x1,y1)) to the mid point x,y
    #calculates the mean x,y,z value of the whole list
    c = len(lst)
    x_mean, y_mean, z_mean = 0, 0, 0
    for elem in lst:
        x,y,z = (elem[0][0] + elem[1][0])/2, (elem[0][1] + elem[1][1])/2, elem[2]
        x_mean += x
        y_mean += y
        z_mean += z
    x_mean, y_mean,z_mean = x_mean/c, y_mean/c, z_mean/c
    return x_mean, y_mean, z_mean
def compare_cluster(lst):
    #takes the list 'lst' with tuples of (k,(mean x,y,z)), which are compareable in terms of
    #mean x,y,z components of the houghline midpoints
    #returns the new key-value, which is closest to the first tuple
    val = lst.pop(0)
    if len(lst) == 1:
        k = lst[0][0]
    elif len(lst) == 0:
        k = -1
    else:
        val_x_m, val_y_m, val_z_m = val
        min_dist = (0,1000)
        for c in lst:
            x_m, y_m, z_m = c[1]
            dist = np.sqrt((val_x_m - x_m)**2 + (val_y_m - y_m)**2 + (val_z_m - z_m)**2)
            if dist < min_dist[1]:
                min_dist = (c[0],dist)
        k = min_dist[0]
    return k
def merge_cluster(dict, t, grad_tol = grad_tol, alpha_tol = alpha_tol):
    #takes the dictionary 'dict' and checks if clusters exist, which are aligned,
    #merges them together and return the new dictionary
    #therefore the algorithm calculates the gradient between the mean x, y, z values
    #of two clusters and compares the gradinent with mean z value of cluster one
    #further it compares the mean z values of both clusters,
    #if every comparison is within the tolerances (grad_tol, aplha_tol) both clusters are aligned
    couples = []
    for k1, v1 in dict.items():
        for k2, v2 in dict.items():
            if (k1 != -1) and (k2 != -1):
                x1,y1,z1 = v1[3]
                x2,y2,z2 = v2[3]
                if (x2-x1) != 0 and (y2-y1) != 0:
                    gk = (x2-x1)
                    ak = (y2-y1)
                    grad = np.degrees(np.arctan(gk/ak))
                    if grad <= (z1+grad_tol) and grad >= (z1-grad_tol)\
                    and z2 <= (z1+alpha_tol) and z2 >= (z1-alpha_tol):
                        if k1 < k2:
                            pair = (k1, k2)
                        else:
                            pair = (k2, k1)
                        if pair not in couples:
                            couples.append(pair)
    print(couples)
    for n_k1, n_k2 in couples[::-1]:
        dict[n_k1][4][0] = dict[n_k1][4][0] + dict[n_k2][4][0]
        dict.pop(n_k2)

    for k, v in dict.items():
        mean =  mean_tuple_lst(v[4][0])
        dict[k][3] = mean
        dict[k][5].append(mean)
        for lst in v[4]:
            lst.append(t)
    return dict
def get_envelopes(pt, w = width):
    #point pt is has x,y,z coordinates, whereas z is the gradient at this point
    #calculates two envelopes around the point pt with the gradient m and the y-intercepts +/- width
    x,y,z = pt
    if z < 90 and z > -90:
        a = np.deg2rad(z)
        m1 = np.tan(a)
        m2 = m1
        b = y - m1*x
        b1 = b + (w/(np.cos(a)))
        b2 = b - (w/(np.cos(a)))
    else:
        print("INF")
        m1 = 30
        m2 = m1
        b1 = 30
        b2 = 30
    return m1, m2, b1, b2
def pt_between_envelopes(pt, p_envelopes):
    #checks if the point pt lies weithin two envelopes charachterised by p_envelopes
    x, y, z = pt
    m1, m2, b1, b2 = p_envelopes
    y1 = x*m1 + b1
    y2 = x*m2 + b2
    return y1 >= y and y2 <= y
def create_final_cluster(dict, lst_pop_cluster, mat_join_cluster):
    #takes the automatically generated dictionary
    #after reviewing the result eventually certain cluster need to be removed or joined
    #lst_pop_cluster conatins k values of the cluster which should be removed
    #mat_join_cluster contains list of k values of cluster which should be joined
    if mat_join_cluster != None:
        for col in mat_join_cluster:
            k = col.pop(0)
            for k_col in col:
                dict[k][4] = dict[k][4] + dict[k_col][4]
                dict.pop(k_col)
    if lst_pop_cluster != None:
        for c in lst_pop_cluster:
            dict.pop(c)
    return dict
def get_labels(loc, val, dtype = int):
    #recalculate a list of labels
    #e.g. form metres to kilometers if loc is in metres and val = 1000
    labels = []
    for c in loc:
        labels.append(dtype(c/val))
    return labels
def hist_mat(mat, bins, weight):
    #creates a histogram matrix with a certain weight
    new_mat = []
    for i in range(len(mat)):
        hist = np.histogram(mat[i], bins = bins, weights = weight[i])
        new_mat.append(hist[0].tolist())
    return new_mat
def new_colormap(r,g,b):
    #creates a color map from white to (rgb)
    c = [(1,1,1), (r,g,b)]
    n_cmap = LSC.from_list('new_cmap', c, 100)
    return n_cmap
def cluster_design(dic, lab = None):
    #creates the dictionary template which is neede for further processing
    random.seed(6)
    marker = '.,*ov^><12348spPhHxXDd'
    if np.any(lab) != None:
        for i in range(-1, max(lab)+1):
            r,g,b = random.uniform(0,1), random.uniform(0,1), random.uniform(0,1)
            m = marker[random.randint(0,len(marker)-1)]
            if i in lab:
                dic[i] = [(r,g,b),m,'Cluster %i'%i,(0,0,0),[[]],[]]
        if -1 in dic.keys():
            r,g,b = 0,0,0
            dic[-1][0] = (r,g,b)
            dic[-1][1] = '+'
            dic[-1][2] = 'Noise'
        return dic
    else:
        #rename the clusters in an increasing order
        dic_new = {}
        if -1 in dic.keys():
            r,g,b = 0,0,0
            dic_new[-1] = dic[-1]
            dic_new[-1][0] = (r,g,b)
            dic_new[-1][1] = '+'
            dic_new[-1][2] = 'Noise'
            dic.pop(-1)
        k_lst = list(dic.keys())
        for i in range(1,len(dic)+1):
            r,g,b = random.uniform(0,1), random.uniform(0,1), random.uniform(0,1)
            m = marker[random.randint(0,len(marker)-1)]
            dic_new[i] = dic[k_lst[i-1]]
            dic_new[i][0] = (r,g,b)
            dic_new[i][1] = m
            dic_new[i][2] = 'Cluster %i'%(i)
        return dic_new
def get_framework(name = FILE_NAME, crop_x = crop_x, crop_y = crop_y):
    #retruns the coordinates of the x,y coordinates of the aspect model

    ######################
    ### Load pvtu data ###
    ######################
    # pvtu loading and gridding largely based on https://stackoverflow.com/questions/23138112/vtk-to-maplotlib-using-numpy
    # load a pvtu file as input
    reader = vtk.vtkXMLPUnstructuredGridReader()
    reader.SetFileName(name%0)
    reader.Update()

    # Get the coordinates of nodes in the mesh
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
    number_of_fields = reader.GetOutput().GetPointData().GetNumberOfArrays()
    for i in range(number_of_fields):
        if reader.GetOutput().GetPointData().GetArrayName(i) == "strain_rate":
            id_strain_rate = i
    #The "strain_rate" field is scalar number 17 in this pvtu file
    #The "plastic_yielding" field is scalar number 21 in this pvtu file
    field_vtk_array = reader.GetOutput().GetPointData().GetArray(id_strain_rate)

    #Get the coordinates of the nodes and their temperatures
    nodes_numpy_array = vtk_to_numpy(nodes_vtk_array)
    x_y_crop = []

    if cropping == False:
        crop_x = (min(nodes_numpy_array[:,0]),max(nodes_numpy_array[:,0]))
        crop_y = (min(nodes_numpy_array[:,1]),max(nodes_numpy_array[:,1]))

    for j,elem in enumerate(nodes_numpy_array):
        if elem[0] >= crop_x[0] and elem[0] <= crop_x[1]\
        and elem[1] >= crop_y[0] and elem[1] <= crop_y[1]:
            x_y_crop.append(elem)
    x_y_crop_arr = np.array(x_y_crop)
    x,y,z= x_y_crop_arr[:,0] , x_y_crop_arr[:,1] , x_y_crop_arr[:,2]
    return x, y
def create_dir(out_dir = OUTPUT_DIRECTORY, mn = M_NAME, folder = None):
    # checks if the output folder exists
    #and creates if it is missing
    if not os.path.isdir(out_dir+mn+'/Results/'):
        os.mkdir(out_dir+mn+'/Results/')
    if not os.path.isdir(out_dir+mn+'/Results/'+folder):
        os.mkdir(out_dir+mn+'/Results/'+folder)
def disjoin(dict, mat):
    #takes the dictionary 'dict' and matrix 'mat'
    #mat has the shape [[(cluster_in, cluster_out),p1,p2,p3,p4],[(cluster_in, cluster_out),p1,p2,p3,p4],...]
    #where cluster_in & cluster_out are key-values of dict
    #it removes the line segments, which lay within a rectangle (p1,p2,p3,p4) in dict[cluster_in]
    #(p1,p2 is below and p3, p4 is above the rearranged parts)
    #and put the line segments into dict[cluster_out]
    if mat != None:
        for elem in mat:
            p1, p2, p3, p4 = elem[1], elem[2], elem[3], elem[4],
            clust_in, clust_out = elem[0]
            m1 = (p2[1]-p1[1])/(p2[0]-p1[0])
            b1 = p1[1]-m1*p1[0]
            m2 = (p4[1]-p3[1])/(p4[0]-p3[0])
            b2 = p3[1]-m2*p3[0]
            lst2 = []
            for i,lines in enumerate(dict[clust_in][4]):
                lst = []
                for j,line in enumerate(lines[0:len(lines)-1]):
                    st_pkt, end_pkt, m = line
                    st = recalculate_coordinates(st_pkt[0],st_pkt[1], m,'grid2real')
                    if pt_between_envelopes(st, (m2,m1,b2,b1)) \
                    and st[0]>p1[0] and st[0]<p2[0]:
                        lst.append((st_pkt, end_pkt, m))
                        lst2.append((i,j))
                lst.append(lines[len(lines)-1])
                if len(lst)>1:
                    dict[clust_out][4].append(lst)
            for elem in lst2[::-1]:
                dict[clust_in][4][elem[0]].pop(elem[1])
    return dict

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib
from matplotlib import (ticker, cm)
from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from skimage.transform import probabilistic_hough_line
from skimage.morphology import (skeletonize)
from sklearn.cluster import DBSCAN
import random
import os
import copy
from build_arr import ba

# check for correct scikit-image library
import pkg_resources
skimageversion = pkg_resources.get_distribution("scikit-image").version     # scikits yerine scikit, ve nokta yerine cizgi
#skimageversion = pkg_resources.get_distribution("scikits.image").version
print('The scikit-image version has to be 0.14.x or above. Yours is: ' + skimageversion)


if import_new_data_set == True:
    if cropping == False:
        crop_x, crop_y = None, None
    ba(num_files,FILE_NAME, ARR_NAME, OUTPUT_DIRECTORY, M_NAME,
    Nx, Ny, res, absolute_min_strain, crop_x, crop_y)

x, y = get_framework()

#############################################################
### Grid the data (images are arrays, so we need to grid) ###
#############################################################
#Draw contours
xmin, xmax = min(x), max(x)
ymin, ymax = min(y), max(y)

x_fac = (xmax-xmin)/Nx
y_fac = (ymax-ymin)/Ny
def recalculate_coordinates(x, y, z, direction, x_fac = x_fac, y_fac = x_fac):
    #takes x,y coordinates and recalculate them from
    #real coordinates to grid coordinates ('real2grid')
    #or from grid coordinates to real coordinates ('grid2real')
    #z is the angle, which stays the same
    if direction == 'grid2real':
        x = x*x_fac+xmin
        y = -(y*y_fac)+ymax
    else:
        x = (x-xmin)/x_fac
        y = -(y-ymax)/y_fac
    return x,y,z

# define grid
xi = np.linspace(xmin, xmax, Nx)
yi = np.linspace(ymin, ymax, Ny)

lvl = []
for i in range(cbar_range[0],(cbar_range[1]+1),1):
    s = '1e%i'%i
    lvl.append(float(s))

cluster = {}
for time in range(num_files):
    print('Timestep %d'%time)
    arr = np.load(ARR_NAME %time)#load reduced data array
    if timestep != None:
        test_time = timestep
    else:
        test_time = time
    if plot_data == True and test_time == time:
        print('saving data image %d/%i'%(time,num_files-1))
        create_dir(folder = 'data')
        fig = plt.figure(dpi=200)
        ax = fig.add_subplot(111)
        plt.contourf(xi,yi,arr,13,locator=ticker.LogLocator(), levels = lvl, cmap=plt.cm.coolwarm)
        ax.set_aspect(1)
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(xmin, xmax)
        xloc, xlab = plt.xticks()
        yloc, ylab = plt.yticks()
        plt.xticks(xloc,get_labels(xloc, 1000, dtype = int))
        plt.yticks(yloc,get_labels(yloc, 1000, dtype = int))
        cbar = plt.colorbar()
        cbar.set_label('Strain rate 1/s')
        ax.set_xlabel('X Axis [km]')
        ax.set_ylabel('Y Axis [km]')
        ax.set_title('Brittle/Ductile divided Data at %dMa'%time)
        plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+'data'+'/%d.png'%time)
        plt.close(fig)
    #######################################################
    ### Apply transformations and create a scikit image ###
    #######################################################

    # stretch contrast to minimum and maximum cutoff strain rate
    mincutoff = min_strain# 1/s
    maxcutoff = max_strain# 1/s
    # shift arr so that it is 0 at mincutoff
    arr = arr - mincutoff
    # scale arr so that it is 1 at maxcutoff
    arr = arr/maxcutoff
    # set all values outside of [0,1] to cutoff limits
    arr[np.isnan(arr)] = 0
    arr[(arr<0)] = 0
    arr[(arr>1)] = 1

    # flip up-down (image convention)
    normTi_flipped = np.flipud(arr)

    # assign image
    img = np.array(normTi_flipped, dtype=float)

    # perform skeletonization
    skeleton = skeletonize(img)

    #########################
    ###HOUGH TRANSFORMATON###
    #########################

    # Line finding using the Probabilistic Hough Transform
    lines = probabilistic_hough_line(skeleton, threshold=threshold, line_length=min_line_length,
                                     line_gap=line_gap, seed = 12)

    ########################
    ### Cluster Analysis ###
    ########################

    # calculate the gradient of the line segments
    segments = []
    for line in lines:
        p0, p1 = line
        temp_mid_point = (((p0[0] + p1[0])//2), ((p0[1] + p1[1])//2))
        if p1[0]-p0[0]!=0:
            gk = (p1[1]-p0[1])
            ak = (p1[0]-p0[0])
            alpha = np.degrees(np.arctan(gk/ak))
        else:
            alpha = 90
        segments.append([stretch_x*temp_mid_point[0],\
        stretch_y*temp_mid_point[1],\
        stretch_z*alpha])

    arr_segments = np.asarray(segments)

    #compute DBSCAN cluster analysis
    if len(lines) != 0:
        db = DBSCAN(eps = max_distance, min_samples = min_samples).fit(arr_segments, y = None, sample_weight = None)
        line_labels = db.labels_

        #differentiate cluster
        cluster_temp = {} #cluster of the certain timestep
        cluster_temp = cluster_design(cluster_temp, line_labels)

        #connect the start, end point of each line segment with the gradient
        #sort it to the corresponding cluster
        for i,p in enumerate(arr_segments):
            p0,p1 = lines[i]
            if p1[0]-p0[0]!=0:
                gk = (p1[1]-p0[1])
                ak = (p1[0]-p0[0])
                alpha = np.degrees(np.arctan(gk/ak))
            else:
                alpha = 90
            cluster_temp[line_labels[i]][4][0].append((lines[i][0], lines[i][1], alpha))

        #calculate the mean x,y,z values for each cluster of line segments
        for k,v in cluster_temp.items():
            x_y_z_mean = mean_tuple_lst(v[4][0])
            cluster_temp[k][3] = x_y_z_mean
            if merge == False:
                v[4][0].append(time)

        if merge == True:
            # merge cluster which seems to belong together within a timestep
            cluster_temp = merge_cluster(cluster_temp, time)

        #pass the temporary cluster analysis for one timestep (cluster_temp)
        #to a the general dictionary 'cluster', which saves the information from all timesteps
        #therefore the algorithm compares the mean x,y,z values per cluster in 'cluster_temp'
        #with the current x,y,z value representing a cluster in 'cluster'
        #if there is a match the cluster from 'cluster_temp' will be saved in the matching cluster in 'cluster'
        #if there is no match a new cluster will be designed in 'cluster'
        if len(cluster) == 0:
            cluster = cluster_temp
            if -1 not in cluster_temp.keys():
                cluster[-1] = [(0,0,0),'+','Noise',(0,0,0),[[]],[]]
        else:
            new_clusters = []
            l = max(cluster.keys())
            m = 1
            for k_temp, v_temp in cluster_temp.items():
                comparable_cluster = [v_temp[3]]
                test = 0
                envelopes = get_envelopes(v_temp[3])
                if k_temp == -1:
                    cluster[-1][4].append(v_temp[4][0])
                    cluster[-1][5].append(v_temp[3])
                    test = 1
                for k,v in cluster.items():
                    if k!=-1 and pt_between_envelopes(v[3], envelopes) and\
                    v[3][2]+alpha_tol_2 >= v_temp[3][2] and v[3][2]-alpha_tol_2 <= v_temp[3][2]:
                        comparable_cluster.append((k,v[3]))
                        test = 1
                new_k = compare_cluster(comparable_cluster)
                cluster[new_k][4].append(v_temp[4][0])
                cluster[new_k][5].append(v_temp[3])
                cluster[new_k][3] = v_temp[3]
                if test == 0:
                    new_clusters.append((l+m,cluster_temp[k_temp]))
                    m += 1
            for elem in new_clusters:
                cluster[elem[0]] = elem[1]
    else:
        print('No lines detected')

    if plot_selected_strain == True and test_time == time:
        print('saving selected strain image %d/%i'%(time,num_files-1))
        create_dir(folder = 'selected_strain')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(img, cmap=plt.cm.gray)
        ax.set_xlim((0, skeleton.shape[1]))
        ax.set_ylim((skeleton.shape[0], 0))
        ax.set_xlabel('X Axis (Gridcoordinates)')
        ax.set_ylabel('Y Axis (Gridcoordinates)')
        ax.set_title('Selected Strain Rates at Time %d'%time)
        fig.tight_layout()
        ax.set_aspect(1)
        plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+'selected_strain'+'/%d.png'%time)
        plt.close(fig)
    if plot_skeleton == True and test_time == time:
        print('saving skeleton image %d/%i'%(time,num_files-1))
        create_dir(folder = 'skeleton')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(skeleton, cmap=plt.cm.gray)
        ax.set_xlim((0, skeleton.shape[1]))
        ax.set_ylim((skeleton.shape[0], 0))
        ax.set_xlabel('X Axis (Gridcoordinates)')
        ax.set_ylabel('Y Axis (Gridcoordinates)')
        ax.set_title('Skeleton at Time %d'%time)
        fig.tight_layout()
        ax.set_aspect(1)
        plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+'skeleton'+'/%d.png'%time)
        plt.close(fig)
    if plot_hough_line_segments == True and test_time == time:
        print('saving hough line segments image %d/%i'%(time,num_files-1))
        create_dir(folder = 'hough_lines')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(skeleton, cmap=plt.cm.gray)
        for line in lines:
            p0, p1 = line
            ax.plot((p0[0], p1[0]), (p0[1], p1[1]))
        ax.set_xlim((0, skeleton.shape[1]))
        ax.set_ylim((skeleton.shape[0], 0))
        ax.set_xlabel('X Axis (Gridcoordinates)')
        ax.set_ylabel('Y Axis (Gridcoordinates)')
        ax.set_title('Probabilistic Hough Line Transformation at Time %d'%time)
        fig.tight_layout()
        ax.set_aspect(1)
        plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+'hough_lines'+'/%d.png'%time)
        plt.close(fig)

##################
###Plot Results###
##################
#final polishing of the clusters with removing, joining or rearrange some
cluster = cluster_design(cluster)
if finalizing_cluster == True:
    final_cluster = copy.deepcopy(cluster)
    final_cluster = disjoin(final_cluster, squares)
    final_cluster = create_final_cluster(final_cluster,pop,join)
    final_cluster = cluster_design(final_cluster)
else:
    final_cluster = copy.deepcopy(cluster)

print(final_cluster.keys())
if timestep == None:
    t_start = 0
    t_end = num_files
else:
    t_start = timestep
    t_end = timestep + 1
if plot_cluster_3D_points == True:
    create_dir(folder = '3D_points')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    labels = []
    for i,p in enumerate(arr_segments):
        x,y,z = p[0],p[1],p[2]
        sc = ax.scatter(x,y,z, color = cluster[line_labels[i]][0], marker = cluster[line_labels[i]][1])
        if line_labels[i] not in labels:
            sc.set_label(cluster[line_labels[i]][2])
            labels.append(line_labels[i])
    ax.legend()
    ax.set_xlim((0, skeleton.shape[1]))
    ax.set_ylim((skeleton.shape[0], 0))
    ax.set_title('Cluster Analysis Points %05d'%time)
    ax.set_xlabel('X')
    ax.set_aspect(1)
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+'3D_points'+'/%d.png'%time)
    plt.close(fig)
if plot_cluster_2D_points == True:
    create_dir(folder = '2D_points')
    print('saving 2D points...')
    for t_current in range (t_start, t_end):
        print('Image %d/%i'%(t_current,num_files-1))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        labels_points = []
        for k,v in cluster.items():
            for lines in v[4]:
                for line in lines[0:len(lines)-1]:
                    x = (line[0][0]+line[1][0])/2
                    y = (line[0][1]+line[1][1])/2
                    z, t_test = line[2], lines[len(lines)-1]
                    if t_test == t_current:
                        x,y,z = recalculate_coordinates(x,y,z,'grid2real')
                        sc = ax.scatter(x,y, color = cluster[k][0], marker = cluster[k][1])
                        if k not in labels_points:
                            sc.set_label(cluster[k][2])
                            labels_points.append(k)

        ax.set_title('Cluster Analysis Midpoints at Time %dMa'%t_current)
        ax.set_xlabel('X Axis [km]')
        ax.set_aspect(1)
        ax.set_ylabel('Y Axis [km]')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        xloc, xlab = plt.xticks()
        yloc, ylab = plt.yticks()
        plt.xticks(xloc,get_labels(xloc, 1000, dtype = int))
        plt.yticks(yloc,get_labels(yloc, 1000, dtype = int))
        ax.legend(loc = 1)
        plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+'2D_points'+'/%d.png'%t_current)
        plt.close(fig)
if plot_all_cluster_line_segments == True or \
plot_new_cluster_line_segments == True:
    print('saving lines...')
    if plot_all_cluster_line_segments == True and\
    plot_new_cluster_line_segments == True:
        create_dir(folder = 'lines_new')
        create_dir(folder = 'lines_all')
        f = 'lines_all'
        f1 = 'lines_new'
        c = cluster.copy()
        c1 = final_cluster.copy()
        for t_current in range(t_start, t_end):
            print('Image %d/%i'%(t_current,num_files-1))
            labels_line = []
            fig = plt.figure(dpi=200)
            ax = fig.add_subplot(111)
            for k, v in c1.items():
                for lines in v[4]:
                    for line in lines[0:len(lines)-1]:
                        t_test = lines[len(lines)-1]
                        p0,p1,alpha = line
                        if t_test == t_current:
                            x0,y0,z0 = recalculate_coordinates(p0[0], p0[1], alpha,'grid2real')
                            x1,y1,z0 = recalculate_coordinates(p1[0], p1[1], alpha,'grid2real')
                            line, = ax.plot((x0,x1), (y0,y1), color = c1[k][0], alpha = .5)
                            if k not in labels_line:
                                line.set_label(c1[k][2])
                                labels_line.append(k)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.legend()
            xloc, xlab = plt.xticks()
            yloc, ylab = plt.yticks()
            plt.xticks(xloc,get_labels(xloc, 1000, dtype = int))
            plt.yticks(yloc,get_labels(yloc, 1000, dtype = int))
            ax.set_aspect(1)
            ax.set_ylabel('Y Axis [km]')
            ax.set_xlabel('X Axis [km]')
            ax.set_title('Cluster Analysis Line Segments at Time %dMa'%t_current)
            plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+f1+'/%d.png'%t_current)
            plt.close(fig)
    elif plot_all_cluster_line_segments == True:
        f = 'lines_all'
        create_dir(folder = f)
        c = cluster.copy()
    elif plot_new_cluster_line_segments == True:
        f = 'lines_new'
        create_dir(folder = f)
        c = final_cluster.copy()

    for t_current in range(t_start, t_end):
        print('Image %d/%i'%(t_current,num_files-1))
        labels_line = []
        fig = plt.figure(dpi=200)
        ax = fig.add_subplot(111)
        for k, v in c.items():
            for lines in v[4]:
                for line in lines[0:len(lines)-1]:
                    t_test = lines[len(lines)-1]
                    p0,p1,alpha = line
                    if t_test == t_current:
                        x0,y0,z0 = recalculate_coordinates(p0[0], p0[1], alpha,'grid2real')
                        x1,y1,z0 = recalculate_coordinates(p1[0], p1[1], alpha,'grid2real')
                        line, = ax.plot((x0,x1), (y0,y1), color = c[k][0], alpha = .5)
                        if k not in labels_line:
                            line.set_label(c[k][2])
                            labels_line.append(k)

        ax.legend(loc = 1)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        xloc, xlab = plt.xticks()
        yloc, ylab = plt.yticks()
        plt.xticks(xloc,get_labels(xloc, 1000, dtype = int))
        plt.yticks(yloc,get_labels(yloc, 1000, dtype = int))
        ax.set_aspect(1)
        ax.set_ylabel('Y Axis [km]')
        ax.set_xlabel('X Axis [km]')
        ax.set_title('Cluster Analysis Line Segments at Time %dMa'%t_current)
        plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+f+'/%d.png'%t_current)
        plt.close(fig)
if final_results == True:
    create_dir(folder = 'final_result')
    bins = np.arange(-90,91,5)
    angles = []
    lengths = []
    t_ax = np.arange(0,41*10e6,10e6)
    c = 0
    for k,v in final_cluster.items():
        angles_temp = [k]
        lengths_temp = []
        for t in range(41):
            angles_t = []
            lengths_t = []
            for lines in v[4]:
                for line in lines[0:len(lines)-1]:
                    if lines[len(lines)-1] == t:
                        p0, p1, alpha = line
                        x0, x1 = p0[0]*x_fac, p1[0]*x_fac
                        y0, y1 = p0[1]*y_fac, p1[1]*y_fac
                        c2 = (x1-x0)**2 + (y1-y0)**2
                        l = np.sqrt(c2)
                        angles_t.append(alpha)
                        lengths_t.append(l/1000)
            angles_temp.append(angles_t)
            lengths_temp.append(lengths_t)
        angles.append(angles_temp)
        lengths.append(lengths_temp)

    length_mat = []
    for j in range(len(angles)):
        k = angles[j].pop(0)
        length_mat.append((k,hist_mat(angles[j], bins, lengths[j])))

    w = 30
    matplotlib.rcParams.update({'font.size': 20})
    print('saving final results...')
    for time in range(t_start, t_end):
        print('Image %d/%i'%(time,num_files-1))
        arr = np.load(ARR_NAME %time)#load reduced data array
        fig = plt.figure(1,dpi=200, figsize = (w/1.7,w))
        gs = gridspec.GridSpec(len(final_cluster)+3,2)
        ax = plt.subplot2grid((len(final_cluster)+3,2),(0,0),colspan = 2, rowspan = 2)
        plt.contourf(xi,yi,arr,levels = lvl, locator = ticker.LogLocator(),cmap = new_colormap(0,0,0))
        for k, v in final_cluster.items():
            for lines in v[4]:
                for line in lines[0:len(lines)-1]:
                    if lines[len(lines)-1] == time:
                        p0,p1,alpha = line
                        p0 = recalculate_coordinates(p0[0],p0[1], alpha, 'grid2real')
                        p1 = recalculate_coordinates(p1[0],p1[1], alpha, 'grid2real')
                        line, = ax.plot((p0[0],p1[0]),(p0[1],p1[1]),color = v[0], alpha = .5, linewidth = 8)

        ax.set_aspect(1)
        plt.clim(10e-24,10e-13)
        cbar=plt.colorbar()
        cbar.set_label('strain rate s^-1')
        xloc, xlab = plt.xticks()
        yloc, ylab = plt.yticks()
        plt.xticks(xloc,get_labels(xloc, 1000, dtype = int))
        plt.yticks(yloc,get_labels(yloc, 1000, dtype = int))
        ax.set_ylabel('Y Axis')
        ax.set_xlabel('X Axis')
        ax.set_title('Time: %d Ma'%time)
        ar_posi = time*(10e+6)
        for i in range(1,len(length_mat)+1):
            r,g,b = final_cluster[length_mat[i-1][0]][0]
            ax = plt.subplot2grid((len(final_cluster)+3,2),(i+2,0), colspan = 2)
            im = plt.imshow(np.transpose(length_mat[i-1][1])[ ::-1,:], extent = (min(t_ax), max(t_ax), -90, 90), interpolation = 'none', cmap = new_colormap(r,g,b), aspect = 'auto')
            if i == 1:
                ax.annotate('',xy=(ar_posi, 90), xytext=(ar_posi, 91),
                        arrowprops=dict(facecolor='black', headwidth = 15))
            elif i == len(length_mat):
                ax.annotate('',xy=(ar_posi, -90), xytext=(ar_posi, -91),
                arrowprops=dict(facecolor='black', headwidth = 15))
            cbar = plt.colorbar()
            cbar.set_label('Length [km]\n%s'%final_cluster[length_mat[i-1][0]][2])
            ax.set_xticklabels([])
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(40))
            ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(20))
            ax.grid(b = True, axis = 'y', which = 'both', color = 'k', alpha = .3, linewidth = 1)
            plt.plot((ar_posi,ar_posi),(-90,90), color = 'k', alpha = .3)
        xloc, xlab = plt.xticks()
        plt.xticks(xloc,get_labels(xloc, 10e6, dtype = int))
        plt.xlabel('Time [Ma]')
        fig.text(0.05,1-5/9 , 'Angle [deg]', va = 'center', rotation = 'vertical')
        plt.savefig(OUTPUT_DIRECTORY+M_NAME+'/Results/'+'final_result'+'/%d.png'%time)
        plt.close(fig)
if plot_data == True or plot_selected_strain == True or final_results == True or\
plot_skeleton == True or plot_hough_line_segments == True or\
plot_cluster_3D_points == True or plot_cluster_2D_points == True or\
plot_all_cluster_line_segments == True or plot_new_cluster_line_segments == True:
    print('Files saved at %s/Results'%(OUTPUT_DIRECTORY+M_NAME))
