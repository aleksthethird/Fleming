from astropy.io import fits

import os
import parse
import numpy as np
from datetime import datetime

from . import Config


## ==================================
## TODO: Make class or pass in config

#format a number i into a string of length 3
#example: 1 becomes 001, 54 becomes 054 etc
def format_index(i):
        
    string = ""
    for j in range(3 - len(str(i))):
        string += "0"
    
    string += str(i)
    return string

#find difference between two lists
def diff(li1, li2): 

    dif = []
    
    for i in range(len(li1)):
        dif.append(li1[i] - li2[i])
        
    
    return dif

#check if two values are the same within a certain fractional difference
def equals(a, b, fraction):
    if a * (1+fraction) > b and a * (1-fraction) < b:
        return True
    return False



def get_image(config, set_number=0, image_number=0):
    
    #build path to image 
    if config.has_sets:
        image_name = config.image_format_str.format(set_number, image_number)
    else:
        image_name = config.image_format_str.format(set_number, image_number)
    
    data_path = os.path.join(config.image_dir, image_name)
        
    #read in image data
    return fits.open(data_path, ext=0)


## TODO: breaks if not using sets
def list_images(config, directory, reduced=True):
    """
    List image files in directory, as well as associated set and image number.

    Parameters
    ==========

    directory: string
        directory to list
    reduces: bool, optional
        should include the reduced prefix?

    Returns
    =======

    Array of strings contaiing image names, as well as set and image number

    """

    
    #print("[DEBUG] listing dirs in {}".format(directory))
    image_list = []

    if not config.has_sets:
        print("[Utilities] list_images currently does not support has_sets=False, sorry")
        return image_list

    if reduced:
        fmt = config.image_format_str
    else:
        fmt = config.raw_image_format_str

    for f in os.listdir(directory):
        res = parse.parse(fmt, f)
        if res != None:
            set_number = int(res.fixed[0])
            image_number = int(res.fixed[1])
            image_list.append((f, set_number, image_number))

    return image_list


## TODO: has_sets=False
def loop_images(config, reduced=True):
    """
    Return array of image names, set number and image number of each image.
    Assumes the correct n_sets and set_size.
    Should be faster than list_images()

    """

    #print("[DEBUG] looping images")
    image_list = []

    if not config.has_sets:
        print("[Utilities] loop_images currently does not support has_sets=False, sorry")
        return image_list

    if reduced:
        fmt = config.image_format_str
    else:
        fmt = config.raw_image_format_str

    for s in range(1, config.n_sets+1):
        for i in range(1, config.set_size+1):
            fname = fmt.format(s, i)
            image_list.append((fname, s, i))

    return image_list


def list_sources(config, directory=None, adjusted=False):
    """
    List list curve (source) files files in directory, as well as associated id number

    Parameters
    ==========

    directory: string, optional
        directory to list for light curves
    adjusted: bool, optional
        has the light curve been adjusted?

    Returns
    =======

    Array of strings containing file paths, as well as id number

    """

    
    if directory != None:
        #print("[DEBUG] listing light curves in {}".format(directory))
        d = directory
    else:
        if adjusted:
            #print("[DEBUG] listing adjusted sources")
            d = config.adjusted_curve_dir
        else:
            #print("[DEBUG] listing sources")
            d = config.light_curve_dir

    source_list = []

    for f in os.listdir(d):
        res = parse.parse(config.source_format_str, f)
        if res != None:
            id_number = int(res.fixed[0])
            path = os.path.join(d, f)
            source_list.append((path, id_number))

    return source_list

def loop_variables(config, ids, adjusted=False):
    """
    Loop over variables in given id table

    """

    if adjusted:
        #print("[DEBUG] listing adjusted sources")
        d = config.adjusted_curve_dir
    else:
        #print("[DEBUG] listing sources")
        d = config.light_curve_dir
            
    variables_list = []

    for i, s in enumerate(ids):
        path = os.path.join(d, config.source_format_str.format(int(s)))
        variables_list.append((i, path, s))

    return variables_list


def read_catalogue(config):
    """
    Reads catalogue file and gives a numpy 2d arrya 
    with named columns

    """

    cat = np.genfromtxt(config.catalogue_path, dtype=[
        ('id', 'int64'),
        ('xcentroid', 'float64'),
        ('ycentroid', 'float64'),
        ('sharpness', 'float64'),
        ('roundness1', 'float64'),
        ('roundness2', 'float64'),
        ('npix', 'float64'),
        ('sky', 'float64'),
        ('peak', 'float64'),
        ('flux', 'float64'),
        ('mag', 'float64'),
        ('RA', 'float64'),
        ('DEC', 'float64'),
        ],
        skip_header=1)

    return cat

#standard deviation of items in a list
def standard_deviation(a):
   
    i = 0
    mean_val = mean(a)
    
    for val in a:
        i += (val - mean_val)**2
    
    i = i / len(a)
    i = i**0.5
    return i

#mean of items in list
def mean(a):
    
    t = 0
    
    for val in a:
        t += val
    
    t = t / len(a)
    
    return t

#check that the given x and y positions are not closer to the edge than
#the specified limit
def is_within_boundaries(x, y, im_x_dim, im_y_dim, edge_distance_limit):
    
    if x >= edge_distance_limit and x <= im_x_dim - edge_distance_limit and y >= edge_distance_limit and y <= im_y_dim - edge_distance_limit:
        return True
    return False

#check that a data point lies above the specified line by a certain fraction
def is_above_line(x, y, m, c, fraction):
    exp = m * x + c

    if y > (1 + fraction) * exp:
        return True
    return False

def partition(a, start, end, ascending):
    follower = leader = start
    
    while leader < end:
        if (a[0][leader] <= a[0][end] and ascending) or (a[0][leader] >= a[0][end] and not ascending):
            for i in range(len(a)):
                a[i][follower], a[i][leader] = a[i][leader], a[i][follower]
            follower += 1
        leader += 1
    
    for i in range(len(a)):
        a[i][follower], a[i][end] = a[i][end], a[i][follower]


    return follower

def _quicksort(a, start, end, ascending):
    if start >= end:
        return
    p = partition(a, start, end, ascending)
    _quicksort(a, start, p-1, ascending)
    _quicksort(a, p+1, end, ascending)
    
def quicksort(a, ascending):
    _quicksort(a, 0, len(a[0])-1, ascending)
    
def finished_job(job_name, start_time):
    current_time = datetime.now()

    time_elapsed = current_time - start_time

    current_time_str = current_time.strftime("%H:%M:%S")
    print("[Job] Finished job '{}' at {} (took {})"
            .format(job_name, current_time_str, str(time_elapsed)))

    return current_time
        
#make .region file for comparing catalogue to actual image
def make_reg_file(directory, name, table):
    
    fname = "{}.reg".format(name)
    path = os.path.join(directory, fname)
    f = open(path, "w")
    
    xs = table['xcentroid']
    ys = table['ycentroid']
    
    #write out xs and ys of sources in catalogue
    for i in range(len(xs)):
        f.write("point " + str(xs[i]) + " " + str(ys[i]) + " # point=circle 4 \r\n")

def n_to_set_and_n(n):
    set = int((n-1)/Constants.set_size) + 1
    i = n % Constants.set_size
    if i == 0:
        i = Constants.set_size
    return set, i

def get_image_data(image_path, n):
        
            set, i = n_to_set_and_n(n)
            
            file = image_path + Constants.file_name + "_" + str(set) + "_" + format_index(i)

            file += Constants.fits_extension
        
            return fits.getdata(file, ext=0)

