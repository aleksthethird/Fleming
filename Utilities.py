from datetime import datetime
import Constants
from astropy.io import fits


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
    
def print_job(job):
    
    now = datetime.now()

    current_time = now.strftime("%H:%M:%S")
    print("Finished " + job + " at " + current_time)
        
#make .region file for comparing catalogue to actual image
def make_reg_file(directory, name, table):
    
    f = open(directory + name + ".reg", "w")
    
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

