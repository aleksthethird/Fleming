import Constants
import Utilities
from astropy.io.fits import HDUList
from astropy.io.fits import PrimaryHDU
from astropy.io.fits import getheader
from astropy.table import Table
import datetime
import random
import math
import os
import Cataloguer
import FluxFinder
import matplotlib.pyplot as plt



import numpy as np

class StreakFinder:
    
    #folder where initial images are stored
    directory = None
    
    #pixel shifts for each image (total shift from first image)
    x_shifts = None
    y_shifts = None
    
    first_image = []
    second_image = []
    
    streaks = []
    cataloguer = None
    ff = None
    
    
   
    def __init__(self, directory, cataloguer, ff):
        self.directory = directory
        self.cataloguer = cataloguer
        self.ff = ff
        

    #find all streaks in the dataset
    def find_all_streaks(self):
        
        #stores which images have been scanned already
        checked = []
        
        #initialise checked array - no images checked thus far
        for i in range(Constants.n_sets*Constants.set_size):
            checked.append(False)
            self.streaks.append([])
        
        #obtain all shifts from the shift file
        self.get_shifts()
        
        #loop through the entire dataset, scanning every fourth image
        #for i in range(0, Constants.set_size*Constants.n_sets, 4):
        for i in range(189, 190, 4):
            
            Utilities.print_job(str(i))
        
            #prepare image i for being scanned for streaks
            self.prepare_image(i)
            
            #indicate that image i has been scanned
            checked[i] = True
            
            #if streaks have been found, check the surrounding images (i-3 to 
            #i+3) for streaks, taking care not to scan an image that has 
            #already been scanned.
            if self.find_streaks(i):
                
                for j in range(-3, 4):
                    
                    if not checked[i+j] and not i == j:
                        
                        #prepare images i+j for being scanned
                        self.prepare_image(i+j)
                        checked[i+j] = True
    
    #output a table containing the id, xcentroid, ycentroid, RA, Dec, RA/s and Dec/s
    def output_results(self):
        
        #get time of observation for each image
        times_file = self.directory + Constants.working_directory + Constants.time_file
        times = [line.rstrip('\n') for line in open(times_file)]
        
        #make streak folder in results folder
        output_dir = self.directory + Constants.working_directory  + Constants.output_directory + Constants.streak_folder
        if not os.path.exists(output_dir):
            os.mkdir(output_dir) 
        
        #path to which the table fo streaks is written
        streak_file = output_dir + Constants.streak_file
        
        f = open(streak_file, "w")
        f.write("id  time                 xcentroid         ycentroid         ra                dec               RA/s               Dec/s              ")
        
        #stores the number of consecutive images with streaks
        count = 0
        #identifier for streak
        objid = 1
        
        #ra and dec for an individual streak
        ras = []
        decs = []
        
        #time each image was taken for a consecutive set of images which contains
        #a streak
        streak_times = []
        
        #loop through array of streaks found. This array contains an x and y position
        #for the centre of each streak in each image (one element for each image).
        #If no streaks were found in that image, the element is empty. 
        #Here we are checking for a consecutive set of images which contain a streak.

        for i in range(len(self.streaks)):
            
            #if image i contains a streak
            if len(self.streaks[i]) > 0:
                
                    count += 1
                    
                    #if the streak is at least 150 pixels away from the edge of the image,
                    #record its x and y position and the time at which it was at 
                    #that position. For use in calculating object velocity
                    if Utilities.is_within_boundaries(self.streaks[i][0][0], self.streaks[i][0][1], Constants.image_width, Constants.image_height, 150):
                        
                        time_string = times[i]
                        time_units = time_string.split("T")
                        time_units = time_units[0].split("-") + time_units[1].split(":")
                        time = datetime.datetime(int(time_units[0]), int(time_units[1]), int(time_units[2]), int(time_units[3]), int(time_units[4]), int(time_units[5]))
                        
                        #add half the exposure time to the time of observation.
                        #This is done because the positions stored are at the centre
                        #of the streak. The object is at this position halfway
                        #through the exposure .
                        time = time + datetime.timedelta(0,17.5)
        
                        streak_times.append(time)
                        
                        #use world coordinate system from Cataloguer to find 
                        #the RA and Dec of the centre of the streak
                        ra, dec = self.cataloguer.wcs.all_pix2world(self.streaks[i][0][0], self.streaks[i][0][1], 0)

                        ras.append(ra)
                        decs.append(dec)
                
            #if image i does not contain a streak, or the end of the array has been reached
            if len(self.streaks[i]) == 0 or i == len(self.streaks)-1:
                
                #if the minimum number of consecutive images with a streak has
                #been reached
                if count > 3:
                    
                    #find a time associated with one of the images and the 
                    #ra and dec of the centre of the streak in that image  
                    time = streak_times[len(streak_times)-2]
                    ra, dec = self.cataloguer.wcs.all_pix2world(self.streaks[i-2][0][0], self.streaks[i-2][0][1], 0)
                    
                    #prepare string for output
                    output_string = "\r\n" + Utilities.format_index(objid)  + " " + time.strftime("%m/%d/%Y, %H:%M:%S") + " " + str(self.streaks[i-2][0][0]) + " " +  str(self.streaks[i-2][0][1]) + " " + str(ra) + " " + str(dec)
                    
                    
                    avg_ra_velocity = 0
                    avg_dec_velocity = 0
                    
                    #calculate the angular velocity in terms of RA/s and Dec/s
                    #in arcseconds/s by finding the mean change in RA and Dec of 
                    #the centre of the streak in consecutive images and dividing
                    #this by the exposure time
                    for j in range(len(ras)-1):
                        
                        time1 = streak_times[j]
                        time2 = streak_times[j+1]
                        
                        diff = (time2-time1).total_seconds()
                        avg_ra_velocity += (ras[j+1]-ras[j])/diff
                        avg_dec_velocity += (decs[j+1]-decs[j])/diff
                    
                    avg_ra_velocity = avg_ra_velocity / (len(ras) - 1)
                    avg_dec_velocity = avg_dec_velocity / (len(decs) - 1)

                    #add angular velocities to output string (*3600 to convert to arcseconds/s)
                    output_string += " " + str(avg_ra_velocity*3600) + " " + str(avg_dec_velocity*3600)
                    
                    #write output string to file
                    f.write(output_string)

                    #create and save thumnbnail of the streak in the results/streaks folder
                    im = self.ff.get_thumbnail(i-2, self.streaks[i-3][0][0], self.streaks[i-3][0][1], 100, False)
                    plt.axis('off')
                    plt.imshow(im, origin='lower', cmap=plt.cm.inferno)
                    plt.savefig(output_dir + "streak_" + Utilities.format_index(objid) + ".jpg")
                    
                    #increment object id (so that each streak has a unique ID)
                    objid += 1

                #since the current image has no streak, or the end of the streak
                #array has been reached, the number of consecutive images is 
                #reset to zero
                count = 0
                x_pos = []
                y_pos = []
                streak_times = []
        
    
        
                
                
            
    #prepare image n to be scanned for streaks. This involves dividing
    #image n by the preceding image to remove stars from the image and removing
    #as many bright pixels as possible - reason explained in the find_streaks()
    def prepare_image(self, n):
        
        self.second_image = []
        
        #read in image n
        path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix 
        temp_image = Utilities.get_image_data(path, n)
        
        set, i = Utilities.n_to_set_and_n(n)
            
        file = path + Constants.file_name + "_" + str(set) + "_" + Utilities.format_index(i) + Constants.fits_extension
        head=getheader(file ,ignore_missing_end=True)

        head=getheader(path + Constants.file_name + "_" + str(set) + "_" + Utilities.format_index(i) + Constants.fits_extension,ignore_missing_end=True)
        
        #read in the image n-1
        first_image = Utilities.get_image_data(path, n-1)
        
        #calculate the shift between the these two images
        x_shift = self.x_shifts[n-1] - self.x_shifts[n-2]
        y_shift = self.y_shifts[n-1] - self.y_shifts[n-2]
        
        #calculate the mean of the preceding image
        mean = np.mean(first_image)
        
        for i in range(len(temp_image)):
            for j in range(len(temp_image[0])):
                
                #find the expected position of pixel (j, i) from image n in 
                #image n-1. This is done to align the images
                y = int(i - round(y_shift))
                x = int(j - round(x_shift))
                
                #if the expected position is within the bounds of the image
                if x > 0 and y > 0 and x < len(temp_image[0]) and y < len(temp_image):
                    
                    #divide the count stored in pixel (j, i) by the count stored
                    #in the equivalent pixel in the preceding image divided by the 
                    #mean of the preceding image. 
                    temp_image[i][j] = temp_image[i][j]/(first_image[y][x]/mean)
                
        #calculate median of image n
        median = np.median(temp_image)
        
        #loop through each pixel in image n 
        for i in range(len(temp_image)):
            
            #prepare array to store the final version of image n for scanning
            self.second_image.append([])
            
            for j in range(len(temp_image[0])):
                
                #if the count in pixel (j, i) is greater than the median
                if temp_image[i][j] > median:
                    
                    #stores the number of pixels surrounding this pixel
                    #that have a count higher than the median
                    count = 0
                    
                    #find number of surrounding pixels which have a count higher
                    #than the median
                    for k in range(-1, 2):
                        for l in range(-1, 2):
                            
                            if not (k == 0 and l == 0):
                                
                                ik = i + k
                                jl = j + l
                                
                                if ik >= 0 and jl >= 0 and ik < len(temp_image) and jl < len(temp_image[0]):
                                   
                                    if temp_image[ik][jl] > median:
                                        count += 1
                    
                    #if the number of surrounding pixels with a count higher 
                    #than the median is less that 3, append the median rather
                    #than the actual pixel value to the final version of the image
                    if count < 3:
                        self.second_image[i].append(median)
                    else:
                        self.second_image[i].append(temp_image[i][j])

                #pixels with a count lower than the median are automatically
                #inserted into the final version
                else:
                    
                    self.second_image[i].append(temp_image[i][j])
         
        
        #hdu = PrimaryHDU(self.second_image, head)
        #hdul = HDUList([hdu], None)
        #hdul.writeto(self.directory + Constants.working_directory + "testimage.fits", overwrite=True)
        
    
    #find all of the streaks in image n 
    def find_streaks(self, n):
        
        #stores x and y positions of the centres of any streaks found in the 
        #image
        streaks = []
        
        #stores the pixels that make up each streak
        streak_pixels = []
        
        #refer to image n as 'data' (for simplicity)
        data = self.second_image
        
        median = np.median(data)
        
        #stores all pixels with counts greater than the median which have been scanned
        all_pixels = set()
        
        #stores all pixels which form streaks in the image
        all_streaks = set()
        
        #loop through each pixel in the image
        for i in range(len(data)):
            
            for j in range(len(data[0])): 
                
                #if the count in the pixel is greater than 1.04 times the median
                #and it has not been scanned already
                if data[i][j] > median*1.04 and not string in all_pixels:
                    
                    #form string storing x and y position of pixel, to be stored
                    #in the various sets of pixels. Cannot simply use an array 
                    #as these cannot be inserted into a set.
                    string = str(i) + " " + str(j)

                    
                    completed = False
                    to_scan = []
                    pixels = []
                    pixels.append(string)
                    to_scan.append(string)
                    all_pixels.add(string)

                    while not completed:
                        arr = to_scan[0].split(" ")
                        y = int(arr[0])
                        x = int(arr[1])
                        
                        found_pixels = self.find_pixels(x, y, data, median, all_pixels)
                        if len(pixels) > 2000:
                            completed = True
                        pixels = pixels + found_pixels
                        to_scan = to_scan + found_pixels
                        
                        to_scan.pop(0)
                        
                        if len(to_scan) == 0:
                            completed = True
                    
                    if len(pixels) > 60:
                        
                        x_centre = 0
                        y_centre = 0
                        max_dist = 0
                        
                        for i in range(len(pixels)):
                            
                            arr = pixels[i].split(" ")
                            y = int(arr[0])
                            x = int(arr[1])
                            x_centre += x
                            y_centre += y
                            
                            for k in range(i, len(pixels)):
                                if not k == i:
                             
                                     
                                    arr = pixels[k].split(" ")
                                    y2 = int(arr[0])
                                    x2 = int(arr[1])
                                            
                                    dist = ((x2-x)**2 + (y2-y)**2)**0.5
                                    
                                    if dist > max_dist:
                                        max_dist = dist
                                
                        
                        x_centre = int(round(x_centre / len(pixels)))
                        y_centre = int(round(y_centre / len(pixels)))
                        
                        if not Utilities.is_within_boundaries(x_centre, y_centre, len(data[0]), len(data), 30):
                            continue
                        
                        vector = [1, 0]
                        occupancies = []
                        max_occupancy = 0
                        max_vector = 0
                        
                        for angle in range(0, 180):
                            
                            count = 0
                            line_pixels = set()
                            
                            for a in range(-30, 31):
                                for b in range(-30, 31):
                                    y = y_centre + a
                                    x = x_centre + b
                                    
                                    if self.distance_to_line([x_centre, y_centre], vector, [x, y]) < 3:
                                        line_pixels.add(str(y) + " " + str(x))
                            
                            for pix in line_pixels:
                                if pix in pixels:
                                    count = count + 1
                            
                            occupancy = count/len(line_pixels)
                            occupancies.append(occupancy)
                            
                           # if occupancy > max_occupancy:
                            #    max_occupancy = occupancy
                             #   max_vector = vector
                            
                            vector = self.rotate(vector, 1)
                            
                        
                        standard_deviation = Utilities.standard_deviation(occupancies)
                        mean = np.mean(occupancies)
                        
                        
                        if not standard_deviation < 0.4 * mean:
                            
                            streak_pixels.append(pixels)

                            #all_streaks = all_streaks|set(pixels)
        
        if len(streak_pixels) == 0:
            return False
        
        if len(streak_pixels) > 1:
            completed = False
            i = 0
            
            while not completed:
                
                same_object = False
                
                for pixel1 in streak_pixels[i]:
                    
                    if same_object:
                        
                        streak_pixels[i+1] = streak_pixels[i+1] + streak_pixels[i]
                        streak_pixels.pop(i)
                        i -= 1
                        
                        break
                    
                    for pixel2 in streak_pixels[i+1]:
                        
                            arr = pixel1.split(" ")
                            y1 = int(arr[0])
                            x1 = int(arr[1])
                            
                            arr = pixel2.split(" ")
                            y2 = int(arr[0])
                            x2 = int(arr[1])
                            
                            dist = ((x2-x1)**2 + (y2-y1)**2)**0.5
                            
                            if dist < 100:
                                same_object = True
                                break
                i += 1
                
                if i == len(streak_pixels)-1:
                    completed = True
            
        for pixels in streak_pixels:
            
            x_centre = 0
            y_centre = 0
            
            for i in range(len(pixels)):
                
                arr = pixels[i].split(" ")
                y = int(arr[0])
                x = int(arr[1])
                x_centre += x
                y_centre += y
            
            x_centre = x_centre / len(pixels)
            y_centre = y_centre / len(pixels)
            streaks.append([x_centre, y_centre])
            
            
        self.streaks[n-1] = streaks
                
                        
# =============================================================================
#         for string in all_streaks:
#             
#             arr = string.split(" ")
#             y = int(arr[0])
#             x = int(arr[1])
#             
#             self.second_image[y][x] = 30000
#             
#         path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix 
# 
#         image = Utilities.get_image_data(path, 189)
#         
#         set_n, i = Utilities.n_to_set_and_n(189)
#             
#         file = path + Constants.file_name + "_" + str(set_n) + "_" + Utilities.format_index(i) + Constants.fits_extension
#         head=getheader(file ,ignore_missing_end=True)
#         
#         hdu = PrimaryHDU(self.second_image, head)
#         hdul = HDUList([hdu], None)
#         hdul.writeto(self.directory + Constants.working_directory + "testimage.fits", overwrite=True)
#     
# =============================================================================
        if not len(streaks)  == 0:
            return True
        
        return False
        
        
    def find_pixels(self, x, y, image, median, all_pixels):
        
        pixels = []
        for i in range(-1, 2):
            
            for j in range(-1, 2):
                if not (i == 0 and j == 0):
                    
                    iy = y + i 
                    jx = x + j
                    
                    if iy >= 0 and jx >= 0 and iy < len(image) and jx < len(image[0]):
                
                        if image[y+i][x+j] > median*1.04:
                            string = str(y+i) + " " + str(x+j)
                            if not string in all_pixels:
                                pixels.append(string)
                                all_pixels.add(string)
        return pixels
                    
                
            
    #read the shifts file
    def get_shifts(self):
        
        #build shift file path
        shifts_path = self.directory + Constants.working_directory + Constants.shift_file
        
        #read shifts in as table
        t = Table.read(shifts_path, format = Constants.table_format)
        
        
        self.x_shifts = t['xshifts']
        self.y_shifts = t['yshifts']
    
    def rotate(self, vector, angle):
        
        new_vector = []
        angle = math.radians(angle)
        
        new_vector.append(math.cos(angle)*vector[0] - math.sin(angle)*vector[1])
        new_vector.append(math.sin(angle)*vector[0] + math.cos(angle)*vector[1])

        mag = (new_vector[0]**2 + new_vector[1]**2)**0.5
        
        new_vector[0] = new_vector[0]/mag
        new_vector[1] = new_vector[1]/mag

        return new_vector
    
    def distance_to_line(self, p_on_line, vector, P):
        
        distance = abs(vector[0]*(P[0]-p_on_line[0]) + vector[1]*(P[1]-p_on_line[1]))
        distance = distance / (vector[0]**2 + vector[1]**2)**0.5
        
        return distance
        
        
            
            
            
            


