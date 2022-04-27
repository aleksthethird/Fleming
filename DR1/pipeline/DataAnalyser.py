from astropy.table import Table
import Constants   
import os 
import Utilities
import matplotlib.pyplot as plt
import FluxFinder
import numpy as np
from PIL import Image
from astropy.visualization import (ZScaleInterval, LinearStretch,
                                   ImageNormalize)

class DataAnalyser:
    
    
    image_names = None
    filesdir = None
    means = []
    stds = []
    var_stds = []
    var_means = []
    id_map = []
    set_size = None
    n_sets = None
    has_sets = None
    avg_fluxes = []
    times = []
    light_curve_dir = None
    results_table = None
    variable_ids = []
    
    
    def __init__(self, filesdir, image_names, has_sets, set_size, n_sets):
        self.filesdir = filesdir
        self.image_names = image_names
        self.has_sets = has_sets
        self.set_size = set_size
        self.n_sets = n_sets
        self.light_curve_dir = filesdir + Constants.working_directory  + Constants.light_curve_directory

        
    #plot the means and standard deviations of all light curves generated
    def get_means_and_stds(self, adjusted):
        
        #build path of the directory in which the light curves are stored
        self.means = []
        self.stds = []
        self.id_map = []
        
        #cat = Table.read(self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + Constants.standard_file_extension, format=Constants.table_format)
        
        
        if not adjusted:
            light_curve_dir = self.filesdir + Constants.working_directory + Constants.light_curve_directory
        else:
            light_curve_dir = self.filesdir + Constants.working_directory + Constants.adjusted_curves_directory

        #for each file in the light curve directory 
        for file in os.listdir(light_curve_dir):
            
            if file[:len(self.image_names)] == self.image_names:
                
                #print(file)
                #read light curve data from file
                t = Table.read(light_curve_dir + file, format = Constants.table_format)
                
                id = file.split("id")[1].split(".")[0]

                if adjusted:
                    if self.remove_cosmics(t):
                        print(id)
                    
                                
                #only plot data point if at least 5 non-zero counts are recorded
                if len(t['counts']) > self.set_size*self.n_sets / 3:
                    
                    
                    mean = Utilities.mean(t['counts'])
                    

                    
                    std = Utilities.standard_deviation(t['counts'])
                    
                    value = std/mean
                                        
                    
                                                       
                    if value > 0 and value < 2: #and mean > 0.02 and mean < 80:
                        
                        self.stds.append(value)
                        self.means.append(mean)
                
                        self.id_map.append(int(id))
                
                    #if Utilities.is_above_line(std, mean, 17, 0, 0.05) and mean > 50:
                    #if value < 0.01 and mean > 3:
                        #print(file.split("id")[1].split(".")[0],Utilities.mean(t['counts']))
        Utilities.quicksort([self.means, self.stds, self.id_map], True)

        if adjusted:
            print("Floor at " + str(min(self.stds)))
        
    def plot_means_and_stds(self, adjusted):
       
        _ = plt.figure(figsize=(16, 10))

        _ = plt.scatter(self.means, self.stds, marker = '.');
        _ = plt.scatter(self.var_means, self.var_stds, marker = '.', color = 'red');
        _ = plt.xscale('log')
        _ = plt.yscale('log')
        _ = plt.xlabel("mean counts")
        _ = plt.ylabel("standard deviation")
        _ = plt.xlim(self.means[len(self.means)-1] * 1.3, self.means[0] * 0.7)
        
        #ensure plot y axis starts from 0
        _ = plt.gca().set_ylim(bottom=0)
        
        suf = ""
        if adjusted:
            suf = "_adjusted"
        _ = plt.savefig(Constants.folder + Constants.working_directory + Constants.output_directory + "SDvMean" + suf)
    
    #returns a score determining the variability of the star
    def get_variable_score(self, index):
        
        #the next part is just the same as the is_variable method
        #in Cataloguer - need to investigate if there is a reason for this
       
        values = []
        
        
        llim = index - Constants.check_radius
        
        if llim < 0:
            llim = 0
        
        ulim = index + Constants.check_radius + 1
        
        if ulim > len(self.means):
            ulim = len(self.means)
        
        for i in range(llim, ulim):
            if i != index:
                values.append(self.stds[i])
        
        #use median instead of mean
        
        median = np.median(values)
        
        if self.stds[index] > median * (1 + Constants.variability_threshold):
            
            self.var_means.append(self.means[index])
            self.var_stds.append(self.stds[index])
            
        return (self.stds[index] - median)/ median
        
        
    #seems to serve same function as same fn in cataloguer?
    #looks like this one is the one that is actually used
    #print plots of all variable stars in dataset, and stores their 
    #ids, x-y coords and variabilities in a table
    #comments in other
    def get_variables(self, adjusted):
        
        t = 0
        
        self.var_means = []
        self.var_stds = []
        
        cat = Table.read(self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + Constants.standard_file_extension, format=Constants.table_format)

        self.results_table = Table(names = ('id', 'xcentroid', 'ycentroid', 'variability', 'RA', 'DEC'))
        
        ff = FluxFinder.FluxFinder(Constants.folder, Constants.file_name, True, 9, 50)

        for i in range(len(self.means)):
            
            variability = self.get_variable_score(i)
            
            if variability > Constants.variability_threshold:
                
                t+= 1
                id = self.id_map[i]
                self.variable_ids.append(id)
                
                row_index = np.where(cat['id'].data==id)
                if 'RA' in cat.colnames:
                    
                    self.results_table.add_row([int(id), cat['xcentroid'][row_index], cat['ycentroid'][row_index], variability, cat['RA'][row_index], cat['DEC'][row_index]])
                else:
         
                    self.results_table.add_row([int(id), cat['xcentroid'][row_index], cat['ycentroid'][row_index], variability, 0, 0])
                #remove False literal here
                if adjusted:
                    ff.plot_light_curve(self.id_map[i], None, True)
        
            
    
    def create_thumbnails(self, ff):
        
        light_curve_dir = self.filesdir + Constants.working_directory + Constants.adjusted_curves_directory
        
        for i in range(len(self.results_table['id'])):

            file = light_curve_dir + Constants.file_name + "id" + str(int(self.results_table['id'][i])) + ".txt"
            t = Table.read(file, format = Constants.table_format)
            
            i_dim = 0
            i_bright = 0
            c = t['counts']
            
            for j in range(len(c)):
                if c[j] > c[i_bright]:
                    i_bright = j
                if c[j] < c[i_dim]:
                    i_dim = j
            
            i_x = self.results_table['xcentroid'][i]
            i_y = self.results_table['ycentroid'][i]
            
            print(i_x, i_y)
            
            print(self.results_table['id'][i])
            
            dim = ff.get_thumbnail(i_dim+1, i_x, i_y, 20, True)
            bright = ff.get_thumbnail(i_bright+1, i_x, i_y, 20, True)
            
            output_dir = self.filesdir + Constants.working_directory  + Constants.output_directory

        
            fig = plt.figure()
            fig.add_subplot(1, 2, 1)
            plt.axis('off')
            

            dim_norm = ImageNormalize(dim, interval=ZScaleInterval(), stretch=LinearStretch())
            plt.imshow(dim, origin='upper', cmap='gray', norm = dim_norm)

            
            fig.add_subplot(1, 2, 2)
            plt.axis('off')

            bright_norm = ImageNormalize(bright, interval=ZScaleInterval(), stretch=LinearStretch())
            plt.imshow(bright, origin='upper', cmap='gray', norm = bright_norm)
            
            plt.savefig(output_dir + "id_" + str(int(self.results_table['id'][i])) + ".jpg")




            
            
            

            
        
        
    #same function in cataloguer - remove one
    #probably from cataloguer
    #comments in other
    def get_ids_for_avg(self):
        
        ids = []
    
    
        for i in range(len(self.means)-1, int((len(self.means)-1) * 0.9), -1):
            #if self.means[i] > 50 and self.stds[i] < 0.04:
            #if not Utilities.is_above_line(-0.0001, 0.03, self.means[i], self.stds[i], 0.01) and self.means[i] > 5:
            #if not Utilities.is_above_line(self.means[i], self.stds[i], 2.2222*10**-9, 0.05777778, 0.001) and self.means[i] > 10^6:
            if not self.id_map[i] in self.variable_ids:
                light_curve_path = self.filesdir + Constants.working_directory + Constants.light_curve_directory + self.image_names + Constants.identifier + str(self.id_map[i]) + Constants.standard_file_extension

                t = Table.read(light_curve_path, format = Constants.table_format)
                
                if len(t['time']) == self.set_size * self.n_sets:
                    ids.append(self.id_map[i])
        return ids
    
    #duplicate method in ff - remove this?
    #comments in other
    def make_avg_curve(self, ids):
                
        for i in range(len(ids)):
            file = self.light_curve_dir + self.image_names + Constants.identifier + str(ids[i]) + Constants.standard_file_extension
            
            t = Table.read(file, format = Constants.table_format)
            
            fluxes = t['counts']
            
            if i == 0:
                for j in range(len(fluxes)):
                    self.avg_fluxes.append(0)
                self.times = t['time']

            
            mean = Utilities.mean(fluxes)
                        
            for k in range(len(fluxes)):
                self.avg_fluxes[k] = self.avg_fluxes[k] + (fluxes[k] / mean)
            
        for l in range(len(self.avg_fluxes)):
            self.avg_fluxes[l] = self.avg_fluxes[l] / len(ids)
        
        light_curve = Table([self.times, self.avg_fluxes], names = ('time','counts') )

        file = self.directory + Constants.working_directory + self.image_names + "_avg" + Constants.standard_file_extension
        light_curve.write(file, format = Constants.table_format, overwrite=True)

    #duplicate method in ff - remove this?
    #comments in other
    def divide_by_average(self):
        
        adjusted_light_curve_dir = self.filesdir + Constants.working_directory  + Constants.adjusted_curves_directory
        if not os.path.exists(adjusted_light_curve_dir):
            os.mkdir(adjusted_light_curve_dir)        
                    
        for file in os.listdir(self.light_curve_dir):
            if file[:len(self.image_names)] == self.image_names:
                
                t = Table.read(self.light_curve_dir + file, format = Constants.table_format)
                                                
                this_fluxes = t['counts']
                this_times = t['time']
                
                id = file.split("id")[1].split(".")[0]
                
                for i in range(len(this_fluxes)):
                    time = this_times[i]
                    
                    for j in range(len(self.avg_fluxes)):
                        if time == self.times[j]:
                            this_fluxes[i] = this_fluxes[i] / self.avg_fluxes[j]
                            #print(this_fluxes[i])
                
                light_curve = Table([this_times, this_fluxes], names = ('time','counts'))
    
                
                out_file = adjusted_light_curve_dir + self.image_names + Constants.identifier + str(id) + Constants.standard_file_extension
                light_curve.write(out_file, format = Constants.table_format, overwrite=True)

        
    #save table of variable stars (in order of decreasing variability)
    def output_results(self):
        
        output_dir = self.filesdir + Constants.working_directory  + Constants.output_directory
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)  
        
        a = [self.results_table['variability'], self.results_table['id'], self.results_table['xcentroid'], self.results_table['ycentroid'], self.results_table['RA'], self.results_table['DEC']]
        
        
        Utilities.quicksort(a, False)
        self.results_table.write(output_dir + self.image_names + "_results" + Constants.standard_file_extension, format = Constants.table_format, overwrite = True)
        Utilities.make_reg_file(output_dir, self.image_names + "_variables", self.results_table)

       
    def remove_cosmics(self, t):
        
        counts = t['counts']
        cosmic_index = -1
        
        if len(counts) == 0:
            return
        
        std = Utilities.standard_deviation(counts)

        
        for i in range(len(counts)):
            
            m = counts[i]
            
            if i == 0:
                l = counts[i+1]
            else:
                l = counts[i-1]
                
            if i == len(counts) - 1:
                r = counts[i-1]
            else:
                r = counts[i+1]
        
            if m - r > Constants.cosmic_threshold * std and m - l > Constants.cosmic_threshold * std:
                
                if cosmic_index != -1:
                    return False
              
                cosmic_index = i
            
        if cosmic_index == -1:
            return False
            
        print("cosmic detected at approx " + str(35*cosmic_index) + "s in id:")
        
        if i == 0:
                replacement = counts[1]
        else:
            replacement = counts[i-1]
                
        t['counts'][cosmic_index] = replacement
            
        return True
          
    