import Reducer
import ShiftFinder
import FluxFinder
import DataAnalyser
import Constants
import Utilities
import MovingObjectFinder
import StreakFinder
import Cataloguer

from datetime import datetime

now = datetime.now()
wcs = None

current_time = now.strftime("%H:%M:%S")
print("Started at " + current_time)


# =============================================================================
#r = Reducer.Reducer(Constants.folder,"No filter",Constants.file_name, 9, 50)
#r.reduce(True)
# # #set_size, n_sets = r.get_set_info()
 
Utilities.print_job("reducing images")
c = Cataloguer.Cataloguer(Constants.folder, Constants.file_name, Constants.has_sets, Constants.set_size, Constants.n_sets)
c.catalogue()
 
# # #mof = MovingObjectFinder.MovingObjectFinder(Constants.folder)
# # #mof.find_moving_objects()
# 
# 
Utilities.print_job("cataloguing stars")
# 
sf = ShiftFinder.ShiftFinder(Constants.folder, Constants.file_name, Constants.has_sets, Constants.set_size, Constants.n_sets)
sf.get_all_shifts()
# 
Utilities.print_job("finding shifts")
# 
# =============================================================================
ff = FluxFinder.FluxFinder(Constants.folder, Constants.file_name, True, 9, 50)
ff.find_all_fluxes()
ff.make_light_curves()

#strf = StreakFinder.StreakFinder(Constants.folder, c, ff)
#strf.find_all_streaks()

Utilities.print_job("making light curves")

r = None
c = None
sf = None

da = DataAnalyser.DataAnalyser(Constants.folder, Constants.file_name, True, 9, 50)

da.get_means_and_stds(False)
da.get_variables(False)
da.output_results() #this should not be here - make it just create results folder with separate method
da.plot_means_and_stds(False)
ids = da.get_ids_for_avg()

ff.make_avg_curve(ids)
ff.divide_by_average()
ff.plot_light_curve(None, Constants.folder + Constants.working_directory + Constants.file_name + "_avg.txt", True)

Utilities.print_job("adjusting light curves")

da.get_means_and_stds(True)
da.get_variables(True)
da.plot_means_and_stds(True)
da.output_results()

da.create_thumbnails(ff)

Utilities.print_job("everything")




