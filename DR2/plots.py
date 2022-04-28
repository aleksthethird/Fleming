from pipeline import Config, PeriodFinder, Utilities

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord 

import os
import shutil

hand_picked = []
hand_picked_unsure = []

def hand_pick():
    """
    Hand-pick variables after running the pipeline.

    Generate a catalogue of their amplitude, period and coordinates.

    """
    config = Config()

    ## Stars that are definitely variable from their light curves
    hand_picked = [
        ["l135",   1005, 1], ## Long period
        ["l135_5",  590, 0], ## Long period?
        ["l136_5",  580, 1], ## Long
        ["l197",    869, 1], ## Long
        ["l197",   1077, 1], ## Long
        ["l135_5",  148, 0], ## Transit
        ["l138_0", 1384, 0], ## Transit
        ["l135_5",  630, 0], ## Transit
        ["l197",    164, 0], ## Transit
        ["l140_5",  479, 0], ## Sharp transit
        ["l136_5",  576, 1], ## Transit
        ["l196",   1287, 0], ## Transit
        ["l196",   1416, 0], ## Weird transit
        ["l196",   1637, 0], ## Transit
        ["l197",   1844, 0], ## Transit
        ["l198",   2640, 0], ## Transit
        ["l136_5",  599, 2], ## Cool
        ["l141_5",  276, 2], ## Cool
        ["l197",     98, 4], ## Multi-period, appears irregular
        ["l138_0", 1194, 3], ## Multi-period
        ["l141_5",  301, 2], ## Multi-period
        ["l197",    371, 3], ## Multi-period
        ["l197",    456, 2], ## Multi-period
        ["l197",   1499, 2], ## Multi-period ?
        ["l197",   2118, 2], ## Irregular?
        ["l198",   2133, 1], ## Flare? Irregular?
        ["l198",   2268, 0], ## Scoop?
        ["l137_0", 1084, 1], ## Sine
        ["l140_0",  760, 1], ## Sine
        ["l141",    125, 1], ## Sharp bottom
        ["l141_5",  448, 1], ## Sine, same star as above
        ["l196",    507, 4], ## Sine, sharp bottom
        ["l197",   2374, 1], ## Sine, long period
        ["l198",    482, 1], ## Sine
        ["l198",   1593, 1], ## Sine
    ]

    ## Stars that might be variable (this was the appendix for my report)
    hand_picked_unsure = [
		["l135_5", 148, 0],
		["l135_5", 392, 0],
		["l135_5", 398, 0],
		["l135_5", 410, 0],
		["l135_5", 590, 0],
		["l135_5", 630, 0],
		["l135_5", 637, 0],
		["l135_5", 839, 0],
		["l135_5", 904, 0],
		["l135_5", 955, 0],
		["l135_5", 1000, 0],
		["l135_5", 1004, 0],
		["l136_5", 576, 0],
		["l136_5", 580, 0],
		["l137_0", 391, 0],
		["l137_0", 512, 0],
		["l137_0", 1039, 0],
		["l137_0", 686, 0],
		["l137_0", 1084, 0],
		["l137_5", 277, 0],
		["l137_5", 576, 0],
		["l137_5", 1034, 0],
		["l137_5", 1137, 0],
		["l137_5", 1171, 0],
		["l138_0", 1528, 0],
		["l140_0", 762, 0],
		["l140_5", 479, 0],
		["l141", 509, 0],
		["l196", 62, 0],
		["l196", 96, 0],
		["l196", 112, 0],
		["l196", 1080, 0],
		["l196", 1152, 0],
		["l196", 1279, 0],
		["l196", 1288, 0],
		["l196", 1416, 0],
		["l196", 1497, 0],
		["l196", 1637, 0],
		["l196_5", 828, 0],
		["l196_5", 930, 0],
		["l197", 444, 0],
		["l197", 479, 0],
		["l197", 841, 0],
		["l197", 869, 0],
		["l197", 957, 0],
		["l197", 1077, 0],
		["l197", 1118, 0],
		["l197", 1293, 0],
		["l197", 1499, 0],
		["l197", 2118, 0],
		["l197", 2140, 0],
		["l197", 2152, 0],
		["l197", 2334, 0],
		["l197", 2374, 0],
		["l197", 2462, 0],
		["l197", 2816, 0],
		["l197", 2856, 0],
		["l198_5", 1340, 0],
		["l198_5", 1411, 0],
		["l198", 745, 0],
		["l198", 790, 0],
		["l198", 1115, 0],
		["l198", 1589, 0],
		["l198", 2133, 0],
		["l198", 2268, 0],
		["l198", 2640, 0],
		["l199", 512, 0],
		["l199", 900, 0],
    ]

    n_variables = len(hand_picked)
    #n_variables = len(hand_picked_unsure)
    print("[HandPicked] Running for {} variables".format(n_variables))

    fname = "results_dr2.txt"
    final_results = open(os.path.join(config.hand_picked_dir, fname), "w")
    final_results.write(
            "RA [deg]\tDEC [deg]\tField id\tsource id\tperiod [h]\tperiod error [h]\tamplitude [%]\tamplitude error [%]\n")

    for field_id, source_id, n_periods in hand_picked:
    #for field_id, source_id, n_periods in hand_picked_unsure:
        config = Config(
            image_prefix = field_id,
        )
        pf = PeriodFinder(config)

        print("[Main] Hand picking for {}, {}".format(field_id, source_id))

        ## Grab the RA and DEC of the source
        cat = Utilities.read_catalogue(config)
        idx = np.where(cat['id'] == source_id)[0][0]
        ra  = cat['RA'][idx]
        dec = cat['DEC'][idx]

        ## Grab the adjusted light curve
        lc_path = os.path.join(config.adjusted_curve_dir,
                config.source_format_str.format(source_id))
        lc = np.genfromtxt(lc_path, dtype=config.light_curve_dtype)

        shutil.copyfile(lc_path, "./workspace/hand_picked/"+config.source_format_str.format(source_id))

        ## Plot the light curve
        plt.scatter(lc['time']/3600, lc['counts'], marker="x")

        ## Get the offset of the potential sine wave
        ps = pf.period_search_curve(
                source_id, lc['time'], lc['counts'], lc['counts_err'],
                n_samples=2000)
        offset = ps[-1]
        attempted_fit = np.ones(len(lc['counts'])) * offset

        try:
            coord = SkyCoord(ra*u.degree, dec*u.degree)
        except Exception as e:
            print("[Main] Couldn't find coordinates for star {}_id{}".format(field_id, source_id))
            continue

        if n_periods > 0:
            for _i in range(n_periods):
                _id, P, P_err, A, A_err, phi, offset = pf.period_search_curve(
                        source_id, lc['time'], lc['counts'], lc['counts_err'],
                        n_samples=2000)

                print("[Main] Fitting period {:5}s".format(P))

                ## Subtract the fit sine curve
                fit_sine = A*np.sin(2*np.pi/P * lc['time'] + phi)
                attempted_fit += fit_sine
                lc['counts'] -= fit_sine

                ## Write the period to the file
                final_results.write("{:.6f}\t{:.6f}\t{}\t{:04}\t{:1.5f}\t{:1.5f}\t{:02.4f}\t{:02.4f}\n"
                        .format(ra, dec, field_id, source_id, P/3600, P_err/3600, 100*A, 100*A_err))

            plt.title("Light curve and fitted plot for variable at\n{}"
                    .format(coord.to_string("hmsdms", precision=0)))
            plt.plot(lc['time']/3600, attempted_fit, color="red")

        else:
            ## If we don't have a period, just write coords
            final_results.write("{:.6f}\t{:.6f}\t{}\t{:04}\t{:.0f}\t{:4.0f}\t{:02.4f}\t{:02.4f}\n"
                    .format(ra, dec, field_id, source_id, 0, 0, 0, 0))

            plt.title("Light curve for variable at\n{}"
                    .format(coord.to_string("hmsdms", precision=0)))

        plt.xlabel("Time [hours]")
        plt.ylabel("Normalised flux [arb. u.]")

        ## Write the plot to disk
        fname = "final_{}_{}{:04}{}".format(config.image_prefix, config.identifier,
                source_id, config.plot_file_extension)
        plt.savefig(os.path.join(config.hand_picked_dir, fname))
        plt.close()


    final_results.close()


def group_plots():
    """
    Create plots of groups of variables.
    Very much not automated
    """

    ncols = 4
    fig, axs = plt.subplots(nrows=3,ncols=ncols, figsize=(16,9))
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.4, wspace=0.2)
    
    for i, (field_id, source_id, n_periods) in enumerate(hand_picked[36:48]):   
        config = Config(
            image_prefix = field_id,
        )
        
        ## Grab the RA and DEC of the source
        cat = Utilities.read_catalogue(config)
        idx = np.where(cat['id'] == source_id)[0][0]
        ra  = cat['RA'][idx]
        dec = cat['DEC'][idx]
        if ra > 360 or dec > 180:
            continue
        ## Grab the adjusted light curve
        lc_path = os.path.join(config.adjusted_curve_dir,
                config.source_format_str.format(source_id))
        lc = np.genfromtxt(lc_path, dtype=config.light_curve_dtype)

    
        ## Plot the light curve
        axs[i//ncols,i%ncols].scatter(lc['time']/3600, lc['counts'], marker="x", s=2)
    
        axs[i//ncols,i%ncols].set_xlabel("Time [hours]")
        axs[i//ncols,i%ncols].set_ylabel("Normalised flux [arb. u.]")
        coord = SkyCoord(ra*u.degree, dec*u.degree)
        axs[i//ncols,i%ncols].set_title("{}: {}"
                    .format(config.source_format_str.format(source_id)[:-4], 
                            coord.to_string("hmsdms", precision=0)))
    #fig.delaxes(axs[2,3])
    #fig.delaxes(axs[2,2])
    #fig.delaxes(axs[2,1])
    #fig.delaxes(axs[1,3])
    #fig.delaxes(axs[1,2])
    #fig.delaxes(axs[1,1])
    plt.savefig("hand_picked_unsure_4.pdf", bbox_inches="tight")


if __name__ == "__main__":
    hand_pick()
    #group_plots()
