from . import Utilities, Config

import os
import numpy as np
import matplotlib.pyplot as plt

class PeriodFinder:

    def __init__(self, config):
        self.config = config


    def period_search(self, source_id, path,
            period_min = 0.5*3600,  # 0.5 hours
            period_max = 8.0*3600,  # 6 hours
            n_samples  = int(2e3),  # 2000 samples
            ):
        """
        Wrapper function for `period_search_curve`.
        Use if you are using an unmodified light curve.

        Parameters
        ----------

        source_id: int
            ID of the source to be analysed.

        path: string
            File path of the light curve to read.

        period_min, period_max: float
            Minimum and maximum period range to search in.

        n_samples: int
            Number of periods to check between the bounds.

        Returns
        -------

        period_min, period_min_err: float
            Most prominent period present, as well as its error.

        amplitude, amplitude_err: float
            Amplitude of most prominent period and its error.

        phi: float
            Phase offset

        offset: float
            Height offset

        """

        curve = np.genfromtxt(path, dtype = self.config.light_curve_dtype).transpose()
        time = curve['time']
        counts = curve['counts']
        errors = curve['counts_err']

        return self.period_search_curve(source_id, time, counts, errors,
                period_min, period_max, n_samples)



    def period_search_curve(self, source_id,
            time, counts, errors,
            #period_min = 0.5*3600,
            period_min = 20.0*60,
            period_max = 5.0*3600,
            n_samples  = 2000
            ):
        """
        Takes a source light curve and finds the most prominent period
        in the curve.

        Default period search range is half an hour to 6 hours at a sample
        rate of 2000 periods in between.

        Uses a model of 
        B + C cos(omega t) + S sin(omega t)
        which is equivalent to
        B + A sin(2 pi/P t + phi)

        "Most prominent period" means period which minimises chi squared.

        Parameters
        ----------

        source_id: int
            ID of the source to be analysed.

        time: numpy array
            Times of each data point

        counts: numpy array
            Our data point. Relative flux of star at that time.

        errors: numpy array
            The uncertainty in the `counts` measurement.

        period_min, period_max: float
            Minimum and maximum period range to search in.

        n_samples: int
            Number of periods to check between the bounds.


        Returns
        -------

        period_min, period_min_err: float
            Most prominent period present, as well as its error.

        amplitude, amplitude_err: float
            Amplitude of most prominent period and its error.

        phi: float
            Phase offset

        offset: float
            Height offset

        """
        
        
        ## Convert periods to angular frequencies
        periods = np.linspace(period_min, period_max, int(n_samples))
        omegas = 2*np.pi/periods


        ## Initial search, get chi squared landscape
        _params, chi2 = self.search_omegas(counts, time, errors, omegas)
        idx_min = np.argmin(chi2)
        omega_min = omegas[idx_min]
        chi2_min  = chi2[idx_min]

        ## If our minimum is on the edge of our range
        ## ie there is no clear period within the range
        if idx_min == n_samples-1 or idx_min == 0:
            #return period_min, period_min_err, amplitude, amplitude_err, phi, B
            ## Return unphysical values
            #print("[PeriodFinder] Error: Minimum chi2 is initially outside bounds for source {:04}"
            #        .format(source_id))
            self.plot_chi2(chi2, omegas, source_id, 0.0)
            return source_id, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0


        ## What's our changes like near the minimum?
        omega_next = omegas[idx_min+1]
        chi2_next  = chi2[idx_min+1]

        ## Save original landscape to plot later
        chi2_orig = chi2
        omegas_orig = omegas



        ## Iterate until we can resolve the minimum of the landscape well
        ## we're aiming for a total change of chi from min to max ~1
        ## Set a max number of iterations just to be safe.
        delta = 10*np.abs(omega_next - omega_min)
        iteration = 1
        while (np.max(chi2) - chi2_min) > self.config.period_chi2_range \
              and iteration < self.config.period_max_iterations:

            #print("[PeriodFinder] Iteration: {:02} for id {:04}".format(iteration, source_id))

            ## NOTE: Original idea, seems to be the best
            ### Find an approximate gradient to guestimate width of next omega range
            approx_gradient = np.sqrt(chi2_next - chi2_min)/(omega_next - omega_min)

            ## Estimation of width we should use to better resolve the minimum
            approx_halfwidth = 1/approx_gradient

            low_bound  = omega_min - self.config.period_width_adjustment*approx_halfwidth
            high_bound = omega_min + self.config.period_width_adjustment*approx_halfwidth

            ## Create new omega range
            omegas = np.linspace(low_bound, high_bound, int(n_samples))

            ## Next approximation, (hopefully) much closer to the minimum
            ## Need to get close enough to be able to approximate to a parabola
            _params, chi2 = self.search_omegas(counts, time, errors, omegas)
            #self.plot_chi2(chi2, omegas, source_id, iteration)
            idx_min = np.argmin(chi2)

            ## If we break here, we cropped too small
            if idx_min == n_samples-1 or idx_min == 0:
                print("[PeriodFinder] Error: Minimum chi2 moved outside of region for source {:04}"
                    .format(source_id))
                self.plot_chi2(chi2_orig, omegas_orig, source_id, 0.0)
                return source_id, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0

            omega_min = omegas[idx_min]
            chi2_min  = chi2[idx_min]

            omega_next = omegas[idx_min+1]
            chi2_next  = chi2[idx_min+1]

            iteration += 1

        if iteration >= self.config.period_max_iterations:
                print("[PeriodFinder] Error: Hit max iterations for source {:04}"
                    .format(source_id))


        ## Assume we are good enough

        is_forwards = True
        large_chi2 = np.where(chi2[idx_min:] > chi2_min + 1)[0]
        #large_chi2 = np.where(chi2 > chi2_min + 1)[0]

        ## If our delta chi2 never goes above 1
        ## ie we have no good minimum
        ## NOTE: Errors are determined here
        if len(large_chi2) == 0:

            ## TODO: Maybe best to always only look backwards. Need testing.
            ## Need to look backwards if the minimum is uneven and has
            ## a shallower gradient in the +omega direction

            ## If we can't find a forward point, look backwards
            large_chi2 = np.where(chi2[:idx_min] > chi2_min + 1)[0]
            is_forwards = False

            if len(large_chi2) == 0:
                ## Return unphysical values
                print("[PeriodFinder] Debug: chi2 landscape too shallow for source {:04}"
                        .format(source_id))
                self.plot_chi2(chi2_orig, omegas_orig, source_id, 0.0)
                self.plot_chi2(chi2, omegas, source_id, 1.0)
                return source_id, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0

        if is_forwards:
            idx_dchi2 = large_chi2[0]
            idx_dchi2 += idx_min
        else:
            idx_dchi2 = large_chi2[-1]
        
        ## Uncertainty of omega is 1/sqrt(curvature at minimum)
        c = (chi2[idx_dchi2] - chi2_min)/(omegas[idx_dchi2] - omega_min)**2
        omega_min_err = np.sqrt(1/c)

        ## Convert from angular frequency to period
        period_min = 2*np.pi/omega_min
        period_min_err = (omega_min_err/omega_min) * period_min

        ## Amplitude estimation and error
        ## Use inverse Hessian to get errors
        params, _H, H_inv = self.periodogram_fit_func(counts, time, errors, omega_min)
        B, C, S = params
        C_err = np.sqrt(H_inv[1,1])
        S_err = np.sqrt(H_inv[2,2])
        A = np.sqrt(C**2 + S**2)

        ## Amplitude and its error
        amplitude_err = np.sqrt( ((C/A)*C_err)**2 + ((S/A)*S_err)**2 )
        amplitude = A

        ## Phase shift for a sin wave
        phi = np.arctan(C/S)

        ## Correct for trig reasons
        if S < 0:
            phi = phi - np.pi

        ## Plot initial chi2 distribution
        if self.config.plot_fit_chi2:
            self.plot_chi2(chi2_orig, omegas_orig, source_id, period_min)

        if self.config.plot_fit_comparison:
            self.plot_fit(time, counts, A, P, phi, offset)

        return source_id, period_min, period_min_err, amplitude, amplitude_err, phi, B



    ## TODO: Docs
    def plot_fit(self, source_id, time_p, counts, A, P_p, phi, offset,
            plot_dir=None, show=False, title=None, units="seconds"):

        if units == "hours":
            time = time_p/3600
            P = P_p/3600
            fname = "compare_{}_{}{:04}_P{:1.4f}h{}".format(
                    self.config.image_prefix,
                    self.config.identifier,
                    source_id,
                    P,
                    self.config.plot_file_extension
                )

        elif units=="seconds":
            time = time_p
            P = P_p
            fname = "compare_{}_{}{:04}_P{:05.0f}s{}".format(
                    self.config.image_prefix,
                    self.config.identifier,
                    source_id,
                    P,
                    self.config.plot_file_extension
                )
        else:
            print("[PeriodFinder] Error in plotting, units '{}' not recognised".format(units))
            return

        ## Sine curve of found period
        if P > 0:
            attempted_fit = offset + A*np.sin(2*np.pi/P * time + phi)
        else:
            attempted_fit = np.ones(len(time)) * np.median(counts)

        ## Plot light curve and overplot sine wave
        plt.scatter(time, counts, marker="x")
        plt.plot(time, attempted_fit, color="red")


        if title == None:
            if units == "seconds":
                plt.title("Overplot period of {:05.0f}s for source id{:04}".format(P, source_id))
            elif units == "hours":
                plt.title("Overplot period of {:1.4f}hours for source id{:04}".format(P, source_id))
        else:
            plt.title(title)

        plt.xlabel("Time [{}]".format(units))
        plt.ylabel("Normalised brightness [arb. u.]")
        if plot_dir == None:
            plt.savefig(os.path.join(self.config.periods_dir, fname))
        else:
            plt.savefig(os.path.join(plot_dir, fname))

        if show:
            plt.show()
            
        plt.close()

    ## TODO: Docs
    def plot_chi2(self, chi2, omegas, source_id, period_min):
        """
        """

        plt.plot(2*np.pi/omegas, chi2)
        fname = "chi2_{}_{}{:04}_P{:05.0f}{}".format(
                self.config.image_prefix,
                self.config.identifier,
                source_id,
                period_min,
                self.config.plot_file_extension
            )
        plt.title("$\chi^2$ plot for source {:04}, period {:5.0f}s".format(source_id, period_min))
        plt.xlabel("Period [s]")
        plt.ylabel("$\chi^2$ [arb. u.]")
        plt.savefig(os.path.join(self.config.periods_dir, fname))
        plt.close()



    def search_omegas(self, counts, time, errors, omegas):
        """
        Loop over range of angular frequencies and get chi squared
        landscape and best fit parameters for each.

        Parameters
        ----------

        counts: numpy array
            Relative flux of source at each point in time

        time: numpy array
            Time elapsed since initial measurement

        errors: numpy array
            Uncertainty on the relative flux

        omegas: numpy array
            Range of angular frequencies to search

        Returns
        -------

        params: numpy array
            Array of the B, C and S parameters for each fit

        chi2: numpy array
            Chi squared value for each angular frequency

        """
        n_samples = len(omegas)
        chi2 = np.zeros(n_samples)
        params = np.zeros((n_samples, 3))

        for i,omega in enumerate(omegas):
            params[i], _H, _H_inv = self.periodogram_fit_func(counts, time, errors, omega)
            chi2[i] = self.periodogram_chi2_(counts, time, errors, omega, params[i])

        return params, chi2
            
        


    def periodogram_fit_func(self, data, time, error, omega):
        """
        Calculates the best fit parameters for a single 
        angular frequency for a given light curve.

        Uses the Hessian to immediately find minimum chi squared.

        Concept from AS5001 ADA.

        Parameters
        ----------

        data: numpy array
            Relative flux of source at each point in time

        time: numpy array
            Time elapsed since initial measurement

        error: numpy array
            Uncertainty on the relative flux

        omega: float64
            Angular frequency of fit.

        Returns
        -------

        fit_params: tuple
            Tuple of the B, C and S parameters

        H: numpy matrix
            Hessian matrix for fit

        H_inv: numpy matrix
            Inverse Hessian matrix for fit

        """

        s2r = 1/error**2
        s = np.sin(omega*time)
        c = np.cos(omega*time)
        
        H = np.array([
            [np.sum(s2r), np.sum(c*s2r), np.sum(s*s2r)],
            [np.sum(c*s2r), np.sum(c*c*s2r), np.sum(c*s*s2r)],
            [np.sum(s*s2r), np.sum(c*s*s2r), np.sum(s*s*s2r)]
        ])
        
        corr = np.array([
            [np.sum(data*s2r)],
            [np.sum(data*c*s2r)],
            [np.sum(data*s*s2r)]
        ])
        
        H_inv = np.linalg.inv(H)
        fit_params = np.reshape(np.matmul(H_inv, corr), 3)
        return fit_params, H, H_inv


    def periodogram_chi2(self, data, time, error, omega, B, C, S):
        """
        Calculates the chi squared value for a set of fit parameters
        """
        F = lambda t: B + C*np.cos(omega*t) + S*np.sin(omega*t)
        chi2 = ((data-F(time))/error)**2
        chi2 = np.sum(chi2)
        return chi2
    
    
    def periodogram_chi2_(self, data, time, error, omega, params):
        """
        Wrapper function for `periodogram_chi2`
        """
        return self.periodogram_chi2(data, time, error, omega, params[0], params[1], params[2])


