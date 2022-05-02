from . import PeriodFinder, Utilities, Config

import numpy as np

import matplotlib.pyplot as plt

class VariableDetector:
    """
    VariableDetector class.

    Does everything to do with determining variability in a source.

    Usually associated with a DataAnalyser object.

    Attributes
    ----------


    config: Config
        Config object for the field

    n_sources: int
        Number of sources the object knows about

    source_ids: numpy array
        IDs of all sources we know about.

    means: numpy array
        Means of the light curves of each source.

    medians: numpy array
        Medians of the light curves of each source.

    stds: numpy array
        Standard deviations of the light curves of each source.

    n_positive: numpy array
        Number of positive values in each light curve

    adjusted: bool
        Are we an object for adjusted or non-adjusted set?

    variable_scores: numpy array
        Scores for standard deviation search
    
    amplitude_scores: numpy array
        Scores for period search

    period_stats: numpy array
        Table of primary period, amplitude, phase, offset and errors for each source.

    variable_ids_s: numpy array
        IDs of all sources deemed variable through the standard deviation check

    variable_ids_a: numpy array
        IDs of all sources deemed variable through the amplitude search

    """

    config = None

    n_sources = 0
    source_ids = None

    means = None
    stds  = None
    medians = None
    n_positive = None

    variable_scores = None
    amplitude_scores = None
    period_stats = None

    variable_ids_s = None
    variable_ids_a = None

    adjusted = False

    def __init__(self, config, source_ids, means, stds, medians, n_positive, adjusted=True):
        self.config = config

        self.n_sources = len(source_ids)
        self.source_ids = source_ids

        self.means = means
        self.stds = stds
        self.medians = medians
        self.n_positive = n_positive

        self.adjusted = adjusted


    def get_variable_score(self, source_id):
        """
        Returns a score determining the variability of the star.

        Standard deviation search.
        
        Score is the standard deviation of the source compared
        to the median of the surrounding stars (in index-space)

        Parameters
        ----------

        source_id: int
            ID of the source to get the score of

        Returns
        -------

        variable_score: int
            Score of the variable

        """
        
        source_index = np.where(self.source_ids == source_id)[0]

        ## Check if we found our source in the IDs
        if len(source_index) > 0:
            source_index = source_index[0]
        else:
            return -1.0
       
        llim = source_index - self.config.check_radius
        ulim = source_index + self.config.check_radius + 1
        
        if llim < 0:
            llim = 0
        
        if ulim >= self.n_sources:
            ulim = self.n_sources-1


        ## Median of standard deviations
        median_std = np.median(self.stds[llim:ulim])
        variable_score = self.stds[source_index]/median_std - 1
        #print("[DEBUG] Source {:04}; {}/{}".format(source_id, self.stds[source_index], median_std))

        return variable_score
        
        
    ## TODO: move signal to noise calculation out of here
    ## TODO: Make quality checks optional (for average )
    def std_dev_search(self, threshold, threshold_upper, min_snr, min_positives):
        """
        Calculates the variability score of all stars in the catalogue
        then builds a list of source indexes which are deemed variable.

        Parameters
        ----------
        
        threshold: float64
            Threshold score, above which is deemed variable

        Returns
        -------

        self.variable_ids_s: numpy array
            Numpy array of all ids which are deemed to be variable candidates.

        """

        self.variable_scores = np.zeros(self.n_sources)

        ## Loop over sources in the catalogue
        for i, path, source_id in Utilities.loop_variables( \
                self.config, self.source_ids, adjusted=self.adjusted):

            self.variable_scores[i] = self.get_variable_score(source_id)
            #lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype)

        signal_to_noise = self.medians/self.stds
        variable_mask = np.where(
                (self.variable_scores > threshold)
                & (self.variable_scores < threshold_upper)
                & (signal_to_noise > min_snr)
                & (self.n_positive > min_positives)
                )[0]

        self.variable_mask = variable_mask


        if len(variable_mask) > 0:
            variable_ids = np.copy(self.source_ids[self.variable_mask])
        else:
            variable_ids = np.array([], dtype='int')

        self.variable_ids_s = variable_ids
        

        print("[VariableDetector] std search: Found {} variables out of {} sources"
                .format(len(variable_ids), self.n_sources))

        return variable_ids


    def amplitude_search(self, amplitude_score_threshold):
        """
        Finds the most prominent period for each source.
        Reports the period, amplitude, phase, offset and uncertainties.

        Parameters
        ----------

        amplitude_score_threshold: float64
            Threshold score, above which is deemed variable

        Returns
        -------

        period_stats: numpy array
            Table of the statistics for each source

        """

        pf = PeriodFinder(self.config)

        amplitude_score = np.zeros(self.n_sources)

        period_stats = np.zeros(self.n_sources, dtype=[
            ('id', 'int64'),
            ('period', 'float64'),
            ('period_err', 'float64'),
            ('amplitude', 'float64'),
            ('amplitude_err', 'float64'),
            ('phi', 'float64'),
            ('offset', 'float64'),
            ])

        for i, path, source_id in Utilities.loop_variables( \
                self.config, self.source_ids, adjusted=self.adjusted):

            print("[VariableDetector] Finding period for source {}/{}".format(i, self.n_sources))

            lc = np.genfromtxt(path, dtype=self.config.light_curve_dtype).transpose()

            ## TODO: Make period range config variable
            period_stats[i] = pf.period_search_curve(source_id, 
                    lc['time'], lc['counts'], lc['counts_err'], n_samples=self.config.n_sample_periods)
            A = period_stats[i]['amplitude']
            A_err = period_stats[i]['amplitude_err']

            ## TODO: Use signal to noise or amplitude per std?
            #amplitude_score[i] = A/A_err
            amplitude_score[i] = A/self.stds[i]

            ## TODO: Make config variable
            if A > 5:
                amplitude_score[i] = 0.0

            #if period_stats['period'][i] > 0:
            #    main_period = A*np.sin(2*np.pi/period_stats['period'][i] * lc['time']
            #            + period_stats['phi'][i]) + period_stats['offset'][i]
            #    new_std = np.std(lc['counts'] - main_period)
            #    amplitude_score[i] = self.stds[i]/new_std

            #else:
            #    amplitude_score[i] = 0

            #if period_stats[i]['period'] > 0:
            #    pf.plot_fit(source_id,
            #        lc['time'], lc['counts'],
            #        A, period_stats[i]['period'],
            #        period_stats[i]['phi'], period_stats[i]['offset'],
            #        plot_dir = self.config.testing_dir)

        self.period_stats = period_stats
        self.amplitude_scores = amplitude_score

        variable_mask = np.where(amplitude_score > amplitude_score_threshold)[0]

        if len(variable_mask) == 1:
            variable_mask = variable_mask[0]
            variable_ids = np.array([self.source_ids[variable_mask]])
        else:
            variable_ids = np.copy(self.source_ids[variable_mask])

        self.variable_ids_a = variable_ids

        print("[VariableDetector] amplitude search: Found {} variables out of {} sources"
                .format(len(variable_ids), self.n_sources))

        return variable_ids


    def get_variables(self):
        """
        Do all searches to find sources which are variable.
        Recommended to call each search yourself, this is just
        here for convenience.

        Returns
        -------

        variable_ids: numpy array
            IDs of each source considered variable

        """
        ## Get variables judged by different methods
        variable_ids_s = self.std_dev_search(self.config.variability_threshold)
        variable_ids_a = self.amplitude_search()

        ## Stars which apppear in both are returned
        #variable_ids = np.intersect1d(variable_ids_s, variable_ids_a)

        ## Stars which were in one or the other are returned
        ## TODO: May contain duplicate IDs
        variable_ids = np.concatenate(variable_ids_s, variable_ids_a)

        return variable_ids

    def get_period_stats(self, ids=[]):
        """
        Retrieve the period stats for the given IDs

        Parameters
        ----------

        ids: array, optional
            IDs to get stats of. If empty, gives all known stats

        Returns
        -------

        period_stats: numpy array
            Table of the statistics for each source

        """
        if len(ids) == 0:
            return self.period_stats
        else:
            ## Sorts IDs in numerical order
            #intersect, indices, _indices2 = np.intersect1d(self.period_stats['id'], ids,
            intersect, indices, _indices2 = np.intersect1d(self.source_ids, ids,
                return_indices=True, assume_unique=True)

            ## Need to 'unsort' to get bacl to  order of ids
            unsort_indices = np.argsort(indices)
            indices = indices[unsort_indices]

            if len(indices) != len(ids):
                print("[VariableDetector] Error: Couldn't find all IDs in variable detector")
            return self.period_stats[indices]


    def get_scores(self, ids=[]):
        """
        Retrieve the scores of each given source.

        Parameters
        ----------

        ids: array, optional
            IDs to get stats of. If empty, gives all known scores

        Returns
        -------

        period_stats: numpy array
            Table of the statistics for each source

        """

        if len(ids) == 0:
            return self.variable_scores, self.amplitude_scores
        else:
            _intersect, indices, _indices2 = np.intersect1d(self.source_ids, ids,
                return_indices=True, assume_unique=True)

            ## Need to 'unsort' to get bacl to  order of ids
            unsort_indices = np.argsort(indices)
            indices = indices[unsort_indices]

            return self.variable_scores[indices], self.amplitude_scores[indices]

    ## TODO: Docs
    ## TODO: add check for zero-period
    def filter_variables(self, variable_ids):
        """
        Remove variables which all have similar periods, amplitudes and phases.
        This is likely due to field effects or something.
        """

        ## Get period stats in order of variable ids
        period_stats = self.get_period_stats(variable_ids)

        ## TODO: Make config variables
        nbins_period = 5
        nbins_phi = 10
        nbins_amp = 40

        hist_period, edge_period = np.histogram(period_stats['period'], bins=nbins_period)
        largest = np.argmax(hist_period)
        period_lower = edge_period[largest]
        period_upper = edge_period[largest+1]
        
        hist_phi, edge_phi = np.histogram(period_stats['phi'], bins=nbins_phi)
        largest = np.argmax(hist_phi)
        phi_lower = edge_phi[largest]
        phi_upper = edge_phi[largest+1]

        hist_amp, edge_amp = np.histogram(period_stats['amplitude'], bins=nbins_amp)
        largest = np.argmax(hist_amp)
        amp_lower = edge_amp[largest]
        amp_upper = edge_amp[largest+1]

        remove_indices = []
        for i, _path, source_id in Utilities.loop_variables(self.config, variable_ids):
            similar_period = period_stats['period'][i] > period_lower \
                    and period_stats['period'][i] < period_upper
            similar_phi = period_stats['phi'][i] > phi_lower \
                    and period_stats['phi'][i] < phi_upper
            similar_amp = period_stats['amplitude'][i] > amp_lower \
                    and period_stats['amplitude'][i] < amp_upper

            ## Only similar amplitude and phase since 'fit' period can change a lot
            if similar_amp and similar_phi:
            #if similar_amp and similar_phi and similar_period:
                remove_indices.append(i)

        print("[VariableDetector] Filtered out {}/{} similar light curves"
                .format(len(remove_indices), len(variable_ids)))
        filtered_variable_ids = np.delete(variable_ids, remove_indices)
        return filtered_variable_ids


