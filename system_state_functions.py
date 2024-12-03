import numpy as np
from scipy.stats import spearmanr, pearsonr, linregress
from scipy.optimize import curve_fit

# Additional static functions used by system_state methods

def tuple_builder(property_list):
    if len(property_list) == 0:
        sd_tuple = (0.0, 0.0, 0.0)
    else:
        mean_property = np.mean(property_list)
        st_dev_property = np.std(property_list)
        sd_tuple = (mean_property - st_dev_property, mean_property, mean_property + st_dev_property)
    return sd_tuple


def linear_model_report(x_val, y_val, is_record_vectors, model_type_str=None, is_shifted=False):
    # checks for typical causes of error, then if valid conducts a linear regression and returns all results in a
    # dictionary structure
    lm_slope, lm_intercept, lm_r, lm_p, lm_std_err = [0, 0, 0, 0, 0]
    is_lm_success = 0
    if len(x_val) == len(y_val) > 1 and not (x_val == float('inf')).any() and not (x_val == float('-inf')).any() and \
            not (y_val == float('inf')).any() and not (y_val == float('-inf')).any():
        try:
            lm_slope, lm_intercept, lm_r, lm_p, lm_std_err = linregress(x_val, y_val)
            is_lm_success = 1
        except ValueError:
            pass
    linear_model = {
        "is_linear_model": True,  # this is for search algorithms to identify that this dictionary is a LM
        "is_success": is_lm_success,
        "slope": lm_slope,
        "intercept": lm_intercept,
        "r": lm_r,
        "p": lm_p,
        "std_err": lm_std_err,
        "x_lim": np.asarray([np.min(x_val), np.max(x_val)]),
        "is_shifted": is_shifted  # this is merely a record of whether the raw data had to be shifted up/right
        # to avoid zeros before taking logs, for log-lin or log-log models
        #
        # be aware that, to print anything here in the final "format_dictionary_to_JSON_string()" it needs to contain
        # only dictionaries and Numpy arrays, and NOT LISTS
    }
    if is_record_vectors:
        linear_model["x_val"] = x_val
        linear_model["y_val"] = y_val
    if model_type_str is not None:
        # pass in a string indicating underlying model, from ["lin-lin", "lin-log", "log-log", "log-lin"]
        linear_model["model_type_str"] = model_type_str
    return linear_model


def single_power_law(n, alpha, dc):
    return alpha * (n - 1) ** dc

def dual_power_law(n, k, alpha_1, alpha_2, dc_1, dc_2, beta):
    return np.exp(-k*(n-1)) * alpha_1 * (n - 1) ** dc_1 + (1 - np.exp(-k*(n-1)))*(alpha_2 * (n - 1) ** dc_2 + beta)

def get_r_squared(x, y, func, fitted_para):
    # calculate the R-squared statistic
    residuals = y - func(x, *fitted_para)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2.0)
    if ss_tot != 0.0:
        r_squared = 1.0 - (ss_res / ss_tot)
    else:
        r_squared = 0.0
    return r_squared

def complexity_scaling_vector_analysis(x_val, y_val):
    # executed during system_state.complexity_analysis() separately
    # for _binary and _population_weighted vectors to analyse the scaling and dimensions of complexity,
    # and attempt to estimate the natural occurrences of complexity organisation within the system.
    #
    # Include first entry (cluster size of one patch has zero complexity), but because of "n-1" in power laws, we must
    # offset this with a small correction term for curve fitting algorithm.
    # Check where y_val next contains zero - then slice off the invalid part if necessary.
    y_zero_index = np.where(y_val[1:] == 0)[0]
    if len(y_zero_index) == 0:
        testable_complexity = y_val  # [C(1) = 0, C(2), ... ]
        testable_x_val = x_val  # [1, 2, 3, ... ]
    else:
        testable_complexity = y_val[0:y_zero_index[0] + 1]
        testable_x_val = x_val[0:y_zero_index[0] + 1]
    # add the offset
    testable_x_val[0] = testable_x_val[0] + 0.0000001

    # overall analysis
    if len(testable_complexity) > 2:

        # fit a single power law if more than two data points
        dual_initial_guess = [0.1, 1, 1, 2, 2, 1]
        try:
            single_para = curve_fit(single_power_law, testable_x_val, testable_complexity, p0=[1, 2],
                                    bounds=((0.0, 0.0),(np.inf, 10.0)))[0]
            is_single_fit_success = 1
            dual_initial_guess = [0.1, single_para[0], single_para[0], single_para[1], single_para[1], 0.0]
            # complexity dimension is power
            single_dc_fitted = single_para[1]
            single_r_squared = get_r_squared(x=testable_x_val, y=testable_complexity,
                                             func=single_power_law, fitted_para=single_para)
        except (Warning, RuntimeError, TypeError):
            is_single_fit_success = 0
            single_dc_fitted = None
            single_para = None
            single_r_squared = None

        if len(testable_x_val) > 6:
            # fit a dual power law if more than six data points
            try:
                # parameters are: k, alpha_1, alpha_2, dc_1, dc_2, beta
                dual_para = curve_fit(dual_power_law, testable_x_val, testable_complexity, p0=dual_initial_guess,
                                      bounds=((0.0000001, 0.0, 0.0, 0.0, 0.0, -10.0),
                                              (10.0, np.inf, np.inf, 10.0, 10.0, 10.0)))[0]
                is_dual_fit_success = 1
                dual_r_squared = get_r_squared(x=testable_x_val, y=testable_complexity,
                                               func=dual_power_law, fitted_para=dual_para)
            except (Warning, RuntimeError, TypeError):
                is_dual_fit_success = 0
                dual_para = None
                dual_r_squared = None
        else:
            is_dual_fit_success = 0
            dual_para = None
            dual_r_squared = None

        # record spectrum of complexity change per n at each value of n
        delta_complexity_per_n = np.zeros(len(testable_complexity))
        for j in range(len(delta_complexity_per_n)-1):
            # from j=2 to end-1 determine the per patch change over [j, j+1]
            delta_complexity_per_n[j] = (testable_complexity[j+1] - testable_complexity[j])/(j+1)
        # record first location of maximum gain per n (as possible upper bound)
        dc_per_n_max_loc = int(np.where(delta_complexity_per_n == np.max(delta_complexity_per_n))[0][0])

        complexity_max_results = {
            "dc_largest_gain": delta_complexity_per_n[dc_per_n_max_loc],
            "dc_largest_gain_n": dc_per_n_max_loc + 1,
        }

        # fit incremental power laws and determine spectrum of d_c over interval of [1, M] for M in [3, end]
        dc_fit = np.zeros(len(testable_complexity))
        for k in range(3, len(dc_fit)):
            # fit dc for C = alpha*(n - 1)^dc over [1, k+1] so at least 3 points
            try:
                dc_fit_para = curve_fit(single_power_law, testable_x_val[0:k], testable_complexity[0:k], p0=[1, 2],
                                        bounds=((0.0, 0.0),(np.inf, 10.0)))[0]
                dc_fit[k] = dc_fit_para[1]
            except (Warning, RuntimeError, TypeError):
                pass
        # record first location of maximum fitted d_c (as possible upper bound)
        dc_fit_max_loc = int(np.where(dc_fit == np.max(dc_fit))[0][0])
        complexity_max_results["dc_largest_fit"] = dc_fit[dc_fit_max_loc]
        complexity_max_results["dc_largest_fit_n"] = dc_fit_max_loc + 1
        # if there is a point where fitted d_c is first > 2, save that (as possible lower bound)
        dc_fit_two_plus = np.where(dc_fit > 2)[0]
        if len(dc_fit_two_plus) > 0:
            dc_fit_two_plus_loc = int(dc_fit_two_plus[0])
        else:
            dc_fit_two_plus_loc = 0
        complexity_max_results["interval_est"] = [max(0, dc_fit_two_plus_loc + 1),
                                                  min(dc_per_n_max_loc + 1, dc_fit_max_loc + 1)]

        graphical_results ={
            # for producing bespoke plots
            "is_complexity_graphical": True,  # this is for search algorithms to identify this dictionary
            "n_vector": [int(x) for x in testable_x_val],  # need to convert int64 elements to int for final JSON saving
            "complexity_vector": list(testable_complexity),
            "delta_c_per_n": list(delta_complexity_per_n),
            "dc_fit": list(dc_fit),
            "is_single_fit_success": is_single_fit_success,
            "single_para": single_para,
            "single_r_squared": single_r_squared,
            "is_dual_fit_success": is_dual_fit_success,
            "dual_para": dual_para,
            "dual_r_squared": dual_r_squared,
        }
    else:
        single_dc_fitted = None
        graphical_results = {}
        complexity_max_results = {}
    return single_dc_fitted, graphical_results, complexity_max_results

def determine_complexity(sub_network, cluster, is_normalised, num_species):
    # sum the absolute state difference values over all unique patch pairs in the cluster
    total_difference = 0
    if num_species > 0:
        if sub_network["num_patches"] > 1:
            for lower_patch_index in range(len(cluster) - 1):
                lower_patch_num = cluster[lower_patch_index]
                for upper_patch_index in range(lower_patch_index + 1, len(cluster)):
                    upper_patch_num = cluster[upper_patch_index]
                    # compare pair
                    if is_normalised:
                        # for populations - difference should be |x_i - x_j| / max{x}
                        lower_vector = sub_network["normalised_population_arrays"][0][lower_patch_num, :]
                        upper_vector = sub_network["normalised_population_arrays"][0][upper_patch_num, :]
                    else:
                        # for binary (occupancy) comparison - difference should be 1 or 0
                        lower_vector = sub_network["population_arrays"][0][lower_patch_num, :]
                        upper_vector = sub_network["population_arrays"][0][upper_patch_num, :]
                    difference = 0
                    for species_index in range(len(lower_vector)):
                        difference += np.abs(lower_vector[species_index] - upper_vector[species_index])
                    total_difference += difference
            # report the total complexity (rather than the per-patch or per-pair average, as this is then
            # compared with log(cluster size) in the subsequent information dimension calculation)

        # however we will, for niceness, normalise by the number of species (in the entire system - not just
        # this subnetwork). But because we take logs and then only consider the fitted gradient for complexity
        # dimension, this will actually have no impact.
        total_difference = total_difference / num_species
    return total_difference


def rank_abundance(sub_networks, is_record_lm_vectors):
    # for each sub_network, conduct a species rank abundance analysis, fitting species rank (0 - N-1, where N is the
    # number of species with non-zero population in this sub_network) against log(relative abundance). As the y-values
    # are normalised, we expect y-intercept of the fit to be -1, so the y-intercept of the original curve is at (0,1).
    # Since the abundances are arranged in decreasing order, we expect the gradient to be negative.
    rank_abundance_report = {}
    for network_key in sub_networks.keys():
        if sub_networks[network_key]["num_patches"] > 1:
            population_array = sub_networks[network_key]["population_arrays"][0]  # NOT species-normalised!
            num_species = np.shape(population_array)[1]

            abundance_list = []
            for species_index in range(num_species):
                species_abundance = np.sum(population_array[:, species_index])
                abundance_list.append(species_abundance)
            abundance_list.sort(reverse=True)
            abundance_vector = np.asarray(abundance_list)
            x_val = np.arange(num_species)

            # remove zero-abundance species
            is_pass = False
            while len(x_val) > 0 and not is_pass:
                is_pass = True
                for relative_species_index in range(len(x_val)):
                    if abundance_vector[relative_species_index] == 0:
                        is_pass = False
                        # remove this entry from both vectors
                        x_val = np.delete(x_val, relative_species_index, axis=0)
                        abundance_vector = np.delete(abundance_vector, relative_species_index, axis=0)
                        break

            # normalise abundance to the maximum global (in sub_network) abundance across all species, then
            # take the logarithm
            if len(x_val) > 0:
                abundance_vector = np.log(abundance_vector / max(abundance_vector))

            # then conduct the linear regression if valid
            if len(x_val) > 0:
                rank_abundance_report[network_key] = linear_model_report(x_val=x_val, y_val=abundance_vector,
                                                                         is_record_vectors=is_record_lm_vectors,
                                                                         model_type_str="log-lin")
            else:
                rank_abundance_report[network_key] = {"is_success": 0}
        else:
            rank_abundance_report[network_key] = {"is_success": 0}

    return rank_abundance_report


def inter_species_predictions_correlation_coefficients(species_1_vector, species_2_vector):
    # must first test for 'nearly constant' vectors to avoid warnings
    species_1_near_constant = np.var(species_1_vector) < 1e-13 * abs(np.mean(species_1_vector))
    species_2_near_constant = np.var(species_2_vector) < 1e-13 * abs(np.mean(species_2_vector))
    is_cc_auto_fail = (species_1_near_constant or species_2_near_constant or np.var(
        species_1_vector) <= 0.0 or np.var(species_2_vector) <= 0.0 or
                       len(species_1_vector) < 2 or len(species_2_vector) < 2)
    is_cc_success = 0
    pearson_cc, pearson_p, spearman_rho, spearman_p = [0.0, 0.0, 0.0, 0.0]
    # for correlation coefficients to work, we need two vectors of at least two pairs and at least some variance
    if not is_cc_auto_fail:
        try:
            pearson_cc, pearson_p = pearsonr(species_1_vector, species_2_vector)
            spearman_rho, spearman_p = spearmanr(species_1_vector, species_2_vector)
            is_cc_success = 1
            if pearson_p == float('NaN') or spearman_p == float('NaN'):
                is_cc_success = 0
        except ValueError:
            is_cc_success = 0
    if is_cc_success == 0:
        pearson_cc, pearson_p, spearman_rho, spearman_p = [0.0, 0.0, 0.0, 0.0]
    corr_coefficients = {
        'is_success': is_cc_success,
        'pearson_cc': pearson_cc,
        'pearson_p_value': pearson_p,
        'spearman_rho': spearman_rho,
        'spearman_p_value': spearman_p,
    }
    return corr_coefficients
