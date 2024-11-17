import numpy as np
from scipy.stats import linregress, spearmanr, pearsonr
from copy import deepcopy

# Additional static methods used by system_state object

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


def complexity_scaling_vector_analysis(x_val, y_val, is_record_lm_vectors, non_log_x_val):
    # executed during system_state.complexity_analysis() separately
    # for _binary and _population_weighted vectors to analyse the scaling and dimensions of complexity,
    # and attempt to estimate the natural occurrences of complexity organisation within the system.

    lm_complexity = {"is_success": 0}
    natural_cluster_size = 0
    dc_estimate_list = []
    dc_natural_range_est = []
    dc_significance = 0
    dc_significance_un_norm = 0

    # ignore first entry (cluster size of one patch - used for diversity, but no meaning for intra-cluster complexity)
    if not (y_val[1:] == 0).any():
        # shift up to a minimum value of 1
        testable_complexity = y_val[1:]
        bin_comp_y_val = np.log(testable_complexity - np.min(testable_complexity) + 1.0)
        lm_complexity = linear_model_report(x_val=x_val[1:],
                                            y_val=bin_comp_y_val,
                                            is_record_vectors=is_record_lm_vectors,
                                            model_type_str="log-log")
        # complexity dimension is +gradient
        lm_complexity["complexity_dimension"] = lm_complexity["slope"]

        # fit log-log plots across rolling intervals
        interval_length = 5
        # num_increments = int(np.floor_divide(len(x_val) - 2, interval_length))
        max_lm_slope = -float('inf')
        dc_estimate_list = []

        # fit dc over cluster size [2, 2+N], [3, 3+N], [4, 4+N], ...
        for k in range(len(x_val) - 1 - interval_length):
            x_start = k + 1
            x_end = x_start + interval_length
            log_n = x_val[x_start:x_end]
            log_c = bin_comp_y_val[x_start:x_end]
            lm_slope = linregress(x=log_n, y=log_c).slope
            dc_estimate_list.append(lm_slope)
            # we try to identify the interval where the rate of growth is greatest
            if lm_slope > max_lm_slope:
                max_lm_slope = lm_slope
                dc_natural_range_est = [int(x_start + 1), int(x_end)]
        # how significant is this "greatest" dc estimate (i.e. how many s.d. above the mean)?
        if np.sum(dc_estimate_list) != 0.0:
            dc_significance = (max_lm_slope - np.mean(dc_estimate_list))/np.std(dc_estimate_list)
        else:
            dc_significance = 0
        # possibly it should NOT be normalised as dc = 2.0 is the universal approximate behaviour
        dc_significance_un_norm = max_lm_slope - np.mean(dc_estimate_list)

    # Now attempt to identify the natural cluster size in the system:
    for x_index in range(1, len(x_val) - 1):  # because of the [_ + 1]
        test_n = non_log_x_val[x_index]
        c_now = y_val[x_index]
        c_next = y_val[x_index + 1]
        c_test = c_now + (0.5 * (test_n + 1)) ** 2 - (0.5 * test_n) ** 2
        if c_next > c_test:
            natural_cluster_size = int(deepcopy(test_n))
            break

    return (lm_complexity, dc_estimate_list, dc_natural_range_est, dc_significance_un_norm,
            dc_significance, natural_cluster_size)


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
