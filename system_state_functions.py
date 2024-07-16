import numpy as np
from scipy.stats import linregress


# Additional static methods used by system_state object

def tuple_builder(property_list):
    if len(property_list) == 0:
        sd_tuple = (0.0, 0.0, 0.0)
    else:
        mean_property = np.mean(property_list)
        st_dev_property = np.std(property_list)
        sd_tuple = (mean_property - st_dev_property, mean_property, mean_property + st_dev_property)
    return sd_tuple


def linear_model_report(x_val, y_val, is_record_vectors, model_type_str=None):
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
        "is_linear_model":  True,  # this is for search algorithms to identify that this dictionary is a LM
        "is_success": is_lm_success,
        "slope": lm_slope,
        "intercept": lm_intercept,
        "r": lm_r,
        "p": lm_p,
        "std_err": lm_std_err,
        "x_lim": np.asarray([np.min(x_val), np.max(x_val)]),
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


def generate_cluster(sub_network, size):
    num_attempts = 0
    max_attempts = 100
    adjacency_matrix = sub_network["adjacency_array"]
    # try to generate and return a connected cluster of the given size from the provided sub_network

    is_success = False
    cluster = []
    # is it possible in principle?
    if size <= sub_network["num_patches"]:
        while num_attempts < max_attempts and not is_success:
            is_success = True
            num_attempts += 1
            cluster = []  # will store the row-column indices relative to the current sub_network
            eligible_patch_indices = [x for x in range(sub_network["num_patches"])]

            # initial
            draw_num = np.random.choice(eligible_patch_indices)
            cluster.append(draw_num)
            eligible_patch_indices.remove(draw_num)

            # attempt to draw an element connected to existing elements
            if size > 1:
                for num_element in range(size - 1):

                    possible_draw = set()
                    # find the neighbours of the current cluster members
                    for potential_element in eligible_patch_indices:
                        for existing_member in cluster:
                            if adjacency_matrix[existing_member, potential_element] == 1 or \
                                    adjacency_matrix[potential_element, existing_member] == 1:
                                possible_draw.add(potential_element)
                                break
                    draw_list = list(possible_draw)
                    if len(draw_list) > 0:
                        draw_num = np.random.choice(draw_list)
                        cluster.append(draw_num)
                        eligible_patch_indices.remove(draw_num)
                    else:
                        is_success = False
                        cluster = []
                        break
    return cluster, is_success


def determine_complexity(sub_network, cluster, is_normalised):
    # sum the absolute state difference values over all unique patch pairs in the cluster
    total_difference = 0
    if sub_network["num_patches"] > 1:
        for lower_patch_index in range(len(cluster) - 1):
            lower_patch_num = cluster[lower_patch_index]
            for upper_patch_index in range(lower_patch_index + 1, len(cluster)):
                upper_patch_num = cluster[upper_patch_index]
                # compare pair
                if is_normalised:
                    # for populations - difference should be |x_i - x_j| / bar{x}
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
