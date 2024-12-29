import random
import numpy as np


def cluster_next_element(adjacency_matrix, patch_list, current_cluster: list,
                         actual_patch_nums: list, cluster_arch_type: str):
    # This method offers much more control on how to specify the topology of clusters to be generated.
    # It is used in:
    # - Perturbations: called only by cluster_builder()
    # - Sample_spatial_data: called by generate_habitat_type() for cluster-based variants only when specified.
    #
    # Returns the list of patches eligible to be drawn as the next element of the cluster being iteratively constructed
    # in the list current_cluster.

    type_patch_nums = []
    if cluster_arch_type == "random":
        type_patch_nums = actual_patch_nums

    elif cluster_arch_type == "star":
        # try to choose a neighbour of the initial node
        for patch_num in actual_patch_nums:
            if adjacency_matrix[current_cluster[0], patch_num] > 0 \
                    or adjacency_matrix[patch_num, current_cluster[0]] > 0:
                type_patch_nums.append(patch_num)

    elif cluster_arch_type == "chain":
        # try to choose a neighbour of the final node
        for patch_num in actual_patch_nums:
            if adjacency_matrix[current_cluster[-1], patch_num] > 0 \
                    or adjacency_matrix[patch_num, current_cluster[-1]] > 0:
                type_patch_nums.append(patch_num)

    elif cluster_arch_type == "disconnected":
        # try to choose a neighbour of None of the current nodes
        for patch_num in actual_patch_nums:
            for already_chosen_num in current_cluster:
                if adjacency_matrix[already_chosen_num, patch_num] == 0 \
                        and adjacency_matrix[patch_num, already_chosen_num] == 0:
                    type_patch_nums.append(patch_num)

    elif cluster_arch_type == "box":
        # try to choose a neighbour of multiple current nodes
        adjacency_counter = []
        for patch_num in actual_patch_nums:
            num_adjacent = 0
            for already_chosen_num in current_cluster:
                if adjacency_matrix[already_chosen_num, patch_num] > 0 \
                        or adjacency_matrix[patch_num, already_chosen_num] > 0:
                    num_adjacent += 1
            adjacency_counter.append([patch_num, num_adjacent])
        # what was the greatest adjacency?
        if len(adjacency_counter) > 0:
            max_adjacency = max([x[1] for x in adjacency_counter])
            # choose from those with greatest adjacency only
            type_patch_nums = [x[0] for x in adjacency_counter if x[1] == max_adjacency]

    elif cluster_arch_type == "position_box":
        # choose a box consisting of positionally-close patches, regardless of actual connectivity
        distance_counter = []
        for patch_num in actual_patch_nums:
            summative_distance = 0.0
            for already_chosen_num in current_cluster:
                summative_distance += np.linalg.norm(patch_list[patch_num].position
                                                     - patch_list[already_chosen_num].position)
            distance_counter.append([patch_num, summative_distance])
        # what was minimum total positional distance of any acceptable patch to the existing members of the cluster?
        min_distance = min([x[1] for x in distance_counter])
        # choose from those with the least distance only
        type_patch_nums = [x[0] for x in distance_counter if x[1] == min_distance]

    else:
        raise "Cluster arch-type not recognised."
    return type_patch_nums


def generate_fast_cluster(sub_network, size, max_attempts, admissible_elements, num_species, is_box=False,
                          is_uniform=False, all_elements_admissible=False, initial_patch=None,
                          return_undersized_cluster=False, is_normalised=False):
    # This is a more limited, but more efficient method to generate entire  clusters of the specified size.
    # Used during complexity_analysis() when we requiring rapidly drawing many 100's of clusters.

    num_attempts = 0
    neighbour_dict = sub_network["neighbour_dict"]
    if is_uniform:
        if is_normalised:
            target_array = sub_network["normalised_population_arrays"][0]
        else:
            target_array = sub_network["population_arrays"][0]
    else:
        target_array = None
    # try to generate and return a connected box-cluster of the given size from the provided sub_network

    is_success = False
    cluster = []

    # is it possible in principle?
    if size <= sub_network["num_patches"]:
        while num_attempts < max_attempts and not is_success:
            is_success = True
            num_attempts += 1
            cluster = []  # will store the row-column indices relative to the current sub_network
            # initialise the preference array
            favour_array = np.zeros(sub_network["num_patches"])
            base_array = np.zeros(sub_network["num_patches"])  # only needed if BOTH is_box AND is_uniform

            # loop through drawing the 1st to Nth elements of the sample
            for num_element in range(size):

                if np.sum(favour_array) > 0:
                    # for all except the first element, draw according to prepared favour arrays for cluster criteria
                    if is_box and is_uniform:
                        # if both, first need to restrict favour_array to where bass_array was maximal
                        bass_maximiser = np.where(base_array == np.max(base_array))[0]
                        favour_array = favour_array * bass_maximiser
                    short_list = np.where(favour_array == np.max(favour_array))[0]
                    draw_num = np.random.choice(short_list)
                    cluster.append(draw_num)
                    favour_array[draw_num] = 0
                else:
                    if len(cluster) == 0:
                        # draw of initial element
                        if initial_patch is not None:
                            draw_num = initial_patch
                        else:
                            draw_num = random.choice(admissible_elements)
                        cluster.append(draw_num)
                    else:
                        # cluster has failed to attain required size
                        is_success = False
                        if not return_undersized_cluster:
                            # reset failed cluster to empty string if small cluster not acceptable
                            cluster = []
                        break

                # check the neighbours of new member (applies to first element also) and update set of possibilities
                for potential_element in neighbour_dict[draw_num]:
                    if not all_elements_admissible:
                        if potential_element not in admissible_elements:
                            continue

                    if potential_element in admissible_elements and potential_element not in cluster:
                        # don't consider patches already chosen!

                        if is_box:
                            # for box style, we only choose from the most connected
                            for potential_element_neighbour in neighbour_dict[potential_element]:
                                # iterate over neighbours of one element rather than the existing cluster, as we
                                # expect < 10 neighbours per cluster, but the clusters *can* be potentially large
                                if potential_element_neighbour in cluster:
                                    if is_uniform:
                                        # we will choose the least distance of the most-connected elements
                                        base_array[potential_element] += 1
                                        # eligibility is 1/(sum of distance to each member of the cluster)
                                        favour_array[potential_element] = determine_difference(
                                            cluster, target_array, num_species, potential_element)

                                    else:
                                        # eligibility is the number of adjacent current members of the cluster
                                        favour_array[potential_element] += 1

                        else:
                            if is_uniform:
                                # eligibility is 1/(sum of distance to each member of the cluster)
                                favour_array[potential_element] = determine_difference(
                                    cluster, target_array, num_species, potential_element)
                            else:
                                # consider eligible as adjacent to a current member of the cluster
                                favour_array[potential_element] = 1

    return cluster, is_success

def determine_difference(cluster, target_array, num_species, potential_element):
    total_difference = 0.0
    for current_element in cluster:
        for species_index in range(num_species):
            total_difference += np.abs(target_array[current_element,
            species_index] - target_array[potential_element, species_index])
    if total_difference == 0.0:
        return_difference = float('inf')
    else:
        return_difference = 1.0 / total_difference
    return return_difference


def draw_partition(sub_network, size, initial_patch, num_species, is_normalised):
    # As far as possible, cluster the elements of the network into box-clusters of the given size (or smaller)
    #
    # partition consists of a numbered dictionary of cluster lists
    partition = {}

    elements_to_partition = [_ for _ in range(sub_network["num_patches"])]
    partition_lookup = np.zeros(sub_network["num_patches"])  # returns cluster of the patch
    cluster_init_patch = initial_patch
    cluster_num = 0

    while len(elements_to_partition) > 0:
        # generate each box cluster
        cluster = generate_fast_cluster(
            sub_network=sub_network,
            size=size,
            max_attempts=1,
            admissible_elements=elements_to_partition,
            num_species=num_species,
            is_box=False,
            is_uniform=True,
            all_elements_admissible=False,
            initial_patch=cluster_init_patch,
            return_undersized_cluster=True,
            is_normalised=is_normalised,
        )[0]

        for element in cluster:
            elements_to_partition.remove(element)
            partition_lookup[element] = cluster_num
        partition[cluster_num] = cluster

        if len(elements_to_partition) > 0:
            # set up for next cluster
            cluster_num += 1
            # choose next largest element from the smallest in the cluster to start with if possible,
            # otherwise choose the smallest overall from those eligible
            temp_init_upper = float('inf')
            temp_init_lower = elements_to_partition[0]
            min_cluster_element = min(cluster)
            for element in elements_to_partition:
                temp_init_lower = min(temp_init_lower, element)
                if element > min_cluster_element:
                    temp_init_upper = min(temp_init_upper, element)
            if temp_init_upper < float('inf'):
                cluster_init_patch = int(temp_init_upper)
            else:
                cluster_init_patch = temp_init_lower

    return partition, partition_lookup

def partition_analysis(sub_network, partition, partition_lookup, num_species, is_normalised):
    # Determine the per-species patch-mean value in each cluster, and the adjacency relationships between clusters,
    # then calculate the mean difference across neighbouring clusters.

    if is_normalised:
        target_array_str = "normalised_population_arrays"
    else:
        target_array_str = "population_arrays"

    num_clusters = len(partition)
    partition_values = {}
    partition_adjacency = np.zeros((num_clusters, num_clusters))
    for cluster_num in range(num_clusters):
        cluster = partition[cluster_num]
        cluster_species_values = {x: 0.0 for x in range(num_species)}
        for species_index in range(num_species):
            cluster_species_values[species_index] = np.mean([sub_network[target_array_str][0][
                                                                _, species_index] for _ in cluster])
        partition_values[cluster_num] = cluster_species_values

        # determine cluster adjacency in the partition
        for element in cluster:
            for neighbour in sub_network["neighbour_dict"][element]:
                neighbour_cluster = int(partition_lookup[neighbour])
                if neighbour_cluster != cluster_num:
                    partition_adjacency[cluster_num, neighbour_cluster] = 1

    # now determine species-averaged difference across adjacent clusters
    total_difference = 0.0
    num_pairs = 0
    for cluster_1 in range(num_clusters-1):
        for cluster_2 in range(cluster_1+1, num_clusters):
            if partition_adjacency[cluster_1, cluster_2] != 0:
                num_pairs += 1
                pair_difference = 0.0
                for species_index in range(num_species):
                    pair_difference += np.abs(partition_values[cluster_1][species_index
                                              ] - partition_values[cluster_2][species_index])
                total_difference += pair_difference/num_species
    if num_pairs > 0:
        total_difference /= num_pairs
    else:
        total_difference = 0.0
    return total_difference