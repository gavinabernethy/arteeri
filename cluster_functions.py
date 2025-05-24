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
    # This is a more limited, but more efficient method to generate entire clusters of the specified size.
    # Used during complexity_analysis() when we requiring rapidly drawing many 100's of clusters.
    cluster = []
    is_success = False
    num_attempts = 0
    internal_complexity = 0.0
    neighbour_dict = sub_network["neighbour_dict"]
    if is_uniform:
        if is_normalised:
            target_array = sub_network["normalised_population_arrays"][0]
        else:
            target_array = sub_network["population_arrays"][0]
    else:
        target_array = None

    # is it possible in principle?
    if size <= sub_network["num_patches"]:
        while num_attempts < max_attempts and not is_success:
            is_success = True
            num_attempts += 1

            # initialise lists
            cluster = []  # will store the row-column indices relative to the current sub_network
            potential_neighbours = []
            # initialise arrays of length equal to the possible_neighbours list
            adjacency_counter = np.array([])
            difference_sum = np.array([])
            internal_complexity = 0.0

            # loop through drawing the 1st to Nth elements of the sample
            for num_element in range(size):

                if len(potential_neighbours) > 0:
                    # for all except the first element, draw from list according to cluster criteria in arrays
                    #
                    # at the most base level, all neighbours of the cluster are eligible
                    base_score = np.ones(len(potential_neighbours))

                    # Modifiers
                    if is_box:
                        # if applicable, FIRST restrict to only the greatest adjacency
                        adj_maximiser = (adjacency_counter == np.max(adjacency_counter))
                        base_score = base_score * adj_maximiser
                    else:
                        adj_maximiser = np.ones_like(difference_sum)
                    if is_uniform:
                        # to choose best score GIVEN the best adjacency!
                        favour_score = adj_maximiser * (1.0 + np.max(difference_sum) - difference_sum)
                        # if applicable, THEN restrict to the highest score for intra-cluster uniformity
                        base_score = base_score * (favour_score == np.max(favour_score))

                    # after possible modifications, draw one of the best
                    short_list = np.where(base_score == np.max(base_score))[0]
                    draw_index = np.random.choice(short_list)
                    draw_num = potential_neighbours[draw_index]
                    cluster.append(draw_num)
                    potential_neighbours.pop(draw_index)
                    adjacency_counter = np.delete(adjacency_counter, draw_index)
                    if is_uniform:
                        # update internal complexity with difference between the newly-added element and all current
                        internal_complexity += difference_sum[draw_index]
                    difference_sum = np.delete(difference_sum, draw_index)
                else:
                    if len(cluster) == 0:
                        # draw of initial element
                        if initial_patch is not None:
                            if all_elements_admissible or initial_patch in admissible_elements:
                                draw_num = initial_patch
                            else:
                                raise Exception("Initial element not admissible.")
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
                    # don't consider patches already chosen!
                    if potential_element not in cluster:
                        if not all_elements_admissible:
                            if potential_element not in admissible_elements:
                                continue

                        # separate treatment required only for those who are NEW eligible neighbours
                        if potential_element not in potential_neighbours:
                            potential_neighbours.append(potential_element)
                            # increase size of all tracking arrays
                            adjacency_counter = np.pad(adjacency_counter, (0,1), mode='constant', constant_values=0)
                            difference_sum = np.pad(difference_sum, (0, 1), mode='constant', constant_values=0)
                            potential_index = len(potential_neighbours) - 1

                            # calculate difference of new neighbour to all except the newest member of the cluster
                            if is_uniform and len(cluster) > 1:
                                for cluster_member in cluster[0:-1]:
                                    difference_sum[potential_index] += determine_difference(
                                        cluster_member, potential_element, target_array, num_species)
                        else:
                            # find the existing index for this neighbour
                            potential_index = potential_neighbours.index(potential_element)

                        if is_box:
                            # for box style, must update all neighbours of the added patch with their additional
                            # adjacency, which may go from 0 to 1 (for new neighbours) or simply be incremented
                            adjacency_counter[potential_index] += 1

                # now if necessary, update ALL neighbours (including new neighbours) with ADDITIONAL difference to the
                # newly-added member of the cluster
                if is_uniform:
                    for potential_index, potential_neighbour in enumerate(potential_neighbours):
                        difference_sum[potential_index] += determine_difference(
                            cluster[-1], potential_neighbour, target_array, num_species)

    # normalise internal complexity
    if is_uniform:
        internal_complexity = internal_complexity * 4.0 / (num_species * (size ** 2 - np.mod(size, 2)))
    return cluster, is_success, internal_complexity

def determine_difference(patch_1, patch_2, target_array, num_species):
    total_difference = 0.0
    for species_index in range(num_species):
        total_difference += np.abs(target_array[patch_1, species_index] - target_array[patch_2, species_index])
    return total_difference

def draw_partition(sub_network, size, initial_patch, num_species, is_normalised):
    # As far as possible, cluster the elements of the network into highly-uniform clusters of the given size.
    # Note that we opt for box-clusters only as this is less ambiguous in what we 'expect' to see, and it feels like
    # a more natural interpretation of the space than the visually-strange but topologically-admissible patterns that
    # *could* get detected otherwise - i.e. it seeks more obvious approximately-square 2D clusters, rather than
    # recognising two clusters joined by a long thin string as a single 'cluster' even though topologically equivalent.
    #
    # partition consists of a numbered dictionary of cluster lists
    partition = {}

    elements_to_partition = [_ for _ in range(sub_network["num_patches"])]
    partition_lookup = np.zeros(sub_network["num_patches"])  # returns cluster of the patch
    cluster_init_patch = initial_patch
    cluster_num = 0
    total_elements_partitioned = 0
    partition_internal_complexity = []

    while len(elements_to_partition) > 0:
        # generate each box cluster
        cluster, is_success, internal_complexity = generate_fast_cluster(
            sub_network=sub_network,
            size=size,
            max_attempts=1,
            admissible_elements=elements_to_partition,
            num_species=num_species,
            is_box=True,  # less regular systems sensitive to this - depends on the topological restrictions desired
            is_uniform=True,
            all_elements_admissible=False,
            initial_patch=cluster_init_patch,
            return_undersized_cluster=True,
            is_normalised=is_normalised,
        )
        total_elements_partitioned += size * is_success  # how many patches put into full-size clusters so far?
        if is_success:
            # only record internal complexity for the properly partitioned elements - i.e. so that we can give a
            # conservative estimate of how much of the system can be partitioned for a given delta
            partition_internal_complexity.append(internal_complexity)

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

    # did we partition more than half the patches into clusters of the required size?
    is_partition_success = (total_elements_partitioned > np.floor(len(elements_to_partition)/2))
    return partition, partition_lookup, is_partition_success, partition_internal_complexity

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

    # now determine min, mean, max species-averaged difference across adjacent clusters (i.e. partition complexity)
    total_difference = 0.0
    min_difference = float('inf')
    max_difference = 0.0
    num_pairs = 0
    for cluster_1 in range(num_clusters-1):
        for cluster_2 in range(cluster_1+1, num_clusters):
            if partition_adjacency[cluster_1, cluster_2] != 0:
                num_pairs += 1
                pair_difference = 0.0
                for species_index in range(num_species):
                    pair_difference += np.abs(partition_values[cluster_1][species_index
                                              ] - partition_values[cluster_2][species_index])
                spec_ave_pair_difference = pair_difference / num_species
                total_difference += spec_ave_pair_difference
                min_difference = min(min_difference, spec_ave_pair_difference)
                max_difference = max(max_difference, spec_ave_pair_difference)
    if num_pairs > 0:
        mean_difference = total_difference / num_pairs
    else:
        mean_difference = 0.0
    return min_difference, mean_difference, max_difference