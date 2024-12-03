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


def generate_box_cluster(sub_network, size, max_attempts, initial_patch=None):
    # This is a more limited, but more efficient method to generate entire box-style clusters of the
    # specified size. Used during complexity_analysis() when we requiring rapidly drawing many 100's of clusters.

    num_attempts = 0
    neighbour_dict = sub_network["neighbour_dict"]
    # try to generate and return a connected box-cluster of the given size from the provided sub_network

    is_success = False
    cluster = []

    # is it possible in principle?
    if size <= sub_network["num_patches"]:
        while num_attempts < max_attempts and not is_success:
            is_success = True
            num_attempts += 1
            cluster = []  # will store the row-column indices relative to the current sub_network

            # initial element
            if initial_patch is not None:
                draw_num = initial_patch
            else:
                draw_num = random.randrange(sub_network["num_patches"])
            cluster.append(draw_num)

            # attempt to draw an element connected to existing elements
            if size > 1:

                # initialise the probability array
                favour_array = np.zeros(sub_network["num_patches"])
                # start with the neighbours of the first element
                possible_draw = list(neighbour_dict[draw_num])
                # discard any that are not eligible
                for candidate_element in possible_draw:
                    if candidate_element != draw_num:
                        favour_array[candidate_element] = 1

                # loop through drawing the 2nd to Nth elements of the sample
                for num_element in range(size - 1):

                    if np.sum(favour_array) > 0:
                        short_list = np.where(favour_array == np.max(favour_array))[0]
                        draw_num = np.random.choice(short_list)
                        cluster.append(draw_num)
                        favour_array[draw_num] = 0
                    else:
                        is_success = False
                        cluster = []
                        break

                    # check the neighbours of the new member and update the set of possibilities
                    for potential_element in neighbour_dict[draw_num]:
                        if potential_element not in cluster:
                            # don't consider patches already chosen!
                            for potential_element_neighbour in neighbour_dict[potential_element]:
                                # iterate over neighbours of one element rather than the existing cluster, as we
                                # expect < 10 neighbours per cluster, but the clusters *can* be potentially very large
                                if potential_element_neighbour in cluster:
                                    favour_array[potential_element] += 1

    return cluster, is_success