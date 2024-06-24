import numpy as np
from copy import deepcopy
from collections import Counter
from degree_distribution import power_law_curve_fit


def tuple_builder(property_list):
    if len(property_list) == 0:
        sd_tuple = (0.0, 0.0, 0.0)
    else:
        mean_property = np.mean(property_list)
        st_dev_property = np.std(property_list)
        sd_tuple = (mean_property - st_dev_property, mean_property, mean_property + st_dev_property)
    return sd_tuple


class System_state:

    def __init__(self, patch_list, species_set, parameters, step=0,
                 patch_adjacency_matrix=None,
                 habitat_type_dictionary=None,
                 habitat_species_traversal=None,
                 habitat_species_feeding=None,
                 current_patch_list=None,
                 dimensions=None,
                 ):
        self.step = step
        self.time = 0
        self.patch_list = patch_list
        self.initial_patch_list = None
        self.species_set = species_set
        self.current_patch_list = current_patch_list
        self.dimensions = dimensions
        self.reserve_list = []
        self.perturbation_history = {}
        self.perturbation_holding = None

        # this requires that the ordering of the species list and of the columns in the external files are identical!
        self.habitat_type_dictionary = habitat_type_dictionary
        self.habitat_species_feeding = habitat_species_feeding
        self.habitat_species_traversal = habitat_species_traversal
        self.patch_adjacency_matrix = patch_adjacency_matrix
        self.initial_patch_adjacency_matrix = None
        # Update all patches
        self.update_all_patches_habitat_based_properties()

        # histories
        # number of patches and the biodiversity are checked every time-step
        self.current_num_patches_history = []
        self.global_biodiversity_history = []
        self.local_biodiversity_history = []
        # other network properties only have changes recorded at initialisation and after a relevant perturbation
        self.habitat_amounts_history = {_: {} for _ in habitat_type_dictionary}  # dict of dict (of time/values)
        self.habitat_auto_correlation_history = {}
        self.num_perturbations = 0
        self.num_perturbations_history = {0: 0}
        self.patch_centrality_history = {}
        self.patch_degree_history = {}
        self.patch_lcc_history = {_: {} for _ in ['all', 'same', 'different'] + list(habitat_type_dictionary.keys())}
        self.degree_distribution_history = {}
        self.degree_dist_power_law_fit_history = {}
        self.total_connections_history = {}  # how many undirected links between different patches?
        self.patch_quality_history = {}
        self.update_habitat_distributions_history()
        self.update_degree_history()
        self.update_centrality_history(parameters=parameters)
        self.update_quality_history()
        self.record_xy_adjacency()

    def update_biodiversity_history(self):
        # this is called EVERY time-step and updates both the local and global properties
        current_found_species_names = set()
        current_local_biodiversity = []
        for patch_num in self.current_patch_list:
            local_biodiversity = 0
            for local_pop in self.patch_list[patch_num].local_populations.values():
                if local_pop.population >= local_pop.species.minimum_population_size:
                    current_found_species_names.add(local_pop.species.name)
                    local_biodiversity += 1
            current_local_biodiversity.append(local_biodiversity)
        self.global_biodiversity_history.append(len(current_found_species_names))
        self.local_biodiversity_history.append(tuple_builder(current_local_biodiversity))

    def update_quality_history(self):
        current_quality_list = []
        for patch_num in self.current_patch_list:
            current_quality_list.append(self.patch_list[patch_num].quality)
        self.patch_quality_history[self.step] = tuple_builder(current_quality_list)

    def update_degree_history(self):
        # record history of network-level averages and SDs of patch degree, patch clustering, and total number of edges
        current_degree_list, current_lcc_list = self.calculate_all_patches_degree()
        self.patch_degree_history[self.step] = tuple_builder(current_degree_list)
        # for average LCC, we need the average over ALL patches, over patches of the SAME habitat type (in general, and
        # separately for EACH habitat type), and over patches of DIFFERENT habitat types (i.e. at least two habitat
        # types represented in the triple).
        habitat_type_lcc_list = {}
        for key in ["all", "same", "different"]:
            habitat_type_lcc_list[key] = [x[key] for x in current_lcc_list]
        # build the lists for each habitat type
        for habitat_type_num in list(self.habitat_type_dictionary.keys()):
            habitat_type_lcc_list[habitat_type_num] = []
        for patch_index, patch_num in enumerate(self.current_patch_list):
            habitat_type_lcc_list[self.patch_list[patch_num].habitat_type_num].append(
                current_lcc_list[patch_index]["same"])
        # take all means and standard deviations and report
        key_list = ["all", "same", "different"] + list(self.habitat_type_dictionary.keys())
        for key in key_list:
            self.patch_lcc_history[key][self.step] = tuple_builder(habitat_type_lcc_list[key])

        # remaining measures:
        self.total_connections_history[self.step] = int(
            (np.sum(current_degree_list) - len(self.current_patch_list)) / 2)
        # determine degree distribution and store
        degree_distribution = Counter(current_degree_list)
        max_degree = np.max(current_degree_list)
        degree_distribution_list = []
        for degree in range(max_degree + 1):
            degree_distribution_list.append(degree_distribution[degree])
        self.degree_distribution_history[self.step] = degree_distribution_list
        self.degree_dist_power_law_fit_history[self.step] = power_law_curve_fit(
            degree_distribution_list=degree_distribution_list)

    def update_centrality_history(self, parameters):
        # update the history of tuples of the average centrality of CURRENT patches (but for all patches in their
        # individual attributes)
        current_centrality_list = self.calculate_all_patches_centrality(parameters=parameters)
        self.patch_centrality_history[self.step] = tuple_builder(current_centrality_list)

    def increment_num_perturbations(self):
        # find the previous amount and add 1
        self.num_perturbations += 1
        self.num_perturbations_history[self.step] = self.num_perturbations

    def update_habitat_distributions_history(self):
        # count the amount of each habitat and the spatial auto-correlation (essentially the a-posteriori probability
        # that two neighbours have the same habitat type) and store history of each.
        norm_sum = 0.0
        auto_cor_sum = 0.0
        temp_habitat_counts = {}
        num_habitat_types = len(self.habitat_type_dictionary)
        if num_habitat_types == 1:
            temp_habitat_counts[0] = len(self.current_patch_list)
            spatial_auto_correlation = 0.0
        else:
            # more than one habitat type
            for habitat_type_num in self.habitat_type_dictionary:
                temp_habitat_counts[habitat_type_num] = 0
            for starting_index, patch_num_1 in enumerate(self.current_patch_list):
                # only consider patches that are currently present
                temp_habitat_counts[self.patch_list[patch_num_1].habitat_type_num] += 1
                # consider neighbours
                for patch_num_2 in self.current_patch_list[starting_index + 1:]:
                    if self.patch_adjacency_matrix[patch_num_1, patch_num_2] == 1.0:
                        norm_sum += 1.0
                        if self.patch_list[patch_num_1].habitat_type_num == \
                                self.patch_list[patch_num_2].habitat_type_num:
                            auto_cor_sum += 1.0

            if norm_sum == 0.0:
                # if norm_sum is zero (i.e.the graph is fully disconnected)
                spatial_auto_correlation = 0.0
            else:
                # now determine the probabilities and thus the expected and SD of same-habitat links
                prob_square_sum = 0.0
                for habitat_type_num in self.habitat_type_dictionary:
                    prob_square_sum += (temp_habitat_counts[habitat_type_num] / len(self.current_patch_list)) ** 2.0
                expectation = prob_square_sum
                st_dev = np.sqrt(prob_square_sum - prob_square_sum ** 2.0)
                # now the habitat-probability-normalised spatial auto-correlation of habitats
                if st_dev == 0.0:
                    # all patches have the same habitat type despite there being multiple possibilities
                    spatial_auto_correlation = 1.0
                else:
                    spatial_auto_correlation = (auto_cor_sum / norm_sum - expectation) / st_dev
        self.habitat_auto_correlation_history[self.step] = spatial_auto_correlation
        for habitat_type_num in self.habitat_type_dictionary:
            self.habitat_amounts_history[habitat_type_num][self.step] = temp_habitat_counts[habitat_type_num]

    def update_all_patches_habitat_based_properties(self):
        for patch in self.patch_list:
            self.update_patch_habitat_based_properties(patch=patch)

    def update_current_patch_history(self):
        self.current_num_patches_history.append(len(self.current_patch_list))

    def calculate_all_patches_degree(self):
        # calculate the degree of each patch and update the set of adjacent patches, and THEN SUBSEQUENTLY we use that
        # to determine the local clustering coefficient
        degree_list = []
        for patch in self.patch_list:
            degree = self.calculate_patch_degree(patch=patch)
            if patch.number in self.current_patch_list:
                degree_list.append(degree)
        # This must be AFTER updating all calculate_patch_degree() because of reliance of patch.set_of_adjacent_patches
        lcc_list = []
        if len(self.current_patch_list) > 0:
            for patch in self.patch_list:
                lcc = self.calculate_lcc(patch=patch)
                if patch.number in self.current_patch_list:
                    lcc_list.append(lcc)  # this is a (patch-length) list of the [all, same, different] LCC's
        return degree_list, lcc_list

    def calculate_all_patches_centrality(self, parameters):
        # calculate the harmonic centrality of the patch
        # use the full patch list to avoid errors
        num_patches = int(len(self.patch_list))
        maximal_iterations = min(int(parameters["main_para"]["MAX_CENTRALITY_MEASURE"]), num_patches)
        minimum_distance_matrix = np.identity(num_patches)
        patch_centrality_vector = np.zeros([num_patches])
        if maximal_iterations > 1:
            path_matrix = deepcopy(self.patch_adjacency_matrix)
            for path_length in range(maximal_iterations):
                for patch_x in range(num_patches):
                    for patch_y in range(num_patches):
                        # matmul causes number of walks to grow so normalise - only need to know when non-zero
                        if path_matrix[patch_x, patch_y] != 0.0:
                            path_matrix[patch_x, patch_y] = 1.0
                            # all this relies on patch_adjacency_matrix always being undirected so [x,y] suffices.
                            if minimum_distance_matrix[patch_x, patch_y] == 0.0:
                                # have found a new connection between previously unreached paths
                                minimum_distance_matrix[patch_x, patch_y] = path_length + 1

                path_matrix = np.matmul(path_matrix, path_matrix)

            # sum reciprocal of non-zero elements (excluding self)
            for patch_x in range(num_patches):
                for patch_y in range(num_patches):
                    if patch_x != patch_y and minimum_distance_matrix[patch_x, patch_y] != 0.0:
                        patch_centrality_vector[patch_x] += 1.0 / minimum_distance_matrix[patch_x, patch_y]
            patch_centrality_vector *= (1.0 / (len(self.current_patch_list) - 1.0))

        centrality_list = []  # to be used for mean, s.d. for the system level history
        for patch in self.patch_list:
            # make sure to iterate through ALL patches as non-current are all still part of the main patch_list
            patch.centrality = patch_centrality_vector[patch.number]
            patch.centrality_history[self.step] = patch_centrality_vector[patch.number]
            if patch.number in self.current_patch_list:
                centrality_list.append(patch.centrality)
        return centrality_list

    def build_all_patches_species_paths_and_adjacency(self, parameters, specified_patch_list=None):
        for patch in self.patch_list:
            if specified_patch_list is None or patch.number in specified_patch_list:
                # By default, rebuild scores and paths for ALL patches
                self.build_species_paths_and_adjacency(patch=patch, parameters=parameters)

    def calculate_lcc(self, patch):
        # determine the lcc of the given patch. This should only be called after all patches have had
        # their .set_of_adjacent_patches updated by calling calculate_patch_degree() for EACH OF THEM
        list_of_neighbours = list(patch.set_of_adjacent_patches)
        num_triangles = [0, 0, 0]  # all, same, different
        num_closed_triangles = [0, 0, 0]  # all, same, different
        lcc = {"all": 0.0, "same": 0.0, "different": 0.0}
        if len(list_of_neighbours) > 1:
            for n_1_index, neighbour_1 in enumerate(list_of_neighbours):
                if neighbour_1 in self.current_patch_list and neighbour_1 != patch.number:
                    for neighbour_2 in list_of_neighbours[n_1_index + 1:]:  # necessarily n_1 not same as n_2
                        if neighbour_2 in self.current_patch_list and neighbour_2 != patch.number:
                            # count this triangle in 'all'
                            num_triangles[0] += 1
                            # now check 'same' or different' habitats
                            is_same_habitats = (patch.habitat_type_num == self.patch_list[neighbour_1].habitat_type_num
                                                == self.patch_list[neighbour_2].habitat_type_num)
                            if is_same_habitats:
                                num_triangles[1] += 1
                            else:
                                num_triangles[2] += 1
                            # now determine if the triangle is closed
                            if neighbour_2 in self.patch_list[neighbour_1].set_of_adjacent_patches:
                                num_closed_triangles[0] += 1
                                if is_same_habitats:
                                    num_closed_triangles[1] += 1
                                else:
                                    num_closed_triangles[2] += 1
            for key_index, key in enumerate(["all", "same", "different"]):
                if num_triangles[key_index] > 0:
                    lcc[key] = float(num_closed_triangles[key_index]) / float(num_triangles[key_index])
        # update patch records with this list
        patch.local_clustering = lcc
        patch.local_clustering_history[self.step] = lcc
        return lcc

    def calculate_patch_degree(self, patch):
        # calculate the degree centrality of the patch
        degree = 0
        set_of_adjacent_patches = set({})
        for patch_number in range(np.size(self.patch_adjacency_matrix, 0)):
            if self.patch_adjacency_matrix[patch.number, patch_number] != 0.0 \
                    or self.patch_adjacency_matrix[patch_number, patch.number] != 0.0:
                if patch_number in self.current_patch_list:
                    degree += 1
                    set_of_adjacent_patches.add(patch_number)
        patch.degree = degree
        patch.degree_history[self.step] = degree
        patch.set_of_adjacent_patches = set_of_adjacent_patches
        patch.set_of_adjacent_patches_history[self.step] = list(set_of_adjacent_patches)
        return degree

    def record_xy_adjacency(self):
        for patch_1_num, patch_1 in enumerate(self.patch_list):
            for patch_2 in self.patch_list[patch_1_num + 1:]:
                if np.linalg.norm(patch_1.position - patch_2.position) == 1:
                    patch_1.set_of_xy_adjacent_patches.add(patch_2.number)
                    patch_2.set_of_xy_adjacent_patches.add(patch_1_num)

    def update_patch_habitat_based_properties(self, patch):
        # this simply builds/resets the dictionary of the species-specific feeding and traversal scores in this
        # patch given the current habitat. It needs to be called when the patch is initiated,
        # or if the habitat in the patch is changed
        patch.this_habitat_species_feeding = {}
        patch.this_habitat_species_traversal = {}

        if len(self.species_set["list"]) == 1 and np.ndim(self.habitat_species_feeding) == 1 \
                and np.ndim(self.habitat_species_traversal) == 1:
            patch.this_habitat_species_feeding[self.species_set["list"][0].name] = \
                self.habitat_species_feeding[patch.habitat_type_num]
            patch.this_habitat_species_traversal[self.species_set["list"][0].name] = \
                self.habitat_species_traversal[patch.habitat_type_num]
        else:
            for species_number, species in enumerate(self.species_set["list"]):
                patch.this_habitat_species_feeding[species.name] = \
                    self.habitat_species_feeding[patch.habitat_type_num, species_number]
                patch.this_habitat_species_traversal[species.name] = \
                    self.habitat_species_traversal[patch.habitat_type_num, species_number]

    def build_species_paths_and_adjacency(self, patch, parameters):
        # Sets a list of dictionaries containing the shortest path cost for each species to travel from the current
        # patch to each other patch.
        # The value of this is stored in patch.species_movement_scores dictionary, where species name is the key.
        patch.species_movement_scores = {}
        #
        # take the movement_scores dictionary of arrays, and identify an iterable list of the patches that each species
        # could in principle reach - however we do not account for species-specific movement/foraging thresholds or
        # path-length restrictions here.
        # These are accounted for later, in (along with prey preferences and distance foraging strategy efficiency)
        # the "interacting populations" and "actual_dispersal_targets" methods.
        patch.adjacency_lists = {}
        patch.stepping_stone_list = []
        stepping_stone_set = set()
        #
        for species in self.species_set["list"]:
            species_name = species.name

            # use Dijkstra's algorithm for weighted undirected graphs
            # create dictionary of final patch costs for different path step-lengths [from this current patch]
            patch_costs = {x: {"best": (float('inf'), float('inf'), 0.0, [])} for x in range(len(self.patch_list))}

            # zero cost to travel to self (i.e. this patch) for any species
            patch_costs[patch.number]["best"] = (0, 0.0, [])  # 0 steps, 0.0 cost, no intermediate steps
            patch_costs[patch.number][0] = (0.0, [])

            # list of patches whose shortest path has been found
            visited = []
            not_visited = [x for x in range(len(self.patch_list))]

            while len(not_visited) > 0:
                # identify the shortest tentative cost/distance to reach a currently-unvisited (in this cycle) node
                best_tentative_cost = float('inf')
                best_tentative_length = float('inf')
                best_tentative_path = []
                next_vertex_num = 0
                for possible_next_patch_num in not_visited:
                    tentative_length = patch_costs[possible_next_patch_num]["best"][0]
                    tentative_cost = patch_costs[possible_next_patch_num]["best"][1]
                    tentative_path = patch_costs[possible_next_patch_num]["best"][2]
                    if tentative_cost < best_tentative_cost:
                        next_vertex_num = possible_next_patch_num
                        best_tentative_length = tentative_length
                        best_tentative_cost = tentative_cost
                        best_tentative_path = tentative_path
                # if this cost is still infinity, then quit algorithm as the graph is disconnected
                if best_tentative_cost == float('inf'):
                    break
                else:
                    visited.append(next_vertex_num)
                    not_visited.remove(next_vertex_num)
                    # now see if a better score to other patches can be achieved through this one
                    for other_patch_num, other_patch in enumerate(self.patch_list):
                        if next_vertex_num != other_patch_num:
                            # not self!
                            if self.patch_adjacency_matrix[next_vertex_num, other_patch_num] == 0.0 \
                                    or other_patch.this_habitat_species_traversal[species_name] <= 0.0:
                                new_path_cost = float('inf')
                                new_path_length = float('inf')
                                new_path = []
                            else:
                                new_path_length = best_tentative_length + 1
                                new_path = best_tentative_path + [next_vertex_num]
                                # single-path-cost = 1 / ( habitat-species-traversal * adjacency-border * patch-size)
                                new_path_cost = \
                                    best_tentative_cost + (1.0 / other_patch.this_habitat_species_traversal[
                                        species_name]) * (1.0 / self.patch_adjacency_matrix[
                                        next_vertex_num, other_patch_num]) / other_patch.size
                                # Note: patch_adjacency_matrix is currently binary, so if this branch is reached
                                # then part of this function will be 1/1;
                                # However it is included because in the future we may wish to alter this matrix
                                # such that there are non-uniform size of borders between patch pairs (separately
                                # from the role of patch size, which linearly reduces the distance cost to access).

                            # is best overall?
                            if new_path_cost < patch_costs[other_patch_num]["best"][1]:
                                patch_costs[other_patch_num]["best"] = (new_path_length, new_path_cost, new_path)
                            # is best for this length?
                            if new_path_length not in patch_costs[other_patch_num] or \
                                    new_path_cost < patch_costs[other_patch_num][new_path_length][0]:
                                patch_costs[other_patch_num][new_path_length] = (new_path_cost, new_path)

            # save
            patch.species_movement_scores[species_name] = patch_costs
            # now select those that are theoretically reachable for foraging/direct dispersal from this patch
            # for each species - i.e. if they had infinite dispersal mobility in the CURRENT spatial network
            # with the given paths and habitats (and that the species' inherent properties of being able to
            # traverse certain habitats remains unchanged.
            reachable_patch_nums = []
            for patch_num in patch.species_movement_scores[species_name]:
                if patch.species_movement_scores[species_name][patch_num]["best"][1] < float('inf'):
                    reachable_patch_nums.append(patch_num)
            patch.adjacency_lists[species_name] = reachable_patch_nums
            # gather a set of the patches used as stepping stones AND the reachable patches (inc. endpoints) using them
            for target, possible_routes in patch_costs.items():
                if len(possible_routes["best"][-1]) > 0:
                    for route in possible_routes.values():
                        path_list = route[-1]
                        # but do not count arbitrarily long paths
                        if 0 < len(path_list) <= parameters["main_para"]["ASSUMED_MAX_PATH_LENGTH"] + 1:
                            stepping_stone_set = stepping_stone_set.union(set(path_list + [int(target)]))
        patch.stepping_stone_list = list(stepping_stone_set)
        print(f"Paths built for patch {patch.number}/{len(self.patch_list) - 1}")

    # --------------------------- SPECIES / COMMUNITY DISTRIBUTION ANALYSIS ----------------------------------------- #
    def distance_metrics(self):
        num_patches = len(self.current_patch_list)
        num_species = len(self.species_set["list"])

        # gather the data
        community_state_presence_array = np.zeros([num_patches, num_species])
        community_state_population_array = np.zeros([num_patches, num_species])  # normalise each column by total pop
        for species_index, species in enumerate(self.species_set["list"]):
            total_population = 0.0
            for patch_num in self.current_patch_list:
                # presence/absence
                community_state_presence_array[patch_num, species_index] = \
                    self.patch_list[patch_num].local_populations[species.name].occupancy
                # population weighted
                this_population = self.patch_list[patch_num].local_populations[species.name].population
                community_state_population_array[patch_num, species_index] = this_population
                total_population += this_population
            if total_population > 0.0:
                community_state_population_array[:, species_index] = \
                    community_state_population_array[:, species_index] / total_population

        # per species distance metrics - and species presence probabilities (overall and per habitat type)
        for species_index in range(num_species):
            self.network_analysis(patch_value_array=community_state_presence_array[:, species_index],
                                  is_binary=True, is_distribution=False)
            self.network_analysis(patch_value_array=community_state_population_array[:, species_index],
                                  is_binary=False, is_distribution=False)

        # community distance metrics
        self.network_analysis(patch_value_array=community_state_presence_array,
                              is_binary=False, is_distribution=True)
        self.network_analysis(patch_value_array=community_state_population_array,
                              is_binary=False, is_distribution=False)

        # each community state probabilities (overall and per habitat type)
        #
        # for each state determine a unique binary identifier
        community_state_binary = np.zeros(num_patches)
        for row_index, row in enumerate(community_state_presence_array):
            binary_identifier = 0
            for column_index, column_value in enumerate(row):
                binary_identifier += column_value * 2.0 ** column_index
            community_state_binary[row_index] = binary_identifier
        # how many UNIQUE states were identified?
        extant_state_set = set({})
        for state in community_state_binary:
            extant_state_set.add(state)
        # then set the patch-vector by presence-absence just based on presence of each state and analyse
        ordered_state_list = list(extant_state_set)
        ordered_state_list.sort()
        for state in ordered_state_list:
            state_array = np.zeros(num_patches)
            for patch_index, patch_state in enumerate(community_state_binary):
                if patch_state == state:
                    state_array[patch_index] = 1.0
            self.network_analysis(patch_value_array=state_array, is_binary=True, is_distribution=False)

        # Shannon entropy
        self.shannon_entropy(patch_state_array=community_state_presence_array,
                             patch_binary_vector=community_state_binary)

        # species-species predictions
        self.inter_species_predictions()

    def inter_species_predictions(self):
        # 1. Species-species presence/absence in same patch
        # 2. Population-population correlation coefficients in same patch (both Spearman and Pearson)
        # 3. For each species, Mantel test for correlation between matrix of path lengths (compare both adjacency and
        # the species-specific shortest distance traversal score) and matrix of population differences across patches.
        # 4. (Applied for radius 1 (i.e. same patch), radius 2 (i.e. path length 1), and radius 3 (
        #   i.e. path length 2) calculation of expected value of population of one species given the presence of
        #   the unit population of a control species. Note carefully that the interpretation is in terms of 'realised
        #   maximum', and not the theoretical maximum in absence of the control species or other parameters. So it
        #   won't give you much information about the target population if the control population is everywhere, unless
        #   you do some control experiments with only the target species or make a decision about its theoretical
        #   maximum and simply rescale by that. The normalised version will tell you about how the RELATIVE distribution
        #   of the target species is correlated to the presence and size of the control species.
        #
        # In cases 1, 2, and 4, calculate these separately for each habitat type subnetwork in addition to overall.
        pass

    def network_analysis(self, patch_value_array, is_binary, is_distribution):
        # need value, habitat type, and neighbours of each patch for presence, auto_correlation, clustering analysis
        num_patches = len(self.current_patch_list)
        if np.ndim(patch_value_array) == 1:
            max_difference = 1
        else:
            max_difference = np.shape(patch_value_array)[1]  # should be either 1 or num_species
        patch_habitat = []
        patch_neighbours = []

        if len(patch_value_array) != num_patches:
            # how many ROWS in the array? Should match length of current_patch_list
            raise Exception("Incorrect dimensions of value array.")

        for patch_num in self.current_patch_list:
            patch_habitat.append(self.patch_list[patch_num].habitat_type_num)
            patch_neighbours.append(self.patch_list[patch_num].set_of_adjacent_patches)

        # set up the nested dictionary to hold results
        template_presence = {"all": (0.0, 0.0)}  # (mean, standard deviation)
        template_auto_corr = {"all": (0.0, 0.0),  # (matching pairs, eligible pairs)
                              "same": (0.0, 0.0),
                              "different": (0.0, 0.0),
                              }
        difference_distribution = {"all": np.zeros(max_difference + 1),  # includes '0' as a possible difference
                                   "same": np.zeros(max_difference + 1),
                                   "different": np.zeros(max_difference + 1),
                                   }

        # add keys for each habitat, habitat type pair
        habitat_type_nums = list(self.habitat_type_dictionary.keys())
        habitat_type_nums.sort()
        for habitat_type_num_1 in habitat_type_nums:
            template_presence[habitat_type_num_1] = (0.0, 0.0)
            for habitat_type_num_2 in habitat_type_nums[habitat_type_num_1:]:
                template_auto_corr[(habitat_type_num_1, habitat_type_num_2)] = (0.0, 0.0)
                difference_distribution[(habitat_type_num_1, habitat_type_num_2)] = np.zeros(max_difference + 1)

        # presence probability of this configuration
        if is_binary:
            # i.e. skip this for population-weighted and full community states
            template_presence["all"] = (np.mean(patch_value_array), np.std(patch_value_array))
            for habitat_type_num_1 in habitat_type_nums:
                habitat_subnet = [patch_value_array[x] for x in range(
                    num_patches) if patch_habitat[x] == habitat_type_num_1]
                if len(habitat_subnet) > 0:
                    template_presence[habitat_type_num_1] = (np.mean(habitat_subnet), np.std(habitat_subnet))

        # auto-correlation
        for patch_num in range(num_patches):
            this_patch_habitat = patch_habitat[patch_num]
            # note that we do *NOT* also calculate this separately for each community state (0-2^N) or
            # species state (0-1, i.e. we do not restrict to counting only over patches where the species was present,
            # but we should be able to easily obtain this average instead if desired since we also store the probability
            # of species presence, and of each community state).
            for patch_neighbour in patch_neighbours[patch_num]:

                # avoid double counting
                if patch_neighbour > patch_num:

                    neighbouring_patch_habitat = self.patch_list[patch_neighbour].habitat_type_num
                    # taxicab / manhattan norm:
                    difference_vector = patch_value_array[patch_num] - patch_value_array[patch_neighbour]
                    if np.ndim(difference_vector) == 0:
                        difference_vector = [difference_vector]
                    l1_difference = int(np.linalg.norm(difference_vector, ord=1))
                    normalised_difference = float(l1_difference) / max_difference
                    similarity = 1.0 - normalised_difference

                    # community difference distributions
                    if is_distribution:
                        # the taxi cab norm (for presence/absence state values) has the advantage of being a
                        # finite set of possible values
                        difference_modifier = 1
                    else:
                        difference_modifier = 0

                    # all
                    template_auto_corr['all'] += (similarity, 1.0)
                    difference_distribution['all'][l1_difference] += difference_modifier
                    # same
                    if this_patch_habitat == neighbouring_patch_habitat:
                        template_auto_corr['same'] += (similarity, 1.0)
                        difference_distribution['same'][l1_difference] += difference_modifier
                    else:
                        template_auto_corr['different'] += (similarity, 1.0)
                        difference_distribution['different'][l1_difference] += difference_modifier
                    # specific habitat combinations
                    small_habitat = min(this_patch_habitat, neighbouring_patch_habitat)  # order them
                    large_habitat = max(this_patch_habitat, neighbouring_patch_habitat)
                    template_auto_corr[(small_habitat, large_habitat)] += (similarity, 1.0)
                    difference_distribution[(small_habitat, large_habitat)][l1_difference] += difference_modifier

        output_dict = {
            "presence": template_presence,
            "auto_correlation": template_auto_corr,
            "difference_distribution": difference_distribution,
        }
        return output_dict

    def shannon_entropy(self, patch_state_array, patch_binary_vector):
        # Shannon information of community probability distributions:
        #   1. Considering each separate species state as a distinct variable value
        #   2. Just of the plain diversity (amount of species)
        # both of these are calculated over the whole system, and over each single-habitat-type subnetwork
        num_patches = len(self.current_patch_list)
        max_difference = np.shape(patch_state_array)[1]  # should be num_species (binary states for each species)

        if len(patch_state_array) != num_patches or len(patch_binary_vector) != num_patches:
            # how many ROWS in the array? Should match length of current_patch_list
            raise Exception("Incorrect dimensions of value array.")

        patch_habitat = []
        for patch_num in self.current_patch_list:
            patch_habitat.append(self.patch_list[patch_num].habitat_type_num)  # note here that if patches are deleted
            # from the current system state, then the position in this list MIGHT NOT MATCH the patch number.

        # add keys for each habitat
        habitat_type_nums = list(self.habitat_type_dictionary.keys())
        habitat_type_nums.sort()

        # this will hold the frequency distributions
        shannon_array = {'all': {
            'total': 0, 'dist_state': np.zeros(max_difference + 1), 'dist_biodiversity': np.zeros(max_difference + 1)}}
        for habitat_type_num in habitat_type_nums:
            shannon_array[habitat_type_num] = {'total': 0,
                                               'dist_state': np.zeros(max_difference + 1),
                                               'dist_biodiversity': np.zeros(max_difference + 1)}

        # iterate only once over the patches
        for temp_patch_index, patch_row in enumerate(patch_state_array):
            # increment counter for all and for this patch's habitat type
            shannon_array['all']['total'] += 1
            shannon_array[patch_habitat[temp_patch_index]]['total'] += 1
            # identify the state of each patch
            #
            # for bio_diversity
            biodiversity = int(np.sum(patch_row))
            # for full system state
            state_value = int(patch_binary_vector[temp_patch_index])
            shannon_array['all']['dist_state'][state_value] += 1
            shannon_array['all']['dist_biodiversity'][biodiversity] += 1
            shannon_array[patch_habitat[temp_patch_index]]['dist_state'][state_value] += 1
            shannon_array[patch_habitat[temp_patch_index]]['dist_biodiversity'][biodiversity] += 1

        # Now calculate the shannon entropy for each distribution. This is done for 2 x (num_habitats+1) configurations.
        shannon_entropy = {}
        for key, value in shannon_array.items():
            state_entropy_sum = 0.0
            biodiversity_entropy_sum = 0.0
            if value['total'] > 0.0:
                for state_frequency in value['dist_state']:
                    if state_frequency > 0.0:
                        state_entropy_sum += (state_frequency / value['total']
                                              ) * np.log(state_frequency / value['total'])
                for biodiversity_frequency in value['dist_biodiversity']:
                    if biodiversity_frequency > 0.0:
                        biodiversity_entropy_sum += (biodiversity_frequency / value['total']
                                                     ) * np.log(biodiversity_frequency / value['total'])
            shannon_entropy[key] = {'state': -state_entropy_sum, 'biodiversity': -biodiversity_entropy_sum}
        return shannon_entropy
