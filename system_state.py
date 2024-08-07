import numpy as np
from copy import deepcopy
from collections import Counter
from degree_distribution import power_law_curve_fit
from scipy.stats import spearmanr, pearsonr
from data_manager_functions import update_local_population_nets
from system_state_functions import tuple_builder, linear_model_report, generate_cluster, \
    determine_complexity, rank_abundance


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
        self.distance_metrics_store = {}  # will hold a very deep nested dictionary of measures

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
        self.habitat_spatial_auto_correlation_history = {}  # normalised by the expectation given habitat amounts
        self.habitat_regular_auto_correlation_history = {}  # plain ratio of same-habitat : any-habitat links
        self.num_perturbations = 0
        self.num_perturbations_history = {0: 0}
        self.patch_centrality_history = {}
        self.patch_degree_history = {}
        self.patch_lcc_history = {_: {} for _ in ['all', 'same', 'different'] + list(habitat_type_dictionary.keys())}
        self.degree_distribution_history = {}
        self.degree_dist_power_law_fit_history = {}
        self.total_connections_history = {}  # how many undirected links between different patches?
        self.patch_quality_history = {}
        self.patch_size_history = {}
        self.update_habitat_distributions_history()
        self.update_degree_history()
        self.update_centrality_history(parameters=parameters)
        self.update_quality_history()
        self.update_size_history()
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

    def update_size_history(self):
        current_size_list = []
        for patch_num in self.current_patch_list:
            current_size_list.append(self.patch_list[patch_num].size)
        self.patch_size_history[self.step] = tuple_builder(current_size_list)

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

    def update_local_capacity_history(self):
        # this is used for calculating source and sink indices of local populations,
        # which could in principle be altered if we introduced a "change patch size" perturbation in future.
        for patch in self.patch_list:
            for local_population in patch.local_populations.values():
                local_population.record_carrying_capacity_history(step=self.step)

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
            regular_auto_correlation = 0.0
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
                regular_auto_correlation = 0.0
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
                    regular_auto_correlation = 1.0
                else:
                    # regular (non-normalised) habitat auto-correlation:
                    regular_auto_correlation = auto_cor_sum / norm_sum

                    # normalised by expectation given the habitat probability distribution:
                    spatial_auto_correlation = (auto_cor_sum / norm_sum - expectation) / st_dev
        # Record:
        self.habitat_regular_auto_correlation_history[self.step] = regular_auto_correlation
        self.habitat_spatial_auto_correlation_history[self.step] = spatial_auto_correlation
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
        # Sets a list of dictionaries containing the shortest path COST (not SCORE) for each species to travel from the
        # current patch to each other patch.
        # The value of this is stored in patch.species_movement_scores dictionary, where species name is the key.
        patch.species_movement_scores = {}
        #
        # take the movement_scores dictionary of arrays, and identify an iterable list of the patches that each species
        # could in principle reach - however we do not account for species-specific movement/foraging thresholds or
        # path-length restrictions here.
        # These are accounted for later, in (along with prey preferences and distance foraging strategy efficiency)
        # the "interacting populations" and "actual_dispersal_targets" methods.
        #
        # Really what we calculate here, and store in the .species_movement_scores, are a set of potential (subject to
        # maximum path length and kappa - for foraging) TRAVEL COSTS to which search cost will later be added, and
        # then conversion to foraging or dispersal score will take place.
        patch.adjacency_lists = {}
        patch.stepping_stone_list = []
        stepping_stone_set = set()
        #
        for species in self.species_set["list"]:
            species_name = species.name

            # use Dijkstra's algorithm for weighted undirected graphs
            # create dictionary of final patch costs for different path step-lengths [from this current patch]
            routes_template = {"routes": {"best": (float('inf'), float('inf'), 0.0, [])}}
            patch_costs = {x: deepcopy(routes_template) for x in range(len(self.patch_list))}

            # zero cost to travel to self (i.e. this patch) for any species
            patch_costs[patch.number]["routes"]["best"] = (0, 0.0, [])  # 0 steps, 0.0 cost, no intermediate steps
            patch_costs[patch.number]["routes"][0] = (0.0, [])

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
                    tentative_length = patch_costs[possible_next_patch_num]["routes"]["best"][0]
                    tentative_cost = patch_costs[possible_next_patch_num]["routes"]["best"][1]
                    tentative_path = patch_costs[possible_next_patch_num]["routes"]["best"][2]
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
                    next_patch = self.patch_list[next_vertex_num]
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
                                # single-path-cost = patch-size / ( habitat-species-traversal * adjacency-border)
                                new_path_cost = best_tentative_cost + next_patch.size / (
                                        next_patch.this_habitat_species_traversal[
                                            species_name] * self.patch_adjacency_matrix[
                                            next_vertex_num, other_patch_num])
                                # Note: patch_adjacency_matrix is currently binary, so if this branch is reached
                                # then part of this function will be 1/1;
                                # However it is included because in the future we may wish to alter this matrix
                                # such that there are non-uniform size of borders between patch pairs (separately
                                # from the role of patch size).

                            # is best overall?
                            if new_path_cost < patch_costs[other_patch_num]["routes"]["best"][1]:
                                patch_costs[other_patch_num]["routes"]["best"] = (
                                    new_path_length, new_path_cost, new_path)
                            # is best for this length?
                            if new_path_length not in patch_costs[other_patch_num]["routes"] or \
                                    new_path_cost < patch_costs[other_patch_num]["routes"][new_path_length][0]:
                                patch_costs[other_patch_num]["routes"][new_path_length] = (new_path_cost, new_path)
                            # store the target patch size and traversal score for this species
                            patch_costs[other_patch_num]["target_patch_size"] = other_patch.size
                            patch_costs[other_patch_num][
                                "target_patch_traversal"] = other_patch.this_habitat_species_traversal[species_name]

            # save
            patch.species_movement_scores[species_name] = patch_costs
            # now select those that are theoretically reachable for foraging/direct dispersal from this patch
            # for each species - i.e. if they had infinite dispersal mobility in the CURRENT spatial network
            # with the given paths and habitats (and that the species' inherent properties of being able to
            # traverse certain habitats remains unchanged.
            reachable_patch_nums = []
            for patch_num in patch.species_movement_scores[species_name]:
                if patch.species_movement_scores[species_name][patch_num]["routes"]["best"][1] < float('inf'):
                    reachable_patch_nums.append(patch_num)
            patch.adjacency_lists[species_name] = reachable_patch_nums
            # gather a set of the patches used as stepping stones AND the reachable patches (inc. endpoints) using them
            for target, target_costs in patch_costs.items():
                possible_routes = target_costs["routes"]
                if len(possible_routes["best"][-1]) > 0:
                    for route in possible_routes.values():
                        path_list = route[-1]
                        # but do not count arbitrarily long paths
                        if 0 < len(path_list) <= parameters["main_para"]["ASSUMED_MAX_PATH_LENGTH"] + 1:
                            stepping_stone_set = stepping_stone_set.union(set(path_list + [int(target)]))
        patch.stepping_stone_list = list(stepping_stone_set)
        print(f"Paths built for patch {patch.number}/{len(self.patch_list) - 1}")

    # --------------------------- SPECIES / COMMUNITY DISTRIBUTION ANALYSIS ----------------------------------------- #
    def update_distance_metrics(self, parameters):
        # Call this to conduct extensive population distribution and community state spatial distribution analysis.
        [network_analysis_species, network_analysis_community_distance, network_analysis_state_probability,
         shannon_entropy, inter_species_predictions_final, inter_species_predictions_average,
         complexity_final, complexity_average, rank_abundance_final, rank_abundance_average] = self.distance_metrics(
            parameters=parameters)
        self.distance_metrics_store = {
            "network_analysis_species": network_analysis_species,
            "network_analysis_community_distance": network_analysis_community_distance,
            "network_analysis_state_probability": network_analysis_state_probability,
            "shannon_entropy": shannon_entropy,
            "inter_species_predictions_final": inter_species_predictions_final,
            "inter_species_predictions_average": inter_species_predictions_average,
            "complexity_final": complexity_final,
            "complexity_average": complexity_average,
            "rank_abundance_final": rank_abundance_final,
            "rank_abundance_average": rank_abundance_average,
        }

    def distance_metrics(self, parameters):
        # Analysis of species and community distributions:
        #
        #  * inter-species predictions
        #     – presence
        #     – similarity
        #     – prediction
        #     – correlation
        #     – linear model
        # * shannon
        #     – determines the entropy of the biodiversity distribution and of the community state distribution
        # * network analysis state probability
        #     – binary-summed state
        #         * presence
        #             - probability of presence (plus irrelevant SD) in each subnetwork of this community state
        #         * auto-correlation
        #             - for each state, the network is classified simply as 0/1 in terms of sharing that exact state,
        #               then we get the auto-correlation for this state across various subnetworks (whole, only
        #               between patches of the same or of different habitats, and between specific habitat pairs).
        # * network analysis community distance
        #     – considering the presence/absence state of all species simultaneously (i.e. ignoring population)
        #         * auto-correlation
        #             - for each subnetwork and link type, record the ratio (same state / link)
        #         * degree distribution
        #             - for each subnetwork and link type, record distribution of the manhattan distance between states
        # * network analysis species
        #     – considering the presence/absence state of only this species (i.e. ignoring population)
        #         * presence
        #             - probability of presence (plus irrelevant SD) in each subnetwork of this species
        #         * auto-correlation
        #             - for each subnetwork and link type, record the
        #               ratio (same state [i.e. species is present or absent] / link)
        #     – population analytics:
        #         * average (of normalised by max local population anywhere in system) population in each sub-network
        #         * average (of normalised by max local population anywhere in system) difference in population sizes
        #           across links - recorded as an array with the total difference and the counted links (for
        #           all, same, different, and every specific habitat pairing).
        # * complexity analysis
        #       - species-area relation (SAR), as the count of the number of different species (diversity) in clusters
        #       - complexity information dimension, as the scaling of state complexity (both binary and
        #           population-weighted versions) with cluster size, where complexity is the summed difference between
        #           all the patch-pairs in the cluster.
        # * species rank-abundance
        #       - the rate of decrease of abundance (total global population) with species' rank by abundance.

        is_record_lm_vectors = parameters["main_para"]["IS_RECORD_METRICS_LM_VECTORS"]
        update_local_population_nets(system_state=self)  # required to update .occupancy attribute of local_populations
        num_patches = len(self.current_patch_list)
        num_species = len(self.species_set["list"])
        species_list = [x.name for x in self.species_set["list"]]

        patch_habitat = []
        patch_neighbours = []
        # Build convenient lists of the important properties for patches which are currently eligible.
        # We should use these (not self.patch_list) for iterating current patches, but when we have the necessary patch
        # number we CAN then use that to access any additional required attributes from patch_list.
        for patch_num in self.current_patch_list:
            patch_habitat.append(self.patch_list[patch_num].habitat_type_num)  # note here that if patches are deleted
            # from the current system state, then the position in this list MIGHT NOT MATCH the patch number.
            patch_neighbours.append(self.patch_list[patch_num].set_of_adjacent_patches)

        # gather the data
        community_state_presence_array = np.zeros([num_patches, num_species])
        community_state_population_array = np.zeros([num_patches, num_species])
        time_averaged_population_array = np.zeros([num_patches, num_species])  # prediction analysis for ave populations
        for species_index, species_name in enumerate(species_list):
            for patch_num in self.current_patch_list:
                # presence/absence
                community_state_presence_array[patch_num, species_index] = \
                    self.patch_list[patch_num].local_populations[species_name].occupancy
                # final population
                community_state_population_array[patch_num, species_index] = self.patch_list[
                    patch_num].local_populations[species_name].population
                # average population
                time_averaged_population_array[patch_num, species_index] = self.patch_list[
                    patch_num].local_populations[species_name].average_population

        # per species distance metrics - and species presence probabilities (overall and per habitat type)
        network_analysis_species = {}
        for species_index, species_name in enumerate(species_list):

            max_population = max(community_state_population_array[:, species_index])
            if max_population > 0.0:
                norm_species_pop_vector = community_state_population_array[:, species_index] / max_population

                network_analysis_species[species_name] = {
                    "species_presence": self.network_analysis(
                        patch_value_array=community_state_presence_array[:, species_index],
                        patch_habitat=patch_habitat, patch_neighbours=patch_neighbours,
                        is_presence=True, is_distribution=False),

                    "species_population": self.network_analysis(
                        patch_value_array=norm_species_pop_vector,
                        patch_habitat=patch_habitat, patch_neighbours=patch_neighbours,
                        is_presence=True, is_distribution=False),
                }

        # community distance metrics
        network_analysis_community_distance = self.network_analysis(
            patch_value_array=community_state_presence_array,
            patch_habitat=patch_habitat, patch_neighbours=patch_neighbours,
            is_presence=False, is_distribution=True)

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
        network_analysis_state_probability = {}
        for state in [0, 1, 2, 3]:  # ordered_state_list
            state_array = np.zeros(num_patches)
            for patch_index, patch_state in enumerate(community_state_binary):
                if patch_state == state:
                    state_array[patch_index] = 1.0
            network_analysis_state_probability[state] = self.network_analysis(
                patch_value_array=state_array, patch_habitat=patch_habitat,
                patch_neighbours=patch_neighbours, is_presence=True, is_distribution=False)
            state_species_list = []
            for species_index in range(len(species_list)):
                if np.mod(state, int(2.0 ** (species_index + 1))) >= int(2.0 ** species_index):
                    state_species_list.append(species_list[species_index])
            network_analysis_state_probability[state]["state_species_list"] = state_species_list

        # Shannon entropy
        shannon_entropy = self.shannon_entropy(patch_state_array=community_state_presence_array,
                                               patch_habitat=patch_habitat, patch_binary_vector=community_state_binary)

        # species-species predictions and information scaling dimensions
        #

        # Prepare keys for each habitat (need this first for consistency)
        habitat_type_nums = list(self.habitat_type_dictionary.keys())
        habitat_type_nums.sort()

        # generate sub-networks
        community_presence_sub_networks = self.generate_sub_networks(population_array=community_state_presence_array,
                                                                     patch_habitat=patch_habitat,
                                                                     habitat_type_nums=habitat_type_nums)

        community_state_sub_networks = self.generate_sub_networks(population_array=community_state_population_array,
                                                                  patch_habitat=patch_habitat,
                                                                  habitat_type_nums=habitat_type_nums)

        time_averaged_sub_networks = self.generate_sub_networks(population_array=time_averaged_population_array,
                                                                patch_habitat=patch_habitat,
                                                                habitat_type_nums=habitat_type_nums)

        # use the non-normalised final population (returns five dictionaries)
        presence_store, similarity_store, prediction_store, correlation_store, linear_model_store = \
            self.inter_species_predictions(sub_networks=community_state_sub_networks,
                                           habitat_type_nums=habitat_type_nums,
                                           is_record_lm_vectors=is_record_lm_vectors)
        inter_species_predictions_final = {
            "presence_store": presence_store,
            "similarity_store": similarity_store,
            "prediction_store": prediction_store,
            "correlation_store": correlation_store,
            "linear_model_store": linear_model_store,
        }
        # and the time-averaged populations
        presence_store, similarity_store, prediction_store, correlation_store, linear_model_store = \
            self.inter_species_predictions(sub_networks=time_averaged_sub_networks,
                                           habitat_type_nums=habitat_type_nums,
                                           is_record_lm_vectors=is_record_lm_vectors)
        inter_species_predictions_average = {
            "presence_store": presence_store,
            "similarity_store": similarity_store,
            "prediction_store": prediction_store,
            "correlation_store": correlation_store,
            "linear_model_store": linear_model_store,
        }

        # Now use both the final and time-averaged community populations to determine SAR, binary and
        # population-weighted information dimension (community spatial complexity), and species-rank abundance
        complexity_final = self.complexity_analysis(
            sub_networks=community_state_sub_networks,
            corresponding_binary=community_presence_sub_networks,
            is_record_lm_vectors=is_record_lm_vectors)
        complexity_average = self.complexity_analysis(
            sub_networks=time_averaged_sub_networks,
            corresponding_binary=None,
            is_record_lm_vectors=is_record_lm_vectors)

        # Species-rank abundance - fit a log-linear model to the curve of relative global abundance vs. species rank
        #
        # Calculated from final populations:
        rank_abundance_final = rank_abundance(sub_networks=community_state_sub_networks,
                                              is_record_lm_vectors=is_record_lm_vectors)
        # Calculated from time-averaged populations:
        rank_abundance_average = rank_abundance(sub_networks=time_averaged_sub_networks,
                                                is_record_lm_vectors=is_record_lm_vectors)

        return [network_analysis_species, network_analysis_community_distance, network_analysis_state_probability,
                shannon_entropy, inter_species_predictions_final, inter_species_predictions_average,
                complexity_final, complexity_average, rank_abundance_final, rank_abundance_average]

    def network_analysis(self, patch_value_array, patch_habitat, patch_neighbours, is_presence, is_distribution):
        # need value, habitat type, and neighbours of each patch for presence, auto_correlation, clustering analysis
        num_patches = len(self.current_patch_list)
        if np.ndim(patch_value_array) == 1:
            max_difference = 1
        else:
            max_difference = np.shape(patch_value_array)[1]  # should be either 1 (for a single species - present or
            # absent, or equal to num_species for the maximum possible difference in states)

        if len(patch_value_array) != num_patches:
            # how many ROWS in the array? Should match length of current_patch_list
            raise Exception("Incorrect dimensions of value array.")

        # set up the required nested dictionaries to hold results
        template_auto_corr = {"all": np.array([0.0, 0.0]),  # (matching pairs, eligible pairs)
                              "same": np.array([0.0, 0.0]),
                              "different": np.array([0.0, 0.0]),
                              }
        template_presence = {}  # define empty dictionaries to clear pre-allocation warnings
        difference_distribution = {}
        if is_presence:
            template_presence = {"all": np.array([0.0, 0.0])}  # (mean, standard deviation)
        if is_distribution:
            difference_distribution = {"all": np.zeros(max_difference + 1),  # includes '0' as a possible difference
                                       "same": np.zeros(max_difference + 1),
                                       "different": np.zeros(max_difference + 1),
                                       }

        # add keys for each habitat, habitat type pair
        habitat_type_nums = list(self.habitat_type_dictionary.keys())
        habitat_type_nums.sort()
        for habitat_type_num_1 in habitat_type_nums:
            if is_presence:
                template_presence[habitat_type_num_1] = np.array([0.0, 0.0])
            for habitat_type_num_2 in habitat_type_nums[habitat_type_num_1:]:
                template_auto_corr[(habitat_type_num_1, habitat_type_num_2)] = np.array([0.0, 0.0])
                if is_distribution:
                    difference_distribution[(habitat_type_num_1, habitat_type_num_2)] = np.zeros(max_difference + 1)

        # presence probability of this configuration
        if is_presence:
            # i.e. skip this for full community states
            template_presence["all"] = np.array([np.mean(patch_value_array), np.std(patch_value_array)])
            for habitat_type_num_1 in habitat_type_nums:
                habitat_subnet = [patch_value_array[x] for x in range(
                    num_patches) if patch_habitat[x] == habitat_type_num_1]
                if len(habitat_subnet) > 0:
                    template_presence[habitat_type_num_1] = np.array([np.mean(habitat_subnet), np.std(habitat_subnet)])

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
                    l1_difference = np.linalg.norm(difference_vector, ord=1)
                    integer_difference = int(l1_difference)  # only for degree distributions
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
                    template_auto_corr['all'] += [similarity, 1.0]
                    # same
                    if this_patch_habitat == neighbouring_patch_habitat:
                        template_auto_corr['same'] += [similarity, 1.0]
                        if is_distribution:
                            difference_distribution['same'][integer_difference] += difference_modifier
                    else:
                        template_auto_corr['different'] += [similarity, 1.0]
                        if is_distribution:
                            difference_distribution['different'][integer_difference] += difference_modifier
                    # specific habitat combinations
                    small_habitat = min(this_patch_habitat, neighbouring_patch_habitat)  # order them
                    large_habitat = max(this_patch_habitat, neighbouring_patch_habitat)
                    template_auto_corr[(small_habitat, large_habitat)] += [similarity, 1.0]
                    if is_distribution:
                        difference_distribution['all'][integer_difference] += difference_modifier
                        difference_distribution[(small_habitat, large_habitat)][
                            integer_difference] += difference_modifier

        output_dict = {
            "auto_correlation": template_auto_corr,
        }
        if is_presence:
            output_dict["presence"] = template_presence
        if is_distribution:
            output_dict["difference_distribution"] = difference_distribution

        return output_dict

    def shannon_entropy(self, patch_state_array, patch_habitat, patch_binary_vector):
        # Shannon information of community probability distributions:
        #   1. Considering each separate species state as a distinct variable value
        #   2. Just of the plain diversity (amount of species)
        # both of these are calculated over the whole system, and over each single-habitat-type subnetwork

        num_patches = len(self.current_patch_list)
        if len(patch_binary_vector) != num_patches:
            # how many ROWS in the array? Should match length of current_patch_list
            raise Exception("Incorrect dimensions of value array.")

        # need to convert patch_binary_vector to just the simplest index of achieved states (up to num_patches)
        found_states = set({})
        for state in patch_binary_vector:
            found_states.add(state)
        found_state_list = list(found_states)
        found_state_list.sort()
        reduced_patch_binary_vector = np.zeros(num_patches)
        for patch_index in range(num_patches):
            reduced_patch_binary_vector[patch_index] = found_state_list.index(patch_binary_vector[patch_index])
        max_state_difference = int(len(found_state_list))
        max_biodiversity = np.shape(patch_state_array)[1] + 1  # add +1 for 'zero' biodiversity

        # add keys for each habitat
        habitat_type_nums = list(self.habitat_type_dictionary.keys())
        habitat_type_nums.sort()

        # this will hold the frequency distributions
        shannon_array = {'all': {
            'total': 0, 'dist_state': np.zeros(max_state_difference), 'dist_biodiversity': np.zeros(max_biodiversity)}}
        for habitat_type_num in habitat_type_nums:
            shannon_array[habitat_type_num] = {'total': 0,
                                               'dist_state': np.zeros(max_state_difference),
                                               'dist_biodiversity': np.zeros(max_biodiversity)}

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
            state_value = int(reduced_patch_binary_vector[temp_patch_index])
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

    def inter_species_predictions(self, sub_networks, habitat_type_nums, is_record_lm_vectors):
        # For each subnetwork (of non-zero size)
        #   – for each (control) species
        #       - for each (control species) radius
        # 		    - for each (response) species
        # 			    - for each (response species) radius
        # 			        A. determine the presence predictions
        #                   B. similarity
        #                   C. species A -> species B (product and root calculation)
        # 				    D. correlation coefficients (for co-variance of populations)
        # 					E. linear model fit
        #

        species_list = [x.name for x in self.species_set["list"]]

        # Prepare the nested storage structures:
        inner_species_val = {species_name: {radius: 0.0 for radius in range(3)} for species_name in species_list}
        outer_species_val = {species_name: {radius: deepcopy(inner_species_val) for radius in range(3)}
                             for species_name in species_list}
        inner_species_dict = {species_name: {radius: {} for radius in range(3)} for species_name in species_list}
        outer_species_dict = {species_name: {radius: deepcopy(inner_species_dict) for radius in range(3)}
                              for species_name in species_list}

        presence_store = {'all': deepcopy(outer_species_val)}
        similarity_store = {'all': deepcopy(outer_species_val)}
        prediction_store = {'all': deepcopy(outer_species_val)}
        correlation_store = {'all': deepcopy(outer_species_dict)}
        linear_model_store = {'all': deepcopy(outer_species_dict)}

        for habitat_type_num in habitat_type_nums:
            presence_store[habitat_type_num] = deepcopy(outer_species_val)
            similarity_store[habitat_type_num] = deepcopy(outer_species_val)
            prediction_store[habitat_type_num] = deepcopy(outer_species_val)
            correlation_store[habitat_type_num] = deepcopy(outer_species_dict)
            linear_model_store[habitat_type_num] = deepcopy(outer_species_dict)

        # Now for each subnetwork:
        for network_key in sub_networks.keys():

            # Check for validity
            current_num_patches = sub_networks[network_key]["num_patches"]
            if current_num_patches > 0:
                current_norm_pop_arrays = sub_networks[network_key]["normalised_population_arrays"]

                # iterate predictor species
                for species_1_index, species_1_name in enumerate(species_list):

                    # iterate ball radius
                    for species_1_ball_radius in range(3):
                        species_1_pop_vector = current_norm_pop_arrays[species_1_ball_radius][:, species_1_index]

                        # check for non-zero predictor population (given this ball size)
                        if np.sum(species_1_pop_vector) > 0.0:

                            # iterate response species (including self)
                            for species_2_index, species_2_name in enumerate(species_list):

                                # iterate ball radius
                                for species_2_ball_radius in range(3):
                                    species_2_pop_vector = current_norm_pop_arrays[
                                                               species_2_ball_radius][:, species_2_index]

                                    # ANALYSIS:
                                    #

                                    # A. Presence prediction
                                    species_1_pres_vector = np.zeros(current_num_patches)
                                    species_2_pres_vector = np.zeros(current_num_patches)
                                    for patch_index in range(current_num_patches):
                                        species_1_pres_vector[patch_index] = int(
                                            species_1_pop_vector[patch_index] > 0.0)
                                        species_2_pres_vector[patch_index] = int(
                                            species_2_pop_vector[patch_index] > 0.0)
                                    presence = np.dot(
                                        species_1_pres_vector, species_2_pres_vector) / np.sum(species_1_pres_vector)
                                    presence_store[network_key][species_1_name][
                                        species_1_ball_radius][species_2_name][species_2_ball_radius] = presence

                                    # B. Similarity
                                    # sum(sqrt(vector element-wise multiplication))/(sqrt(sum(predictor)*sum(response)))
                                    if np.sum(species_1_pop_vector) == 0.0 or np.sum(species_2_pop_vector) == 0.0:
                                        similarity = 0.0
                                    else:
                                        similarity = np.sum(np.sqrt(np.multiply(
                                            species_1_pop_vector, species_2_pop_vector))) / (
                                                         np.sqrt(np.sum(species_1_pop_vector) * np.sum(
                                                             species_2_pop_vector)))
                                    similarity_store[network_key][species_1_name][
                                        species_1_ball_radius][species_2_name][species_2_ball_radius] = similarity

                                    # C. Prediction
                                    # sum(sqrt(vector element-wise multiplication))/sum(predictor)
                                    prediction = np.sum(np.sqrt(np.multiply(species_1_pop_vector, species_2_pop_vector))
                                                        ) / np.sum(species_1_pop_vector)
                                    prediction_store[network_key][species_1_name][
                                        species_1_ball_radius][species_2_name][species_2_ball_radius] = prediction

                                    # D. (Two) correlation coefficients
                                    # must first test for 'nearly constant' vectors to avoid warnings
                                    species_1_near_constant = np.var(species_1_pop_vector
                                                                     ) < 1e-13 * abs(np.mean(species_1_pop_vector))
                                    species_2_near_constant = np.var(species_2_pop_vector
                                                                     ) < 1e-13 * abs(np.mean(species_2_pop_vector))
                                    is_cc_auto_fail = (species_1_near_constant or species_2_near_constant or np.var(
                                        species_1_pop_vector) <= 0.0 or np.var(species_2_pop_vector) <= 0.0 or
                                                       len(species_1_pop_vector) < 2 or len(species_2_pop_vector) < 2)
                                    is_cc_success = 0
                                    pearson_cc, pearson_p, spearman_rho, spearman_p = [0.0, 0.0, 0.0, 0.0]
                                    # for correlation coefficients to work, we need two vectors of at least two pairs
                                    # and at least some variance
                                    if not is_cc_auto_fail:
                                        try:
                                            pearson_cc, pearson_p = pearsonr(species_1_pop_vector, species_2_pop_vector)
                                            spearman_rho, spearman_p = spearmanr(
                                                species_1_pop_vector, species_2_pop_vector)
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
                                    correlation_store[network_key][species_1_name][
                                        species_1_ball_radius][species_2_name][
                                        species_2_ball_radius] = corr_coefficients

                                    # E. Linear model
                                    # note that we could use curve_fit() for a more general model specification here
                                    linear_model = linear_model_report(
                                        x_val=species_1_pop_vector, y_val=species_2_pop_vector,
                                        is_record_vectors=is_record_lm_vectors, model_type_str="lin-lin")
                                    linear_model_store[network_key][species_1_name][
                                        species_1_ball_radius][species_2_name][species_2_ball_radius] = linear_model

        # output five nested dictionaries of results
        return presence_store, similarity_store, prediction_store, correlation_store, linear_model_store

    def complexity_analysis(self, sub_networks, corresponding_binary, is_record_lm_vectors):
        # This function takes a set of sub_network partitions and, if it is possible to do so with at least three data
        # points after random sampling (multiple attempts for each) of several random clusters of connected patches
        # within the sub_network, conducts a linear regression to analyse:
        # - the SAR (species abundance relation, i.e. how does species diversity increase with patch size)
        # - possibly two forms of complexity/information dimension:
        #       - a binary version based on community states, only calculated if a corresponding_binary matrix is
        #           passed in, as we assume sub_networks to contain population sizes rather than binary occupancies.
        #           The expectation is that this is conducted only for the final meta-community snapshot, and not for
        #            time-averaged versions for which we would also conduct population (rather than binary) analysis.
        #       - a population-weighted version using the sub_networks' populations, normalised for each sub_network.
        complexity_report = {}

        # need to treat each sub_network entirely separately
        for network_key in sub_networks.keys():
            current_num_patches = sub_networks[network_key]["num_patches"]
            max_delta = int(min(current_num_patches / 4, np.sqrt(current_num_patches)))

            # do we have at least three data points for reasonable size of clusters IN THIS SUB_NETWORK?
            if max_delta > 2:

                # prepare holding arrays
                successful_clusters = np.zeros(max_delta)
                species_diversity = np.zeros(max_delta)
                binary_complexity = np.zeros(max_delta)
                population_weighted_complexity = np.zeros(max_delta)

                # iterate over cluster radius
                for delta in range(1, max_delta + 1):
                    num_clusters = 10

                    # iterate over clusters of this radius
                    for cluster_attempt in range(num_clusters):
                        cluster, is_success = generate_cluster(sub_network=sub_networks[network_key], size=delta)

                        if is_success:
                            successful_clusters[delta - 1] += 1

                            # determine total diversity in cluster
                            species_diversity[delta - 1] += self.count_diversity(
                                sub_network=sub_networks[network_key],
                                cluster=cluster)

                            # ----- Complexity (information dimension) ----- #
                            #
                            # For this we compare all unique pairs of patches within the cluster, and sum the total of
                            # their differences in state (either binary or weighted relative to the
                            # sub-network-species-normalised populations).
                            # Then the complexity dimension examines how the rate of complexity increase compares with
                            # the rate of the increase in cluster size.

                            # determine binary complexity within cluster
                            if corresponding_binary is not None:
                                # NOTE: we ONLY conduct this analysis for the final and not the time-averaged version,
                                # passing in the corresponding final-occupancy-based sub_network
                                binary_complexity[delta - 1] += determine_complexity(
                                    sub_network=corresponding_binary[network_key],
                                    cluster=cluster,
                                    is_normalised=False)

                            # determine population-weighted complexity within cluster
                            population_weighted_complexity[delta - 1] += determine_complexity(
                                sub_network=sub_networks[network_key],
                                cluster=cluster,
                                is_normalised=True)

                    # normalise output arrays
                    if successful_clusters[delta - 1] > 0:
                        species_diversity[delta - 1] = species_diversity[delta - 1] / successful_clusters[delta - 1]
                        binary_complexity[delta - 1] = binary_complexity[delta - 1] / successful_clusters[delta - 1]
                        population_weighted_complexity[delta - 1] = population_weighted_complexity[
                                                                        delta - 1] / successful_clusters[delta - 1]

                # prepare x-dimension array
                x_val = np.log(np.arange(1, max_delta + 1))

                # exclude any zero-size elements
                is_pass = False
                while len(x_val) > 0 and not is_pass:
                    is_pass = True
                    for delta_index in range(len(x_val)):
                        if successful_clusters[delta_index] == 0:
                            is_pass = False
                            # remove this entry from all vectors
                            x_val = np.delete(x_val, delta_index, axis=0)
                            successful_clusters = np.delete(successful_clusters, delta_index, axis=0)
                            species_diversity = np.delete(species_diversity, delta_index, axis=0)
                            binary_complexity = np.delete(binary_complexity, delta_index, axis=0)
                            population_weighted_complexity = np.delete(population_weighted_complexity,
                                                                       delta_index, axis=0)
                            break
            else:
                x_val = np.array([])
                species_diversity = np.array([])
                binary_complexity = np.array([])
                population_weighted_complexity = np.array([])

            # check arrays still have sufficiently-many entries
            if max_delta > 2 and len(x_val) > 2:
                # linear regressions
                #
                # ensure only non-zero values of diversity and complexity
                if not (species_diversity == 0).any():
                    sar_y_val = np.log(species_diversity)
                    lm_sar = linear_model_report(x_val=x_val, y_val=sar_y_val,
                                                 is_record_vectors=is_record_lm_vectors, model_type_str="log-log")
                    # for the SAR, fitted relationship is diversity = exp(intercept) * size ^ (gradient)
                else:
                    lm_sar = {"is_success": 0}

                # ignore first entry (patch size 0 - useful for diversity, but meaningless for intra-cluster complexity)
                if not (binary_complexity[1:] == 0).any():
                    bin_comp_y_val = np.log(binary_complexity[1:])
                    lm_binary_complexity = linear_model_report(x_val=x_val[1:], y_val=bin_comp_y_val,
                                                               is_record_vectors=is_record_lm_vectors,
                                                               model_type_str="log-log")
                    # complexity dimension is -gradient
                    lm_binary_complexity["complexity_dimension"] = -lm_binary_complexity["slope"]
                else:
                    lm_binary_complexity = {"is_success": 0}

                if not (population_weighted_complexity[1:] == 0).any():
                    pop_weight_y_val = np.log(population_weighted_complexity[1:])
                    lm_pop_weighted_complexity = linear_model_report(x_val=x_val[1:], y_val=pop_weight_y_val,
                                                                     is_record_vectors=is_record_lm_vectors,
                                                                     model_type_str="log-log")
                    # complexity dimension is -gradient
                    lm_pop_weighted_complexity["complexity_dimension"] = -lm_pop_weighted_complexity["slope"]
                else:
                    lm_pop_weighted_complexity = {"is_success": 0}

                # record
                complexity_report[network_key] = {
                    "is_cluster_success": 1,
                    "lm_sar": lm_sar,
                    "lm_binary_complexity": lm_binary_complexity,
                    "lm_pop_weighted_complexity": lm_pop_weighted_complexity,
                }
            else:
                complexity_report[network_key] = {
                    "is_cluster_success": 0,
                    "lm_sar": {},
                    "lm_binary_complexity": {},
                    "lm_pop_weighted_complexity": {},
                }

        return complexity_report

    def count_diversity(self, sub_network, cluster):
        # count the number of species present in the population array of this sub_network, whose relatively-nth patches
        # are indexed by the list "cluster" - cluster does NOT contain inherent patch numbers (unless the sub_network
        # is 'all' and zero patches have been deleted.)
        found_species_set = set()
        for patch_index in cluster:
            population_vector = sub_network["population_arrays"][0][patch_index, :]
            for possible_species_index in range(len(population_vector)):
                this_species_min = self.species_set["list"][possible_species_index].minimum_population_size
                if population_vector[possible_species_index] > this_species_min:
                    found_species_set.add(possible_species_index)
        diversity = len(found_species_set)
        return diversity

    def generate_sub_networks(self, population_array, patch_habitat, habitat_type_nums):
        # For each subnetwork
        # 	– Determine the adjacency matrix
        #   – Determine the vectors of radius-1 and radius-2 population sizes
        # 	– Determine the network-normalised (by radius-specific local maximum) population vectors

        # produce habitat subnetworks first - including adjacency arrays (store them in dict with same network keys)
        sub_networks = {}
        sub_network_list = [x for x in habitat_type_nums]
        sub_network_list.append('all')
        for network_key in sub_network_list:
            # need adjacency matrix and population array
            temp_num_patches = len(self.current_patch_list)
            temp_adjacency = deepcopy(self.patch_adjacency_matrix)
            temp_population = deepcopy(population_array)
            temp_patch_habitat = deepcopy(patch_habitat)

            # if restricting to a single-habitat sub-network, iterate to eliminate unwanted rows and columns
            if network_key != 'all':
                is_pass = False
                while temp_num_patches > 0 and not is_pass:
                    is_pass = True
                    for temp_patch_index in range(temp_num_patches):
                        if temp_patch_habitat[temp_patch_index] != network_key:
                            # remove
                            temp_patch_habitat.pop(temp_patch_index)
                            temp_population = np.delete(temp_population, temp_patch_index, axis=0)
                            temp_adjacency = np.delete(temp_adjacency, temp_patch_index, axis=0)
                            temp_adjacency = np.delete(temp_adjacency, temp_patch_index, axis=1)
                            # failed to cycle uninterrupted
                            is_pass = False
                            temp_num_patches -= 1
                            break

            # check for non-zero size of sub-network:
            if temp_num_patches > 0:
                # now generate the radius-averaged population vectors for each species in this habitat sub-network
                radius_one_population_array = np.zeros(np.shape(temp_population))
                radius_two_population_array = np.zeros(np.shape(temp_population))
                for temp_patch_index in range(temp_num_patches):
                    for ball_radius in range(1, 3):
                        if ball_radius == 1:
                            # ball (radius 1)
                            ball_patch_indices = list(np.where(temp_adjacency[temp_patch_index, :] != 0)[0])
                            ball_size = len(ball_patch_indices)
                        elif ball_radius == 2:
                            # ball (radius 2)
                            composite_adjacency = np.matmul(temp_adjacency, temp_adjacency)
                            ball_patch_indices = list(np.where(composite_adjacency[temp_patch_index, :] != 0)[0])
                            ball_size = len(ball_patch_indices)
                        else:
                            raise Exception("Invalid ball radius.")
                        # now iterate over each species and take the average population over the elements of the ball
                        for species_index in range(np.shape(temp_population)[1]):
                            ball_population = 0.0
                            for ball_index in ball_patch_indices:
                                ball_population += temp_population[ball_index, species_index]
                            if ball_radius == 1:
                                radius_one_population_array[
                                    temp_patch_index, species_index] = ball_population / ball_size
                            elif ball_radius == 2:
                                radius_two_population_array[
                                    temp_patch_index, species_index] = ball_population / ball_size
                            else:
                                raise Exception("Invalid ball radius.")

                # For each radius, identify max local population for each species and create normalised pop. matrix:
                population_array_dict = {
                    0: temp_population,
                    1: radius_one_population_array,
                    2: radius_two_population_array,
                }
                normalised_pop_array_dict = {}
                for ball_radius in range(3):
                    normalised_pop_array = np.zeros(np.shape(population_array_dict[ball_radius]))
                    for species_index in range(np.shape(population_array_dict[ball_radius])[1]):
                        species_max_population = np.max(population_array_dict[ball_radius][:, species_index])
                        if species_max_population > 0.0:
                            normalised_pop_array[:, species_index] = population_array_dict[ball_radius][
                                                                     :, species_index] / species_max_population
                    normalised_pop_array_dict[ball_radius] = deepcopy(normalised_pop_array)

                # now store the single-habitat subnetwork
                sub_networks[network_key] = deepcopy({
                    "num_patches": temp_num_patches,
                    "population_arrays": population_array_dict,
                    "normalised_population_arrays": normalised_pop_array_dict,
                    "adjacency_array": temp_adjacency,
                })
            else:
                sub_networks[network_key] = {"num_patches": 0}
        return sub_networks
