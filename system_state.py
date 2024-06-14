import numpy as np
from copy import deepcopy


def tuple_builder(property_list):
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
        self.patch_lcc_history = {}
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
        self.patch_lcc_history[self.step] = tuple_builder(current_lcc_list)
        self.total_connections_history[self.step] = int(
            (np.sum(current_degree_list) - len(self.current_patch_list)) / 2)

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
                    lcc_list.append(lcc)
        print(lcc_list)
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
        num_triangles = 0
        num_closed_triangles = 0
        list_of_neighbours = list(patch.set_of_adjacent_patches)
        lcc = 0.0
        if len(list_of_neighbours) > 1:
            for n_1_index, neighbour_1 in enumerate(list_of_neighbours):
                if neighbour_1 in self.current_patch_list and neighbour_1 != patch.number:
                    for neighbour_2 in list_of_neighbours[n_1_index+1:]:  # necessarily n_1 not same as n_2
                        if neighbour_2 in self.current_patch_list and neighbour_2 != patch.number:
                            num_triangles += 1
                            if neighbour_2 in self.patch_list[neighbour_1].set_of_adjacent_patches:
                                num_closed_triangles += 1
            if num_triangles > 0:
                lcc = float(num_closed_triangles)/float(num_triangles)
        # update patch records
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
