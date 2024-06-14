# generates sample spatial networks and habitat data in .CSV files for testing
import random
import shutil
import numpy as np
import os
import networkx  # https://networkx.org/documentation/stable/reference/generators.html


# ----------------------------- FOLDER PREPARATION ----------------------- #

def save_array(file, print_array):
    with open(file, mode='w') as f:
        np.savetxt(f, print_array, delimiter=', ', newline='\n', fmt='%.20f')


def check_and_create_directory(test_set, dir_path, can_overwrite_existing_dataset):
    if os.path.exists(dir_path):
        if not can_overwrite_existing_dataset:
            raise Exception("This spatial dataset already exists.")
        else:
            # overwrite: delete the existing directory then re-create a fresh empty directory
            print(f"Test set {test_set} already existed - deleting previous files and clearing directory.")
            shutil.rmtree(dir_path)
            os.makedirs(dir_path)
    else:
        os.makedirs(dir_path)


def create_description_file(desc, dir_path):
    with open(file=f'{dir_path}/description.txt', mode='w') as f:
        f.write(desc)


# ----------------------------- CONSTRUCTING THE SPATIAL NETWORK ----------------------- #

# THIS HAPPENS FIRST, AS SOME OTHER PROPERTIES MAY DEPEND ON POSITION AND ADJACENCY
def generate_patch_position_adjacency(num_patches, graph_para):
    # Automatically place them in a rectangular grid
    num_rows = int(np.ceil(np.sqrt(num_patches)))
    num_columns = int(np.ceil(num_patches / num_rows))
    position_array = np.zeros([num_patches, 2])
    for patch in range(num_patches):
        x = np.mod(patch, num_columns)
        y = np.floor(patch / num_columns)
        position_array[patch, 0] = x
        position_array[patch, 1] = y

        # Now adjacency
    graph_type = graph_para["GRAPH_TYPE"]
    adjacency_array = np.zeros([num_patches, num_patches])
    adjacency_spec = graph_para["ADJACENCY_MANUAL_SPEC"]
    if graph_type == "manual":
        # the "ADJACENCY_MANUAL_SPEC" should be a list (length = num_patches) of lists (length = num_patches)
        if adjacency_spec is not None and type(adjacency_spec) == list and len(adjacency_spec) == num_patches:
            # check dimensions and values and that the result is a symmetric matrix
            for x in range(num_patches):
                if type(adjacency_spec[x]) == list and len(adjacency_spec[x]) == num_patches:
                    for y in range(num_patches):
                        if adjacency_spec[x][y] not in [0, 1]:
                            raise Exception("Error in graph_para['ADJACENCY_MANUAL_SPEC']. "
                                            "Values should only be '0' or '1'.")
                        if adjacency_spec[x][y] != adjacency_spec[y][x]:
                            raise Exception("Error in graph_para['ADJACENCY_MANUAL_SPEC']. Matrix is not symmetric.")
                else:
                    raise Exception(f"Error in graph_para['ADJACENCY_MANUAL_SPEC']. Row {x} is not a list with "
                                    f"the correct number of columns.")
                if adjacency_spec[x][x] != 1:
                    raise Exception("Error in graph_para['ADJACENCY_MANUAL_SPEC']. Diagonals should all be 1 unless the"
                                    " patch is removed.")
            # convert list of lists from the parameters to an array
            adjacency_array = np.asarray(adjacency_spec)
        else:
            raise Exception("Error in graph_para['ADJACENCY_MANUAL_SPEC']. Incorrect number of rows.")
    else:
        # check that user was not attempting to manually specify adjacency and forgot to set the graph type
        if adjacency_spec is not None and type(adjacency_spec) == list and len(adjacency_spec) != 0:
            raise Exception("Check that graph type is `manual' or clear the manual adjacency specification.")

        if graph_type == "lattice":
            for x in range(num_patches):
                for y in range(num_patches):
                    is_suitable = False

                    # include diagonals?
                    if graph_para["IS_LATTICE_INCLUDE_DIAGONALS"]:
                        test_distance = 1.999
                    else:
                        test_distance = 1.001
                    # check regular within-grid adjacency
                    if np.linalg.norm(np.array([position_array[x, 0] - position_array[y, 0],
                                                position_array[x, 1] - position_array[y, 1]])) < test_distance:
                        is_suitable = True

                    # wrap around?
                    if graph_para["IS_LATTICE_WRAPPED"]:
                        # we need to check the minimum distance across all the types of edges of the lattice
                        # to ensure we catch all cases (esp. when diagonals are possible), we go look at both patches
                        # and determine if we will need to check their 'modulated/wrapped' x and y values.
                        if position_array[x, 0] == num_columns - 1:
                            patch_x_x_possible = [position_array[x, 0], -1]
                        else:
                            patch_x_x_possible = [position_array[x, 0]]
                        if position_array[y, 0] == num_columns - 1:
                            patch_y_x_possible = [position_array[y, 0], -1]
                        else:
                            patch_y_x_possible = [position_array[y, 0]]

                        if position_array[x, 1] == num_rows - 1:
                            patch_x_y_possible = [position_array[x, 1], -1]
                        else:
                            patch_x_y_possible = [position_array[x, 1]]
                        if position_array[y, 1] == num_rows - 1:
                            patch_y_y_possible = [position_array[y, 1], -1]
                        else:
                            patch_y_y_possible = [position_array[y, 1]]
                        # then we look at all combinations to determine if there is a shortcut
                        min_x_dist = num_columns
                        for patch_x_x in patch_x_x_possible:
                            for patch_y_x in patch_y_x_possible:
                                min_x_dist = min(min_x_dist, np.abs(patch_x_x - patch_y_x))
                        min_y_dist = num_rows
                        for patch_x_y in patch_x_y_possible:
                            for patch_y_y in patch_y_y_possible:
                                min_y_dist = min(min_y_dist, np.abs(patch_x_y - patch_y_y))
                        if np.linalg.norm(np.array([min_x_dist, min_y_dist])) < test_distance:
                            is_suitable = True

                    if is_suitable:
                        draw = np.random.binomial(n=1, p=graph_para["LATTICE_GRAPH_CONNECTIVITY"])
                        adjacency_array[x, y] = draw
                        adjacency_array[y, x] = draw
        elif graph_type == "line":
            for x in range(num_patches):
                if x > 0:
                    adjacency_array[x - 1, x] = 1
                if x < num_patches - 1:
                    adjacency_array[x, x + 1] = 1
        elif graph_type == "star":
            # all patches only adjacent to the first patch
            for x in range(num_patches):
                adjacency_array[x, 0] = 1
                adjacency_array[0, x] = 1
        elif graph_type == "random":
            for x in range(num_patches):
                for y in range(x):
                    draw = np.random.binomial(n=1, p=graph_para["RANDOM_GRAPH_CONNECTIVITY"])
                    adjacency_array[x, y] = draw
                    adjacency_array[y, x] = draw
        elif graph_type == "small_world":
            # https://networkx.org/documentation/stable/reference/generated/
            # networkx.generators.random_graphs.watts_strogatz_graph.html
            graph = networkx.watts_strogatz_graph(n=num_patches,
                                                  k=graph_para["SMALL_WORLD_NUM_NEIGHBOURS"],
                                                  p=graph_para["SMALL_WORLD_SHORTCUT_PROBABILITY"])
            adjacency_array = networkx.to_numpy_array(graph)
            if len(graph.nodes) == 0:
                raise Exception("Small World graph failed to generate - probably unsuitable number "
                                "of neighbours for this number of patches.")

        elif graph_type == "scale_free":
            # https://networkx.org/documentation/stable/reference/generated/
            # networkx.generators.directed.scale_free_graph.html
            graph = networkx.scale_free_graph(n=num_patches)  # directed
            adjacency_array = networkx.to_numpy_array(networkx.to_undirected(graph))
            adjacency_array[adjacency_array > 1] = 1
        elif graph_type == "cluster":
            # https://networkx.org/documentation/stable/reference/generated/
            # networkx.generators.random_graphs.powerlaw_cluster_graph.html
            graph = networkx.powerlaw_cluster_graph(n=num_patches,
                                                    m=graph_para["CLUSTER_NUM_NEIGHBOURS"],
                                                    p=graph_para["CLUSTER_PROBABILITY"], )
            adjacency_array = networkx.to_numpy_array(graph)
        else:
            raise Exception("Which type of graph is the spatial network?")

    # check that adjacency array has correctly generated:
    if adjacency_array.shape[0] != num_patches or adjacency_array.shape[1] != num_patches:
        raise Exception("Graph generation process has failed to produce adjacency array of size NxN "
                        "where N is the number of patches specified.")
    else:
        # ensure every patch is always considered adjacent to itself
        for x in range(num_patches):
            adjacency_array[x, x] = 1

        # ensure that the adjacency graphs are always symmetric
        for x in range(num_patches):
            for y in range(num_patches):
                if adjacency_array[x, y] == 1:
                    adjacency_array[y, x] = 1

    return adjacency_array, position_array


def generate_patch_quality(num_patches, adjacency_array, position_array, graph_para):
    quality_type = graph_para["QUALITY_TYPE"]
    max_quality = graph_para["MAX_QUALITY"]
    min_quality = graph_para["MIN_QUALITY"]
    if max_quality < min_quality or max_quality > 1.0 or min_quality < 0.0:
        raise Exception("Max and min quality parameters not as expected in [0, 1].")
    #
    # quality types are {'manual', 'random', 'auto_correlation', 'gradient'}
    if quality_type == "manual":
        # the "QUALITY_MANUAL_SPEC" should be a list of length = num_patches
        quality_spec = graph_para["QUALITY_MANUAL_SPEC"]
        if quality_spec is not None and type(quality_spec) == list and len(quality_spec) == num_patches:
            quality_array = np.zeros(shape=(num_patches, 1))
            for patch_num in range(num_patches):
                patch_quality = quality_spec[patch_num]
                # check valid
                if 0.0 <= patch_quality <= 1.0:
                    quality_array[patch_num] = patch_quality
                else:
                    raise Exception("Unsuitable patch quality is specified.")
        else:
            raise Exception("Error in the expected manual specification of patch quality (QUALITY_MANUAL_SPEC)")
    elif quality_type == "random":
        quality_array = min_quality + np.random.rand(num_patches, 1) * (max_quality - min_quality)
    elif quality_type == "auto_correlation":
        auto_correlation = graph_para["QUALITY_SPATIAL_AUTO_CORRELATION"]
        quality_array = np.zeros(shape=(num_patches, 1))
        # initial value:
        quality_array[0, 0] = np.random.rand()
        for patch in range(1, num_patches):
            # construct mean
            num_neighbours = 0
            quality_sum = 0.0
            # iterate over those neighbors who have already been assigned their quality
            for other_patch in range(patch):
                if adjacency_array[patch, other_patch] == 1:
                    num_neighbours += 1
                    quality_sum += quality_array[other_patch, 0]
            draw = np.random.rand()
            if num_neighbours > 0:
                quality_mean = quality_sum / float(num_neighbours)
                temp_value = draw + auto_correlation * (quality_mean - draw)
                # auto_corr = 0 => draw
                # auto_corr = 1 => mean
                # auto_corr = -1 => 2 draw - mean (i.e. same distance on other side of draw)
            else:
                temp_value = draw
            # check against max and min
            quality_array[patch, 0] = min(max(min_quality, temp_value), max_quality)
    elif quality_type == "gradient":
        if num_patches == 1:
            quality_array = min_quality + np.random.rand(1, 1) * (max_quality - min_quality)
        else:
            fluctuation = graph_para["QUALITY_FLUCTUATION"]
            axis = graph_para["QUALITY_AXIS"]  # x, y, x+y
            if axis == "x":
                value_vector = position_array[:, 0]
            elif axis == "y":
                value_vector = position_array[:, 1]
            elif axis == "x+y":
                value_vector = position_array[:, 0] + position_array[:, 1]
            else:
                raise Exception("Axis not chosen correctly.")
            min_pos = np.amin(value_vector)
            max_pos = np.amax(value_vector)
            if max_pos == min_pos:
                raise Exception("No variation along the axis specified.")
            else:
                quality_array = np.zeros(shape=(num_patches, 1))
                for patch in range(0, num_patches):
                    quality_array[patch, 0] = min(1.0, max(0.0, min_quality + (value_vector[patch] - min_pos) *
                                                           (max_quality - min_quality) / (max_pos - min_pos) +
                                                           fluctuation * np.random.rand()))
    else:
        raise Exception("Which type of scheme is used for patch quality in the spatial network?")
    return quality_array


def generate_patch_size(num_patches, min_patch_size, max_patch_size, graph_para):
    specified_size_list = graph_para["PATCH_SIZE_MANUAL_SPEC"]
    if specified_size_list is not None and len(specified_size_list) == num_patches:
        # in this case we have manually specified the patch sizes as a list
        patch_size_array = np.zeros(shape=(num_patches, 1))
        for patch_num in range(num_patches):
            # check is valid:
            patch_size = specified_size_list[patch_num]
            if 0.0 <= patch_size <= 1.0:
                patch_size_array[patch_num, 0] = patch_size
            else:
                raise Exception("Unsuitable patch size is specified.")
    elif specified_size_list is None:
        # otherwise generate sizes randomly
        patch_size_array = min_patch_size + (max_patch_size - min_patch_size) * np.random.rand(num_patches, 1)
    else:
        raise Exception("The graph_para option 'PATCH_SIZE_MANUAL_SPEC' should be either None or a list of length"
                        " equal to the number of patches.")
    return patch_size_array


def generate_habitat_type(generated_habitat_set, num_patches, generated_habitat_probabilities,
                          adjacency_array, graph_para):
    actual_habitat_list = list(generated_habitat_set)  # needs to be ordered for weighting the probability vector
    actual_num_habitats = len(actual_habitat_list)
    auto_correlation = graph_para["HABITAT_SPATIAL_AUTO_CORRELATION"]
    habitat_array = np.zeros(shape=(num_patches, 1))

    specified_habitat_list = graph_para["HABITAT_TYPE_MANUAL_SPEC"]
    if specified_habitat_list is not None and len(specified_habitat_list) == num_patches:
        # in this case we have manually specified the habitat types as a list
        for patch in range(num_patches):
            # check is valid:
            habitat_type_num = specified_habitat_list[patch]
            if habitat_type_num in actual_habitat_list:
                habitat_array[patch, 0] = habitat_type_num
            else:
                raise Exception("Unsuitable habitat type num is specified.")
    elif specified_habitat_list is None:
        # otherwise generate habitats probabilistically
        # initial value:
        habitat_array[0, 0] = int(random.choice(actual_habitat_list))
        # base probability vector
        if generated_habitat_probabilities is not None:
            base_probability = np.zeros(shape=(actual_num_habitats, 1))
            for habitat_type_num in actual_habitat_list:
                base_probability[habitat_type_num, 0] = generated_habitat_probabilities[habitat_type_num]
        else:
            # uniform
            base_probability = np.ones(shape=(actual_num_habitats, 1))
        # normalise
        normalised_base_probability = base_probability / np.sum(base_probability)
        for patch in range(1, num_patches):
            # construct probability array
            habitat_probability = np.zeros(shape=(actual_num_habitats, 1))
            # iterate over those neighbors who have already been assigned their habitats
            for other_patch in range(patch):
                if adjacency_array[patch, other_patch] == 1:
                    habitat_probability[actual_habitat_list.index(int(habitat_array[other_patch, 0])), 0] += 1

            # normalise the previous-patch-weighted distribution (so that auto-correlation is independent of the number
            # of patches that have already been assigned their habitat types)
            if np.sum(habitat_probability) != 0.0:
                normalised_habitat_probability = habitat_probability / np.sum(habitat_probability)
            else:
                normalised_habitat_probability = np.zeros(shape=(actual_num_habitats, 1))
            # combine both parts and weight by auto-correlation and complement respectively
            combined_habitat_probability = auto_correlation * normalised_habitat_probability + (
                    1.0 - auto_correlation) * normalised_base_probability
            # check if negative entries
            norm_combined_habitat_probability = np.zeros(shape=(actual_num_habitats, 1))
            for habitat_type in range(actual_num_habitats):
                norm_combined_habitat_probability[habitat_type] = max(0.0, combined_habitat_probability[habitat_type])
            # then re-normalise
            if np.sum(norm_combined_habitat_probability) == 0.0:
                norm_combined_habitat_probability = normalised_base_probability
            else:
                norm_combined_habitat_probability = norm_combined_habitat_probability / \
                                                    np.sum(norm_combined_habitat_probability)
            # then draw
            habitat_array[patch, 0] = actual_habitat_list[np.random.choice(
                actual_num_habitats, p=np.transpose(norm_combined_habitat_probability)[0])]
    else:
        raise Exception("The graph_para option 'HABITAT_TYPE_MANUAL_SPEC' should be either None or a list of length"
                        " equal to the number of patches.")
    return habitat_array


def generate_habitat_species_scores(num_species, num_habitats, generated_spec, score_type):
    if generated_spec[score_type]["IS_SPECIES_SCORES_SPECIFIED"]:
        # The feeding or traversal scores are manually specified for each (habitat, species) pair.
        array = np.zeros([num_habitats, num_species])
        data = generated_spec[score_type]["HABITAT_SCORES"]
        if list(data.keys()) != [x for x in range(num_habitats)]:
            raise Exception(f'Error with HABITAT_SCORES: {score_type}')
        for habitat in data:
            if len(data[habitat]) != num_species:
                raise Exception(f'Error with HABITAT_SCORES: {score_type}')
            else:
                array[habitat, :] = data[habitat]
    else:
        # Otherwise generate them from the provided interval.
        min_score = generated_spec[score_type]["MIN_SCORE"]
        max_score = generated_spec[score_type]["MAX_SCORE"]
        if 0 <= min_score < max_score <= 1.0:
            array = min_score + np.random.rand(num_habitats, num_species) * (max_score - min_score)
        else:
            raise Exception(f"Maximum and minimum values of {score_type} scores not correctly ordered in [0,1]")

    return array


def generate_all_spatial_settings(is_output_files, desc, dir_path, test_set, can_overwrite_existing_dataset,
                                  num_species, num_patches, num_habitats, graph_para, generated_habitat_set,
                                  generated_habitat_probabilities, generated_spec):
    if is_output_files:
        check_and_create_directory(test_set=test_set, dir_path=dir_path,
                                   can_overwrite_existing_dataset=can_overwrite_existing_dataset)
        create_description_file(desc, dir_path=dir_path)
    adjacency_array, position_array = generate_patch_position_adjacency(num_patches=num_patches, graph_para=graph_para)
    patch_size_array = generate_patch_size(num_patches=num_patches, min_patch_size=graph_para["MIN_SIZE"],
                                           max_patch_size=graph_para["MAX_SIZE"], graph_para=graph_para)
    patch_quality_array = generate_patch_quality(num_patches=num_patches, adjacency_array=adjacency_array,
                                                 position_array=position_array, graph_para=graph_para)
    patch_habitat_array = generate_habitat_type(generated_habitat_set=generated_habitat_set, num_patches=num_patches,
                                                generated_habitat_probabilities=generated_habitat_probabilities,
                                                adjacency_array=adjacency_array, graph_para=graph_para)
    scores_dict = {}
    for score_type in ["FEEDING", "TRAVERSAL"]:
        scores_dict[score_type.lower()] = generate_habitat_species_scores(num_species=num_species,
                                                                          num_habitats=num_habitats,
                                                                          generated_spec=generated_spec,
                                                                          score_type=score_type)
    x = 1
    if is_output_files:
        save_array(f'{dir_path}/patch_position.csv', position_array)
        save_array(f'{dir_path}/patch_adjacency.csv', adjacency_array)
        save_array(f'{dir_path}/patch_size.csv', patch_size_array)
        save_array(f'{dir_path}/patch_quality.csv', patch_quality_array)
        save_array(f'{dir_path}/patch_habitat_type.csv', patch_habitat_array)
        save_array(f'{dir_path}/habitat_species_feeding.csv', scores_dict["feeding"])
        save_array(f'{dir_path}/habitat_species_traversal.csv', scores_dict["traversal"])

    return position_array, adjacency_array, patch_size_array, patch_quality_array, \
        patch_habitat_array, scores_dict["feeding"], scores_dict["traversal"]


# ---------------------- EXECUTE ---------------------- #

def run_sample_spatial_data(parameters, is_output_files=False):
    test_set = parameters["graph_para"]["SPATIAL_TEST_SET"]
    description = parameters["graph_para"]["SPATIAL_DESCRIPTION"]
    can_overwrite_existing_dataset = True
    dir_path = f'spatial_data_files/test_{test_set}'

    print(f"Beginning generation of spatial network {test_set}.")
    # check the habitats
    master_habitat_types_set = parameters["main_para"]["HABITAT_TYPES"]
    for habitat_num in parameters["main_para"]["INITIAL_HABITAT_SET"]:
        if habitat_num not in master_habitat_types_set:
            raise Exception(f'Habitat type {habitat_num} to be used in generation but is not part of the global set.')
    for habitat_num in master_habitat_types_set:
        if habitat_num > len(master_habitat_types_set) or type(habitat_num) is not int:
            raise Exception(f'Habitat type {habitat_num} is not a suitable number.')
    for check_num in range(len(master_habitat_types_set)):
        if check_num not in master_habitat_types_set:
            raise Exception(f'{check_num} is missing from the global set of habitat types.')

    position_array, adjacency_array, patch_size_array, patch_quality_array, patch_habitat_array, \
        habitat_feeding_scores, habitat_traversal_scores \
        = generate_all_spatial_settings(
            is_output_files=is_output_files,
            desc=description,
            dir_path=dir_path,
            test_set=test_set,
            can_overwrite_existing_dataset=can_overwrite_existing_dataset,
            num_species=len(parameters["main_para"]["SPECIES_TYPES"]),
            num_patches=parameters["main_para"]["NUM_PATCHES"],
            num_habitats=len(master_habitat_types_set),
            graph_para=parameters["graph_para"],
            generated_habitat_set=parameters["main_para"]["INITIAL_HABITAT_SET"],
            generated_habitat_probabilities=parameters["main_para"]["INITIAL_HABITAT_BASE_PROBABILITIES"],
            generated_spec=parameters["main_para"]["GENERATED_SPEC"],
        )
    print(f"Test set {test_set} generation complete.")
    return position_array, adjacency_array, patch_size_array, patch_quality_array, \
        patch_habitat_array, habitat_feeding_scores, habitat_traversal_scores


# # so that this method is called when the script is executed
if __name__ == '__main__':
    import sys
    import importlib
    if len(sys.argv) > 1:
        # optionally pass in an argument specifying the particular parameters_???.py file to use
        parameters_file = importlib.import_module(sys.argv[1])
    else:
        # otherwise use "parameters.py" as the default
        parameters_file = importlib.import_module("parameters")
    master_para = getattr(parameters_file, "master_para")
    run_sample_spatial_data(parameters=master_para, is_output_files=True)
