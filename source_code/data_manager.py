from source_code.data_save_functions import *
import shutil
import os


# ----------------------------- FUNCTIONS USED IN SYSTEM INITIALISATION ----------------------- #

def generate_simulation_number(minimum=99, save_data=True, is_sub_folders=False, sub_folder_capacity=100):
    # finds the next simulation number, and generates the required folder structure

    sim_path = '../results'  # over-written default to avoid "possibly unassigned" warning
    sim_number = 100  # over-written default to avoid "possibly unassigned" warning
    if is_sub_folders:
        current_parent_folder_num = 0
        next_parent_folder_num = 1
        is_next_parent_folder_exists = True

        # check if the first parent folder par_0/ exists already
        if os.path.exists(f'results/par_0'):

            # does the next iteration of parent folder exist?
            while is_next_parent_folder_exists:
                next_parent_folder_num = current_parent_folder_num + 1
                next_parent_folder_path = f'results/par_{next_parent_folder_num}'
                if not os.path.exists(next_parent_folder_path):
                    is_next_parent_folder_exists = False
                else:
                    current_parent_folder_num = next_parent_folder_num

            # now identify the highest folder of this highest sub_folder
            folder_list = os.listdir(f'results/par_{current_parent_folder_num}')
            folder_no_hidden = [f for f in folder_list if not f.startswith('.')]  # remove hidden files
            folder_int_list = [int(folder_str) for folder_str in folder_no_hidden]
            sim_number = max(folder_int_list) + 1

            # now check if it has met (or exceeded) capacity
            if len(folder_int_list) >= sub_folder_capacity:
                # if so, need to create the path to the next (not currently-existing) sub_folder
                sim_path = f'results/par_{next_parent_folder_num}/{sim_number}'
            else:
                # there is capacity, so just put it in this current sub_folder
                sim_path = f'results/par_{current_parent_folder_num}/{sim_number}'

        else:
            # if no parent folders exist, just default to par_0/100
            sim_path = f'results/par_0/{minimum+1}'

    else:
        # No sub_folder system, just the folder located straight in the root result/ directory
        sim_number = minimum
        is_folder_used = True
        while is_folder_used:
            sim_number += 1
            sim_path = f'results/{sim_number}'
            if not os.path.exists(sim_path):
                is_folder_used = False

    if save_data:
        os.makedirs(sim_path)

    return sim_number, sim_path


def write_initial_files(parameters, metadata, sim_path, parameters_filename):
    parameters_file = f"{sim_path}/parameters.json"
    dump_json(data=parameters, filename=parameters_file)
    metadata_file = f"{sim_path}/metadata.json"
    metadata["write_time"] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dump_json(data=metadata, filename=metadata_file)
    # now copy the entire parameter scripts to ensure that we really have everything for reproducibility
    output_directory = sim_path
    shutil.copy(parameters_filename, output_directory)
    shutil.copy("parameters_species_repository.py", output_directory)


def save_adj_variables(patch_list, spatial_set_number):
    # Save patch.stepping_stone_list, patch.species_movement_scores, patch.adjacency_lists in the spatial
    # network description folder
    dir_path = f'spatial_data_files/test_{spatial_set_number}/adj_variables/'
    for patch in patch_list:
        sms_file = dir_path + f'species_movement_scores/sms_{patch.number}.json'
        dump_json(data=patch.species_movement_scores, filename=sms_file)
        ssl_file = dir_path + f'stepping_stone_list/ssl_{patch.number}.json'
        dump_json(data=patch.stepping_stone_list, filename=ssl_file)
        al_file = dir_path + f'adjacency_list/al_{patch.number}.json'
        dump_json(data=patch.adjacency_lists, filename=al_file)


def load_adj_variables(patch_list, spatial_set_number):
    # Load initial patch.stepping_stone_list, patch.species_movement_scores, patch.adjacency_lists from the spatial
    # network description folder
    dir_path = f'spatial_data_files/test_{spatial_set_number}/adj_variables/'
    for patch in patch_list:
        sms_file = dir_path + f'species_movement_scores/sms_{patch.number}.json'
        patch.species_movement_scores = load_json(input_file=sms_file)
        ssl_file = dir_path + f'stepping_stone_list/ssl_{patch.number}.json'
        patch.stepping_stone_list = load_json(input_file=ssl_file)
        al_file = dir_path + f'adjacency_list/al_{patch.number}.json'
        patch.adjacency_lists = load_json(input_file=al_file)


def save_reserve_list(reserve_list, spatial_set_number):
    # Save reserve list in the spatial network description folder
    file_path = f'spatial_data_files/test_{spatial_set_number}/reserve_list/reserve_list.json'
    dump_json(data=reserve_list, filename=file_path)


def load_reserve_list(spatial_set_number):
    # Load reserve list from the spatial network description folder
    file_path = f'spatial_data_files/test_{spatial_set_number}/reserve_list/reserve_list.json'
    reserve_list = load_json(input_file=file_path)
    return reserve_list


# ----------------------------- SAVING DATA AT THE END OF THE SIMULATION ----------------------- #

def print_key_outputs_to_console(simulation_obj):
    # prints the final and averaged-final
    #   - presence/absence
    #   - local population size, standard deviation, and periodicity
    #   - enter/leave and source and sink value
    # for each local population of each species in each patch
    #
    # update present net_enter, occupancy etc.
    update_local_population_nets(system_state=simulation_obj.system_state)
    # Increase the limit to print full large arrays (necessary for adjacency matrices).
    np.set_printoptions(threshold=sys.maxsize)

    print("\n******************** SIMULATION OUTPUTS: BEGIN ********************\n")
    print("Numpy seed:")
    print(f"{simulation_obj.metadata['numpy_seed']}")
    print("Random seed:")
    print(f"{simulation_obj.metadata['random_seed']}")
    print("\n********** SPATIAL NETWORK DESCRIPTION **********\n")
    print(f"{len(simulation_obj.system_state.patch_list)}, "
          f"{len(simulation_obj.system_state.habitat_type_dictionary)}")
    print("Habitat species traversal:")
    print(simulation_obj.system_state.habitat_species_traversal)
    print("Habitat species feeding:")
    print(simulation_obj.system_state.habitat_species_feeding)
    print("Habitat amounts (final):")
    for habitat_num, habitat_history in simulation_obj.system_state.habitat_amounts_history.items():
        print(f"{habitat_num}: {habitat_history[max(x for x in habitat_history)]}")
    print("Posterior habitat spatial auto-correlation (final):")
    print("Regular:")
    print(simulation_obj.system_state.habitat_regular_auto_correlation_history[
              max(x for x in simulation_obj.system_state.habitat_regular_auto_correlation_history)])
    print("Normalised:")
    print(simulation_obj.system_state.habitat_spatial_auto_correlation_history[
              max(x for x in simulation_obj.system_state.habitat_spatial_auto_correlation_history)])
    print("LCCs (final):")
    print("All:", simulation_obj.system_state.patch_lcc_history['all'][
        max(x for x in simulation_obj.system_state.patch_lcc_history['all'])])
    print("Same:", simulation_obj.system_state.patch_lcc_history['same'][
        max(x for x in simulation_obj.system_state.patch_lcc_history['same'])])
    print("Different:", simulation_obj.system_state.patch_lcc_history['different'][
        max(x for x in simulation_obj.system_state.patch_lcc_history['different'])])
    print("Degree (final):")
    print(simulation_obj.system_state.patch_degree_history[
              max(x for x in simulation_obj.system_state.patch_degree_history)])
    print("Size:")
    print(simulation_obj.system_state.patch_size_history[
              max(x for x in simulation_obj.system_state.patch_size_history)])
    print("Quality:")
    print(simulation_obj.system_state.patch_quality_history[
              max(x for x in simulation_obj.system_state.patch_quality_history)])
    print("Centrality:")
    print(simulation_obj.system_state.patch_centrality_history[
              max(x for x in simulation_obj.system_state.patch_centrality_history)])
    print("Patch adjacency matrix (final):")
    print(simulation_obj.system_state.patch_adjacency_matrix)
    print("\n********** LOCAL POPULATION OUTPUTS **********\n")
    for patch in simulation_obj.system_state.patch_list:
        for local_population in patch.local_populations.values():
            output_str = f"{simulation_obj.sim_number}, {patch.number}, {patch.habitat_type_num}, " \
                         f"{local_population.name}, {local_population.occupancy}, {local_population.population}, " \
                         f"{local_population.internal_change}, {local_population.population_enter}, " \
                         f"{local_population.population_leave}, {local_population.source}, {local_population.sink}, " \
                         f"{local_population.maximum_foraging_distance}, " \
                         f"{local_population.weighted_foraging_distance}, " \
                         f"{local_population.average_population}, {local_population.average_internal_change}, " \
                         f"{local_population.average_population_enter}, {local_population.average_population_leave}, " \
                         f"{local_population.average_source}, {local_population.average_sink}, " \
                         f"{local_population.st_dev_population}, {local_population.max_abs_population}, " \
                         f"{local_population.population_period_weak}, {local_population.population_period_med}, " \
                         f"{local_population.population_period_strong}, " \
                         f"{local_population.recent_occupancy_change_frequency};"
            print(output_str)
    if simulation_obj.parameters["plot_save_para"]["IS_PRINT_DISTANCE_METRICS_TO_CONSOLE"]:
        print("\n********** SPECIES AND DISTANCE METRIC OUTPUTS **********\n")
        print('{')
        for counter, (key, value) in enumerate(simulation_obj.system_state.distance_metrics_store.items()):
            output_str = f'"METRIC_KEY_{key}": ' \
                         f'{json.dumps(value, ensure_ascii=True, default=set_default, skipkeys=True)}'
            if counter == len(simulation_obj.system_state.distance_metrics_store) - 1:
                is_final_item = True
            else:
                is_final_item = False
            print(format_dictionary_to_JSON_string(output_str, is_final_item=is_final_item, is_indenting=False))
        print('}')

    print("\n******************** SIMULATION OUTPUTS: END ********************\n")
    # Restore default configuration.
    np.set_printoptions(threshold=1000)


def save_all_data(simulation_obj):
    from source_code.data_plot_functions import global_species_time_series_properties
    # Access required properties:
    metadata = simulation_obj.metadata
    species_set = simulation_obj.system_state.species_set
    patch_list = simulation_obj.system_state.patch_list
    parameters = simulation_obj.parameters
    sim_path = simulation_obj.sim_path
    step = simulation_obj.system_state.step
    perturbation_history = simulation_obj.system_state.perturbation_history
    current_num_patches_history = simulation_obj.system_state.current_num_patches_history

    # update present net_enter, occupancy etc.
    update_local_population_nets(system_state=simulation_obj.system_state)

    # ---- Then enact the saving procedure: ---- #
    #
    # Core files:
    write_parameters_file(parameters=parameters, sim_path=sim_path, step=step)
    write_metadata_file(metadata=metadata, sim_path=sim_path, step=step)
    # Low-detail output on species behaviours:
    write_average_population_data(patch_list=patch_list, sim_path=sim_path, step=step)
    write_perturbation_history_data(perturbation_history=perturbation_history, sim_path=sim_path, step=step)
    global_species_time_series_properties(patch_list=patch_list, species_set=species_set, parameters=parameters,
                                          sim_path=sim_path, step=step,
                                          current_num_patches_history=current_num_patches_history,
                                          is_save_plots=False, is_save_data=True)
    # All species full local population size, internal change, and dispersal in .CSVs for each local_pop object:
    if simulation_obj.parameters["plot_save_para"]["IS_SAVE_LOCAL_POP_HISTORY_CSV"]:
        write_population_history_data(patch_list=patch_list, sim_path=sim_path, step=step)
    # JSON file of system_state with distance-metrics and histories of network-average properties (e.g. local
    # biodiversity, patch size and quality) and perturbation history:
    if simulation_obj.parameters["plot_save_para"]["IS_SAVE_SYSTEM_STATE_DATA"]:
        write_system_state(system_state=simulation_obj.system_state, sim_path=sim_path, step=step)
    # JSON files with the full history of each patch and additional JSON per local population:
    if simulation_obj.parameters["plot_save_para"]["IS_SAVE_PATCH_DATA"]:
        write_patch_list_local_populations(patch_list=patch_list, sim_path=sim_path, step=step,
                                           is_save_local_populations=simulation_obj.parameters[
                                               "plot_save_para"]["IS_SAVE_PATCH_LOCAL_POP_DATA"])
    if simulation_obj.parameters["plot_save_para"]["IS_PICKLE_SAVE"]:
        # Pickle save of the Python objects
        pickle_save(simulation_obj=simulation_obj, sim_path=sim_path, step=step)
    if simulation_obj.parameters["plot_save_para"]["IS_SAVE_DISTANCE_METRICS"]:
        # produce JSON of species and community distribution analysis
        distance_metrics_save(simulation_obj=simulation_obj, sim_path=sim_path, step=step)
    if simulation_obj.parameters["plot_save_para"]["IS_SAVE_CURRENT_MOVE_SCORES"]:
        # save a JSON file per patch containing the SMS of each species to each other patch
        write_current_species_movement_scores(patch_list=patch_list, sim_path=sim_path, step=step)


# ----------------------------- SAVING PLOTS AT THE END OF THE SIMULATION ----------------------- #

def all_plots(simulation_obj):
    from source_code.data_plot_functions import (plot_network_properties, create_time_series_plot,
                                     plot_current_local_population_attribute, plot_accessible_sub_graphs,
                                     global_species_time_series_properties)
    #
    # This should in principle be callable at any step in the simulation, and not just at the end.
    #
    # Note that we may add capabilities to plot global properties/histories, stored in the sim_obj and not in sys_state
    #
    # Access required properties:
    species_set = simulation_obj.system_state.species_set
    patch_list = simulation_obj.system_state.patch_list
    parameters = simulation_obj.parameters
    sim_path = simulation_obj.sim_path
    step = simulation_obj.system_state.step
    current_num_patches_history = simulation_obj.system_state.current_num_patches_history

    # update present net_enter, occupancy etc.
    update_local_population_nets(system_state=simulation_obj.system_state)

    # ---- TYPE I: Patch plot (heat maps) of non-species-specific properties of the spatial network ---- #
    adjacency_path_list = create_adjacency_path_list(
        patch_list=patch_list,
        patch_adjacency_matrix=simulation_obj.system_state.patch_adjacency_matrix)
    plot_network_properties(patch_list=patch_list, sim_path=sim_path, step=step,
                            adjacency_path_list=adjacency_path_list, is_biodiversity=True,
                            is_reserves=True, is_partition=True, is_label_habitat_patches=True, is_retro=False)

    # ---- TYPE II: Time-series of non-species-specific properties of the spatial network ---- #
    #
    # (note that most patch properties should be relative to the currently-present patches only)
    # - every TS: global biodiversity
    # - every TS: mean and s.d. local biodiversity
    # - every TS: number of patches present
    every_step_dict = {
        "global_biodiversity": {"y_label": "Global Biodiversity",
                                "data": [simulation_obj.system_state.global_biodiversity_history],
                                "legend": None,
                                "shading": False},
        "local_biodiversity": {"y_label": "Local Biodiversity",
                               "data": [[x[0] for x in simulation_obj.system_state.local_biodiversity_history],
                                        [x[1] for x in simulation_obj.system_state.local_biodiversity_history],
                                        [x[2] for x in simulation_obj.system_state.local_biodiversity_history]],
                               "legend": None,
                               "shading": True},
        "num_patches": {"y_label": "Number of current patches",
                        "data": [simulation_obj.system_state.current_num_patches_history],
                        "legend": None,
                        "shading": False},
    }
    for prop_name, prop in every_step_dict.items():
        file_path = f"{sim_path}/{step}/figures/network_time_series/ts_{prop_name}.png"
        create_time_series_plot(data=prop["data"], parameters=parameters,
                                file_path=file_path, y_label=prop["y_label"], is_shading=prop["shading"])
    # The following properties are only recorded when they are changed, due to a perturbation. This therefore requires
    # some conversion from the form {step_0: value_0, step_n: value_n} to a list with values for every time-step.
    # Additional unpacking is required in cases of lists of dictionaries recorded for different habitat types, or
    # where we have recorded the tuples [{step_0: (mean-SD, mean, mean+SD), step_n: ...}].
    # - number and auto-correlation of each type of habitat
    # - total number of network connections
    # - total perturbations
    # - mean and s.d. of patch degree
    # - mean and s.d. of patch quality
    # - mean and s.d. of patch centrality
    change_step_dict = {
        "habitat_amounts": {"y_label": "Amounts of each habitat present",
                            "data": [x for x in simulation_obj.system_state.habitat_amounts_history.values()],
                            "legend": [x for x in simulation_obj.system_state.habitat_amounts_history.keys()],
                            "shading": False,
                            "is_triple": False},
        "habitat_auto_correlation": {"y_label": "Spatial auto-correlation of habitat types",
                                     "data": [simulation_obj.system_state.habitat_regular_auto_correlation_history,
                                              simulation_obj.system_state.habitat_spatial_auto_correlation_history],
                                     "legend": ["Regular", "Normalised against expectation"],
                                     "shading": False,
                                     "is_triple": False},
        "total_network_connections": {"y_label": "Total edges in the spatial network",
                                      "data": [simulation_obj.system_state.total_connections_history],
                                      "legend": None,
                                      "shading": False,
                                      "is_triple": False},
        "total_perturbations": {"y_label": "Total perturbation events",
                                "data": [simulation_obj.system_state.num_perturbations_history],
                                "legend": None,
                                "shading": False,
                                "is_triple": False},
        "patch_degree": {"y_label": "Patch degree",
                         "data": [simulation_obj.system_state.patch_degree_history],
                         "legend": None,
                         "shading": True,
                         "is_triple": True},
        "patch_lcc": {"y_label": "Patch LCC",
                      "data": [x for x in simulation_obj.system_state.patch_lcc_history.values()],
                      "legend": [x for x in simulation_obj.system_state.patch_lcc_history.keys()],
                      "shading": True,
                      "is_triple": True},
        "patch_quality": {"y_label": "Patch quality",
                          "data": [simulation_obj.system_state.patch_quality_history],
                          "legend": None,
                          "shading": True,
                          "is_triple": False},
        "patch_centrality": {"y_label": "Patch centrality",
                             "data": [simulation_obj.system_state.patch_centrality_history],
                             "legend": None,
                             "shading": True,
                             "is_triple": True},
    }

    # iterate through the different properties
    for prop_name, prop in change_step_dict.items():
        file_path = f"{sim_path}/{step}/figures/network_time_series/ts_{prop_name}.png"
        time_series_data = []
        # iterate over the dictionaries contained in a list, will be multiple if multiple series to go on same plot
        for series in prop["data"]:
            new_series = []
            # for a given dictionary of the form {step_0: value_0, step_1: value_1, ... , step_n: value_n},
            # get the keys and convert to a sorted list to ENSURE correct lowest-to-highest order
            sorted_list_of_steps = list(series.keys())
            sorted_list_of_steps.sort()
            # add the end point for the final section (need to add one more to count 'including this step')
            sorted_list_of_steps.append(step + 1)
            # now iterate over the steps where initial values or change recorded (omit the added endpoint)
            for step_num, start_step in enumerate(sorted_list_of_steps[: -1]):
                # then append the running list with the current value for the length of steps until the next change
                new_series += [series[start_step] for _ in range(sorted_list_of_steps[step_num + 1] - start_step)]
            time_series_data.append(new_series)
        create_time_series_plot(data=time_series_data, parameters=parameters, legend_list=prop["legend"],
                                file_path=file_path, y_label=prop["y_label"], is_shading=prop["shading"],
                                is_triple=prop["is_triple"])

    # Confirmed that the following retrospective plots match the pre-plots generated at the initial time-step
    # retrospective_network_plots(
    #     initial_patch_list=simulation_obj.system_state.initial_patch_list,
    #     actual_patch_list=simulation_obj.system_state.patch_list,
    #     initial_patch_adjacency_matrix=simulation_obj.system_state.initial_patch_adjacency_matrix,
    #     sim=sim, step=28)

    # ---- Type III: Patch plots (heat maps) of different species-specific properties ---- #
    for species in species_set["list"]:
        # can iterate through list of any local_population attributes that you wish to create patch plots of:
        attribute_to_plot = ["occupancy", "population", "internal_change", "net_internal", "population_enter",
                             "population_leave", "net_enter", "source", "sink",
                             "maximum_foraging_distance", "weighted_foraging_distance",
                             "average_population", "average_internal_change", "average_net_internal",
                             "average_population_enter", "average_population_leave", "average_net_enter",
                             "average_source", "average_sink",
                             "population_period_weak", "population_period_med", "population_period_strong",
                             "recent_occupancy_change_frequency"]
        for attr in attribute_to_plot:
            plot_current_local_population_attribute(species=species, patch_list=patch_list, sim_path=sim_path,
                                                    adjacency_path_list=adjacency_path_list,
                                                    attribute_name=attr, step=step)
        if parameters["plot_save_para"]["IS_PLOT_ACCESSIBLE_SUB_GRAPHS"]:
            # patch plots showing the fully-connected network sub-graphs from the POV of each species
            plot_accessible_sub_graphs(patch_list=patch_list, parameters=parameters,
                                       species=species, sim_path=sim_path, step=step)

    # ---- Type IV: Species-specific time-series ---- #
    global_species_time_series_properties(patch_list=patch_list, species_set=species_set, parameters=parameters,
                                          sim_path=sim_path, step=step,
                                          current_num_patches_history=current_num_patches_history,
                                          is_save_plots=True, is_save_data=False)
    is_local_plots = parameters["plot_save_para"]["LOCAL_PLOTS"]  # individual plot files per patch?
    if parameters["plot_save_para"]["IS_PLOT_LOCAL_TIME_SERIES"]:
        from source_code.data_plot_functions import plot_local_time_series
        plot_local_time_series(patch_list=patch_list, species_set=species_set, parameters=parameters,
                               sim_path=sim_path, step=step, is_local_plots=is_local_plots)

    # ---- Type V: Optional extras ---- #
    if parameters["plot_save_para"]["IS_PLOT_UNRESTRICTED_PATHS"]:
        # plots all reachable foraging/direct dispersal paths per species WITHOUT threshold and path-length restriction
        from source_code.data_plot_functions import plot_unrestricted_shortest_paths
        plot_unrestricted_shortest_paths(patch_list=patch_list, species_set=species_set, sim_path=sim_path, step=step)
    if parameters["plot_save_para"]["IS_BIODIVERSITY_ANALYSIS"]:
        # species-area curves (SAR analysis)
        from source_code.data_plot_functions import biodiversity_analysis
        biodiversity_analysis(patch_list=patch_list, species_set=species_set, parameters=parameters,
                              sim_path=sim_path, step=step)
    if parameters["plot_save_para"]["IS_PLOT_ADJACENCY_SUB_GRAPHS"]:
        # plots the undirected adjacency-based sub-graphs, regardless of any species ability to traverse them
        from source_code.data_plot_functions import plot_adjacency_sub_graphs
        plot_adjacency_sub_graphs(system_state=simulation_obj.system_state, sim_path=sim_path)
    if parameters["plot_save_para"]["IS_PLOT_INTERACTIONS"]:
        from source_code.data_plot_functions import plot_interactions
        plot_interactions(patch_list=patch_list, adjacency_path_list=adjacency_path_list, sim_path=sim_path, step=step)
    if parameters["plot_save_para"]["IS_PLOT_DEGREE_DISTRIBUTION"]:
        from source_code.data_plot_functions import plot_degree_distribution
        plot_degree_distribution(
            degree_distribution_history=simulation_obj.system_state.degree_distribution_history,
            degree_dist_power_law_fit_history=simulation_obj.system_state.degree_dist_power_law_fit_history,
            sim_path=sim_path, step=step)
    if parameters["plot_save_para"]["IS_PLOT_DISTANCE_METRICS_LM"]:
        from source_code.data_plot_functions import (plot_distance_metrics_lm,
                                                     complexity_plotting, partition_spectrum_plotting)
        plot_distance_metrics_lm(distance_metrics_store=simulation_obj.system_state.distance_metrics_store,
                                 sim_path=sim_path, step=step)
        complexity_plotting(distance_metrics_store=simulation_obj.system_state.distance_metrics_store,
                            sim_path=sim_path, step=step)
        partition_spectrum_plotting(distance_metrics_store=simulation_obj.system_state.distance_metrics_store,
                            sim_path=sim_path, step=step)
    print("Completed all_plots().")


# ---------------------------------- SNAPSHOTS (FOR BEFORE/AFTER PERTURBATIONS) ---------------------------------- #

def population_snapshot(system_state, sim_path, update_stored, output_figures):
    update_local_population_nets(system_state=system_state)
    if update_stored:
        update_stored_populations(patch_list=system_state.patch_list, step=system_state.step)

    # population recording
    save_current_local_population_attribute(patch_list=system_state.patch_list, sim_path=sim_path,
                                            attribute_name="population", step=system_state.step)
    if output_figures:
        from source_code.data_plot_functions import plot_current_local_population_attribute
        adjacency_path_list = create_adjacency_path_list(
            patch_list=system_state.patch_list,
            patch_adjacency_matrix=system_state.patch_adjacency_matrix)
        # population visualisation
        for species in system_state.species_set["list"]:
            plot_current_local_population_attribute(
                patch_list=system_state.patch_list,
                adjacency_path_list=adjacency_path_list,
                sim_path=sim_path,
                attribute_name="population",
                step=system_state.step,
                species=species
            )


def change_snapshot(system_state, sim_path, output_figures):
    update_local_population_nets(system_state=system_state)
    for patch in system_state.patch_list:
        for local_pop in patch.local_populations.values():
            local_pop.stored_change["absolute"] = local_pop.population - \
                                                  local_pop.stored_pop_values[system_state.step - 1]
            if local_pop.stored_pop_values[system_state.step - 1] != 0.0:
                local_pop.stored_change["relative"] = (local_pop.population -
                                                       local_pop.stored_pop_values[system_state.step - 1]) / \
                                                      local_pop.stored_pop_values[system_state.step - 1]
            else:
                local_pop.stored_change["relative"] = 0.0

    # change recording
    save_current_local_population_attribute(patch_list=system_state.patch_list,
                                            sim_path=sim_path,
                                            attribute_name="stored_change",
                                            step=system_state.step,
                                            sub_attr="absolute"
                                            )
    save_current_local_population_attribute(patch_list=system_state.patch_list,
                                            sim_path=sim_path,
                                            attribute_name="stored_change",
                                            step=system_state.step,
                                            sub_attr="relative",
                                            )
    if output_figures:
        # change visualisation
        from source_code.data_plot_functions import plot_current_local_population_attribute
        adjacency_path_list = create_adjacency_path_list(
            patch_list=system_state.patch_list,
            patch_adjacency_matrix=system_state.patch_adjacency_matrix)
        for species in system_state.species_set["list"]:
            plot_current_local_population_attribute(patch_list=system_state.patch_list,
                                                    adjacency_path_list=adjacency_path_list,
                                                    sim_path=sim_path,
                                                    attribute_name="stored_change",
                                                    step=system_state.step,
                                                    species=species,
                                                    sub_attr="absolute"
                                                    )
            plot_current_local_population_attribute(patch_list=system_state.patch_list,
                                                    adjacency_path_list=adjacency_path_list,
                                                    sim_path=sim_path,
                                                    attribute_name="stored_change",
                                                    step=system_state.step,
                                                    species=species,
                                                    sub_attr="relative"
                                                    )
