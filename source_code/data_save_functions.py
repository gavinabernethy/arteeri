from source_code.data_core_functions import *
from datetime import datetime


# -------------------------------- FUNCTIONS FOR SAVING AND LOADING THE ENTIRE SYSTEM -------------------------------- #

def write_parameters_file(parameters, sim_path, step):
    parameters_file = f"{sim_path}/{step}/parameters.json"
    dump_json(data=parameters, filename=parameters_file)


def write_metadata_file(metadata, sim_path, step):
    metadata_file = f"{sim_path}/{step}/metadata.json"
    metadata["write_time"] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dump_json(data=metadata, filename=metadata_file)


def write_perturbation_history_data(perturbation_history, sim_path, step):
    pert_history_file = f"{sim_path}/{step}/perturbation_history.json"
    dump_json(data=perturbation_history, filename=pert_history_file)


def write_current_species_movement_scores(patch_list, sim_path, step):
    for patch in patch_list:
        species_movement_scores_file = f"{sim_path}/{step}/data/species_movement_scores/sms_{patch.number}.json"
        dump_json(data=patch.species_movement_scores, filename=species_movement_scores_file)
    # recall that in previous versions these were stored as an array, and thus had to be converted to a
    # nested dictionary first


def write_average_population_data(patch_list, sim_path, step):
    print("Saving text and JSON file of all patches averaged pop.s, internal change, dispersal for each local_pop.")
    # build a dictionary
    average_population = {}
    # Text file is more readable but less easy to analyse:
    txt_file_name = f"{sim_path}/{step}/data/average_populations.txt"
    with safe_open_w(txt_file_name) as f:
        for patch in patch_list:
            this_patch = {}
            for local_pop in patch.local_populations.values():
                this_patch[local_pop.name] = [local_pop.average_population, local_pop.average_internal_change,
                                              local_pop.average_population_enter, local_pop.average_population_leave,
                                              local_pop.average_source, local_pop.average_sink,
                                              local_pop.st_dev_population, local_pop.max_abs_population,
                                              local_pop.population_period_weak, local_pop.population_period_med,
                                              local_pop.population_period_strong,
                                              local_pop.recent_occupancy_change_frequency]
                f.write(f"{patch.number}, {local_pop.name}, {local_pop.average_population}, "
                        f"{local_pop.average_internal_change}, {local_pop.average_population_enter}, "
                        f"{local_pop.average_population_leave}, {local_pop.average_source}, {local_pop.average_sink}, "
                        f"{local_pop.st_dev_population}, {local_pop.max_abs_population}, "
                        f"{local_pop.population_period_weak}, {local_pop.population_period_med}, "
                        f"{local_pop.population_period_strong}, "
                        f"{local_pop.recent_occupancy_change_frequency};\n")
            average_population[patch.number] = this_patch
    # then also write the dictionary to a JSON
    json_file_name = f"{sim_path}/{step}/data/average_populations.json"
    dump_json(data=average_population, filename=json_file_name)


def write_population_history_data(patch_list, sim_path, step):
    print("Saving full local_population history (pop. size, internal change, dispersal) in individual .csv files.")
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            file_name = f"{sim_path}/{step}/data/local_pop_csv/patch_{patch.number}_{local_pop.name}.csv"
            with safe_open_w(file_name) as f:
                pop_history_array = np.asarray(local_pop.population_history)
                internal_change_array = np.asarray(local_pop.internal_change_history)
                pop_enter_array = np.asarray(local_pop.population_enter_history)
                pop_leave_array = np.asarray(local_pop.population_leave_history)
                combined_array = np.transpose(np.vstack([pop_history_array, pop_enter_array,
                                                         pop_leave_array, internal_change_array]))
                # noinspection PyTypeChecker
                np.savetxt(f, combined_array, newline='\n', fmt='%.20f')


def write_system_state(system_state, sim_path, step):
    print("Saving system_state object in JSON file.")
    json_file_name = f"{sim_path}/{step}/data/system_state.json"
    dump_json(data=system_state.__dict__, filename=json_file_name)


def write_patch_list_local_populations(patch_list, sim_path, step, is_save_local_populations):
    print("Saving patch objects in JSON files.")
    for patch in patch_list:
        json_file_name = f"{sim_path}/{step}/data/patch_data/patch_{patch.number}.json"
        dump_json(data=patch.__dict__, filename=json_file_name)
        if is_save_local_populations:
            for species_name, local_pop in patch.local_populations.items():
                json_local_file_name = f"{sim_path}/{step}/data/local_pop_json/patch" \
                                       f"_{patch.number}_{species_name}.json"
                dump_json(data=local_pop.__dict__, filename=json_local_file_name)


def distance_metrics_save(simulation_obj, sim_path, step):
    print("Saving distance metrics in JSON file.")
    json__file_name = f"{sim_path}/{step}/data/distance_metrics.json"
    dump_json(data=simulation_obj.system_state.distance_metrics_store, filename=json__file_name)


def pickle_save(simulation_obj, sim_path, step):
    pickle_file_name = f"{sim_path}/{step}/data/simulation_obj.pkl"
    save_object(simulation_obj, pickle_file_name)


def pickle_load(sim_path, step):
    pickle_file_name = f"{sim_path}/{step}/data/simulation_obj.pkl"
    file = open(pickle_file_name, "rb")
    simulation_obj = pickle.load(file)
    file.close()
    return simulation_obj


# ------------------------ UPDATING AND SAVING CURRENT VALUES OF LOCAL POPULATION ATTRIBUTES ------------------------ #

def update_local_population_nets(system_state):
    # current net_values and occupancy - no point calculating these all the time
    for patch in system_state.patch_list:
        for local_pop in patch.local_populations.values():
            local_pop.update_local_nets()


def update_stored_populations(patch_list, step):
    # store all current local pop values
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            local_pop.stored_pop_values[step] = local_pop.population
    # these may then be printed for any species at any time using
    # plot_current_local_population_attribute(patch_list=patch_list, sim=sim, attribute_name="stored_pop_values",
    #                                         step=step, species=species, sub_attr=STEP_TO_PRINT)


def save_current_local_population_attribute(patch_list, sim_path, attribute_name, step, sub_attr=None):
    attribute_dictionary = {}
    for patch in patch_list:
        patch_sub_dictionary = {}
        for local_pop in patch.local_populations.values():
            value = getattr(local_pop, attribute_name)
            if sub_attr is not None:
                value = value[sub_attr]
            patch_sub_dictionary[local_pop.name] = value
        attribute_dictionary[patch.number] = patch_sub_dictionary
    if sub_attr is not None:
        json_file_name = f"{sim_path}/{step}/data/{attribute_name}_{sub_attr}.json"
    else:
        json_file_name = f"{sim_path}/{step}/data/{attribute_name}.json"
    dump_json(data=attribute_dictionary, filename=json_file_name)


# ----------------------------------- RECORDING BASIC SPATIAL NETWORK PROPERTIES ----------------------------------- #

def save_network_properties(system_state, sim_path, step):
    # save .txt file with the essential spatial network properties - useful to use at the beginning of the program
    # alongside plot_network_properties() to describe the abiotic system
    txt_file_name = f"{sim_path}/{step}/data/network_properties.txt"
    with safe_open_w(txt_file_name) as f:
        f.write("Habitat amounts:\n")
        for habitat_num, habitat_history in system_state.habitat_amounts_history.items():
            f.write(f"{habitat_num}: {habitat_history[max(x for x in habitat_history)]}\n")
        f.write("\nPosterior habitat spatial auto-correlation:\n")
        f.write("Regular:\n")
        f.write(str(system_state.habitat_regular_auto_correlation_history[
                    max(x for x in system_state.habitat_regular_auto_correlation_history)])+'\n')
        f.write("Normalised:\n")
        f.write(str(system_state.habitat_spatial_auto_correlation_history[
                    max(x for x in system_state.habitat_spatial_auto_correlation_history)])+'\n')
        f.write("\nLCCs:\n")
        f.write("All:\n")
        f.write(str(system_state.patch_lcc_history['all'][max(x for x in system_state.patch_lcc_history['all'])])+'\n')
        f.write("Same:\n")
        f.write(str(system_state.patch_lcc_history['same'][
                        max(x for x in system_state.patch_lcc_history['same'])])+'\n')
        f.write("Different:\n")
        f.write(str(system_state.patch_lcc_history['different'][
                    max(x for x in system_state.patch_lcc_history['different'])])+'\n')
        f.write("\nDegree:\n")
        f.write(str(system_state.patch_degree_history[max(x for x in system_state.patch_degree_history)])+'\n')
        f.write("\nSize:\n")
        f.write(str(system_state.patch_size_history[max(x for x in system_state.patch_size_history)])+'\n')
        f.write("\nQuality:\n")
        f.write(str(system_state.patch_quality_history[max(x for x in system_state.patch_quality_history)])+'\n')
        f.write("\nCentrality:\n")
        f.write(str(system_state.patch_centrality_history[max(x for x in system_state.patch_centrality_history)])+'\n')
