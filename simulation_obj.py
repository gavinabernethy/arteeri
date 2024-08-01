from data_manager import save_all_data, generate_simulation_number, all_plots, population_snapshot, \
    change_snapshot, write_initial_files, save_adj_variables, load_adj_variables, load_reserve_list, \
    save_reserve_list, print_key_outputs_to_console
from data_manager_functions import plot_network_properties, create_adjacency_path_list
from sample_spatial_data import run_sample_spatial_data
import os
from patch import Patch
from local_population import Local_population
from species import Species
from datetime import datetime
from population_dynamics import *
from system_state import System_state
from perturbation import *
import json.decoder


def load_dataset(filename, force_dimension=None):
    # import an array from file
    if force_dimension is not None:
        data_array = np.loadtxt(filename, dtype='float', delimiter=',', ndmin=force_dimension)
    else:
        data_array = np.genfromtxt(filename, dtype='float', delimiter=',', autostrip=True)
    return data_array


######################################################################################################

class Simulation_obj:
    def __init__(self, parameters, metadata, parameters_filename):
        self.parameters = parameters
        self.metadata = metadata
        self.is_allow_file_creation = parameters["plot_save_para"]["IS_ALLOW_FILE_CREATION"]
        self.is_save = parameters["plot_save_para"]["IS_SAVE"]
        self.is_plot = parameters["plot_save_para"]["IS_PLOT"]
        self.is_print_key_outputs_to_console = parameters["plot_save_para"]["IS_PRINT_KEY_OUTPUTS_TO_CONSOLE"]
        self.total_steps = self.parameters["main_para"]["NUM_TRANSIENT_STEPS"] + self.parameters[
            "main_para"]["NUM_RECORD_STEPS"]
        self.sim_number, self.sim_path = generate_simulation_number(
            save_data=self.is_allow_file_creation,
            is_sub_folders=parameters["plot_save_para"]["IS_SUB_FOLDERS"],
            sub_folder_capacity=parameters["plot_save_para"]["SUB_FOLDER_CAPACITY"]
        )
        self.system_state = self.construction()
        self.parameters_filename = parameters_filename
        if self.is_allow_file_creation:
            write_initial_files(parameters=self.parameters, metadata=self.metadata, sim_path=self.sim_path,
                                parameters_filename=self.parameters_filename)
        print(f"Constructed simulation object for number {self.sim_number}.")

    def construction(self):
        # import all habitat and spatial network data from files
        test_set = self.parameters["graph_para"]["SPATIAL_TEST_SET"]
        dir_path = f'spatial_data_files/test_{test_set}/'
        # first, check that all .CSV files are present and strip any trailing commas!
        habitat_type_dictionary = self.parameters["main_para"]["HABITAT_TYPES"]
        try:
            for filename in os.listdir(dir_path):
                filepath = dir_path + filename
                if filepath.endswith('.csv'):
                    with open(filepath) as f:
                        content = f.read().rstrip()
                    if self.is_allow_file_creation:
                        if content[-1] == ',':
                            new_filename = filepath + '.tmp'
                            with open(new_filename, 'w') as f:
                                f.write(content[:-1])
                            os.rename(new_filename, filepath)
            habitat_species_traversal_array = load_dataset(
                f'{dir_path}/habitat_species_traversal.csv', force_dimension=2)
            habitat_species_feeding_array = load_dataset(f'{dir_path}/habitat_species_feeding.csv', force_dimension=2)
            patch_habitat_type_array = load_dataset(f'{dir_path}/patch_habitat_type.csv', force_dimension=1)
            patch_quality_array = load_dataset(f'{dir_path}/patch_quality.csv', force_dimension=1)
            patch_size_array = load_dataset(f'{dir_path}/patch_size.csv', force_dimension=1)
            patch_position_array = load_dataset(f'{dir_path}/patch_position.csv')
            patch_adjacency_array = load_dataset(f'{dir_path}/patch_adjacency.csv')

            # A test set has been successfully loaded - but we must check that it is suitable for this parameter setup.
            # Otherwise - halt. It will not work and the user probably forgot to check their spatial test setting.
            if patch_quality_array.shape[0] != self.parameters["main_para"][
                "NUM_PATCHES"] or habitat_species_traversal_array.shape[0] != \
                    len(habitat_type_dictionary) or habitat_species_traversal_array.shape[1] != \
                    len(self.parameters["main_para"]["SPECIES_TYPES"]):
                raise Exception(f"Existing spatial test set {test_set} not commensurate with specified number"
                                f" of patches, habitats and/or species. Change these parameters or allow creation of"
                                f" a replacement spatial test set.")
        except FileNotFoundError:
            # if (any) of the spatial network files do not exist, we MUST now create them here
            patch_position_array, patch_adjacency_array, patch_size_array, patch_quality_array, \
            patch_habitat_type_array, habitat_species_feeding_array, habitat_species_traversal_array = \
                run_sample_spatial_data(parameters=self.parameters,
                                        is_output_files=self.is_allow_file_creation)
            # flatten the Nx1 arrays to N vectors
            patch_quality_array = np.ndarray.flatten(patch_quality_array)
            patch_size_array = np.ndarray.flatten(patch_size_array)
            patch_habitat_type_array = np.ndarray.flatten(patch_habitat_type_array)

        # If there is only one patch and/or habitat, extend the arrays so they can be correctly sliced in the call
        if np.ndim(patch_position_array) == 1:
            patch_position_array = np.array([patch_position_array, [0.0, 0.0]])
        if np.ndim(patch_quality_array) == 0:
            patch_quality_array = np.asarray([patch_quality_array])
        if np.ndim(patch_size_array) == 0:
            patch_size_array = np.asarray([patch_size_array])
        if np.ndim(patch_habitat_type_array) == 0:
            patch_habitat_type_array = np.asarray([patch_habitat_type_array])
        if np.ndim(habitat_species_feeding_array) == 0:
            habitat_species_feeding_array = np.array([[habitat_species_feeding_array, 0.0], [0.0, 0.0]])
        elif np.ndim(habitat_species_feeding_array) == 1:
            habitat_species_feeding_array = np.array([habitat_species_feeding_array, [0.0, 0.0]])
        if np.ndim(habitat_species_traversal_array) == 0:
            habitat_species_traversal_array = np.array([[habitat_species_traversal_array, 0.0], [0.0, 0.0]])
        elif np.ndim(habitat_species_traversal_array) == 1:
            habitat_species_traversal_array = np.array([habitat_species_traversal_array, [0.0, 0.0]])
        if np.ndim(patch_adjacency_array) < 2:
            patch_adjacency_array = np.array([[patch_adjacency_array, 0.0], [0.0, 0.0]])
        # Record the maximum x and y values in the spatial network:
        dimensions = np.max(patch_position_array, axis=0)

        # Check the species inputs.
        initial_species_set = self.parameters["main_para"]["INITIAL_SPECIES_SET"]
        species_types = self.parameters["main_para"]["SPECIES_TYPES"]
        # Confirm keys in initialisation declaration are all found as species keys:
        for spec_num in initial_species_set:
            if spec_num not in species_types:
                raise Exception(f'Species type {spec_num} to be used in sim initiation but not part of the global set.')
        # Confirm that species keys are suitable integers in [0, num_spec - 1]:
        for spec_num in species_types:
            if spec_num > len(species_types) or type(spec_num) is not int:
                raise Exception(f'Species type {spec_num} is not a suitable number.')
            # Confirm that species names are unique:
            for spec_num_2 in species_types:
                if spec_num_2 != spec_num and species_types[spec_num_2] == species_types[spec_num]:
                    raise Exception(f"Species {species_types[spec_num_2]} assigned keys {spec_num} and {spec_num_2}.")
        # Ensure all integers in [0, num_spec - 1] are used as species keys:
        for check_num in range(len(species_types)):
            if check_num not in species_types:
                raise Exception(f'{check_num} is missing from the global set of species type keys.')
        # Confirm that the number of unique assigned keys matches the final size of the species repository:
        if len(species_types) != len(self.parameters["species_para"]):
            raise Exception("Number of expected species types does not match the parsed species parameter dictionary.")

        # generate ordered list of the species from names list in the parameters
        species_list = []
        species_dictionary = {}
        for species_num in self.parameters["main_para"]["INITIAL_SPECIES_SET"]:
            species_name = self.parameters["main_para"]["SPECIES_TYPES"][species_num]
            new_species = Species(
                name=species_name,
                species_num=species_num,
                lifespan=self.parameters["species_para"][species_name]["LIFESPAN"],
                minimum_population_size=self.parameters["species_para"][species_name]["MINIMUM_POPULATION_SIZE"],
                predator_list=self.parameters["species_para"][species_name]["PREDATOR_LIST"],
                resource_usage_conversion=self.parameters["species_para"][species_name]["GROWTH_PARA"][
                    "RESOURCE_USAGE_CONVERSION"],
                is_dispersal=self.parameters["species_para"][species_name]["DISPERSAL_PARA"]["IS_DISPERSAL"],
                dispersal_para=self.parameters["species_para"][species_name]["DISPERSAL_PARA"],
                is_dispersal_path_restricted=self.parameters["species_para"][species_name]["DISPERSAL_PARA"][
                    "IS_DISPERSAL_PATH_RESTRICTED"],
                always_move_with_minimum=self.parameters["species_para"][species_name]["DISPERSAL_PARA"][
                    "ALWAYS_MOVE_WITH_MINIMUM"],
                dispersal_penalty=max(self.parameters["species_para"][species_name]["DISPERSAL_PARA"][
                                          "SS_DISPERSAL_PENALTY"],
                                      self.parameters["pop_dyn_para"]["GENERAL_DISPERSAL_PENALTY"]),
                initial_population_mechanism=self.parameters["species_para"][species_name]["INITIAL_POPULATION_PARA"][
                    "INITIAL_POPULATION_MECHANISM"],
                initial_population_para=self.parameters["species_para"][species_name]["INITIAL_POPULATION_PARA"],
                growth_function=self.parameters["species_para"][species_name]["GROWTH_PARA"]["GROWTH_FUNCTION"],
                growth_para=self.parameters["species_para"][species_name]["GROWTH_PARA"],
                seasonal_period=self.parameters["species_para"][species_name]["SEASONAL_PERIOD"],
                is_growth_offset=self.parameters["species_para"][species_name]["GROWTH_PARA"]["ANNUAL_OFFSET"][
                    "IS_GROWTH_OFFSET"],
                growth_annual_duration=self.parameters["species_para"][species_name]["GROWTH_PARA"]["ANNUAL_OFFSET"][
                    "ANNUAL_DURATION"],
                growth_offset_species=self.parameters["species_para"][species_name]["GROWTH_PARA"]["ANNUAL_OFFSET"][
                    "GROWTH_OFFSET_SPECIES"],
                is_growth_offset_local=self.parameters["species_para"][species_name]["GROWTH_PARA"]["ANNUAL_OFFSET"][
                    "IS_GROWTH_OFFSET_LOCAL"],
                growth_offset_local=self.parameters["species_para"][species_name]["GROWTH_PARA"]["ANNUAL_OFFSET"][
                    "GROWTH_OFFSET_LOCAL"],
                predation_para=self.parameters["species_para"][species_name]["PREDATION_PARA"],
                is_predation_only_prevents_death=self.parameters["species_para"][species_name]["PREDATION_PARA"][
                    "IS_PREDATION_ONLY_PREVENTS_DEATH"],
                is_nonlocal_foraging=self.parameters["species_para"][species_name]["PREDATION_PARA"][
                    "IS_NONLOCAL_FORAGING"],
                is_foraging_path_restricted=self.parameters["species_para"][species_name]["PREDATION_PARA"][
                    "IS_NONLOCAL_FORAGING_PATH_RESTRICTED"],
                is_pure_direct_impact=self.parameters["species_para"][species_name]["IS_PURE_DIRECT_IMPACT"],
                pure_direct_impact_para=self.parameters["species_para"][species_name]["PURE_DIRECT_IMPACT_PARA"],
                is_direct_offset=self.parameters["species_para"][species_name]["PURE_DIRECT_IMPACT_PARA"][
                    "ANNUAL_OFFSET"]["IS_DIRECT_OFFSET"],
                direct_annual_duration=self.parameters["species_para"][species_name]["PURE_DIRECT_IMPACT_PARA"][
                    "ANNUAL_OFFSET"]["ANNUAL_DURATION"],
                direct_offset_species=self.parameters["species_para"][species_name]["PURE_DIRECT_IMPACT_PARA"][
                    "ANNUAL_OFFSET"]["DIRECT_OFFSET_SPECIES"],
                is_direct_offset_local=self.parameters["species_para"][species_name]["PURE_DIRECT_IMPACT_PARA"][
                    "ANNUAL_OFFSET"]["IS_DIRECT_OFFSET_LOCAL"],
                direct_offset_local=self.parameters["species_para"][species_name]["PURE_DIRECT_IMPACT_PARA"][
                    "ANNUAL_OFFSET"]["DIRECT_OFFSET_LOCAL"],
                direct_impact_on_me=self.parameters["species_para"][species_name]["DIRECT_IMPACT_ON_ME"],
                is_perturbs_environment=self.parameters["species_para"][species_name]["IS_PERTURBS_ENVIRONMENT"],
                perturbation_para=self.parameters["species_para"][species_name]["PERTURBATION_PARA"],
            )
            species_list.append(new_species)
            species_dictionary[species_name] = new_species
        species_set = {
            "dict": species_dictionary,
            "list": species_list,
        }

        # generate patches and place in a list
        patch_list = []
        for patch_num in range(self.parameters["main_para"]["NUM_PATCHES"]):
            new_patch = Patch(
                position=patch_position_array[patch_num, :],
                patch_number=patch_num,
                patch_quality=patch_quality_array[patch_num],
                patch_size=patch_size_array[patch_num],
                habitat_type_num=int(patch_habitat_type_array[patch_num]),
                habitat_type=habitat_type_dictionary[int(patch_habitat_type_array[patch_num])],
            )
            patch_list.append(new_patch)

        current_patch_list = [x for x in range(self.parameters["main_para"]["NUM_PATCHES"])]

        system_state = System_state(patch_list=patch_list,
                                    species_set=species_set,
                                    step=0,
                                    patch_adjacency_matrix=patch_adjacency_array,
                                    habitat_type_dictionary=habitat_type_dictionary,
                                    habitat_species_traversal=habitat_species_traversal_array,
                                    habitat_species_feeding=habitat_species_feeding_array,
                                    current_patch_list=current_patch_list,
                                    parameters=self.parameters,
                                    dimensions=dimensions,
                                    )

        return system_state

    ######################################################################################################

    def species_pathing(self):
        # generate shortest path cost for each species using Dijkstra's algorithm, and list of reachable patches
        is_generate_fresh = True
        if self.parameters["main_para"]["IS_LOAD_ADJ_VARIABLES"]:
            print("Attempting to load pre-existing adjacency variables.")
            try:
                load_adj_variables(patch_list=self.system_state.patch_list,
                                   spatial_set_number=self.parameters["graph_para"]["SPATIAL_TEST_SET"])
                is_generate_fresh = False
                print("Successfully loaded pre-existing adjacency variables.")
            except (FileNotFoundError, json.decoder.JSONDecodeError):
                is_generate_fresh = True
                print("Need to generate new adjacency variables as none to load or load unsuccessful.")

        if is_generate_fresh:
            # build fresh
            print("Generating new adjacency variables.")
            self.system_state.build_all_patches_species_paths_and_adjacency(parameters=self.parameters)
            print("Adjacency variables successfully generated.\n")
        if self.is_allow_file_creation and self.parameters["main_para"]["IS_SAVE_ADJ_VARIABLES"]:
            save_adj_variables(patch_list=self.system_state.patch_list,
                               spatial_set_number=self.parameters["graph_para"]["SPATIAL_TEST_SET"])
            print("Adjacency variables saved.\n")

    ######################################################################################################

    def full_simulation(self):
        print(f"Initialising simulation number {self.sim_number}.\n")
        self.species_pathing()
        self.simulation()
        self.metadata["simulation_end_time"] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if self.is_allow_file_creation:
            if self.is_save:
                print(f"\n{self.sim_number}: Beginning data saves.")
                save_all_data(simulation_obj=self)
                print(f"{self.sim_number}: Completed data saves.")
            if self.is_plot:
                print(f"{self.sim_number}: Beginning plot exports.")
                all_plots(simulation_obj=self)
                print(f"{self.sim_number}: Completed plot exports.")
        if self.is_print_key_outputs_to_console:
            print_key_outputs_to_console(simulation_obj=self)
        print(f"Completed simulation number {self.sim_number}.\n")

    ######################################################################################################

    def simulation(self):

        # construct reserves if necessary and flag them in patch properties
        is_generate_fresh = True
        reserve_clusters = [[]]
        test_set = self.parameters["graph_para"]["SPATIAL_TEST_SET"]
        if self.parameters["perturbation_para"]["IS_RESERVE_LOAD"]:
            try:
                reserve_clusters = load_reserve_list(spatial_set_number=test_set)
                is_generate_fresh = False
            except (FileNotFoundError, json.decoder.JSONDecodeError):
                print("Need to generate new reserve list as none to load or load unsuccessful.")
                is_generate_fresh = True
        if is_generate_fresh:
            if self.parameters["perturbation_para"]["IS_RESERVE_GENERATION"]:
                reserve_clusters = reserve_construction(system_state=self.system_state,
                                                        cluster_spec=self.parameters["perturbation_para"][
                                                            "RESERVE_CLUSTERS"],
                                                        clusters_must_be_separated=self.parameters["perturbation_para"][
                                                            "CLUSTERS_MUST_BE_SEPARATED"])
            else:
                reserve_clusters = self.parameters["perturbation_para"]["RESERVE_PATCH_CLUSTERS"]
        self.system_state.reserve_list = reserve_clusters
        for cluster_num, cluster in enumerate(reserve_clusters):
            for reserve_order, reserve_patch_num in enumerate(cluster):
                self.system_state.patch_list[reserve_patch_num].is_reserve = 1
                self.system_state.patch_list[reserve_patch_num].latest_perturbation_code = 0.0  # 0.1 otherwise
                self.system_state.patch_list[reserve_patch_num].reserve_order = [cluster_num, reserve_order]
        if self.is_allow_file_creation and self.parameters["perturbation_para"]["IS_RESERVE_SAVE"]:
            save_reserve_list(reserve_list=self.system_state.reserve_list, spatial_set_number=test_set)

        # initial populations
        for patch in self.system_state.patch_list:
            patch.local_populations = {}
            for species in self.system_state.species_set["list"]:
                # all patches should contain a local_population object for each species, which is persistent and
                # is NOT deleted if the local population should go extinct (this may be temporary and it is referenced)
                patch.local_populations[species.name] \
                    = Local_population(species=species,
                                       patch=patch,
                                       parameters=self.parameters,
                                       current_patch_list=self.system_state.current_patch_list,
                                       )
        is_nonlocal_foraging = self.parameters["pop_dyn_para"]["IS_NONLOCAL_FORAGING_PERMITTED"]
        is_local_foraging_ensured = self.parameters["pop_dyn_para"]["IS_LOCAL_FORAGING_ENSURED"]
        build_interacting_populations_list(
            patch_list=self.system_state.patch_list,
            species_list=self.system_state.species_set["list"],
            is_nonlocal_foraging=is_nonlocal_foraging,
            is_local_foraging_ensured=is_local_foraging_ensured,
            time=0,
        )
        is_dispersal = self.parameters["pop_dyn_para"]["IS_DISPERSAL_PERMITTED"]
        build_actual_dispersal_targets(
            patch_list=self.system_state.patch_list,
            species_list=self.system_state.species_set["list"],
            is_dispersal=is_dispersal,
            time=0,
        )

        # remove any known not to exist
        pass

        # now save a deepcopy of the patch list and the adjacency matrix - this is the initial abiotic system state!
        self.system_state.initial_patch_list = self.build_minimal_initial_patch_list()
        initial_patch_adjacency_matrix = deepcopy(self.system_state.patch_adjacency_matrix)
        self.system_state.initial_patch_adjacency_matrix = initial_patch_adjacency_matrix

        # do we require a visualisation of the system after initialisation BEFORE the first (0th) time-step has run?
        if self.is_allow_file_creation:
            if -1 in self.parameters["plot_save_para"]["MANUAL_SPATIAL_NETWORK_SAVE_STEPS"]:
                adjacency_path_list = create_adjacency_path_list(
                    patch_list=self.system_state.patch_list,
                    patch_adjacency_matrix=self.system_state.patch_adjacency_matrix)
                plot_network_properties(patch_list=self.system_state.patch_list, sim_path=self.sim_path, step=-1,
                                        adjacency_path_list=adjacency_path_list, is_biodiversity=True,
                                        is_reserves=True, is_retro=False)

        # MAIN LOOP - conduct simulation
        for step in range(self.total_steps):
            time = int(step / self.parameters["main_para"]["STEPS_TO_DAYS"])  # "int" truncates, round down as input>0
            self.system_state.time = time
            self.system_state.step = step

            # check for (system-level, not species-induced) perturbation time
            if step in self.parameters["perturbation_para"]["PERT_STEP_DICTIONARY"]:
                pert_archetype = self.parameters["perturbation_para"]["PERT_STEP_DICTIONARY"][step]
                pert_paras = self.parameters["perturbation_para"]["PERT_ARCHETYPE_DICTIONARY"][pert_archetype]
                if self.is_allow_file_creation and self.parameters["perturbation_para"]["IS_OUTPUT_DATAFILES"]:
                    save_all_data(simulation_obj=self)
                    population_snapshot(system_state=self.system_state, sim_path=self.sim_path, update_stored=True,
                                        output_figures=self.parameters["perturbation_para"]["IS_PLOTS"])
                # then enact the perturbation and reset the lists
                perturbation(system_state=self.system_state, parameters=self.parameters, pert_paras=pert_paras)

            # immediately after the perturbation and one iteration:
            if step - 1 in self.parameters["perturbation_para"]["PERT_STEP_DICTIONARY"]:
                population_snapshot(system_state=self.system_state, sim_path=self.sim_path, update_stored=True,
                                    output_figures=self.parameters["perturbation_para"]["IS_PLOTS"])
                if self.is_allow_file_creation and self.parameters["perturbation_para"]["IS_OUTPUT_DATAFILES"]:
                    save_all_data(simulation_obj=self)
                    # call a function to calculate the resulting change due to the perturbation
                    change_snapshot(system_state=self.system_state, sim_path=self.sim_path,
                                    output_figures=self.parameters["perturbation_para"]["IS_PLOTS"])

            # ---- Call the main update of all local populations ---- #
            update_populations(patch_list=self.system_state.patch_list,
                               species_list=self.system_state.species_set["list"],
                               parameters=self.parameters,
                               time=time,
                               step=step,
                               current_patch_list=self.system_state.current_patch_list,
                               is_ode_recordings=self.parameters["plot_save_para"]["IS_ODE_RECORDINGS"],
                               )

            # ---- Check for species-induced perturbations ---- #
            if self.parameters["pop_dyn_para"]["IS_SPECIES_PERTURBS_ENVIRONMENT"]:
                self.species_induced_perturbations()

            # ---- Check if the spatial environment has a natural tendency to change/restore ---- #
            if self.parameters["graph_para"]["IS_ENVIRONMENT_NATURAL_RESTORATION"]:
                self.environment_natural_restoration()

            # ---- Update the history of the number of available patches and biodiversity every time-step ---- #
            self.system_state.update_current_patch_history()
            self.system_state.update_biodiversity_history()

            # --------------------------------------------------------#

            # print the spatial network of the system during the simulation at the END OF STEP if manually specified
            if self.is_allow_file_creation and \
                    step in self.parameters["plot_save_para"]["MANUAL_SPATIAL_NETWORK_SAVE_STEPS"]:
                adjacency_path_list = create_adjacency_path_list(
                    patch_list=self.system_state.patch_list,
                    patch_adjacency_matrix=self.system_state.patch_adjacency_matrix)
                plot_network_properties(patch_list=self.system_state.patch_list, sim_path=self.sim_path, step=step,
                                        adjacency_path_list=adjacency_path_list, is_biodiversity=True,
                                        is_reserves=True, is_retro=False)

            # print current step
            if np.mod(step, 100) == 99:
                print(f'{self.sim_number}: Completed step: {step}/{self.total_steps - 1}')

        # -----------------------------------------------------------------------------------------------------------#
        # FINAL CALCULATIONS FOR THE LOCAL_POPULATION OBJECTS
        #
        # normalise average populations
        for patch in self.system_state.patch_list:
            for local_population in patch.local_populations.values():

                # build averages from recent histories - note that if there are M = M1 + M2 total steps, then the
                # population history indexes are from 0 to M-1, thus the "final step" (in terms of
                # history list indices) should be M-1:
                local_population.build_recent_time_averages(current_step=self.total_steps - 1,
                                                            back_steps=self.parameters["main_para"]["NUM_RECORD_STEPS"])

                if self.parameters["main_para"]["IS_CALCULATE_HURST"]:
                    import hurst
                    import warnings
                    warnings.filterwarnings(action='error',
                                            category=RuntimeWarning)
                    # this is used because of Hurst Exp throwing warnings, even for parts of the code not being
                    # executed because of the failure (constant time-series not suitable for Hurst).
                    # warnings.resetwarnings()  - you will need to run this if Warnings break the program!
                    # Calculate Hurst Exponent of each local population history time-series:
                    try:
                        local_population.population_history_hurst_exponent = hurst.compute_Hc(
                            series=local_population.population_history,
                            kind="random_walk", simplified=True)
                    except (RuntimeWarning, ValueError):
                        local_population.population_history_hurst_exponent = None

                # ??? Correlation dimension here ???

        # ??? Calculate Hurst Exponent of the global and average local diversity time-series:

        # -----------------------------------------------------------------------------------------------------------#

        # Analyse the final and averaged population distributions:
        self.system_state.update_distance_metrics(parameters=self.parameters)

        # return the resulting spatial distribution of species and average population sizes
        results = []
        for patch in self.system_state.patch_list:
            local_results = {
                "habitat_type": patch.habitat_type,
                "habitat_quality": patch.quality,
            }
            for local_pop in patch.local_populations.values():
                local_results[local_pop.name] = local_pop.average_population
            results.append(local_results)

    def environment_natural_restoration(self):
        # allows the environment to naturally restore towards a default habitat and quality
        is_quality_change = self.parameters["graph_para"]["RESTORATION_PARA"]["IS_QUALITY_CHANGE"]
        quality_change_probability = self.parameters["graph_para"]["RESTORATION_PARA"]["QUALITY_CHANGE_PROBABILITY"]
        quality_change_scale = self.parameters["graph_para"]["RESTORATION_PARA"]["QUALITY_CHANGE_SCALE"]
        quality_desired = self.parameters["graph_para"]["RESTORATION_PARA"]["QUALITY_DESIRED"]
        is_habitat_change = self.parameters["graph_para"]["RESTORATION_PARA"]["IS_HABITAT_CHANGE"]
        habitat_change_probability = self.parameters["graph_para"]["RESTORATION_PARA"]["HABITAT_CHANGE_PROBABILITY"]
        habitat_type_num_desired = self.parameters["graph_para"]["RESTORATION_PARA"]["HABITAT_TYPE_NUM_DESIRED"]
        # empty lists - we will build these with the patches and perturbations to implement all at once later
        patch_quality_change_list = []
        actual_quality_change_list = []
        patch_habitat_change_list = []
        actual_habitat_change_list = []
        for patch in self.system_state.patch_list:
            if is_quality_change:
                # if a quality change in principle, check if this patch differs and if so test if we will move
                if patch.quality != quality_desired and np.random.binomial(1, quality_change_probability):
                    patch_quality_change_list.append(patch.number)
                    # the 'relative add' to current quality = scale * difference from target
                    actual_quality_change_list.append(quality_change_scale * (quality_desired - patch.quality))
            if is_habitat_change:
                # if a habitat change in principle, check if this patch differs and if so test if we will move
                if patch.habitat_type_num != habitat_type_num_desired and np.random.binomial(
                        1, habitat_change_probability):
                    patch_habitat_change_list.append(patch.number)
                    actual_habitat_change_list.append(habitat_type_num_desired)
        # now enact both perturbations if necessary
        if len(patch_quality_change_list) > 0:
            perturbation(system_state=self.system_state, parameters=self.parameters,
                         pert_paras={
                             "perturbation_type": "patch_perturbation",
                             "perturbation_subtype": "change_quality",
                             "patch_list_overwrite": patch_quality_change_list,
                             "quality_change": actual_quality_change_list,
                             "is_absolute": False,
                             "is_relative_add": True,
                         },
                         )
        if len(patch_habitat_change_list) > 0:
            perturbation(system_state=self.system_state, parameters=self.parameters,
                         pert_paras={
                             "perturbation_type": "patch_perturbation",
                             "perturbation_subtype": "change_habitat",
                             "patch_list_overwrite": patch_habitat_change_list,
                             "habitat_nums_to_change_to": actual_habitat_change_list,
                         },
                         )

    def species_induced_perturbations(self):
        perturbation_has_occurred = False
        # reset
        self.system_state.perturbation_holding = {"remove": set({}),  # set of removed patches
                                                  "habitat_change": {},  # dictionary of patch: list habitat changes
                                                  "adjacency_change": {},  # dictionary of pair: list of adj changes
                                                  "quality_change": {},  # dict of patch: sum value quality changes
                                                  }
        for species in self.system_state.species_set["list"]:
            if species.is_perturbs_environment:
                perturbation_para = species.perturbation_para["PERTURBATION"]
                is_adjacency_pert = perturbation_para["IS_ADJACENCY_CHANGE"]
                # if there needs to be a possibility of making new connections, we draw against all (x,y)-neighbours.
                is_inc_adjacency_pert = (perturbation_para["ABSOLUTE_ADJACENCY_CHANGE"] > 0.0)
                is_patch_pert = perturbation_para["IS_REMOVAL"] or perturbation_para[
                    "IS_HABITAT_TYPE_CHANGE"] or perturbation_para["IS_QUALITY_CHANGE"]

                for patch in self.system_state.patch_list:
                    for local_pop in patch.local_populations.values():
                        if local_pop.species == species:
                            # enact for this local population
                            if "same" in species.perturbation_para["TO_IMPACT"]:
                                patch_num = patch.number
                                if is_patch_pert and draw_perturbation(local_pop=local_pop, patch_relation="SAME"):
                                    # draw once for this patch for within-patch properties
                                    perturbation_has_occurred = True
                                    self.record_perturbation(patch_num=patch_num, perturbation_para=perturbation_para)
                                if is_adjacency_pert:
                                    if is_inc_adjacency_pert:
                                        list_to_check = patch.set_of_xy_adjacent_patches
                                    else:
                                        list_to_check = patch.set_of_adjacent_patches
                                    # draw separately for every individual link between this patch and neighbours
                                    for other_patch_num in list_to_check:
                                        if other_patch_num != patch_num:
                                            if draw_perturbation(local_pop=local_pop, patch_relation="SAME"):
                                                perturbation_has_occurred = True
                                                self.record_adj_perturbation(patch_nums=[patch_num, other_patch_num],
                                                                             perturbation_para=perturbation_para)
                            # enact for any directly connected patches to this one
                            if "adjacent" in species.perturbation_para["TO_IMPACT"]:
                                for patch_num in patch.set_of_adjacent_patches:
                                    if is_patch_pert and draw_perturbation(local_pop=local_pop,
                                                                           patch_relation="ADJACENT"):
                                        perturbation_has_occurred = True
                                        self.record_perturbation(patch_num=patch_num,
                                                                 perturbation_para=perturbation_para)
                                    if is_adjacency_pert:
                                        if is_inc_adjacency_pert:
                                            list_to_check = self.system_state.patch_list[
                                                patch_num].set_of_xy_adjacent_patches
                                        else:
                                            list_to_check = self.system_state.patch_list[
                                                patch_num].set_of_adjacent_patches
                                        for other_patch_num in list_to_check:
                                            if other_patch_num != patch_num:
                                                # now look at patches connected to this neighbour (inc the original)
                                                if draw_perturbation(local_pop=local_pop, patch_relation="ADJACENT"):
                                                    perturbation_has_occurred = True
                                                    self.record_adj_perturbation(
                                                        patch_nums=[patch_num, other_patch_num],
                                                        perturbation_para=perturbation_para)
                            # enact for any patches with distance exactly 1.0 to this patch
                            if "xy-adjacent" in species.perturbation_para["TO_IMPACT"]:
                                for patch_num in patch.set_of_xy_adjacent_patches:
                                    if is_patch_pert and draw_perturbation(local_pop=local_pop,
                                                                           patch_relation="XY_ADJACENT"):
                                        perturbation_has_occurred = True
                                        self.record_perturbation(patch_num=patch_num,
                                                                 perturbation_para=perturbation_para)
                                    if is_adjacency_pert:
                                        if is_inc_adjacency_pert:
                                            list_to_check = self.system_state.patch_list[
                                                patch_num].set_of_xy_adjacent_patches
                                        else:
                                            list_to_check = self.system_state.patch_list[
                                                patch_num].set_of_adjacent_patches
                                        for other_patch_num in list_to_check:
                                            if other_patch_num != patch_num:
                                                # now look at patches directly connected to this xy-neighbour
                                                if draw_perturbation(local_pop=local_pop, patch_relation="XY_ADJACENT"):
                                                    perturbation_has_occurred = True
                                                    self.record_adj_perturbation(
                                                        patch_nums=[patch_num, other_patch_num],
                                                        perturbation_para=perturbation_para)
        if perturbation_has_occurred:
            # resolve final outcomes and implement them with the (up to four) patch perturbations
            if len(self.system_state.perturbation_holding["remove"]) != 0:
                # removal first as a single perturbation
                perturbation(system_state=self.system_state, parameters=self.parameters,
                             pert_paras={
                                 "perturbation_type": "patch_perturbation",
                                 "perturbation_subtype": "patch_removal",
                                 "patch_list_overwrite": self.system_state.perturbation_holding["remove"],
                             },
                             )
            # for habitat change
            if len(self.system_state.perturbation_holding["habitat_change"]) != 0:
                final_patch_num_list = []  # build this simultaneously to ensure ordering preserved
                final_habitat_change_list = []
                for patch_num in self.system_state.perturbation_holding["habitat_change"]:
                    final_patch_num_list.append(patch_num)
                    final_habitat_change_list.append(random.choice(
                        self.system_state.perturbation_holding["habitat_change"][patch_num]))
                perturbation(system_state=self.system_state, parameters=self.parameters,
                             pert_paras={
                                 "perturbation_type": "patch_perturbation",
                                 "perturbation_subtype": "change_habitat",
                                 "patch_list_overwrite": final_patch_num_list,
                                 "habitat_nums_to_change_to": final_habitat_change_list,
                             },
                             )
            # for quality change
            if len(self.system_state.perturbation_holding["quality_change"]) != 0:
                final_patch_num_list = []  # build this simultaneously to ensure ordering preserved
                final_quality_change_list = []
                for patch_num in self.system_state.perturbation_holding["quality_change"]:
                    final_patch_num_list.append(patch_num)
                    final_quality_change_list.append(self.system_state.perturbation_holding[
                                                         "quality_change"][patch_num])
                perturbation(system_state=self.system_state, parameters=self.parameters,
                             pert_paras={
                                 "perturbation_type": "patch_perturbation",
                                 "perturbation_subtype": "change_quality",
                                 "patch_list_overwrite": final_patch_num_list,
                                 "quality_change": final_quality_change_list,
                                 "is_absolute": False,
                                 "is_relative_add": True,
                             },
                             )
            # for adjacency change
            if len(self.system_state.perturbation_holding["adjacency_change"]) != 0:
                final_pair_list = []
                final_adjacency_change_list = []
                for pair in self.system_state.perturbation_holding["adjacency_change"]:
                    mean_change = np.mean(self.system_state.perturbation_holding["adjacency_change"][pair])
                    if mean_change >= 0.5:
                        actual_change = 1
                    else:
                        actual_change = 0
                    final_pair_list.append(pair)
                    final_adjacency_change_list.append(actual_change)
                perturbation(system_state=self.system_state, parameters=self.parameters,
                             pert_paras={
                                 "perturbation_type": "patch_perturbation",
                                 "perturbation_subtype": "change_adjacency",
                                 "patch_list_overwrite": final_pair_list,
                                 "is_pairs": True,
                                 "adjacency_change": final_adjacency_change_list,
                             },
                             )

    def record_perturbation(self, patch_num, perturbation_para):
        # for a given patch and the species-specific perturbation parameters, update the running lists of how this patch
        # is being impacted in this time-step (excluding adjacency changes, which require a dedicated function)
        if perturbation_para["IS_REMOVAL"]:
            # removal overwrites all other effects, simply record a set of all patches to remove
            self.system_state.perturbation_holding["remove"].add(patch_num)
        else:
            if perturbation_para["IS_HABITAT_TYPE_CHANGE"]:
                # record a list of all habitat types to change this patch to (will choose one at random)
                if patch_num in self.system_state.perturbation_holding["habitat_change"]:
                    self.system_state.perturbation_holding["habitat_change"][patch_num].append(
                        perturbation_para["HABITAT_TYPE_NUM_TO_CHANGE_TO"])
                else:
                    self.system_state.perturbation_holding["habitat_change"][patch_num] = [
                        perturbation_para["HABITAT_TYPE_NUM_TO_CHANGE_TO"]]
            if perturbation_para["IS_QUALITY_CHANGE"]:
                # sum all the relative quality changes for this patch (will take the mean)
                if patch_num in self.system_state.perturbation_holding["quality_change"]:
                    self.system_state.perturbation_holding["quality_change"][patch_num] += \
                        perturbation_para["RELATIVE_QUALITY_CHANGE"]
                else:
                    self.system_state.perturbation_holding["quality_change"][patch_num] = \
                        perturbation_para["RELATIVE_QUALITY_CHANGE"]

    def record_adj_perturbation(self, patch_nums, perturbation_para):
        # update list of {0, 1} adj changes for this patch pair (will check if mean >= 0.5)
        pair = (min(patch_nums), max(patch_nums))  # consistent order
        if pair in self.system_state.perturbation_holding["adjacency_change"]:
            # the (,) pair can indeed be used as a dictionary key
            self.system_state.perturbation_holding["adjacency_change"][pair].append(
                perturbation_para["ABSOLUTE_ADJACENCY_CHANGE"])
        else:
            self.system_state.perturbation_holding["adjacency_change"][pair] = [
                perturbation_para["ABSOLUTE_ADJACENCY_CHANGE"]]

    def build_minimal_initial_patch_list(self):
        # this creates a simplified copy of the initial patch list with a minimal record of the essential physical
        # properties (creating a deepcopy of the entire patch list after local population objects added leads to an
        # error with too deep a recursion tree).
        # We thus here record all properties that we may wish to visualise during the retrospective_network_plots()
        initial_patch_list = []
        for patch in self.system_state.patch_list:
            new_patch = Patch(
                position=patch.position,
                patch_number=patch.number,
                patch_quality=patch.quality,
                patch_size=patch.size,
                habitat_type_num=patch.habitat_type_num,
            )
            new_patch.degree = deepcopy(patch.degree)
            new_patch.centrality = deepcopy(patch.centrality)
            new_patch.num_times_perturbed = deepcopy(patch.num_times_perturbed)
            new_patch.is_reserve = deepcopy(patch.is_reserve)
            new_patch.reserve_order = deepcopy(patch.reserve_order)
            new_patch.latest_perturbation_code = deepcopy(patch.latest_perturbation_code)
            initial_patch_list.append(new_patch)
        return initial_patch_list


def draw_perturbation(local_pop, patch_relation):
    # for a given local population and affected patch, determine if a species-induced perturbation probabilistically
    # occurs based on the population density and the species-specific coefficient parameters
    coefficients = local_pop.species.perturbation_para[
        "IMPLEMENTATION_PROBABILITY_COEFFICIENTS"][patch_relation]
    density = local_pop.population / local_pop.carrying_capacity
    # build the probability function
    probability = coefficients[0] * np.heaviside(density, 0.0) + coefficients[
        1] * density + coefficients[2] * density ** 2.0 + coefficients[3] * density ** 3.0
    # draw value and return it
    return np.random.binomial(1, probability)
