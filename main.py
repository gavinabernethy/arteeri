#!/usr/bin/env python3

from simulation_obj import Simulation_obj
from data_save_functions import save_network_properties
from data_core_functions import create_adjacency_path_list
import random
from datetime import datetime
from data_manager import load_json
import numpy as np
from sample_spatial_data import run_sample_spatial_data
import importlib
import sys


# --------------------------------- MAIN PROGRAMS --------------------------------- #

def new_program(master_para, parameters_basename):
    print(f"\nBeginning a fresh simulation.")
    np_seed = np.random.randint(4294967296)
    random_seed = np.random.randint(4294967296)
    np.random.seed(np_seed)
    random.seed(random_seed)
    metadata = {
        "numpy_seed": np_seed,
        "random_seed": random_seed,
        "program_start_time": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    }
    call_program(parameters=master_para, metadata=metadata, parameters_basename=parameters_basename)


def repeat_program(parameters_basename, sim_number: int, sim_path: str):
    print(f"\nRepeating simulation: {sim_number}.")
    loaded_parameters = load_json(f"{sim_path}/parameters.json")
    # NOTE: should this be replaced with loading the actual parameter.py and species para scripts? Would need to
    # be careful not to overwrite the existing scripts in the main directory from which the function is being called.
    # Probably best not too, if we can actually access everything needed as variables from the JSON load.
    #
    # make changes as desired to the program
    loaded_metadata = load_json(f"{sim_path}/metadata.json")
    np.random.seed(loaded_metadata["numpy_seed"])
    random.seed(loaded_metadata["random_seed"])
    loaded_metadata["Copy of simulation"] = sim_number
    call_program(parameters=loaded_parameters, metadata=loaded_metadata, parameters_basename=parameters_basename)


def call_program(parameters, metadata, parameters_basename):
    simulation_obj = Simulation_obj(parameters=parameters, metadata=metadata,
                                    parameters_filename=parameters_basename + ".py")
    if parameters["plot_save_para"]["IS_ALLOW_FILE_CREATION"] and parameters["plot_save_para"]["PLOT_INIT_NETWORK"]:
        # print figures of the abiotic system (patches, habitats, centrality, reserves, etc.)
        adjacency_path_list = create_adjacency_path_list(
            patch_list=simulation_obj.system_state.patch_list,
            patch_adjacency_matrix=simulation_obj.system_state.patch_adjacency_matrix)
        from data_plot_functions import plot_network_properties
        plot_network_properties(patch_list=simulation_obj.system_state.patch_list,
                                sim_path=simulation_obj.sim_path,
                                step="initial_network",  # name of sub-folder
                                adjacency_path_list=adjacency_path_list,
                                is_biodiversity=False,
                                is_reserves=True,
                                is_partition=False,
                                is_label_habitat_patches=False,
                                is_retro=False)
        if parameters["plot_save_para"]["IS_SAVE"]:
            # save a .txt file with clustering and auto-correlation amounts
            save_network_properties(system_state=simulation_obj.system_state,
                                    sim_path=simulation_obj.sim_path,
                                    step="initial_network",  # name of sub-folder
                                    )
    if parameters["main_para"]["IS_SIMULATION"]:
        simulation_obj.full_simulation()


#
# --------------------------------- EXECUTE --------------------------------- #
#

def execution():
    if len(sys.argv) > 1:
        # optionally pass in an argument specifying the particular parameters_???.py file to use
        parameters_basename = sys.argv[1]
    else:
        # otherwise use "parameters.py" as the default
        parameters_basename = "parameters"
    parameters_file = importlib.import_module(parameters_basename)
    importlib.reload(parameters_file)  # needed if the execution() function is being called again after the parameters
    # file has been changed.

    master_para = getattr(parameters_file, "master_para")
    meta_para = getattr(parameters_file, "meta_para")
    if meta_para["IS_RUN_SAMPLE_SPATIAL_DATA_FIRST"]:
        run_sample_spatial_data(parameters=master_para,
                                is_output_files=master_para["plot_save_para"]["IS_ALLOW_FILE_CREATION"])
    num_repeats = meta_para["NUM_REPEATS"]

    for simulation in range(num_repeats):
        if meta_para["IS_NEW_PROGRAM"]:
            new_program(master_para=master_para, parameters_basename=parameters_basename)
        else:
            repeat_program(parameters_basename=parameters_basename,
                           sim_number=meta_para["REPEAT_PROGRAM_CODE"],
                           sim_path=meta_para["REPEAT_PROGRAM_PATH"])


if __name__ == '__main__':
    execution()
