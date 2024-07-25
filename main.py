#!/usr/bin/env python3

from simulation_obj import Simulation_obj
import random
from datetime import datetime
from data_manager import load_json
import numpy as np
from sample_spatial_data import run_sample_spatial_data
import importlib
import sys


# --------------------------------- MAIN PROGRAMS --------------------------------- #

def new_program(master_para, parameters_filename_base):
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
    simulation_obj = Simulation_obj(parameters=master_para, metadata=metadata,
                                    parameters_filename=parameters_filename_base + ".py")
    simulation_obj.full_simulation()


def repeat_program(parameters_filename_base, sim_number: int, sim_path: str):
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
    simulation_obj = Simulation_obj(parameters=loaded_parameters, metadata=loaded_metadata,
                                    parameters_filename=parameters_filename_base + ".py")
    simulation_obj.full_simulation()


#
# --------------------------------- EXECUTE --------------------------------- #
#

def execution():
    if len(sys.argv) > 1:
        # optionally pass in an argument specifying the particular parameters_???.py file to use
        parameters_filename_base = sys.argv[1]
    else:
        # otherwise use "parameters.py" as the default
        parameters_filename_base = "parameters"
    parameters_file = importlib.import_module(parameters_filename_base)
    importlib.reload(parameters_file)  # needed if the execution() function is being called again after the parameters
    # file has been changed.

    master_para = getattr(parameters_file, "master_para")
    meta_para = getattr(parameters_file, "meta_para")
    if meta_para["IS_RUN_SAMPLE_SPATIAL_DATA_FIRST"]:
        run_sample_spatial_data(parameters=master_para, is_output_files=True)
    num_repeats = meta_para["NUM_REPEATS"]

    for simulation in range(num_repeats):
        if meta_para["IS_NEW_PROGRAM"]:
            new_program(master_para=master_para, parameters_filename_base=parameters_filename_base)
        else:
            repeat_program(parameters_filename_base=parameters_filename_base,
                           sim_number=meta_para["REPEAT_PROGRAM_CODE"], sim_path=meta_para["REPEAT_PROGRAM_PATH"])


if __name__ == '__main__':
    execution()
