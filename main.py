from simulation_obj import Simulation_obj
from parameters import master_para, meta_para
import random
from datetime import datetime
from data_manager import load_json
import numpy as np
from sample_spatial_data import run_sample_spatial_data


# --------------------------------- MAIN PROGRAMS --------------------------------- #

def new_program():
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
    simulation_obj = Simulation_obj(parameters=master_para, metadata=metadata)
    simulation_obj.full_simulation()


def repeat_program(sim: int):
    print(f"\nRepeating simulation: {sim}.")
    loaded_parameters = load_json(f"results/{sim}/parameters.json")
    # NOTE: should this be replaced with loading the actual parameter.py and species para scripts? Would need to
    # be careful not to overwrite the existing scripts in the main directory from which the function is being called.
    # Probably best not too, if we can actually access everything needed as variables from the JSON load.
    #
    # make changes as desired to the program
    loaded_metadata = load_json(f"results/{sim}/metadata.json")
    np.random.seed(loaded_metadata["numpy_seed"])
    random.seed(loaded_metadata["random_seed"])
    loaded_metadata["Copy of simulation"] = sim
    simulation_obj = Simulation_obj(parameters=loaded_parameters, metadata=loaded_metadata)
    simulation_obj.full_simulation()


#
# --------------------------------- EXECUTE --------------------------------- #
#
if meta_para["IS_RUN_SAMPLE_SPATIAL_DATA_FIRST"]:
    run_sample_spatial_data(is_output_files=True)
num_repeats = meta_para["NUM_REPEATS"]
for simulation in range(num_repeats):
    if meta_para["IS_NEW_PROGRAM"]:
        new_program()
    else:
        repeat_program(sim=meta_para["REPEAT_PROGRAM_CODE"])
