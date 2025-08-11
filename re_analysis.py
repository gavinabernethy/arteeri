#!/usr/bin/env python3

# this script loads all overview data from a specified complete simulation,
# then conducts whatever additional analysis is required
#
# In particular, you can also request any specific data point from a local population recording in the JSON files,
# and (very helpful) you can append data from the "ode_recording" nested dictionary, to obtain the time-series of
# only a specific property (e.g. local_growth or r_mod) to easily inspect them, rather than having to spend time
# formatting the large JSONs to be able to actually read these data streams.
import numpy as np
from source_code.data_core_functions import load_json
from source_code.data_save_functions import pickle_load
from source_code.data_plot_functions import retrospective_network_plots
import os.path
import functools

# ---------------------- REQUEST ---------------------- #
# Set the simulation number and the time-step of the data you want to reload:
SIM_NUMBER = 111
TIME = 1999

# Load specific nested values from a JSON file. To combine values from a JSON (e.g. produce the time-series of r_mod,
# or local_growth for a local population) specify the local population JSON file that you want (i.e. patch and species),
# then two paths in a list - specifying where the root of the iterable (list of dictionary) is, and then the local path
# of what you want to retrieve for the nest of *each item* in this iterable object.
SPECIFIC_DATA = [
    {"patch_number": 1,
     "species_name": "predator",
     "property_path": [("ode_recording",), ("new_population",)], # e.g. [("name", )], [("ode_recording",), ("local_growth",)]
     }
]


# ----------------------------------------------------- #

def load_overview_data(sim, time):
    parameters = load_json(f"results/{sim}/{time}/parameters.json")
    metadata = load_json(f"results/{sim}/{time}/metadata.json")
    average_population = load_json(f"results/{sim}/{time}/data/average_populations.json")
    species_dictionary = parameters["main_para"]["SPECIES_TYPES"]
    num_patches = parameters["main_para"]["NUM_PATCHES"]
    population_history_dictionary = load_population_history_data(
        sim=sim,
        time=time,
        num_patches=num_patches,
        species_dictionary=species_dictionary
    )
    overview_data_int = {
        "Parameters": parameters,
        "Metadata": metadata,
        "average_population": average_population,
        "population_history_dictionary": population_history_dictionary,
    }
    return overview_data_int

def load_data_stream(sim, time, patch_number, species_name, property_path):
    file_name = f"results/{sim}/{time}/data/local_pop_json/patch_{patch_number}_{species_name}.json"
    if os.path.exists(file_name):
        json_file = load_json(file_name)
        root_path = property_path[0]
        root_property = functools.reduce(dict.get, root_path, json_file)
        if len(property_path) == 1:
            return root_property
        elif len(property_path) == 2:
            relative_path = property_path[1]
            combined_data_stream = []
            if type(root_property) is list:
                iterable_object = root_property
            elif type(root_property) is dict:
                iterable_object = root_property.values()
            else:
                raise 'Error - not a list or dictionary to iterate through.'

            for iter_item in iterable_object:
                if type(iter_item) is dict:
                    combined_data_stream.append(functools.reduce(dict.get, relative_path, iter_item))
                else:
                    raise 'Error - not a dictionary to continue searching.'
            return combined_data_stream
        else:
            raise Exception('Data stream not suitable.')
    else:
        print (f"File patch_{patch_number}_{species_name}.json not found. Did you enable patch JSON storage "
                        f"for this simulation?")

def load_population_history_data(sim, time, num_patches, species_dictionary):
    population_history_dictionary = {}
    for patch_number in range(num_patches):
        patch_dictionary = {}
        for species_name in species_dictionary.values():
            file_name = f"results/{sim}/{time}/data/patch_{patch_number}_{species_name}.csv"
            if os.path.exists(file_name):
                patch_dictionary[species_name] = np.genfromtxt(file_name, dtype='float', delimiter=',', autostrip=True)
        population_history_dictionary[patch_number] = patch_dictionary
    return population_history_dictionary


# ---------------------- EXECUTE ---------------------- #
overview_data = load_overview_data(sim=SIM_NUMBER, time=TIME)
if len(SPECIFIC_DATA) > 0:
    combined_data = []
    for element in SPECIFIC_DATA:
        combined_data.append(
            load_data_stream(
                sim=SIM_NUMBER,
                time=TIME,
                patch_number=element["patch_number"],
                species_name=element["species_name"],
                property_path=element["property_path"],
            )
        )
simulation_obj = pickle_load(sim_path=f'results/{SIM_NUMBER}', step=TIME)
# want to produce some spatial network plots at an arbitrary time-step?
if simulation_obj is not None:  # i.e. if there WAS a pickle file to find...
    retrospective_network_plots(
        initial_patch_list=simulation_obj.system_state.initial_patch_list,
        actual_patch_list=simulation_obj.system_state.patch_list,
        initial_patch_adjacency_matrix=simulation_obj.system_state.initial_patch_adjacency_matrix,
        sim_path=f'results/{SIM_NUMBER}', step=TIME)

# break here to show data in the console when run in debugger
x=1
