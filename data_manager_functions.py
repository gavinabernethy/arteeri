import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as pe
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import os.path
import numpy as np
import random
from datetime import datetime
import json
import pickle
import sys
from copy import deepcopy


# ----------------------------- AUXILIARY FUNCTIONS FOR FILE SAVING AND OBJECT HANDLING ----------------------------- #

def safe_open_w(path):
    # Open "path" for writing, creating any parent directories as needed.
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, 'w')


def load_json(input_file):
    with open(input_file, mode='r') as f:
        output_array = json.load(f)
        if type(output_array) == dict:
            output_array = convert_keys_to_int(output_array)
    return output_array


def convert_keys_to_int(d: dict):
    # This is needed because I sometimes use integers as keys for the species and habitat master sets
    # See:
    # https://stackoverflow.com/questions/1450957/pythons-json-module-converts-int-dictionary-keys-to-strings
    new_dict = {}
    for k, v in d.items():
        try:
            new_key = int(k)
        except ValueError:
            new_key = k
        if type(v) == dict:
            v = convert_keys_to_int(v)
        new_dict[new_key] = v
    return new_dict


def dump_json(data, filename):
    with safe_open_w(filename) as f:
        json.dump(data, f, ensure_ascii=False, default=set_default, skipkeys=True)


def save_object(obj, filename):
    sys.setrecursionlimit(3000)
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        try:
            pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
        except RecursionError:
            print(f"Object {obj} was too large to pickle.")


def set_default(obj):
    # use this option to convert any nested sets (e.g. in the parameters) to a list for the JSON serialising
    if isinstance(obj, set):
        return list(obj)
    # convert numpy arrays to nest lists
    if isinstance(obj, np.ndarray):
        return obj.tolist()


# -------------------------------- FUNCTIONS FOR SAVING AND LOADING THE ENTIRE SYSTEM -------------------------------- #

def write_parameters_file(parameters, sim, step):
    parameters_file = f"results/{sim}/{step}/parameters.json"
    dump_json(data=parameters, filename=parameters_file)


def write_metadata_file(metadata, sim, step):
    metadata_file = f"results/{sim}/{step}/metadata.json"
    metadata["write_time"] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dump_json(data=metadata, filename=metadata_file)


def write_perturbation_history_data(perturbation_history, sim, step):
    pert_history_file = f"results/{sim}/{step}/perturbation_history.json"
    dump_json(data=perturbation_history, filename=pert_history_file)


def write_current_species_movement_scores(patch_list, sim, step):
    for patch in patch_list:
        species_movement_scores_file = f"results/{sim}/{step}/data/species_movement_scores/sms_{patch.number}.json"
        dump_json(data=patch.species_movement_scores, filename=species_movement_scores_file)
    # recall that in previous versions these were stored as an array, and thus had to be converted to a
    # nested dictionary first


def write_average_population_data(patch_list, sim, step):
    print("Saving text and JSON file of all patches averaged pop.s, internal change, dispersal for each local_pop.")
    # build a dictionary
    average_population = {}
    # Text file is more readable but less easy to analyse:
    txt_file_name = f"results/{sim}/{step}/data/average_populations.txt"
    with safe_open_w(txt_file_name) as f:
        for patch in patch_list:
            this_patch = {}
            for local_pop in patch.local_populations.values():
                this_patch[local_pop.name] = [local_pop.average_population, local_pop.average_internal_change,
                                              local_pop.average_population_enter, local_pop.average_population_leave,
                                              local_pop.average_source, local_pop.average_sink,
                                              local_pop.st_dev_population, local_pop.max_abs_population,
                                              local_pop.population_period]
                f.write(f"{patch.number}, {local_pop.name}, {local_pop.average_population}, "
                        f"{local_pop.average_internal_change}, {local_pop.average_population_enter}, "
                        f"{local_pop.average_population_leave}, {local_pop.average_source}, {local_pop.average_sink}, "
                        f"{local_pop.st_dev_population}, {local_pop.max_abs_population}, "
                        f"{local_pop.population_period}, {local_pop.recent_occupancy_change_frequency};\n")
            average_population[patch.number] = this_patch
    # then also write the dictionary to a JSON
    json_file_name = f"results/{sim}/{step}/data/average_populations.json"
    dump_json(data=average_population, filename=json_file_name)


def write_population_history_data(patch_list, sim, step):
    print("Saving full local_population history (pop. size, internal change, dispersal) in individual .csv files.")
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            file_name = f"results/{sim}/{step}/data/local_pop_csv/patch_{patch.number}_{local_pop.name}.csv"
            with safe_open_w(file_name) as f:
                pop_history_array = np.asarray(local_pop.population_history)
                internal_change_array = np.asarray(local_pop.internal_change_history)
                pop_enter_array = np.asarray(local_pop.population_enter_history)
                pop_leave_array = np.asarray(local_pop.population_leave_history)
                combined_array = np.transpose(np.vstack([pop_history_array, pop_enter_array,
                                                         pop_leave_array, internal_change_array]))
                # noinspection PyTypeChecker
                np.savetxt(f, combined_array, newline='\n', fmt='%f')


def write_patch_list_local_populations(patch_list, sim, step, is_save_local_populations):
    print("Saving patch objects in JSON files.")
    for patch in patch_list:
        json_file_name = f"results/{sim}/{step}/data/patch_data/patch_{patch.number}.json"
        dump_json(data=patch.__dict__, filename=json_file_name)
        if is_save_local_populations:
            for species_name, local_pop in patch.local_populations.items():
                json_local_file_name = f"results/{sim}/{step}/data/local_pop_json/patch" \
                                       f"_{patch.number}_{species_name}.json"
                dump_json(data=local_pop.__dict__, filename=json_local_file_name)


def pickle_save(simulation_obj, sim, step):
    pickle_file_name = f"results/{sim}/{step}/data/simulation_obj.pkl"
    save_object(simulation_obj, pickle_file_name)


def pickle_load(sim, step):
    pickle_file_name = f"results/{sim}/{step}/data/simulation_obj.pkl"
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


def save_current_local_population_attribute(patch_list, sim, attribute_name, step, sub_attr=None):
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
        json_file_name = f"results/{sim}/{step}/data/{attribute_name}_{sub_attr}.json"
    else:
        json_file_name = f"results/{sim}/{step}/data/{attribute_name}.json"
    dump_json(data=attribute_dictionary, filename=json_file_name)


# ---------------------------------------- BASE FUNCTIONS FOR CREATING PLOTS ---------------------------------------- #

def create_time_series_plot(data, parameters, file_path, y_label='Population', legend_list=None,
                            start_time=1, end_time=None, is_shading=False):
    if end_time is None:
        end_time = parameters["main_para"]["NUM_TRANSIENT_STEPS"] + parameters["main_para"]["NUM_RECORD_STEPS"]
    x_data = np.linspace(start_time, end_time, num=1 + end_time - start_time)
    fig = plt.figure()
    for series_number, data_time_series in enumerate(data):
        plt.plot(x_data, data_time_series)
        if is_shading and series_number > 0:
            # shade between this series and the previous one
            plt.fill_between(x_data, data[series_number - 1], data_time_series, alpha=0.2)
    if legend_list is not None:
        if len(legend_list) < 20:
            # avoid trying to draw a legend if too many objects
            plt.legend(legend_list)
    plt.xlabel('Time step')
    plt.ylabel(y_label)
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(file_path, dpi=400)
    plt.close(fig)


def create_patches_plot(patch_list, color_property, file_path, path_list=None, patch_label_attr=None,
                        label_paths=False, path_color=None, label_patches=False, patch_label_color=None,
                        use_colors=False, use_color_bar=False):
    # create the basemap image of the patches and their relative positions
    if path_list is None:
        path_list = []
    fig, ax = plt.subplots()
    min_val = np.zeros([2])
    max_val = np.zeros([2])
    min_color_value = np.min(color_property)
    max_color_value = np.max(color_property)
    # Some set up for the colorings
    color_offset = 0.0
    is_one_color = False
    if np.max(color_property) - np.min(color_property) == 0.0:
        is_one_color = True
    elif use_colors:
        color_offset = 0.15  # for cyclic color maps (e.g. hsv) [0] and [1] are the same, so need an offset
    # Now the primary iteration of each patch to determine its position, value and color
    for patch_num, patch in enumerate(patch_list):
        length = np.sqrt(patch.size)
        # space them out slightly so that size=1 patches can still be distinguished
        position = (patch.position[0]*1.1, patch.position[1]*1.1)
        min_val = np.asarray([min(min_val[0], position[0] - length), min(min_val[1], position[1] - length)])
        max_val = np.asarray([max(max_val[0], position[0] + 2 * length), max(max_val[1], position[1] + 2 * length)])

        if is_one_color:
            # only one entry
            if np.max(color_property) == 0.0:
                max_color_value = 1.0
                c = [0.0, 0.0, 0.0]
            else:
                min_color_value = 0.0
                c = [1.0, 1.0, 1.0]
        else:
            if use_colors:
                color_map = cm.get_cmap('hsv')
                c = color_map((1.0 - color_offset) * float((color_property[patch_num] - np.min(color_property)) /
                                                           (np.max(color_property) - np.min(color_property)))
                              + color_offset)
            else:
                c = [(1.0 - color_offset) * float((color_property[patch_num] - np.min(color_property)) /
                                                  (np.max(color_property) - np.min(color_property))) + color_offset
                     for _ in range(3)]

        # check if coloring is max or min for possible color bar
        # if float(color_property[patch_num]) in [min_color_value, max_color_value]:
        #     if len(custom_map) == 0:
        #         custom_map.append((float(color_property[patch_num]), c))
        #     elif float(color_property[patch_num]) not in [x[0] for x in custom_map]:
        #         custom_map.append((float(color_property[patch_num]), c))

        rectangle = patches.Rectangle(position, length, length, linewidth=1, edgecolor='black', facecolor=c)
        ax.add_patch(rectangle)
        if label_patches:
            # annotate the centre of the rectangle with value
            rx, ry = rectangle.get_xy()
            if patch_label_attr is not None:
                data_val = getattr(patch_list[patch_num], patch_label_attr)
                if type(data_val) == list and len(data_val) == 0:
                    # do not clutter figure by printing empty values
                    data_val = None
            else:
                # label patches with the colour matrix unless an attribute specified
                data_val = color_property[patch_num, 0]
            if data_val is not None:
                # format the patch labels
                if type(data_val) in [int, float]:
                    # if a single number, set format
                    data_str = "{0:.3g}".format(data_val)
                else:
                    # cast to string and remove all whitespace
                    data_str = str(data_val).replace(" ", "")
                if patch_label_color is None:
                    color = 'red'
                else:
                    color = patch_label_color
                ax.annotate(data_str, (rx + length / 2, ry + length / 2),
                            color=color, fontsize=5, ha='center', va='center')

    # draw paths (showing traversable links) if required:
    if path_list is not None and len(path_list) > 0:

        min_x = min([element.position[0] for element in patch_list])
        max_x = max([element.position[0] for element in patch_list])
        min_y = min([element.position[1] for element in patch_list])
        max_y = max([element.position[1] for element in patch_list])

        for path_tuple in path_list:

            # check for wrap
            is_wrap = False
            if min_x != max_x and min_y != max_y:
                if patch_list[path_tuple[0]].position[0] == patch_list[path_tuple[1]].position[0]:
                    if {min_y, max_y} <= {patch_list[path_tuple[0]].position[1], patch_list[path_tuple[1]].position[1]}:
                        is_wrap = True
                if patch_list[path_tuple[0]].position[1] == patch_list[path_tuple[1]].position[1]:
                    if {min_x, max_x} <= {patch_list[path_tuple[0]].position[0], patch_list[path_tuple[1]].position[0]}:
                        is_wrap = True

            if is_wrap:
                patch_one = [patch_list[path_tuple[0]].position[0] + 0.5 * np.sqrt(patch_list[path_tuple[0]].size),
                             patch_list[path_tuple[0]].position[1] + 0.5 * np.sqrt(patch_list[path_tuple[1]].size)]
                patch_two = [patch_list[path_tuple[1]].position[0] + 0.5 * np.sqrt(patch_list[path_tuple[0]].size),
                             patch_list[path_tuple[1]].position[1] + 0.5 * np.sqrt(patch_list[path_tuple[1]].size)]
                both_patches = [patch_one, patch_two]
                for line_num in range(2):
                    x_start = both_patches[line_num][0]
                    x_end = both_patches[line_num][0] - (both_patches[line_num-1][0] - both_patches[line_num][0]) / (
                            np.sqrt(len(patch_list)))
                    y_start = both_patches[line_num][1]
                    y_end = both_patches[line_num][1] - (both_patches[line_num - 1][1] - both_patches[line_num][1]) / (
                            np.sqrt(len(patch_list)))
                    if type(path_tuple[3]) is list and len(path_tuple[3]) == 3:
                        ax.plot([x_start, x_end], [y_start, y_end], color=path_tuple[3])
                    elif path_color is None:
                        ax.plot([x_start, x_end], [y_start, y_end])
                    else:
                        ax.plot([x_start, x_end], [y_start, y_end], color=path_color)

            else:
                # not a wrp - just a direct path
                start_x = patch_list[path_tuple[0]].position[0] + 0.5 * np.sqrt(patch_list[path_tuple[0]].size)
                end_x = patch_list[path_tuple[1]].position[0] + 0.5 * np.sqrt(patch_list[path_tuple[1]].size)
                start_y = patch_list[path_tuple[0]].position[1] + 0.5 * np.sqrt(patch_list[path_tuple[0]].size)
                end_y = patch_list[path_tuple[1]].position[1] + 0.5 * np.sqrt(patch_list[path_tuple[1]].size)
                mid_x = (start_x + end_x) / 2
                mid_y = (start_y + end_y) / 2

                if type(path_tuple[3]) is list and len(path_tuple[3]) == 3:
                    # is color specified for this path?
                    line = ax.plot([start_x, end_x], [start_y, end_y], color=path_tuple[3])
                else:
                    # otherwise the general coloring schemes
                    if path_color is None:
                        line = ax.plot([start_x, end_x], [start_y, end_y])
                    else:
                        line = ax.plot([start_x, end_x], [start_y, end_y], color=path_color)

                if label_paths:
                    # annotate with link strength
                    line_color = line[-1].get_color()
                    if path_tuple[2] is not None:
                        if type(path_tuple[2]) == str:
                            data_str = path_tuple[2]
                        elif type(path_tuple[2]) in [float, int, np.float64]:
                            data_str = "{0:.3g}".format(path_tuple[2])
                        else:
                            raise Exception("Don't recognise TYPE of path label.")
                        ax.annotate(data_str, (mid_x, mid_y), color=line_color, fontsize=6, ha='center', va='center',
                                    path_effects=[pe.withStroke(linewidth=1, foreground="black")])
    plt.xlim([min_val[0], max_val[0]])
    plt.ylim([min_val[1], max_val[1]])
    if use_color_bar:
        # Draw a basic color map from 0 - 1
        if use_colors and not is_one_color:
            custom_color_map = cm.get_cmap('hsv')
        else:
            custom_map = [(0.0, [0.0, 0.0, 0.0]), (1.0, [1.0, 1.0, 1.0])]
            custom_color_map = LinearSegmentedColormap.from_list("", custom_map)
        c_mappable = ScalarMappable(cmap=custom_color_map)
        c_bar = plt.colorbar(c_mappable)
        # Then change the legend of the color bar to match the actual data
        c_bar.ax.get_yaxis().set_ticks([])
        c_bar_min_pos = 3.4
        c_bar_max_pos = 3.4
        if min_color_value < 0.0:
            c_bar_min_pos = 3.8
        if max_color_value < 0.0:
            c_bar_max_pos = 3.8
        c_bar.ax.text(c_bar_min_pos, color_offset, f'${"{:.4g}".format(min_color_value)}$', ha='center', va='center')
        c_bar.ax.text(c_bar_max_pos, 1.0, f'${"{:.4g}".format(max_color_value)}$', ha='center', va='center')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(file_path, dpi=400)
    plt.close(fig)


def data_extractor(dictionary_object):
    # this is a small helper function that unpacks the patch_{quality/degree/centrality} change-step dictionaries into
    # a form suitable for the time-series plotting
    extracted_data = []
    for _ in range(3):
        extracted_data.append({x: y[_] for x, y in dictionary_object.items()})
    return extracted_data


# ------------------------------------- MAIN FUNCTIONS FOR PLOTTING TIME-SERIES ------------------------------------- #

def global_species_time_series_properties(
        patch_list, species_set, parameters, sim, step, current_num_patches_history, is_save_data, is_save_plots):
    # This is used to prepare global time-series of species-specific data streams. Depending on the arguments passed,
    # it may produce the outputs as .csv, or as plots, or both. These options enabling separate handling by the
    # all_plots() and save_all_data() functions, although it will be inefficient if we are using both.
    num_time_steps = step + 1
    for species in species_set["list"]:
        species_global_population, species_global_population_change, species_global_patches_occupied, \
            species_global_patches_occupied_change, species_internal_population_change, species_dispersal_out, \
            species_dispersal_in, species_patches_colonised, species_patches_extinct, species_source_average, \
            species_sink_average = (np.zeros([num_time_steps, 1]) for _ in range(11))
        species_source_index = np.zeros([num_time_steps, 3])
        species_sink_index = np.zeros([num_time_steps, 3])

        for patch in patch_list:
            # locate the species
            for local_pop in patch.local_populations.values():
                if local_pop.species == species:
                    # access the local population time-series
                    for _ in range(num_time_steps):
                        species_global_population[_] += float(local_pop.population_history[_])
                        species_internal_population_change[_] += float(local_pop.internal_change_history[_])
                        species_dispersal_out[_] += float(local_pop.population_leave_history[_])
                        species_dispersal_in[_] += float(local_pop.population_enter_history[_])
                        sink = 0.0
                        source = 0.0
                        if float(local_pop.population_history[_]) >= local_pop.species.minimum_population_size:
                            species_global_patches_occupied[_] += 1
                            sink = local_pop.population_enter_history[_] / local_pop.population_history[_]
                            if _ > 0 and float(local_pop.population_history[_ - 1]) <\
                                    local_pop.species.minimum_population_size:
                                species_patches_colonised[_] += 1
                                species_global_patches_occupied_change[_] += 1
                        elif _ > 0 and float(local_pop.population_history[_ - 1]) >=\
                                local_pop.species.minimum_population_size:
                            species_patches_extinct[_] += 1
                            species_global_patches_occupied_change[_] += 1
                        if _ > 0 and local_pop.population_history[_ - 1] > 0.0:
                            source = local_pop.population_leave_history[_] / local_pop.population_history[_ - 1]

                        # update running totals for global source and sink calculations at that step
                        species_sink_average[_] += sink
                        species_source_average[_] += source
                        if sink >= 0.2:
                            species_sink_index[_, 0] += 1
                            if sink >= 0.5:
                                species_sink_index[_, 1] += 1
                                if sink >= 0.8:
                                    species_sink_index[_, 2] += 1
                        if source >= 0.2:
                            species_source_index[_, 0] += 1
                            if source >= 0.5:
                                species_source_index[_, 1] += 1
                                if source >= 0.8:
                                    species_source_index[_, 2] += 1
                    break

        # normalising
        for _ in range(num_time_steps):
            species_sink_average[_] /= current_num_patches_history[_]
            species_source_average[_] /= current_num_patches_history[_]
            species_sink_index[_, :] /= current_num_patches_history[_]
            species_source_index[_, :] /= current_num_patches_history[_]
            if _ > 0:
                species_global_population_change[_] = species_global_population[_] - species_global_population[_ - 1]

        global_property_dict = {
            "population": {
                "data": species_global_population,
                "y_label": "Global population",
                "legend_list": None,
            },
            "population_change": {
                "data": species_global_population_change,
                "y_label": "Global population change",
                "legend_list": None,
            },
            "internal_population_change": {
                "data": species_internal_population_change,
                "y_label": "Global internal population change",
                "legend_list": None,
            },
            "dispersal_out": {
                "data": species_dispersal_out,
                "y_label": "Global dispersal (out)",
                "legend_list": None,
            },
            "dispersal_in": {
                "data": species_dispersal_in,
                "y_label": "Global dispersal (in)",
                "legend_list": None,
            },
            "patches_occupied": {
                "data": species_global_patches_occupied,
                "y_label": "Number of patches occupied",
                "legend_list": None,
            },
            "patches_occupied_status_change": {
                "data": species_global_patches_occupied_change,
                "y_label": "Number of patches whose status changed \n (colonisation or extinction)",
                "legend_list": None,
            },
            "patches_colonised": {
                "data": species_patches_colonised,
                "y_label": "Number of new patches colonised",
                "legend_list": None,
            },
            "patches_extinct": {
                "data": species_patches_extinct,
                "y_label": "Number of patch extinctions of local populations",
                "legend_list": None,
            },
            "sink_average": {
                "data": species_sink_average,
                "y_label": "Patch-average sink value",
                "legend_list": None,
            },
            "source_average": {
                "data": species_source_average,
                "y_label": "Patch-average source value",
                "legend_list": None,
            },
            "sink_index": {
                "data": species_sink_index,
                "y_label": "Sink indexes",
                "legend_list": ["0.2", "0.5", "0.8"],
            },
            "source_index": {
                "data": species_source_index,
                "y_label": "Source indexes",
                "legend_list": ["0.2", "0.5", "0.8"],
            },
        }
        for _ in global_property_dict:
            if is_save_plots:
                file_path = f"results/{sim}/{step}/figures/species_global_ts_{species.name}_" + _ + ".png"
                create_time_series_plot(data=[global_property_dict[_]["data"]],
                                        parameters=parameters,
                                        y_label=global_property_dict[_]["y_label"],
                                        legend_list=global_property_dict[_]["legend_list"],
                                        end_time=step + 1,
                                        file_path=file_path)
            if is_save_data:
                file_path = f"results/{sim}/{step}/data/species_global_ts_{species.name}_" + _ + ".csv"
                with safe_open_w(file_path) as f:
                    # noinspection PyTypeChecker
                    np.savetxt(f, global_property_dict[_]["data"], delimiter=', ', newline='\n', fmt='%f')


def plot_local_time_series(patch_list, species_set, parameters, sim, step, is_local_plots):
    local_properties_to_plot = [
        {"attribute_y_label": "Population",
         "attribute_filename": "population",
         "attribute_id": "population_history",
         "rolling_average_steps": 1},
        {"attribute_y_label": "Population",
         "attribute_filename": "population_RA",
         "attribute_id": "population_history",
         "rolling_average_steps": 100},
    ]

    for prop_dict in local_properties_to_plot:
        num_time_steps = step + 1
        ra_end_time = num_time_steps - prop_dict["rolling_average_steps"] + 1
        if ra_end_time > 0:
            time_series_array = np.zeros([len(patch_list), len(species_set["list"]), ra_end_time])
            for species_number, species in enumerate(species_set["list"]):
                for patch in patch_list:
                    for local_population in patch.local_populations.values():
                        if local_population.species == species:
                            # access the local population time-series of the desired attribute
                            local_attribute_base = np.array(getattr(local_population, prop_dict["attribute_id"])[1:])

                            # calculate rolling average
                            for _ in range(ra_end_time):
                                for step_ahead in range(prop_dict["rolling_average_steps"]):
                                    time_series_array[patch.number, species_number, _] += local_attribute_base[
                                        _ + step_ahead]
                                time_series_array[
                                    patch.number, species_number, _] /= prop_dict["rolling_average_steps"]

                            # plot
                            if is_local_plots:
                                file_path = f"results/{sim}/{step}/figures/local_time_series/species_specific/" \
                                            f"species_local_ts_{species.name}_{patch.number}_" \
                                            f"{prop_dict['attribute_filename']}.png"
                                create_time_series_plot(data=[time_series_array[patch.number, species_number, :]],
                                                        parameters=parameters, y_label=prop_dict["attribute_y_label"],
                                                        end_time=ra_end_time, file_path=file_path)

            # plot entire-species:
            for _ in range(len(species_set["list"])):
                legend_list = [patch.number for patch in patch_list]
                file_path = f"results/{sim}/{step}/figures/species_local_ts_{species_set['list'][_].name}_" \
                            f"{prop_dict['attribute_filename']}_all.png"
                create_time_series_plot(data=time_series_array[:, _], parameters=parameters,
                                        y_label=prop_dict["attribute_y_label"],
                                        legend_list=legend_list, end_time=ra_end_time, file_path=file_path)
            # plot entire-patch:
            if is_local_plots:
                for _, patch in enumerate(patch_list):
                    legend_list = [species.name for species in species_set["list"]]
                    file_path = f"results/{sim}/{step}/figures/local_time_series/all_species/" \
                                f"species_local_ts_all_species_{patch.number}_{prop_dict['attribute_filename']}.png"
                    create_time_series_plot(data=time_series_array[_, :], parameters=parameters,
                                            y_label=prop_dict["attribute_y_label"], legend_list=legend_list,
                                            end_time=ra_end_time, file_path=file_path)
        else:
            print(f"Could not create a rolling average plot of {prop_dict['attribute_y_label']} as"
                  f" ra_end_time = {ra_end_time}.")


def plot_current_local_population_attribute(patch_list, sim, attribute_name, step, species, sub_attr=None):
    attribute_matrix = np.zeros([len(patch_list), 1])
    for patch in patch_list:
        # locate the species
        for local_population in patch.local_populations.values():
            if local_population.species == species:
                value = getattr(local_population, attribute_name)
                if sub_attr is not None:
                    value = value[sub_attr]
                attribute_matrix[patch.number, 0] = value
                break
    if sub_attr is not None:
        file_path = f"results/{sim}/{step}/figures/species_{attribute_name}_{sub_attr}_{species.name}.png"
    else:
        file_path = f"results/{sim}/{step}/figures/species_{attribute_name}_{species.name}.png"
    create_patches_plot(patch_list=patch_list, color_property=attribute_matrix,
                        file_path=file_path, use_color_bar=True)


# ---------------------------------------------- CREATING PATCH PLOTS ---------------------------------------------- #

def plot_network_properties(patch_list, sim, step, adjacency_path_list, is_biodiversity, is_reserves, is_retro=False):
    properties = {
        'habitat_type':
            {'attribute_id': 'habitat_type_num',
             'use_color_bar': True,
             'label_patches': True,
             'patch_label_attr': None,
             'path_list': None,
             'path_color': None,
             'use_colors': False,
             'patch_label_color': None,
             },
        'patch_quality':
            {'attribute_id': 'quality',
             'use_color_bar': True,
             'label_patches': False,
             'patch_label_attr': None,
             'path_list': None,
             'path_color': None,
             'use_colors': False,
             'patch_label_color': None,
             },
        'patch_adjacency':
            {'attribute_id': 'degree',
             'use_color_bar': True,
             'label_patches': False,
             'patch_label_attr': None,
             'path_list': adjacency_path_list,
             'path_color': [0.7, 0.7, 0.7],
             'use_colors': False,
             'patch_label_color': None,
             },
        'centrality':
            {'attribute_id': 'centrality',
             'use_color_bar': True,
             'label_patches': False,
             'patch_label_attr': None,
             'path_list': adjacency_path_list,
             'path_color': [0.7, 0.7, 0.7],
             'use_colors': False,
             'patch_label_color': None,
             },
        'perturbations':
            {'attribute_id': 'latest_perturbation_code',  # 0 = reserve, 0.1 = normal, 0.2 - 1.0 is pert. cluster
             'use_color_bar': True,
             'label_patches': True,
             'patch_label_attr': 'num_times_perturbed',
             'path_list': adjacency_path_list,
             'path_color': [0.01, 0.2, 0.2],
             'use_colors': True,  # easier to see differences in perturbation clusters and reserves.
             'patch_label_color': 'white',
             },
        'lcc':
            {'attribute_id': 'local_clustering',
             'use_color_bar': True,
             'label_patches': False,
             'patch_label_attr': None,
             'path_list': adjacency_path_list,
             'path_color': [0.7, 0.7, 0.8],
             'use_colors': False,
             'patch_label_color': None,
             },
    }

    if is_biodiversity:
        # count biodiversity - we don't include this in a retrospective spatial network plot as that only
        # rebuilds physical properties of the network as they would have been at that time step
        for patch in patch_list:
            patch.update_biodiversity()
        properties['biodiversity'] = {
            'attribute_id': 'biodiversity',
            'use_color_bar': True,
            'label_patches': False,
            'patch_label_attr': None,
            'path_list': None,
            'path_color': None,
            'use_colors': False,
            'patch_label_color': None,
        }

    if is_reserves:
        # we do not currently have the ability to change reserves, so until that is implemented there is no benefit
        # to re-drawing them from the perspective of a given time-step - hence this option will be False when calling
        # from the retrospective function
        properties['reserves'] = {
            'attribute_id': 'is_reserve',
            'use_color_bar': True,
            'label_patches': True,
            'patch_label_attr': 'reserve_order',
            'path_list': adjacency_path_list,
            'path_color': [0.7, 0.7, 0.7],
            'use_colors': False,
            'patch_label_color': None,
        }

    # then iterate over list of patch properties to plot according to their options
    for prop in properties:
        color_matrix = np.zeros([len(patch_list), 1])
        for patch in patch_list:
            color_matrix[patch.number, 0] = getattr(patch, properties[prop]["attribute_id"])
        if is_retro:
            # modify the file name so that retrospective and 'live' plots are distinguishable
            file_path = f"results/{sim}/{step}/figures/{prop}_retro.png"
        else:
            file_path = f"results/{sim}/{step}/figures/{prop}.png"
        create_patches_plot(patch_list=patch_list, color_property=color_matrix, file_path=file_path,
                            use_color_bar=properties[prop]["use_color_bar"],
                            label_patches=properties[prop]["label_patches"],
                            patch_label_attr=properties[prop]["patch_label_attr"],
                            path_color=properties[prop]["path_color"],
                            path_list=properties[prop]["path_list"],
                            use_colors=properties[prop]["use_colors"],
                            patch_label_color=properties[prop]["patch_label_color"],
                            )


def retrospective_network_plots(initial_patch_list, actual_patch_list, initial_patch_adjacency_matrix, sim, step):
    # this function is currently unused, but the idea is that it could be called after a reload or re-analysis
    # and then used to plot the abiotic properties of the spatial network at any given time (i.e. after various
    # perturbations) without having to re-run the simulation, since we save both the initial state of the network
    # (as the initial patch list) and record for each patch all of the changes that have been made due to perturbation!
    modified_patch_list = deepcopy(initial_patch_list)  # don't change the initial list!
    mod_patch_adjacency_matrix = deepcopy(initial_patch_adjacency_matrix)  # don't change the initial matrix!
    #
    # Now work through the time-recorded patch perturbations and modify the network as required!
    for patch in actual_patch_list:
        scalar_histories = {
            "habitat_history": "habitat_type_num",
            "quality_history": "quality",
            "size_history": "size",
            "degree_history": "degree",
            "centrality_history": "centrality",
            "latest_perturbation_code_history": "latest_perturbation_code",
        }
        for history_dict, change_quantity in scalar_histories.items():
            if len(getattr(patch, history_dict)) > 0:
                for change_time, change_value in getattr(patch, history_dict).items():
                    if change_time <= step:
                        # implement the change in the INITIAL patch list (the deepcopy)
                        setattr(modified_patch_list[patch.number], change_quantity, change_value)
                    else:
                        # Too far - have passed the requested time
                        break
        # the adjacency_history_list is a bit more complicated:
        for adj_change in patch.adjacency_history_list:
            # check the step
            if adj_change[0] <= step:
                # again we assume adjacency is undirected with a symmetric matrix
                mod_patch_adjacency_matrix[patch.number, adj_change[1]] = adj_change[2]
                mod_patch_adjacency_matrix[adj_change[1], patch.number] = adj_change[2]
        # count the total number of (potentially species-specific) perturbations to date in this patch
        num_perturbations_increase = getattr(modified_patch_list[patch.number], "num_times_perturbed")
        for perturbation in patch.perturbation_history_list:
            if perturbation[0] <= step:
                num_perturbations_increase += 1
        setattr(modified_patch_list[patch.number], "num_times_perturbed", num_perturbations_increase)
    # now build the list of adjacency paths to draw on, based on this matrix
    adjacency_path_list = create_adjacency_path_list(patch_list=modified_patch_list,
                                                     patch_adjacency_matrix=mod_patch_adjacency_matrix)
    plot_network_properties(patch_list=modified_patch_list, sim=sim, step=step,
                            adjacency_path_list=adjacency_path_list,
                            is_biodiversity=False, is_reserves=False, is_retro=True)


# --------------------------------------- SPECIALIST PATCH PLOTTING FUNCTIONS --------------------------------------- #

def plot_interactions(patch_list, adjacency_path_list, sim, step):
    # This is a tool to allow you to actually see the range and spatial feeding relationships impacting a given local
    # population, and the dispersal that is actually occurring at this particular time-step.
    # You can't see the amounts by which these individually contribute to prey-gain (see JSON files for net of this),
    # but you can at least observe which local populations are accessible and are actively being predated upon to some
    # amount.
    habitat_matrix = np.zeros([len(patch_list), 1])
    for patch in patch_list:
        # What is the habitat type number?
        habitat_matrix[patch.number, 0] = patch.habitat_type_num

    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            interacting_populations = local_pop.interacting_populations
            leaving_array = local_pop.leaving_array

            # prey sources
            prey_path_list = []
            for other_population in interacting_populations:
                # skip same patch populations and non-prey populations (e.g. same species) as this
                # can cause overlapping paths that are not visible on the image
                if other_population["object"].name in local_pop.species.prey_dict and \
                        not other_population["is_same_patch"]:
                    # color according to the feeding relationship:
                    if other_population["object"].population > 0.0:
                        # REALLY ACTUALLY FEEDING NOW
                        if local_pop.population > 0.0:
                            color = [0.08, 0.93, 0.1]  # green
                        else:
                            # would really feed if MY own population WAS alive
                            color = [0.93, 0.08, 0.1]  # red
                    else:
                        # Would feed if THEIR population was alive (and if mine was alive, which may or may not be)
                        color = [0.08, 0.13, 0.93]  # blue
                    prey_path_list.append((patch.number, other_population["object"].patch_num, None, color))

            # general interacting species list (possible predators and prey, and reachable same species, regardless of
            # any parameter options that my restrict and prevent the actualisation of some of these relationships)
            interacting_pop_path_list = []
            for other_population in interacting_populations:
                if not other_population["is_same_patch"]:
                    color = [0.3, 0.3, 0.3]
                    interacting_pop_path_list.append((patch.number, other_population["object"].patch_num, None, color))

            # dispersal
            dispersal_path_list = []
            for destination in range(len(patch_list)):
                if leaving_array[destination] != 0.0:
                    # this shows the ACTUAL final destinations, so no checks required
                    dispersal_path_list.append((patch.number, destination,
                                                leaving_array[destination][0], [0.3, 0.3, 0.3]))

            path_lists = {
                "prey": {
                    "data": prey_path_list,
                    "name": "prey",
                },
                "interacting_populations": {
                    "data": interacting_pop_path_list,
                    "name": "ip",
                },
                "dispersal": {
                    "data": dispersal_path_list,
                    "name": "dl",
                }
            }

            for path_list_dict in path_lists.values():
                if path_list_dict["data"] is not None and len(path_list_dict["data"]) > 0:
                    file_path = f"results/{sim}/{step}/figures/interactions/" \
                                f"local_population_{patch.number}_{local_pop.name}_{path_list_dict['name']}.png"
                    create_patches_plot(patch_list=patch_list,
                                        file_path=file_path,
                                        color_property=habitat_matrix,
                                        path_list=adjacency_path_list + path_list_dict["data"],
                                        label_paths=True,
                                        use_color_bar=True,
                                        path_color=[0.7, 0.7, 0.7]
                                        )


def plot_unrestricted_shortest_paths(patch_list, species_set, sim, step):
    # this plots all of the reachable foraging/direct migration paths for each species WITHOUT the restrictions on
    # threshold and path length (these actual interactions are shown in the individual "interaction" plots)
    #
    # Hence this is not likely to be of much use anymore, unless you had a very sparse, disconnected graph and wanted
    # to see the traversable paths through the network - however you would probably see this more easily (apart from
    # the movement scores) with the the species-specific network maps that color the network by disconnections as
    # perceived by the species (although it may not be possible, even in principle, to travel between any two patches
    # in a connected subgraph in one time-step - which is what these shortest paths are supposed to illustrate).
    #
    # patch color is the degree of the patch
    habitat_matrix = np.zeros([len(patch_list), 1])
    for patch in patch_list:
        habitat_matrix[patch.number, 0] = patch.degree
    for species in species_set["list"]:
        # generate list of [ (start, stop, annotation) ] tuples for paths between patches
        path_list = []
        for patch in patch_list:
            adjacency_list = getattr(patch, "adjacency_lists")
            for reachable_patch_num in adjacency_list[species.name]:
                if reachable_patch_num != patch.number:
                    reachable_patch_score = patch.species_movement_scores[
                        species.name][reachable_patch_num]["best"][2]
                    path_list.append((patch.number, reachable_patch_num, reachable_patch_score, []))
        file_path = f"results/{sim}/{step}/figures/unrestricted_best_patch_paths_{species.name}.png"
        create_patches_plot(patch_list=patch_list, color_property=habitat_matrix,
                            file_path=file_path, path_list=path_list, label_paths=True, use_color_bar=False)


def plot_adjacency_sub_graphs(system_state, sim):
    # This visualises the undirected adjacency-based sub-graphs, regardless of any species actual ability to
    # traverse them. Since adjacency can, in principle, be one-way, it is not necessarily the case that
    # all sub-graphs produced are fully-connected. However, the only alternative is to bias by some rule (such as
    # favouring larger sub-graphs as in "plot_accessible_sub_graphs()") or to produce a unique plot per patch.
    #
    # Iterate through the patches
    max_patches = len(system_state.patch_list)
    unclassified_patches = [x for x in range(max_patches)]
    color_classification = -1 * np.ones([max_patches, 1])
    next_classifier = -1
    patch_adjacency_matrix = system_state.patch_adjacency_matrix

    while len(unclassified_patches) > 0:
        # next color
        next_classifier += 1
        # next unvisited patch
        starting_patch = unclassified_patches[0]

        # for this starting patch, find all accessible patches from (OR TO) this patch:
        current_list = []
        not_checked = [starting_patch]
        while len(not_checked) > 0:
            for patch_from in not_checked:
                for patch_to in range(max_patches):
                    if patch_from != patch_to and (patch_adjacency_matrix[patch_from, patch_to] != 0.0 or
                                                   patch_adjacency_matrix[patch_to, patch_from] != 0.0) and \
                            patch_to not in current_list and patch_to not in not_checked:
                        not_checked.append(patch_to)
                not_checked.remove(patch_from)
                current_list.append(patch_from)
        # remove classified patches from list and color the array
        for patch_num in current_list:
            color_classification[patch_num] = next_classifier
            unclassified_patches.remove(patch_num)

    file_path = f"results/{sim}/{system_state.step}/figures/adjacency_subgraph.png"
    create_patches_plot(patch_list=system_state.patch_list, color_property=color_classification,
                        file_path=file_path, use_colors=True, label_patches=True, patch_label_color='black')


def plot_accessible_sub_graphs(patch_list, parameters, species, sim, step):
    # Visualising the subsets of all patches that are accessible to the given species;
    # Separate colors for each disconnected subgraph in the spatial network that they can eventually fully-explore.
    #
    # Iterate through the patches
    max_patches = parameters["main_para"]["NUM_PATCHES"]
    unclassified_patches = [x for x in range(max_patches)]
    color_classification = -1 * np.ones([max_patches, 1])
    subgraph_size = np.zeros([max_patches, 1])
    next_classifier = -1

    while len(unclassified_patches) > 0:
        # next color
        next_classifier += 1
        # next unvisited patch
        starting_patch = unclassified_patches[0]

        # for this starting patch, find all accessible patches
        current_list = []
        not_checked = [starting_patch]
        while len(not_checked) > 0:
            for patch_from in not_checked:
                for patch_to in patch_list[patch_from].adjacency_lists[species.name]:
                    if patch_to not in current_list and patch_to not in not_checked:
                        not_checked.append(patch_to)
                not_checked.remove(patch_from)
                current_list.append(patch_from)
        # count the total number of included nodes:
        subgraph_size[next_classifier] = len(current_list)
        # remove classified patches from list and color the array
        for patch_num in current_list:
            if color_classification[patch_num] == -1:
                color_classification[patch_num] = next_classifier
            else:
                # if already classified, prioritise the largest subgraph for coloring
                # (This is still not rigorous as it biases those constructed earlier. In theory this later subgraph
                # could be larger but repeatedly precluded from reaching its full size and overtaking.
                # However it should give an idea of domains for the species.)
                if subgraph_size[int(color_classification[patch_num])] < len(current_list):
                    color_classification[patch_num] = next_classifier
            if patch_num in unclassified_patches:
                # note that as the adjacency graph is directed, there may be some patches that can be reached from
                # multiple other sub-graphs
                unclassified_patches.remove(patch_num)

    # now print
    file_path = f"results/{sim}/{step}/figures/species_{species.name}_subgraph.png"
    create_patches_plot(patch_list=patch_list, color_property=color_classification,
                        file_path=file_path, use_colors=True)


# --------------------- AUXILIARY FUNCTIONS REQUIRED BY THE SPECIALIST PATCH PLOTTING FUNCTIONS --------------------- #

def create_adjacency_path_list(patch_list, patch_adjacency_matrix):
    path_list = []
    for patch_1 in range(len(patch_list)):
        for patch_2 in range(len(patch_list)):
            if patch_1 != patch_2 and patch_adjacency_matrix[patch_1, patch_2] != 0.0:
                path_list.append((patch_1, patch_2, None, []))
    return path_list


def find_connected_sets(patch_list, starting_node, scale, list_of_lists):
    this_set = {starting_node}
    reachable_patches = patch_list[starting_node].set_of_adjacent_patches
    if starting_node in reachable_patches:
        reachable_patches.remove(starting_node)  # should be unnecessary
    for next_node in range(1, scale):
        # choose the next one and add
        if len(reachable_patches) == 0:
            break
        else:
            node = random.choice(list(reachable_patches))
            this_set.add(node)
            reachable_patches.remove(node)
            for x in patch_list[node].set_of_adjacent_patches:
                if x not in this_set:
                    reachable_patches.add(x)
    if len(this_set) == scale:
        is_new = True
        # check it is not a duplicate
        if len(list_of_lists) > 0:
            for prev_list in list_of_lists:
                if set(prev_list) == this_set:
                    is_new = False
                    break
        if is_new:
            temp_list = list(this_set)
            temp_list.sort()
            list_of_lists.append(temp_list)  # list of sorted lists of patches
    return list_of_lists


# --------------------------------------- PRODUCING SPECIES-AREA CURVES (SAR) --------------------------------------- #

def biodiversity_analysis(patch_list, species_set, parameters, sim, step):
    # note that this can take a long time if there is a large number (>200) of patches
    max_patches = parameters["main_para"]["NUM_PATCHES"]
    minimum_tries = parameters["plot_save_para"]["MIN_BIODIVERSITY_ATTEMPTS"]
    biodiversity_output = np.zeros([max_patches, 3])

    # consider biodiversity over a sample of areas with size from 1-N patches (where N is the total number of patches)
    for scale in range(1, max_patches + 1):
        # Choosing every patch at least once as the starting patch:
        # - Build a list of reachable patches which is initially just the adjacent patches to the starting patch
        # - Randomly choose one and add it to the set of visited patches, remove it from the reachable patches, and
        #       add any of its adjacent patches to the reachable patches(not including repeats or already visited)
        # - Repeat until the size is reached or there are no reachable patches
        # - Finally order the sets of patches, and check that the same set of patches are not repeated.

        list_of_lists = []
        # try once starting from every patch
        for starting_node in range(max_patches):
            list_of_lists = find_connected_sets(patch_list=patch_list,
                                                starting_node=starting_node,
                                                scale=scale,
                                                list_of_lists=list_of_lists)
        # if don't have as many as desired, try some more random attempts
        if len(list_of_lists) < minimum_tries:
            for attempt in range(minimum_tries):
                starting_node = random.randint(0, max_patches - 1)
                list_of_lists = find_connected_sets(patch_list=patch_list,
                                                    starting_node=starting_node,
                                                    scale=scale,
                                                    list_of_lists=list_of_lists)
        if len(list_of_lists) > 0:
            # count the unique species
            ave_diversity = 0.0
            for current_patch_list in list_of_lists:

                # Count the biodiversity in a given connected subgraph:
                unique_species_set = set({})  # of species objects (not names)
                break_loop = False
                for patch_num in current_patch_list:
                    for local_pop in patch_list[patch_num].local_populations.values():
                        if local_pop.species not in unique_species_set and local_pop.population \
                                > local_pop.species.minimum_population_size:
                            unique_species_set.add(local_pop.species)
                            # stop if we have now found all species in the system
                            if len(unique_species_set) == len(species_set["list"]):
                                break_loop = True
                                break
                    if break_loop:
                        break

                actual_diversity = len(unique_species_set)

                ave_diversity += float(actual_diversity)
            ave_diversity = ave_diversity / float(len(list_of_lists))
        else:
            # set to zero by default
            ave_diversity = 0.0
        # record results including how many unique sets of patches were tested
        biodiversity_output[scale - 1, 0] = scale
        biodiversity_output[scale - 1, 1] = ave_diversity
        biodiversity_output[scale - 1, 2] = len(list_of_lists)

    # plot the species-area curve
    x_data = biodiversity_output[:, 0]
    y_data = biodiversity_output[:, 1]
    z_data = biodiversity_output[:, 2]
    fig, ax1 = plt.subplots()
    ax1_color = 'blue'
    ax1.plot(x_data, z_data, color=ax1_color)
    ax1.set_xlabel('Scale')
    ax1.set_ylabel('Sets of patches tested', color=ax1_color)
    ax2 = ax1.twinx()
    ax2_color = 'red'
    ax2.plot(x_data, y_data, color=ax2_color)
    ax2.set_ylabel('Average biodiversity', color=ax2_color)
    file_path = f"results/{sim}/{step}/figures/species_area_curve.png"
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(file_path, dpi=400)
    plt.close(fig)
