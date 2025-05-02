import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as pe
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import os.path
from copy import deepcopy

from scipy.stats import alpha

from data_core_functions import *


# ----------------------------------------------- FIGURE MANAGEMENT ----------------------------------------------- #

def print_and_close(fig, file_path):
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(file_path, dpi=400)
    plt.close(fig)

# ---------------------------------------- BASE FUNCTIONS FOR CREATING PLOTS ---------------------------------------- #

def create_time_series_plot(data, parameters, file_path, y_label='Population', legend_list=None,
                            start_time=1, end_time=None, is_shading=False, is_triple=False):
    if end_time is None:
        end_time = parameters["main_para"]["NUM_TRANSIENT_STEPS"] + parameters["main_para"]["NUM_RECORD_STEPS"]
    x_data = np.linspace(start_time, end_time, num=1 + end_time - start_time)
    fig = plt.figure()
    for series_number, data_time_series in enumerate(data):
        if is_triple:
            # if we are passing lists of triple tuples, then this is [(mean-1SD, mean, mean+1SD), ...] time-series
            #
            # plot the mean and ±1SD
            data_lower = [y[0] for y in data_time_series]
            data_mid = [y[1] for y in data_time_series]
            data_upper = [y[2] for y in data_time_series]
            line, = plt.plot(x_data, data_mid)  # mean
            current_colour = line.get_color()  # set the surrounding tube to be the same colour
            plt.plot(x_data, data_lower, label='_nolegend_', c=current_colour, alpha=0.3)  # mean - 1SD
            plt.plot(x_data, data_upper, label='_nolegend_', c=current_colour, alpha=0.3)  # mean + 1SD
            if is_shading:
                plt.fill_between(x_data, data_lower, data_mid, alpha=0.2, label='_nolegend_', color=current_colour)
                plt.fill_between(x_data, data_mid, data_upper, alpha=0.2, label='_nolegend_', color=current_colour)
        else:
            # regular single-series data streams are being passed to the function
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

    # space them out slightly so that size=1 patches can still be distinguished
    patch_scaling_factor = 1.5

    # Some set up for the colorings
    color_offset = 0.0
    is_one_color = False
    if np.max(color_property) - np.min(color_property) == 0.0:
        is_one_color = True
    elif use_colors:
        color_offset = 0.15  # for cyclic color maps (e.g. hsv) [0] and [1] are the same, so need an offset

    # Now the primary iteration of each patch to determine its position, value and color
    line_width = min(1.0, 20.0 / np.sqrt(len(patch_list)))  # both borders and paths
    for patch_num, patch in enumerate(patch_list):
        length = np.sqrt(patch.size)

        position = (patch.position[0] * patch_scaling_factor, patch.position[1] * patch_scaling_factor)
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

        rectangle = patches.Rectangle(position, length, length, linewidth=line_width, edgecolor='black',
                                      facecolor=c, zorder=1.5)
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
                            color=color, fontsize=5, ha='center', va='center', zorder=1.9)  # highest z-order for labels

    # draw paths (showing traversable links) if required (with a lower z-order than the patches):
    #
    if path_list is not None and len(path_list) > 0:

        min_x = min([element.position[0] for element in patch_list])
        max_x = max([element.position[0] for element in patch_list])
        min_y = min([element.position[1] for element in patch_list])
        max_y = max([element.position[1] for element in patch_list])

        for path_tuple in path_list:

            # check for wrap
            is_x_wrap = False
            is_y_wrap = False
            if min_x + 1 < max_x and min_y + 1 < max_y:
                if patch_list[path_tuple[0]].position[0] == patch_list[path_tuple[1]].position[0]:
                    if {min_y, max_y} <= {patch_list[path_tuple[0]].position[1], patch_list[path_tuple[1]].position[1]}:
                        is_y_wrap = True
                if patch_list[path_tuple[0]].position[1] == patch_list[path_tuple[1]].position[1]:
                    if {min_x, max_x} <= {patch_list[path_tuple[0]].position[0], patch_list[path_tuple[1]].position[0]}:
                        is_x_wrap = True

            start_patch_radius = 0.5 * np.sqrt(patch_list[path_tuple[0]].size)
            end_patch_radius = 0.5 * np.sqrt(patch_list[path_tuple[1]].size)
            if is_x_wrap or is_y_wrap:
                patch_one = [patch_list[path_tuple[0]].position[0] * patch_scaling_factor + start_patch_radius,
                             patch_list[path_tuple[0]].position[1] * patch_scaling_factor + start_patch_radius]
                patch_two = [patch_list[path_tuple[1]].position[0] * patch_scaling_factor + end_patch_radius,
                             patch_list[path_tuple[1]].position[1] * patch_scaling_factor + end_patch_radius]
                both_patches = [patch_one, patch_two]
                for line_num in range(2):
                    x_start = both_patches[line_num][0]
                    if is_x_wrap:
                        x_end = both_patches[line_num][0] - (both_patches[line_num - 1][0] -
                                                             both_patches[line_num][0]) / (np.sqrt(len(patch_list)))
                    else:
                        x_end = both_patches[line_num][0] + (both_patches[line_num - 1][0] -
                                                             both_patches[line_num][0]) / 2
                    y_start = both_patches[line_num][1]
                    if is_y_wrap:
                        y_end = both_patches[line_num][1] - (both_patches[line_num - 1][1] -
                                                             both_patches[line_num][1]) / (np.sqrt(len(patch_list)))
                    else:
                        y_end = both_patches[line_num][1] + (both_patches[line_num - 1][1] -
                                                             both_patches[line_num][1]) / 2
                    if type(path_tuple[3]) is list and len(path_tuple[3]) == 3:
                        ax.plot([x_start, x_end], [y_start, y_end], color=path_tuple[3],
                                linewidth=line_width, zorder=1.1)
                    elif path_color is None:
                        ax.plot([x_start, x_end], [y_start, y_end], linewidth=line_width, zorder=1.1)
                    else:
                        ax.plot([x_start, x_end], [y_start, y_end], color=path_color, linewidth=line_width, zorder=1.1)

            else:
                # not a wrap - just a direct path
                start_x = patch_list[path_tuple[0]].position[0] * patch_scaling_factor + start_patch_radius
                end_x = patch_list[path_tuple[1]].position[0] * patch_scaling_factor + end_patch_radius
                start_y = patch_list[path_tuple[0]].position[1] * patch_scaling_factor + start_patch_radius
                end_y = patch_list[path_tuple[1]].position[1] * patch_scaling_factor + end_patch_radius
                mid_x = (start_x + end_x) / 2
                mid_y = (start_y + end_y) / 2

                if type(path_tuple[3]) is list and len(path_tuple[3]) == 3:
                    # is color specified for this path?
                    line = ax.plot([start_x, end_x], [start_y, end_y], color=path_tuple[3],
                                   linewidth=line_width, zorder=1.1)
                else:
                    # otherwise the general coloring schemes
                    if path_color is None:
                        line = ax.plot([start_x, end_x], [start_y, end_y], linewidth=line_width, zorder=1.1)
                    else:
                        line = ax.plot([start_x, end_x], [start_y, end_y], color=path_color,
                                       linewidth=line_width, zorder=1.1)

                if label_paths:
                    # annotate with link strength - labels should have the highest z-order
                    line_color = line[-1].get_color()
                    if path_tuple[2] is not None:
                        if type(path_tuple[2]) == str:
                            data_str = path_tuple[2]
                        elif type(path_tuple[2]) in [float, int, np.float64]:
                            data_str = "{0:.3g}".format(path_tuple[2])
                        else:
                            raise Exception("Don't recognise TYPE of path label.")
                        ax.annotate(data_str, (mid_x, mid_y), color=line_color, fontsize=6, ha='center', va='center',
                                    path_effects=[pe.withStroke(linewidth=1, foreground="black")], zorder=1.9)
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
    print_and_close(fig, file_path)


# ------------------------------------- MAIN FUNCTIONS FOR PLOTTING TIME-SERIES ------------------------------------- #

def global_species_time_series_properties(
        patch_list, species_set, parameters, sim_path, step, current_num_patches_history, is_save_data, is_save_plots):
    # This is used to prepare global time-series of species-specific data streams. Depending on the arguments passed,
    # it may produce the outputs as .csv, or as plots, or both. These options enabling separate handling by the
    # all_plots() and save_all_data() functions, although it will be inefficient if we are using both.
    num_time_steps = step + 1
    for species in species_set["list"]:
        species_global_population, species_global_population_change, species_global_patches_occupied, \
        species_global_patches_occupied_change, species_internal_population_change, species_dispersal_out, \
        species_dispersal_in, species_patches_colonised, species_patches_extinct, species_source_average, \
        species_sink_average = (np.zeros([num_time_steps, 1]) for _ in range(11))
        species_source_index = np.zeros([num_time_steps, 5])
        species_sink_index = np.zeros([num_time_steps, 5])

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

                        # source and sink
                        sink = 0.0
                        source = 0.0
                        positive_change = max(0.0, local_pop.population_enter_history[
                            _] - local_pop.population_leave_history[_]) + max(0.0, local_pop.internal_change_history[_])
                        if positive_change > 0:
                            sink = max(0.0, local_pop.population_enter_history[_] -
                                       local_pop.population_leave_history[_]) / positive_change
                        if local_pop.potential_dispersal_history[_] > 0.0:
                            source = max(0.0, local_pop.population_leave_history[
                                _] - local_pop.population_enter_history[_]) / local_pop.potential_dispersal_history[_]

                        # occupancy and change in status
                        if float(local_pop.population_history[_]) >= local_pop.species.minimum_population_size:
                            species_global_patches_occupied[_] += 1

                            if _ > 0 and float(local_pop.population_history[_ - 1]) < \
                                    local_pop.species.minimum_population_size:
                                species_patches_colonised[_] += 1
                                species_global_patches_occupied_change[_] += 1

                        elif _ > 0 and float(local_pop.population_history[_ - 1]) >= \
                                local_pop.species.minimum_population_size:
                            species_patches_extinct[_] += 1
                            species_global_patches_occupied_change[_] += 1

                        # update running totals for global source and sink calculations at that step
                        species_sink_average[_] += sink
                        species_source_average[_] += source
                        if sink >= 0.01:
                            species_sink_index[_, 0] += 1
                            if sink >= 0.05:
                                species_sink_index[_, 1] += 1
                                if sink >= 0.1:
                                    species_sink_index[_, 2] += 1
                                    if sink >= 0.4:
                                        species_sink_index[_, 3] += 1
                                        if sink >= 0.8:
                                            species_sink_index[_, 4] += 1
                        if source >= 0.01:
                            species_source_index[_, 0] += 1
                            if source >= 0.05:
                                species_source_index[_, 1] += 1
                                if source >= 0.1:
                                    species_source_index[_, 2] += 1
                                    if sink >= 0.4:
                                        species_source_index[_, 3] += 1
                                        if sink >= 0.8:
                                            species_source_index[_, 4] += 1
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
                "legend_list": ["0.01", "0.05", "0.1", "0.4", "0.8"],
            },
            "source_index": {
                "data": species_source_index,
                "y_label": "Source indexes",
                "legend_list": ["0.01", "0.05", "0.1", "0.4", "0.8"],
            },
        }
        for _ in global_property_dict:
            if is_save_plots:
                file_path = f"{sim_path}/{step}/figures/species_global_time_series/" \
                            f"species_global_ts_{species.name}_" + _ + ".png"
                create_time_series_plot(data=[global_property_dict[_]["data"]],
                                        parameters=parameters,
                                        y_label=global_property_dict[_]["y_label"],
                                        legend_list=global_property_dict[_]["legend_list"],
                                        end_time=step + 1,
                                        file_path=file_path)
            if is_save_data:
                file_path = f"{sim_path}/{step}/data/species_global_ts_{species.name}_" + _ + ".csv"
                with safe_open_w(file_path) as f:
                    # noinspection PyTypeChecker
                    np.savetxt(f, global_property_dict[_]["data"], delimiter=', ', newline='\n', fmt='%.20f')


def plot_local_time_series(patch_list, species_set, parameters, sim_path, step, is_local_plots):
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
                                file_path = f"{sim_path}/{step}/figures/local_time_series/species_specific/" \
                                            f"species_local_ts_{species.name}_{patch.number}_" \
                                            f"{prop_dict['attribute_filename']}.png"
                                create_time_series_plot(data=[time_series_array[patch.number, species_number, :]],
                                                        parameters=parameters, y_label=prop_dict["attribute_y_label"],
                                                        end_time=ra_end_time, file_path=file_path)

            # plot entire-species:
            for _ in range(len(species_set["list"])):
                legend_list = [patch.number for patch in patch_list]
                file_path = f"{sim_path}/{step}/figures/species_time_series/species_local" \
                            f"_ts_{species_set['list'][_].name}_{prop_dict['attribute_filename']}_all.png"
                create_time_series_plot(data=time_series_array[:, _], parameters=parameters,
                                        y_label=prop_dict["attribute_y_label"],
                                        legend_list=legend_list, end_time=ra_end_time, file_path=file_path)
            # plot entire-patch:
            if is_local_plots:
                for _, patch in enumerate(patch_list):
                    legend_list = [species.name for species in species_set["list"]]
                    file_path = f"{sim_path}/{step}/figures/local_time_series/all_species/" \
                                f"species_local_ts_all_species_{patch.number}_{prop_dict['attribute_filename']}.png"
                    create_time_series_plot(data=time_series_array[_, :], parameters=parameters,
                                            y_label=prop_dict["attribute_y_label"], legend_list=legend_list,
                                            end_time=ra_end_time, file_path=file_path)
        else:
            print(f"Could not create a rolling average plot of {prop_dict['attribute_y_label']} as"
                  f" ra_end_time = {ra_end_time}.")


def plot_current_local_population_attribute(patch_list, sim_path, attribute_name, step, species,
                                            adjacency_path_list, sub_attr=None):
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
        file_path = f"{sim_path}/{step}/figures/species_attributes/" \
                    f"species_{attribute_name}_{sub_attr}_{species.name}.png"
    else:
        file_path = f"{sim_path}/{step}/figures/species_attributes/species_{attribute_name}_{species.name}.png"
    create_patches_plot(patch_list=patch_list, color_property=attribute_matrix, path_list=adjacency_path_list,
                        path_color=[0.7, 0.7, 0.7], file_path=file_path, use_color_bar=True)


# ---------------------------------------------- CREATING PATCH PLOTS ---------------------------------------------- #

def plot_network_properties(patch_list, sim_path, step, adjacency_path_list,
                            is_biodiversity, is_reserves, is_partition, is_label_habitat_patches=True, is_retro=False):
    properties = {
        'habitat_type_numbered':
            {'attribute_id': 'habitat_type_num',
             "sub_attribute_list": [None],  # this should be a list containing None
             'use_color_bar': True,
             'label_patches': True,
             'patch_label_attr': 'number',
             'path_list': adjacency_path_list,
             'path_color': [0.7, 0.7, 0.7],
             'use_colors': False,
             'patch_label_color': None,
             },
        'habitat_type':
            {'attribute_id': 'habitat_type_num',
             "sub_attribute_list": [None],  # this should be a list containing None
             'use_color_bar': True,
             'label_patches': is_label_habitat_patches,
             'patch_label_attr': None,
             'path_list': adjacency_path_list,
             'path_color': [0.7, 0.7, 0.7],
             'use_colors': False,
             'patch_label_color': None,
             },
        'patch_quality':
            {'attribute_id': 'quality',
             "sub_attribute_list": [None],  # this should be a list containing None
             'use_color_bar': True,
             'label_patches': False,
             'patch_label_attr': None,
             'path_list': adjacency_path_list,
             'path_color': [0.7, 0.7, 0.7],
             'use_colors': False,
             'patch_label_color': None,
             },
        'patch_adjacency':
            {'attribute_id': 'degree',
             "sub_attribute_list": [None],  # this should be a list containing None
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
             "sub_attribute_list": [None],  # this should be a list containing None
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
             "sub_attribute_list": [None],  # this should be a list containing None
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
             "sub_attribute_list": ['all', 'same', 'different'],
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
            'sub_attribute_list': [None],  # this should be a list containing None
            'use_color_bar': True,
            'label_patches': False,
            'patch_label_attr': None,
            'path_list': adjacency_path_list,
            'path_color': [0.7, 0.7, 0.7],
            'use_colors': False,
            'patch_label_color': None,
        }

    if is_reserves:
        # we do not currently have the ability to change reserves, so until that is implemented there is no benefit
        # to re-drawing them from the perspective of a given time-step - hence this option will be False when calling
        # from the retrospective function
        properties['reserves'] = {
            'attribute_id': 'is_reserve',
            'sub_attribute_list': [None],  # this should be a list containing None
            'use_color_bar': True,
            'label_patches': True,
            'patch_label_attr': 'reserve_order',
            'path_list': adjacency_path_list,
            'path_color': [0.7, 0.7, 0.7],
            'use_colors': False,
            'patch_label_color': None,
        }

    if is_partition:
        properties['final_all_binary_partition'] = {
             'attribute_id': 'partition_code',
             "sub_attribute_list": [None],  # this should be a list containing None
             'use_color_bar': True,
             'label_patches': True,
             'patch_label_attr': 'partition_code',
             'path_list': adjacency_path_list,
             'path_color': [0.7, 0.7, 0.7],
             'use_colors': True,  # easier to visualise differences
             'patch_label_color': None,
        }

    file_pre_suffix = ""
    if is_retro:
        # modify the file name so that retrospective and 'live' plots are distinguishable
        file_pre_suffix += "_retro"

    # then iterate over list of patch properties to plot according to their options
    for prop in properties:
        for sub_attr in properties[prop]["sub_attribute_list"]:
            if sub_attr is None:
                file_suffix = file_pre_suffix
            else:
                file_suffix = file_pre_suffix + '_' + str(sub_attr)
            color_matrix = np.zeros([len(patch_list), 1])
            for patch in patch_list:
                if sub_attr is None:
                    patch_attribute = getattr(patch, properties[prop]["attribute_id"])
                    file_suffix = file_pre_suffix
                else:
                    patch_attribute = getattr(patch, properties[prop]["attribute_id"])[sub_attr]
                color_matrix[patch.number, 0] = patch_attribute

            file_path = f"{sim_path}/{step}/figures/network_properties/{prop}{file_suffix}.png"
            create_patches_plot(patch_list=patch_list, color_property=color_matrix, file_path=file_path,
                                use_color_bar=properties[prop]["use_color_bar"],
                                label_patches=properties[prop]["label_patches"],
                                patch_label_attr=properties[prop]["patch_label_attr"],
                                path_color=properties[prop]["path_color"],
                                path_list=properties[prop]["path_list"],
                                use_colors=properties[prop]["use_colors"],
                                patch_label_color=properties[prop]["patch_label_color"],
                                )


def retrospective_network_plots(initial_patch_list, actual_patch_list, initial_patch_adjacency_matrix, sim_path, step):
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
    plot_network_properties(patch_list=modified_patch_list, sim_path=sim_path, step=step,
                            adjacency_path_list=adjacency_path_list,
                            is_biodiversity=False, is_reserves=False, is_partition=False, is_retro=True)


# --------------------------------------- SPECIALIST PATCH PLOTTING FUNCTIONS --------------------------------------- #

def plot_interactions(patch_list, adjacency_path_list, sim_path, step):
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
                    file_path = f"{sim_path}/{step}/figures/interactions/" \
                                f"local_population_{patch.number}_{local_pop.name}_{path_list_dict['name']}.png"
                    create_patches_plot(patch_list=patch_list,
                                        file_path=file_path,
                                        color_property=habitat_matrix,
                                        path_list=adjacency_path_list + path_list_dict["data"],
                                        label_paths=True,
                                        use_color_bar=True,
                                        path_color=[0.7, 0.7, 0.7]
                                        )


def plot_unrestricted_shortest_paths(patch_list, species_set, sim_path, step):
    # this plots all the reachable foraging/direct migration paths for each species WITHOUT the restrictions on
    # threshold and path length (these actual interactions are shown in the individual "interaction" plots)
    #
    # Hence this is not likely to be of much use anymore, unless you had a very sparse, disconnected graph and wanted
    # to see the traversable paths through the network - however you would probably see this more easily (apart from
    # the movement scores) with the species-specific network maps that color the network by disconnections as
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
        file_path = f"{sim_path}/{step}/figures/unrestricted_best_patch_paths_{species.name}.png"
        create_patches_plot(patch_list=patch_list, color_property=habitat_matrix,
                            file_path=file_path, path_list=path_list, label_paths=True, use_color_bar=False)


def plot_adjacency_sub_graphs(system_state, sim_path):
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
            for patch_from in list(not_checked):  # iterate over DUMMY COPY of list
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

    file_path = f"{sim_path}/{system_state.step}/figures/adjacency_subgraph.png"
    create_patches_plot(patch_list=system_state.patch_list, color_property=color_classification,
                        file_path=file_path, use_colors=True, label_patches=True, patch_label_color='black')


def plot_accessible_sub_graphs(patch_list, parameters, species, sim_path, step):
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
            for patch_from in list(not_checked):  # iterate over DUMMY COPY of the list
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
    file_path = f"{sim_path}/{step}/figures/species_subgraph/species_{species.name}_subgraph.png"
    create_patches_plot(patch_list=patch_list, color_property=color_classification,
                        file_path=file_path, use_colors=True)


# ------------------------- PRODUCING DISTANCE, NETWORK, COMPLEXITY LINEAR REGRESSION PLOTS ------------------------- #

def partition_spectrum_plotting(distance_metrics_store, sim_path, step):
    # Plots the spectrum of max and min possible complexity of delta-partitions
    list_of_key_paths = recursive_dict_search(nested_dict=distance_metrics_store,
                                              seek_contains_key="is_partition_graphical",
                                              upper_key_path=None)
    if len(list_of_key_paths) > 0:
        # results found
        for path in list_of_key_paths:
            key_list = path.split("|")
            target_dict = return_dict_from_address(key_list, distance_metrics_store)
            # now target_dict is the results dictionary for this sub_network
            num_delta = len(target_dict["pw_sup_spectrum"])
            n_values = np.linspace(1, num_delta, num_delta)
            # pop-weighted and binary (if applicable)
            if len(target_dict["binary_sup_spectrum"]) > 0:
                y_key_list = ["pw", "binary"]
            else:
                y_key_list = ["pw"]
            for y_key in y_key_list:
                y_values_inf = np.asarray(target_dict[f"{y_key}_inf_spectrum"])
                y_values_mean = np.asarray(target_dict[f"{y_key}_mean_spectrum"])
                y_values_sup = np.asarray(target_dict[f"{y_key}_sup_spectrum"])
                natural_delta = target_dict[f"{y_key}_sup_minmax_delta"]
                target_height = target_dict[f"{y_key}_sup_target"]
                fig = plt.figure()
                plt.plot(n_values, y_values_sup, c='k', markersize=5, marker='o', mfc='white', mec='k')
                plt.plot(n_values, y_values_mean, c=[0.3, 0.3, 0.3], linewidth=1,
                         markersize=3, marker='o', mfc='white', mec='k', linestyle=':')
                plt.axhline(y=target_height, color=[0.2, 0.2, 0.2], linestyle=':', linewidth=1.5)
                plt.fill_between(n_values, y_values_inf, y_values_sup, alpha=0.2, color='grey')
                plt.axvline(x=natural_delta, color='k', linestyle='--')
                plt.xlabel("Partition delta")
                plt.ylabel("Partition complexity")
                plt.legend(["Supremum complexity over partitions", "Mean complexity over partitions",
                            "Target supremum complexity"], framealpha=1.0)
                plt.ylim([0, 1])
                print_name = path.replace('|', '_')
                file_path = f"{sim_path}/{step}/figures/complexity/{print_name}_{y_key}.png"
                print_and_close(fig, file_path)
    else:
        # If IS_PARTITION_ANALYSIS was false but IS_PLOT_DISTANCE_METRICS_LM was true, this routine will fail to find
        # anything, as the default partition outputs will be empty dictionaries without the identifying keys.
        print("Complexity analysis - no partition data found for printing spectra.")
        pass

def complexity_plotting(distance_metrics_store, sim_path, step):
    # While plot_distance_metrics_lm() will plot SAR and complexity dimension across the full range, we need this
    # additional function to plot the dc spectrum and Delta_C_per_N.
    list_of_key_paths = recursive_dict_search(nested_dict=distance_metrics_store,
                                              seek_contains_key="is_complexity_graphical",
                                              upper_key_path=None)
    if len(list_of_key_paths) > 0:
        # results found

        for path in list_of_key_paths:
            key_list = path.split("|")
            target_dict = return_dict_from_address(key_list, distance_metrics_store)
            # now target_dict is the graphical_results dictionary
            n_values = np.asarray(target_dict["n_vector"])

            # plot the single and dual fitted power laws for complexity
            plot_y = {}
            if target_dict["is_single_fit_success"] == 1:
                plot_y["single"] = target_dict["single_para"][0] * (n_values - 1) ** target_dict["single_para"][1]
            if target_dict["is_dual_fit_success"] == 1:
                p = target_dict["dual_para"]
                plot_y["dual"] = np.exp(-p[0] * (n_values - 1)) * p[1] * (n_values - 1) ** p[3] + (
                        1 - np.exp(-p[0] * (n_values - 1))) * (p[2] * (n_values - 1) ** p[4] + p[5])

            for law_type in ["single", "dual"]:
                if target_dict[f"is_{law_type}_fit_success"]:
                    fig = plt.figure()
                    plt.plot(n_values, plot_y[law_type], c='b')
                    plt.scatter(n_values, np.asarray(target_dict["complexity_vector"]), c='r', s=5, edgecolors='k')
                    plt.xlabel(r"$n$")
                    plt.ylabel(r"$C(n)$")
                    r_str = "{0:.5g}".format(target_dict[f"{law_type}_r_squared"])
                    plt.legend((f'Fit: $R^{2}$ = {r_str}', 'Actual values'))
                    print_name = path.replace('|', '_')
                    file_path = f"{sim_path}/{step}/figures/complexity/{print_name}_{law_type}.png"
                    print_and_close(fig, file_path)

            # plot the spectral properties
            y_label = [r"$\frac{\Delta C(n)}{\Delta n}$", r"$d_{c}$"]
            s_e_val = [[0, len(n_values)-1], [3, len(n_values)]]  # range of acceptable values to plot
            for _, y_key in enumerate(["delta_c_per_n", "dc_fit"]):
                y_values = np.asarray(target_dict[y_key])
                fig = plt.figure()
                plt.plot(n_values[s_e_val[_][0]:s_e_val[_][1]], y_values[s_e_val[_][0]:s_e_val[_][1]],
                         c='b', markersize=5, marker='o', mfc='r', mec='k')
                plt.xlabel(r"$n$")
                plt.ylabel(y_label[_])
                print_name = path.replace('|', '_')
                file_path = f"{sim_path}/{step}/figures/complexity/{print_name}_{y_key}.png"
                print_and_close(fig, file_path)


def plot_distance_metrics_lm(distance_metrics_store, sim_path, step):
    # produce plots of the linear models associated with distance metrics analysis, together with scatter plots of the
    # data they were fitted to, if this was collected.

    # recursively search the tree of the dictionary
    list_of_key_paths = recursive_dict_search(nested_dict=distance_metrics_store,
                                              seek_contains_key="is_linear_model",
                                              upper_key_path=None)
    if len(list_of_key_paths) > 0:
        # results found
        for path in list_of_key_paths:
            key_list = path.split("|")
            target_dict = return_dict_from_address(key_list, distance_metrics_store)
            # now target_dict is a linear model output dictionary in standard format
            if target_dict["is_success"] == 1:
                model_type = target_dict["model_type_str"]  # identifies if log-log, lin-log, lin-lin, or log-lin
                slope = target_dict["slope"]
                intercept = target_dict["intercept"]
                is_shifted = bool(target_dict["is_shifted"])
                # plot the linear model, reconstructed to the appropriate base model
                x_start = target_dict["x_lim"][0]
                x_end = target_dict["x_lim"][1]
                num_plot_values = 1000
                if model_type in ["lin-log", "log-log"]:
                    # transform fitting x-values back if x-axis was on log scale
                    x_start = np.exp(x_start)
                    x_end = np.exp(x_end)
                plot_x_values = np.linspace(x_start, x_end, num_plot_values)
                plot_y_values = np.zeros(num_plot_values)
                for test_val_index in range(num_plot_values):
                    plot_y_values[test_val_index] = linear_model_function(model_type=model_type,
                                                                          slope=slope, intercept=intercept,
                                                                          x_val=plot_x_values[test_val_index])
                fig = plt.figure()
                plt.plot(plot_x_values, plot_y_values, c='r')

                # check if x_val and y_val were printed
                try:
                    x_values = target_dict["x_val"]
                    y_values = target_dict["y_val"]
                    # transform both x- and y- empirical values back to original domains if necessary
                    if model_type in ["lin-log", "log-log"]:
                        x_values = np.exp(x_values)
                    if model_type in ["log-lin", "log-log"]:
                        y_values = np.exp(y_values)
                except KeyError:
                    pass
                else:
                    # if so calculate the r-squared
                    ss_res = 0.0
                    ss_tot = 0.0
                    y_actual_mean = np.mean(y_values)
                    for x_index, x_val in enumerate(x_values):
                        # observed y_val
                        y_val = y_values[x_index]
                        # then what is the fitted y_val for this?
                        y_plot_val = linear_model_function(model_type, slope, intercept, x_val)
                        ss_res += (y_val - y_plot_val) ** 2.0
                        ss_tot += (y_val - y_actual_mean) ** 2.0
                    if ss_tot == 0.0:
                        # default to zero if no variation in output - even if this means a perfect y=c line can be fit!
                        r_str = "0.0"
                    else:
                        r_squared = 1.0 - ss_res / ss_tot
                        r_str = "{0:.3g}".format(r_squared)
                    # scatter the empirical data points on the same plot
                    plt.scatter(x_values, y_values)
                    plt.legend((f'Fit: $R^{2}$ = {r_str}', 'Actual values'))
                    plt.xlabel('Input')
                    plt.ylabel('Output')
                    print_name = path.replace('|', '_')
                    if is_shifted:
                        # if the log- / -log data had to be shifted prior to fitting, that is noted here in the suffix.
                        print_name += "_shift"
                    file_path = f"{sim_path}/{step}/figures/linear_models/{print_name}.png"
                    print_and_close(fig, file_path)


def linear_model_function(model_type, slope, intercept, x_val):
    if model_type == "lin-lin":
        # linear model
        y_val = slope * x_val + intercept
    elif model_type == "log-log":
        # power law
        y_val = np.exp(intercept) * x_val ** slope
    elif model_type == "log-lin":  # (y, x)
        # exponential model
        y_val = np.exp(slope * x_val + intercept)
    elif model_type == "lin-log":  # (y, x)
        # logarithmic model
        y_val = slope * np.log(x_val) + intercept
    else:
        raise Exception("Invalid type of model underlying linear model fit.")
    return y_val


def recursive_dict_search(nested_dict, seek_contains_key, upper_key_path):
    # recursively search a nested dictionary for any sub-dictionaries containing the provided "seek key"
    paths_to_target_dict = []
    for key, value in nested_dict.items():
        if upper_key_path is not None:
            key_path = str(upper_key_path) + '|' + str(key)
        else:
            key_path = str(key)
        if key == seek_contains_key:
            paths_to_target_dict.append(upper_key_path)
        elif isinstance(value, dict):
            next_results = recursive_dict_search(value, seek_contains_key, key_path)
            for result in next_results:
                paths_to_target_dict.append(result)
    return paths_to_target_dict


def return_dict_from_address(key_list, target_dict):
    # after recursive_dict_search() provides the list of addresses, iterate over them and call this function to
    # return the ACTUAL dictionary object that you were looking for (which directly CONTAINS the target key identifier).
    for key in key_list:
        # iterate through the keys in a given path
        try:
            target_dict = target_dict[key]
        except KeyError:
            # possibly the dictionary key is an integer (habitat type number) and not a string
            try:
                target_dict = target_dict[int(key)]
            except KeyError:
                raise Exception("Key failure.")
    return target_dict

# --------------------------------------- PRODUCING DEGREE DISTRIBUTION --------------------------------------- #

def plot_degree_distribution(degree_distribution_history, degree_dist_power_law_fit_history, sim_path, step):
    # scatter the most recent degree distribution and plot the overlaid power law fit at the present step
    #
    # default to the initial system
    this_distribution = degree_distribution_history[0]
    this_power_law = degree_dist_power_law_fit_history[0]
    # sort the times and find the lowest after the given time step
    possible_steps = list(degree_distribution_history.keys())
    possible_steps.sort()
    for test_step in possible_steps:
        if test_step > step:
            this_distribution = degree_distribution_history[test_step]  # list of frequencies starting at k=0
            this_power_law = degree_dist_power_law_fit_history[test_step]  # [fit_success, a, b, start_x, covariance]
            break
    x_dist_data = np.linspace(0, len(this_distribution) - 1, len(this_distribution))
    y_dist_data = np.asarray(this_distribution)
    # scatter plot
    fig = plt.figure()
    plt.scatter(x_dist_data, y_dist_data)
    # do we plot the fitted power law?
    if this_power_law[0] == 1:  # was success
        a = this_power_law[1]
        b = this_power_law[2]
        start_x = this_power_law[3]
        min_y = np.min(this_distribution)
        non_zero_distribution = this_distribution[start_x:]
        x_plot_data = np.linspace(start_x, start_x + len(non_zero_distribution) - 1, 100000)
        y_plot_data = a * (x_plot_data - start_x + 0.1) ** b - np.abs(min_y) - 0.1
        plt.plot(x_plot_data, y_plot_data, c='r')
        # determine r-squared
        ss_res = 0.0
        ss_tot = 0.0
        y_actual_mean = np.mean(non_zero_distribution)
        for y_index, y_val in enumerate(non_zero_distribution):
            x_val = x_dist_data[start_x + y_index]  # what x_value does this true y_val correspond to?
            # then what is the fitted y_val for this?
            y_plot_val = a * (x_val - start_x + 0.1) ** b - np.abs(min_y) - 0.1
            ss_res += (y_val - y_plot_val) ** 2.0
            ss_tot += (y_val - y_actual_mean) ** 2.0
        if ss_tot == 0.0:
            r_str = "0.0"
        else:
            r_squared = 1.0 - ss_res / ss_tot
            r_str = "{0:.3g}".format(r_squared)
        plt.legend(('Actual values', f'Fit: $R^{2}$ = {r_str}'))
    plt.xlabel('Degree')
    plt.ylabel('Frequency')
    file_path = f"{sim_path}/{step}/figures/network_degree_distribution.png"
    print_and_close(fig, file_path)


# --------------------------------------- PRODUCING SPECIES-AREA CURVES (SAR) --------------------------------------- #

def biodiversity_analysis(patch_list, species_set, parameters, sim_path, step):
    #
    # This conducts a slightly more comprehensive SAR investigation, for plotting, than the version included as part of
    # the network complexity analysis in system_state.complexity_analysis() (that is,
    # species_diversity += system_state.count_diversity()).

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
    file_path = f"{sim_path}/{step}/figures/species_area_curve.png"
    print_and_close(fig, file_path)
