import os.path
import numpy as np
import json
import pickle
import sys
import random

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


def format_dictionary_to_JSON_string(input_string, is_final_item, is_indenting):
    # This takes a string obtained from passing a nested dictionary to json.dumps and recreates the indented structure
    # of the original dictionary in the string, suitable for printing to screen or file as a human-readable JSON that
    # can also then be read back from the file by passing the raw text to json.loads() to convert to a dictionary.
    #
    # The purpose of this is to:
    # i) print highly-nested dictionary structures better than permitted by pprint()
    # ii) allow writing and reading of dictionaries as JSON data WITHOUT EVER PRINTING TO AN ACTUAL .JSON FILE, which
    #       makes it easier to include in the single STOUT > .txt printing in Condor.
    #
    modified_string = input_string.replace('},', '}')
    modified_string = modified_string.replace('} ', '}')
    open_list = modified_string.split('{')
    indent_count = 0
    output_string = ''
    if is_indenting:
        indent_str = '    '
    else:
        indent_str = ''

    # first split and iterate by { which increases the nesting
    for element_index, element in enumerate(open_list):
        base_string = ''
        indent_count += 1
        if '}' in element:
            closed_sublist = element.split('}')

            # within them, split and iterate by } which decreases nesting
            for sub_index, sub_element in enumerate(closed_sublist[0:len(closed_sublist) - 1]):
                indent_count -= 1

                if len(closed_sublist[sub_index]) > 0 and len(closed_sublist[sub_index + 1]) > 0:
                    base_string += sub_element + '\n' + indent_str * indent_count + '},\n' \
                                   + indent_str * (indent_count - 1)
                elif len(closed_sublist[sub_index]) > 0 and len(closed_sublist[sub_index + 1]) == 0:
                    base_string += sub_element + '\n' + indent_str * indent_count + '}\n' \
                                   + indent_str * (indent_count - 1)
                elif len(closed_sublist[sub_index]) == 0 and len(closed_sublist[sub_index + 1]) > 0:
                    base_string += sub_element + '},\n' + indent_str * (indent_count - 1)
                elif len(closed_sublist[sub_index]) == 0 and len(closed_sublist[sub_index + 1]) == 0:
                    if element_index == len(open_list) - 1 and sub_index == len(closed_sublist) - 2:
                        # no new-line at the very end of the dictionary
                        if is_final_item:
                            # only add a final comma if required
                            base_string += sub_element + '}'
                        else:
                            base_string += sub_element + '},'
                    else:
                        base_string += sub_element + '}\n' + indent_str * (indent_count - 1)

            if element_index == len(open_list) - 1:
                # for the final element of the whole, do not add an additional { after the final sub-element
                base_string += closed_sublist[-1]
            else:
                base_string += closed_sublist[-1] + '{\n' + indent_str * indent_count
        else:
            # only one "sub-element" i.e. { ... { with no } in between
            if element_index == len(open_list) - 1:
                # for the final element, do not add an additional {
                base_string = element
            else:
                base_string = element + '{\n' + indent_str * indent_count
        output_string += base_string
    return output_string


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
