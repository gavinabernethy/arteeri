import numpy as np
from copy import deepcopy
from population_dynamics import build_interacting_populations_list, build_actual_dispersal_targets, \
    reset_dispersal_values, pre_dispersal_of_local_population


# VISUALISATION: be aware that perturbations will not show up in the species dispersal time series, and that a full
# iteration will also occur before the recording of data, so that isolated impact of the perturbation alone will never
# be visible. You will instead observe the normally-indexed time-series, and note that part of the TS change was due
# to perturbation during the time-step, in addition to subsequent reproduction, dispersal, etc.
#
# To include:
#
#   - Program the recording, analysis and visualisation of the impact of perturbations of all kinds.
#       - save the simulation in a reloadable state prior to the perturbation
#       - [X] visualisations before and after the perturbation
#       - [X] visualisation of the change
#       - how does the change compare to average change per iteration prior to the perturbation?
#
#   - Spatial network connectivity:
#       - Adding wildlife corridors or removing existing connections.
#   - Patch quality (affecting all patches globally):
#       - Seasonal variation.
#       - Stochastic variation due to bad winters.
#   - Habitat types (pass in a list of patches to alter):
#       - Specified changes in the nature recovery plan.
#       - Climate change and gradual gradient change to habitat types.
#       - Periodic variation due to crop rotation (periodic and random timing; fixed or random cycle of habitat types).
#       - Flooding (sudden, rapid change to habitat types of set duration).
# 	- Cold snap impacting reproductive rates, mortality, carrying capacity and patch quality (e.g. could reduce all
# 	     patch quality, but freeze waterways making them traversable to terrestrial species, for a set number of days)
# 	     - This would be implemented as two sets of two perturbations, the first time reducing patch quality and
# 	     changing to the ``frozen'' habitat type which would need to already exist in the set with appropriate traversal
# 	     scores for each species. The second would then reverse both of these effects to the original quality and
# 	     habitat types.
#   - Note that change_patch_parameter() can apply to either quality (thus, reproductive rate) or size for implementing
#        either patch degradation or restoration (impacting carrying capacity) as parts of the same connected patch
#        become less or more usable by residents.

#
# ---------------------------------------------- RESERVE CONSTRUCTION ---------------------------------------------- #
#

def reserve_construction(system_state, cluster_spec, clusters_must_be_separated):
    # run this at the beginning of the simulation if reserves are to be created but are not specified
    eligible_patch_nums = list(set(system_state.current_patch_list))
    reserve_list = cluster_builder(system_state=system_state,
                                   cluster_spec=cluster_spec,
                                   eligible_patch_nums=eligible_patch_nums,
                                   clusters_must_be_separated=clusters_must_be_separated)
    return reserve_list


#
# ------------------------------------------------ MAIN PERTURBATION ------------------------------------------------ #
#

# Implementing the perturbation
def perturbation(system_state, parameters, pert_paras):
    print(f" ...Step {system_state.step} - implementing "
          f"{pert_paras['perturbation_type']}: {pert_paras['perturbation_subtype']}.")

    # prepare defaults to pass in
    try:
        patch_list_overwrite = pert_paras["patch_list_overwrite"]
    except KeyError:
        patch_list_overwrite = None
    try:
        is_pairs = pert_paras["is_pairs"]
    except KeyError:
        is_pairs = False
    try:
        habitat_nums_to_change_to = pert_paras["habitat_nums_to_change_to"]
    except KeyError:
        habitat_nums_to_change_to = None
    try:
        parameter_change = pert_paras["parameter_change"]
    except KeyError:
        parameter_change = None
    try:
        parameter_change_type = pert_paras["parameter_change_type"]
    except KeyError:
        parameter_change_type = None
    try:
        parameter_change_attr = pert_paras["parameter_change_attr"]
    except KeyError:
        parameter_change_attr = None
    try:
        adjacency_change = pert_paras["adjacency_change"]
    except KeyError:
        adjacency_change = None
    try:
        is_reserves_overwrite = pert_paras["is_reserves_overwrite"]  # if True then ignore reserves
    except KeyError:
        is_reserves_overwrite = False
    try:
        clusters_must_be_separated = pert_paras["clusters_must_be_separated"]
    except KeyError:
        clusters_must_be_separated = False
    try:
        proximity_to_previous = pert_paras["proximity_to_previous"]
    except KeyError:
        proximity_to_previous = None
    try:
        prev_weighting = pert_paras["prev_weighting"]
    except KeyError:
        prev_weighting = None
    try:
        all_weighting = pert_paras["all_weighting"]
    except KeyError:
        all_weighting = None
    try:
        rebuild_all_patches = pert_paras["rebuild_all_patches"]
    except KeyError:
        rebuild_all_patches = False

    # population perturbation only
    try:
        probability = pert_paras["probability"]
    except KeyError:
        probability = None
    try:
        fraction_of_population = pert_paras["fraction_of_population"]
    except KeyError:
        fraction_of_population = None
    try:
        species_affected = pert_paras["species_affected"]  # if None (list of names of species) then default is all
    except KeyError:
        species_affected = None
    try:
        all_patches_habitats_affected = pert_paras["all_patches_habitats_affected"]
    except KeyError:
        all_patches_habitats_affected = False
    try:
        patches_affected = pert_paras["patches_affected"]  # if None (list of patch numbers) default all except reserve
    except KeyError:
        patches_affected = None

    if pert_paras["perturbation_type"] == "patch_perturbation":

        patch_perturbation(system_state=system_state,
                           parameters=parameters,
                           perturbation_subtype=pert_paras["perturbation_subtype"],
                           patch_list_overwrite=pert_paras["patch_list_overwrite"],
                           patches_affected=patches_affected,
                           is_pairs=is_pairs,
                           habitat_nums_to_change_to=habitat_nums_to_change_to,
                           parameter_change=parameter_change,
                           parameter_change_type=parameter_change_type,
                           parameter_change_attr=parameter_change_attr,
                           adjacency_change=adjacency_change,
                           is_reserves_overwrite=is_reserves_overwrite,
                           clusters_must_be_separated=clusters_must_be_separated,
                           proximity_to_previous=proximity_to_previous,
                           prev_weighting=prev_weighting,
                           all_weighting=all_weighting,
                           rebuild_all_patches=rebuild_all_patches,
                           )

    elif pert_paras["perturbation_type"] == "population_perturbation":

        population_perturbation(system_state=system_state,
                                parameters=parameters,
                                perturbation_subtype=pert_paras["perturbation_subtype"],
                                # 'dispersal' or 'extinction' or 'displacement'
                                probability=probability,
                                fraction_of_population=fraction_of_population,
                                patch_list_overwrite=patch_list_overwrite,
                                is_reserves_overwrite=is_reserves_overwrite,
                                species_affected=species_affected,
                                all_patches_habitats_affected=all_patches_habitats_affected,
                                patches_affected=patches_affected,
                                clusters_must_be_separated=clusters_must_be_separated,
                                proximity_to_previous=proximity_to_previous,
                                prev_weighting=prev_weighting,
                                all_weighting=all_weighting,
                                )
    else:
        # Note that habitat types (particularly species-specific feeding and traversal scores) should NEVER be altered.
        raise 'Perturbation main type not recognised. Expected "population_perturbation" or "patch_perturbation".'

    system_state.increment_num_perturbations()


#
# ----------------------------------------- PERTURBATION TYPES: POPULATION ----------------------------------------- #
#
# Perturb local populations of species.
# CAREFUL - recall that making any changes e.g. to species properties of a local population object, will change that
# species reference and thus both the species in principle and all local population references to this species!
def population_perturbation(
        system_state,
        parameters,
        perturbation_subtype,  # 'dispersal' or 'extinction' or 'displacement'
        probability,
        fraction_of_population,
        patch_list_overwrite=None,  # if a nom-empty list is passed, we ignore everything else!
        is_reserves_overwrite=False,  # if True then ignore reserves
        species_affected=None,  # if None (list of names of species) then default is all
        patches_affected=None,  # list of dictionaries of patch clusters
        all_patches_habitats_affected=None,
        clusters_must_be_separated=False,
        proximity_to_previous: int = 0,
        prev_weighting=None,
        all_weighting=None,
):
    if perturbation_subtype in ['dispersal', 'displacement']:
        reset_dispersal_values(patch_list=system_state.patch_list)
    if perturbation_subtype not in ['dispersal', 'extinction', 'displacement']:
        raise 'Unexpected local population perturbation type. Expected "dispersal" or "extinction" or "displacement".'

    if species_affected is None:
        # list of the names of all species
        species_affected = list(system_state.species_set["dict"].keys())

    # if a list of specific patch numbers (rather than a list of dictionaries) is specified as a priority, we enact
    # them and ignore all other checks and features.
    # (Note that we can also pass a list of patches as one cluster ("patch_numbers") in the regular sense)
    if patch_list_overwrite is not None and len(patch_list_overwrite) > 0:
        # specific patches are specified
        perturbation_list = [patch_list_overwrite]  # list of a single 'cluster' list
    elif patches_affected is not None and len(patches_affected) > 0:
        # cluster building process
        perturbation_list = perturbation_patch_selection(system_state=system_state,
                                                         is_reserves_overwrite=is_reserves_overwrite,
                                                         patches_affected=patches_affected,
                                                         clusters_must_be_separated=clusters_must_be_separated,
                                                         proximity_to_previous=proximity_to_previous,
                                                         prev_weighting=prev_weighting,
                                                         all_weighting=all_weighting)
    else:
        # Can choose all but then remove reserves and unwanted habitat types.
        # Note that all_patches_habitat_types only applies to this case where no patches are explicitly specified (in
        # which case habitat type is not considered) and no clusters are generated (since habitat restrictions are
        # specified PER CLUSTER in that case).
        #
        # Reserves first:
        if not is_reserves_overwrite:
            # remove the reserves from final consideration
            set_of_reserve_patches = set([reserve for cluster in system_state.reserve_list for reserve in cluster])
            eligible_patch_nums = list(set(deepcopy(system_state.current_patch_list)) - set_of_reserve_patches)
        else:
            eligible_patch_nums = deepcopy(system_state.current_patch_list)
        # Then habitats:
        if all_patches_habitats_affected is None:
            perturbation_list = [eligible_patch_nums]
        else:
            single_cluster = []
            for patch_num in eligible_patch_nums:
                if system_state.patch_list[patch_num].habitat_type_num in all_patches_habitats_affected:
                    single_cluster.append(patch_num)
            perturbation_list = [single_cluster]

    # record history of the perturbations enacted
    system_state.perturbation_history[system_state.step] = perturbation_list

    # now flatten list for implementation
    flat_list_to_perturb = [patch for cluster in perturbation_list for patch in cluster]

    # ---------------- IMPLEMENT PERTURBATION ---------------- #
    #
    for patch_num in flat_list_to_perturb:
        for local_pop in system_state.patch_list[patch_num].local_populations.values():
            if local_pop.name in species_affected:
                # probabilistic implementation
                if np.random.binomial(1, max(0.0, min(1.0, probability))):
                    # implement!
                    system_state.patch_list[patch_num].increment_perturbation_count()  # count (per species) perturb.
                    system_state.patch_list[patch_num].perturbation_history_list.append(
                        [system_state.step, perturbation_subtype, local_pop.name])
                    if perturbation_subtype == 'extinction':
                        old_population = deepcopy(local_pop.population)
                        temp_pop = local_pop.population * (1.0 - fraction_of_population)
                        if temp_pop < local_pop.species.minimum_population_size:
                            local_pop.population = 0.0
                        else:
                            local_pop.population = temp_pop
                        if old_population != local_pop.population:
                            system_state.patch_list[patch_num].increment_meaningful_perturbation_count()
                    elif perturbation_subtype == 'dispersal':
                        # uses (and may be limited by) CURRENT species-specific dispersal scores and mechanisms
                        pre_dispersal_of_local_population(
                            patch_list=system_state.patch_list,
                            parameters=parameters,
                            local_pop=local_pop,
                            specified_fraction=fraction_of_population,
                            is_dispersal_override=True,
                        )
                    elif perturbation_subtype == 'displacement':
                        # force displacement to adjacent neighbours regardless of species dispersal parameters
                        local_pop.population_leave = fraction_of_population * local_pop.population
                        for patch_to_num in system_state.patch_list[patch_num].set_of_adjacent_patches:
                            if patch_to_num != patch_num:
                                species_find = system_state.patch_list[
                                    patch_to_num].local_populations[local_pop.name]
                                species_find.population_enter += local_pop.population_leave / (
                                        system_state.patch_list[patch_num].degree - 1)
                    else:
                        raise "Error."
    if perturbation_subtype in ['dispersal', 'displacement']:
        # For population perturbation types involving emigration, the movements are conducted simultaneously.
        for patch in system_state.patch_list:
            for local_pop in patch.local_populations.values():
                old_population = deepcopy(local_pop.population)
                temp_pop = local_pop.population + local_pop.population_enter - local_pop.population_leave
                if temp_pop < local_pop.species.minimum_population_size:
                    temp_pop = 0.0
                local_pop.population = max(0.0, float(temp_pop))
                if old_population != local_pop.population:
                    # for dispersal / displacement, count every time any population has a NET change - including from
                    # a neighbouring patch being perturbed, while perturbations that cancel each other exactly will
                    # NOT be counted in the "meaningful" metric
                    patch.increment_meaningful_perturbation_count()


#
# -------------------------------------------- PERTURBATION TYPES: PATCH -------------------------------------------- #
#

def patch_perturbation(
        system_state,
        parameters,
        perturbation_subtype,
        patch_list_overwrite=None,
        patches_affected=None,
        is_pairs=False,
        habitat_nums_to_change_to=None,
        parameter_change=None,
        parameter_change_type=None,
        parameter_change_attr=None,
        adjacency_change=None,
        is_reserves_overwrite=False,
        clusters_must_be_separated=False,
        proximity_to_previous=None,
        prev_weighting=None,
        all_weighting=None,
        rebuild_all_patches=False,
):
    # Note that habitat properties (in particular species-specific habitat feeding and traversal scores)
    # SHOULD NEVER BE CHANGED even in a perturbation.
    # Instead declare an additional habitat type and change patch to that (rather than altering the habitat itself).
    #
    # establish the dictionary of functions
    function_dictionary = {
        "change_habitat": change_patch_to_habitat,
        "change_parameter": change_patch_parameter,
        "remove_patch": remove_patch,
        "change_adjacency": change_patch_adjacency,
    }

    # Is the patch selection process to be manually overwritten with a list of patches or pairs of patches?
    if patch_list_overwrite is not None and len(patch_list_overwrite) > 0:
        if type(patch_list_overwrite[0]) == int and not is_pairs:
            # list of a single 'cluster' list
            patches_to_alter = patch_list_overwrite
            patch_pairs_to_change = None
            perturbation_list = [patch_list_overwrite]
        elif type(patch_list_overwrite[0]) in [list, tuple] and is_pairs:
            # list of patch pairs
            patches_to_alter = None
            patch_pairs_to_change = patch_list_overwrite
            perturbation_list = [[x for y in patch_pairs_to_change for x in y]]
        else:
            raise "Error"
    elif patches_affected is not None and len(patches_affected) > 0:
        # go through the full selection process
        perturbation_list = perturbation_patch_selection(system_state=system_state,
                                                         is_reserves_overwrite=is_reserves_overwrite,
                                                         patches_affected=patches_affected,
                                                         clusters_must_be_separated=clusters_must_be_separated,
                                                         proximity_to_previous=proximity_to_previous,
                                                         prev_weighting=prev_weighting,
                                                         all_weighting=all_weighting)

        # then assign these selected patches or pairs of patches
        #
        # perturbation_list should be a list of lists of clusters. The clustering is for how they are chosen and
        # recorded, but is not actually involved in the implementation. Thus we move to the following forms before
        # passing to the perturbation functions:
        # - patches_to_alter should be a single list of the values of all affected patches (i.e. from all clusters)
        # - patch_pairs_to_change should be a list of lists of pairs
        if is_pairs:
            patches_to_alter = None
            patch_pairs_to_change = [(perturbation_list[2 * i], perturbation_list[2 * i + 1]) for i in
                                     range(int(len(perturbation_list) / 2))]
        else:
            patches_to_alter = [patch for cluster in perturbation_list for patch in cluster]
            patch_pairs_to_change = None
    else:
        raise "No patches provided to patch perturbation."

    # record history of the perturbations enacted
    system_state.perturbation_history[system_state.step] = perturbation_list
    # assign a colour code to all patches in each cluster
    for cluster in perturbation_list:
        draw_cluster_code = np.random.uniform(low=0.2, high=1)
        for patch_num in cluster:
            system_state.patch_list[patch_num].latest_perturbation_code = draw_cluster_code
            system_state.patch_list[patch_num].latest_perturbation_code_history[system_state.step] = draw_cluster_code

    # replace any random habitats or patch adjacency required:
    if patches_to_alter is not None and len(patches_to_alter) > 0 and not is_pairs:
        habitat_nums_to_change_to = set_random_choices(list_to_check=habitat_nums_to_change_to,
                                                       replacement_type="list",
                                                       replacement_possibilities=[x for x in parameters[
                                                           "main_para"]["HABITAT_TYPES"]])

    # Check specified patch pairs and randomly draw if required subject to reserves
    if patch_pairs_to_change is not None and len(patch_pairs_to_change) > 0 and is_pairs:
        adjacency_change = set_random_choices(list_to_check=adjacency_change,
                                              replacement_type="distribution",
                                              replacement_possibilities=[0, 1])

    # ------------- IMPLEMENTING THE PERTURBATION ------------- #
    #
    # call the function
    function_dictionary[perturbation_subtype](
        system_state=system_state,
        parameters=parameters,
        patches_to_change=patches_to_alter,
        habitat_nums_to_change_to=habitat_nums_to_change_to,
        parameter_change=parameter_change,
        parameter_change_type=parameter_change_type,
        parameter_change_attr=parameter_change_attr,
        patch_pairs_to_change=patch_pairs_to_change,
        adjacency_change=adjacency_change,
    )

    # build the list of all patches that were directly perturbed
    altered_patch_numbers = []
    if patches_to_alter is not None and len(patches_to_alter) > 0:
        altered_patch_numbers = list(set(altered_patch_numbers + patches_to_alter))
    if patch_pairs_to_change is not None and len(patch_pairs_to_change) > 0:
        altered_patch_numbers = list(set(altered_patch_numbers + [y for x in patch_pairs_to_change for y in x]))

    # count the patch perturbations (whether or not they resulted in a meaningful change)
    for patch_num in altered_patch_numbers:
        system_state.patch_list[patch_num].increment_perturbation_count()
        system_state.patch_list[patch_num].perturbation_history_list.append([system_state.step, perturbation_subtype])
    # pass only the non-duplicated, numbers of patches that have been changed, to rebuild purely internal properties
    reset_local_population_attributes(patch_list=system_state.patch_list, altered_patch_numbers=altered_patch_numbers)

    # re-determining the species movement scores for both foraging and dispersal:
    if rebuild_all_patches:
        # (slow for large networks) rebuild for all patches rather than just the estimate of those closely impacted
        system_state.build_all_patches_species_paths_and_adjacency(parameters=parameters)
    else:
        # start with the patches literally changed
        likely_affected_patches = altered_patch_numbers
        # what patches had routes going through them?
        for patch_to_check in system_state.patch_list:
            for changed_patch_num in altered_patch_numbers:
                if changed_patch_num in patch_to_check.stepping_stone_list:
                    likely_affected_patches = list(set(likely_affected_patches + [patch_to_check.number]))
        # now check those who are NOW (after the perturbation) the N-th degree neighbours (note that
        # because the diagonal of the patch_adjacency_matrix are all 1, then an entry being non-zero indicates that
        # we are identifying UP TO Nth Degree Neighbours (as we could also have an M < N degree neighbour, plus
        # revisiting the same final node N-M times. So this does not tell us shortest paths but DOES identify the
        # existence of any connections/reachability, as desired here.
        #
        # The point is that this catches anyone (anywhere) who MAY now have a route that uses this patch.
        path_matrix = deepcopy(system_state.patch_adjacency_matrix)
        for _ in range(parameters["main_para"]["ASSUMED_MAX_PATH_LENGTH"] - 1):
            # matrix multiplication
            path_matrix = np.matmul(path_matrix, path_matrix)
            # this causes number of walks to grow exponentially so normalise - we only need to know when non-zero
            for patch_x in range(len(system_state.patch_list)):
                for patch_y in range(len(system_state.patch_list)):
                    if path_matrix[patch_x, patch_y] != 0.0:
                        path_matrix[patch_x, patch_y] = 1.0
        # no point checking patch numbers we already consider affected
        remaining_to_check = list(set([x for x in range(len(system_state.patch_list))]) - set(likely_affected_patches))
        for patch_num_to_check in remaining_to_check:
            for changed_patch_num in altered_patch_numbers:
                if path_matrix[patch_num_to_check, changed_patch_num] > 0.0 \
                        or path_matrix[changed_patch_num, patch_num_to_check] > 0.0:
                    likely_affected_patches = list(set(likely_affected_patches + [patch_num_to_check]))
        # pass list to rebuild function
        system_state.build_all_patches_species_paths_and_adjacency(parameters=parameters,
                                                                   specified_patch_list=likely_affected_patches)

    is_nonlocal_foraging = parameters["pop_dyn_para"]["IS_NONLOCAL_FORAGING_PERMITTED"]
    is_local_foraging_ensured = parameters["pop_dyn_para"]["IS_LOCAL_FORAGING_ENSURED"]
    build_interacting_populations_list(
        patch_list=system_state.patch_list,
        species_list=system_state.species_set["list"],
        is_nonlocal_foraging=is_nonlocal_foraging,
        is_local_foraging_ensured=is_local_foraging_ensured,
        time=system_state.time,
    )
    is_dispersal = parameters["pop_dyn_para"]["IS_DISPERSAL_PERMITTED"]
    build_actual_dispersal_targets(
        patch_list=system_state.patch_list,
        species_list=system_state.species_set["list"],
        is_dispersal=is_dispersal,
        time=system_state.time,
    )


# ------------------------------------------- PATCH PERTURBATION SUBTYPES ------------------------------------------- #

# Change habitat type
def change_patch_to_habitat(system_state,
                            parameters,
                            patches_to_change,
                            habitat_nums_to_change_to,
                            parameter_change=None,
                            parameter_change_type=None,
                            parameter_change_attr=None,
                            patch_pairs_to_change=None,
                            adjacency_change=None,
                            ):
    # pass in a list of the numbers of the patches to change, and the habitats (as a name string) to change them to
    # remember that the patch numbers start at 0!
    habitat_types = parameters["main_para"]["HABITAT_TYPES"]
    if len(patches_to_change) != len(habitat_nums_to_change_to):
        raise "List of patches to alter is not the same length as the list of new habitats."
    else:
        for num_in_temp_list, patch_number in enumerate(patches_to_change):
            new_habitat_type_num = habitat_nums_to_change_to[num_in_temp_list]
            # count if it actually makes a difference
            if system_state.patch_list[patch_number].habitat_type_num != new_habitat_type_num:
                system_state.patch_list[patch_number].increment_meaningful_perturbation_count()
                system_state.patch_list[patch_number].habitat_type_num = new_habitat_type_num
                system_state.patch_list[patch_number].habitat_type = habitat_types[new_habitat_type_num]
                system_state.update_patch_habitat_based_properties(system_state.patch_list[patch_number])
                # record history
                system_state.patch_list[patch_number].habitat_history[system_state.step] = new_habitat_type_num
    # record changes at a system_state level in the history
    system_state.update_habitat_distributions_history()


# Change patch floating point parameter value in range [0, 1] (quality or size)
def change_patch_parameter(system_state,
                           parameters,
                           patches_to_change,
                           parameter_change,
                           parameter_change_attr,
                           parameter_change_type,
                           habitat_nums_to_change_to=None,
                           patch_pairs_to_change=None,
                           adjacency_change=None,
                           ):
    if parameter_change_attr not in ["size", "quality"]:
        # ensure that either "size" or "quality" is specified:
        raise Exception(f"Perturbation parameter {parameter_change_attr} should be either 'size' or 'quality'.")
    elif type(parameter_change) is list and len(patches_to_change) != len(parameter_change):
        # pass in a list of the numbers of the patches to change, and the change required to their quality/size
        # either as a scaling factor (is_relative) or as a new specified value to adopt (is_absolute)
        raise Exception("List of patches to alter is not the same length as the list of specified alterations.")
    elif parameter_change_type not in ["absolute", "relative_multiply", "relative_add"]:
        raise Exception("Exactly one of these must be the provided key.")
    else:
        for num_in_temp_list, patch_number in enumerate(patches_to_change):
            old_parameter_value = deepcopy(getattr(system_state.patch_list[patch_number], parameter_change_attr))
            # determine the value of the change - either a single value for all, or a specified patch-length list
            if type(parameter_change) is list:
                this_parameter_change = parameter_change[num_in_temp_list]
            elif type(parameter_change) in [float, int]:
                this_parameter_change = parameter_change
            else:
                raise "Value change needs to be either a single value or a list of equal length to" \
                      " the number of patches to change."
            # what type of change is it?
            if parameter_change_type == "relative_multiply":
                new_parameter_value = getattr(system_state.patch_list[patch_number],
                                              parameter_change_attr) * this_parameter_change
            elif parameter_change_type == "relative_add":
                new_parameter_value = getattr(system_state.patch_list[patch_number],
                                              parameter_change_attr) + this_parameter_change
            elif parameter_change_type == "absolute":
                new_parameter_value = this_parameter_change
            else:
                raise "Parameter value error."
            # now set the new value
            new_parameter_value = min(max(new_parameter_value, 0.0), 1.0)
            setattr(system_state.patch_list[patch_number], parameter_change_attr, new_parameter_value)
            # record history for local patches
            if parameter_change_attr == "size":
                system_state.patch_list[patch_number].size_history[system_state.step] = new_parameter_value
            elif parameter_change_attr == "quality":
                system_state.patch_list[patch_number].quality_history[system_state.step] = new_parameter_value
            # record if real change
            if new_parameter_value != old_parameter_value:
                system_state.patch_list[patch_number].increment_meaningful_perturbation_count()
    # record network-level values (mean and s.d. of values for all extant patches)
    if parameter_change_attr == "size":
        system_state.update_size_history()
    elif parameter_change_attr == "quality":
        system_state.update_quality_history()


# Remove a patch
def remove_patch(system_state,
                 parameters,
                 patches_to_change,
                 patch_pairs_to_change=None,
                 adjacency_change=None,
                 parameter_change=None,
                 parameter_change_type=None,
                 parameter_change_attr=None,
                 habitat_nums_to_change_to=None,
                 ):
    patches_to_remove = patches_to_change
    # Set populations and adjacency scores to zero for all species but do not delete it from the patch list so that:
    # (a) it can be re-connected later,
    # (b) no need to resize the arrays,
    # (c) the patch numbers will remain constant in the patch_list. THIS IS ESSENTIAL!
    for patch_number in patches_to_remove:
        # check not already "removed"
        if patch_number not in system_state.current_patch_list:
            # count the change and proceed
            system_state.patch_list[patch_number].increment_meaningful_perturbation_count()
            # set all local populations of this patch to zero
            for local_population in system_state.patch_list[patch_number].local_populations.values():
                local_population.population = 0.0
            for species_name in system_state.patch_list[patch_number].species_movement_scores:
                system_state.patch_list[patch_number].species_movement_scores[species_name] = {
                    x: {"best": (float('inf'), float('inf'), 0.0)} for x in range(len(system_state.patch_list))}

            # set degree to zero and record new values
            system_state.patch_list[patch_number].degree = 0
            system_state.patch_list[patch_number].degree_history[system_state.step] = 0

            # set self-adjacency to zero, record history of set of adjacent patches to zero
            system_state.patch_adjacency_matrix[patch_number, patch_number] = 0
            system_state.patch_list[patch_number].set_of_adjacent_patches = set({})
            system_state.patch_list[patch_number].set_of_adjacent_patches_history[system_state.step] = []

            # now look at other patches - need to set the corresponding row and column in the adjacency matrix to zero
            # and remove the patch number from the list of currently accessible patches:
            for patch in system_state.patch_list:
                if patch.number != patch_number:
                    # Now look at all connected patches (but not the same patch again)
                    if system_state.patch_adjacency_matrix[patch_number, patch.number] != 0.0 or \
                            system_state.patch_adjacency_matrix[patch.number, patch_number] != 0.0:
                        # record changes to degree
                        patch.degree = max(0, patch.degree - 1)
                        patch.degree_history[system_state.step] = patch.degree

                        # record changes to sets of adjacent patches
                        patch.set_of_adjacent_patches.remove(patch_number)
                        patch.set_of_adjacent_patches_history[system_state.step] = list(patch.set_of_adjacent_patches)

                        # record adjacency change and zero the corresponding patch adjacency matrix entries
                        patch.adjacency_history_list.append([system_state.step, patch_number, 0])
                        system_state.patch_list[patch_number].adjacency_history_list.append(
                            [system_state.step, patch.number, 0])
                        system_state.patch_adjacency_matrix[patch_number, patch.number] = 0
                        system_state.patch_adjacency_matrix[patch.number, patch_number] = 0

            # only do this after reducing the degree of connected patches
            system_state.current_patch_list.remove(patch_number)
            system_state.patch_list[patch_number].removal_history[system_state.step] = 'removed'
    # after removing a patch, update all patch centrality and record history of average properties for CURRENT patches
    system_state.update_centrality_history(parameters=parameters)
    system_state.update_habitat_distributions_history()
    system_state.update_quality_history()
    system_state.update_degree_history()


# Add or remove corridors between specified patches
def change_patch_adjacency(system_state,
                           parameters,
                           patch_pairs_to_change,
                           adjacency_change,
                           patches_to_change=None,
                           parameter_change=None,
                           parameter_change_type=None,
                           parameter_change_attr=None,
                           habitat_nums_to_change_to=None,
                           ):
    # patch_pairs_to_change should be a N-length list of 2x1 tuples containing the pairs of patches
    # adjacency_change should be a length N list of new adjacency values from the set {0,1}
    if len(patch_pairs_to_change) != len(adjacency_change):
        raise "List of patch pairs to alter is not the same length as the list of adjacency alterations."
    else:
        for pair_num, pair in enumerate(patch_pairs_to_change):
            # store previous values
            old_adjacency = [deepcopy(system_state.patch_adjacency_matrix[pair[0], pair[1]]),
                             deepcopy(system_state.patch_adjacency_matrix[pair[1], pair[0]])]
            # update degree if necessary
            if adjacency_change[pair_num] == 0.0:
                if system_state.patch_adjacency_matrix[pair[0], pair[1]] != 0.0 or \
                        system_state.patch_adjacency_matrix[pair[1], pair[0]] != 0.0:
                    # were adjacent (in at least one direction), now not
                    system_state.patch_list[pair[0]].degree = max(0, system_state.patch_list[pair[0]].degree - 1)
                    system_state.patch_list[pair[1]].degree = max(0, system_state.patch_list[pair[1]].degree - 1)
                    system_state.patch_list[pair[0]].set_of_adjacent_patches.remove(pair[1])
                    system_state.patch_list[pair[1]].set_of_adjacent_patches.remove(pair[0])
            elif adjacency_change[pair_num] == 1.0:
                if system_state.patch_adjacency_matrix[pair[0], pair[1]] == 0.0 and \
                        system_state.patch_adjacency_matrix[pair[1], pair[0]] == 0.0:
                    # were NOT adjacent (in any direction), now are (in both directions)
                    system_state.patch_list[pair[0]].degree = min(system_state.patch_list[pair[0]].degree + 1,
                                                                  len(system_state.patch_list))
                    system_state.patch_list[pair[1]].degree = min(system_state.patch_list[pair[1]].degree + 1,
                                                                  len(system_state.patch_list))
                    system_state.patch_list[pair[0]].set_of_adjacent_patches.add(pair[1])
                    system_state.patch_list[pair[1]].set_of_adjacent_patches.add(pair[0])
            else:
                raise Exception('Adjacency matrix should only contain (or change to) values of either 0 or 1.')
            # finally replace the corresponding entry in the adjacency matrix with the new adjacency value
            system_state.patch_adjacency_matrix[pair[0], pair[1]] = adjacency_change[pair_num]
            system_state.patch_adjacency_matrix[pair[1], pair[0]] = adjacency_change[pair_num]
            # record change histories
            system_state.patch_list[pair[0]].adjacency_history_list.append(
                [system_state.step, pair[1], adjacency_change[pair_num]])
            system_state.patch_list[pair[1]].adjacency_history_list.append(
                [system_state.step, pair[0], adjacency_change[pair_num]])
            # record new history of degree, adjacency
            system_state.patch_list[pair[0]].degree_history[system_state.step] = system_state.patch_list[pair[0]].degree
            system_state.patch_list[pair[1]].degree_history[system_state.step] = system_state.patch_list[pair[1]].degree
            system_state.patch_list[pair[0]].set_of_adjacent_patches_history[system_state.step] = list(
                system_state.patch_list[pair[0]].set_of_adjacent_patches)
            system_state.patch_list[pair[1]].set_of_adjacent_patches_history[system_state.step] = list(
                system_state.patch_list[pair[1]].set_of_adjacent_patches)

            # check actual change has occurred
            if old_adjacency[0] != system_state.patch_adjacency_matrix[pair[0], pair[1]] or \
                    old_adjacency[1] != system_state.patch_adjacency_matrix[pair[1], pair[0]]:
                system_state.patch_list[pair[0]].increment_meaningful_perturbation_count()
                system_state.patch_list[pair[1]].increment_meaningful_perturbation_count()

    # after an adjacency change, update all patch centrality and habitat spatial auto-correlation
    system_state.calculate_all_patches_centrality(parameters=parameters)
    system_state.update_centrality_history(parameters=parameters)
    system_state.update_habitat_distributions_history()
    system_state.update_degree_history()


#
# ----------------------------------------------- AUXILIARY FUNCTIONS ----------------------------------------------- #
#

def set_random_choices(list_to_check, replacement_type, replacement_possibilities):
    if list_to_check is not None:
        for number, item in enumerate(list_to_check):
            if item == "r":
                if replacement_type == "list":
                    list_to_check[number] = int(np.random.choice(replacement_possibilities))
                elif replacement_type == "distribution":
                    list_to_check[number] = float(np.random.uniform(replacement_possibilities[0],
                                                                    replacement_possibilities[1]))
                else:
                    raise "Invalid type of replacement for random values."
    return list_to_check


def reset_local_population_attributes(patch_list, altered_patch_numbers):
    # This is to be certain that altered patches have their properties (built from references) updated.
    # May be unnecessary, but good to be sure!
    for patch_num in altered_patch_numbers:
        patch = patch_list[patch_num]
        for local_pop in patch.local_populations.values():
            local_pop.rebuild_patch_dependent_properties(patch=patch)


#
# ----------------------------- CLUSTER CHOICE FUNCTIONS FOR RESERVES AND PERTURBATIONS ----------------------------- #
#

def perturbation_patch_selection(system_state,
                                 is_reserves_overwrite,
                                 patches_affected,
                                 clusters_must_be_separated,
                                 proximity_to_previous,
                                 prev_weighting,
                                 all_weighting):
    # this function generates the list of patches to perturb in either a patch or population perturbation,
    # unless a specific priority set of patch numbers was passed in which case this function will not have been called.

    eligible_patch_nums = system_state.current_patch_list

    if not is_reserves_overwrite:
        # we DO remove the reserves from final consideration
        set_of_reserve_patches = set([reserve for cluster in system_state.reserve_list for reserve in cluster])
        eligible_patch_nums = list(set(eligible_patch_nums) - set_of_reserve_patches)

    # {0, 1, 2, 3} - (hard-permitted) PROXIMITY TO PREVIOUS PERTURBATION:
    # [0 = can choose same, 1 = can be adj, 2 = not adj, 3 = cannot be adj to a patch adj to a previous choice]
    #
    # remove patches adjacent to prior clusters if necessary
    # Note that this does not apply in cases of removed patches - as they are now considered no longer adjacent to any!
    if len(system_state.perturbation_history) > 0 and proximity_to_previous > 0:
        eligible_patch_nums_copy = deepcopy(eligible_patch_nums)
        # need to use a copy as removing from the same list while iterating over it does not behave properly!
        most_recently_perturbed_patch_nums = list(set([
            patch_num for cluster in system_state.perturbation_history[
                max(system_state.perturbation_history)] for patch_num in cluster]))
        for patch_num in eligible_patch_nums_copy:
            keep_checking = True
            for prev_patch in most_recently_perturbed_patch_nums:
                if not keep_checking:
                    # this is needed for the case of proximity_to_previous > 2 as there is a nested loop to escape
                    break
                if patch_num == prev_patch:
                    # remove patches that were in previous choices
                    eligible_patch_nums.remove(patch_num)
                    break
                # Check if ONE stepping stone required and applies here:
                if proximity_to_previous > 1:
                    # remove patches adjacent to those present in previous waves
                    if system_state.patch_adjacency_matrix[prev_patch, patch_num] \
                            > 0 or system_state.patch_adjacency_matrix[patch_num, prev_patch] > 0:
                        eligible_patch_nums.remove(patch_num)
                        break
                # THEN (not "else") check if TWO stepping stones are required and applies here:
                if proximity_to_previous > 2:
                    # remove other patches adjacent to those adjacent to those present in previous waves
                    for intermediate_patch_num in system_state.current_patch_list:
                        if intermediate_patch_num != patch_num:
                            if (system_state.patch_adjacency_matrix[prev_patch, intermediate_patch_num] > 0 or
                                system_state.patch_adjacency_matrix[intermediate_patch_num, prev_patch] > 0) and \
                                    (system_state.patch_adjacency_matrix[patch_num, intermediate_patch_num] > 0 or
                                     system_state.patch_adjacency_matrix[intermediate_patch_num, patch_num] > 0):
                                eligible_patch_nums.remove(patch_num)
                                keep_checking = False
                                break
        del eligible_patch_nums_copy  # clear this variable as no longer required

    # Weighting the probability distribution of the eligible choices:
    if len(system_state.perturbation_history) > 0 and (prev_weighting is not None or all_weighting is not None):
        most_recently_perturbed_patch_nums = list(set([
            patch_num for cluster in system_state.perturbation_history[
                max(system_state.perturbation_history)] for patch_num in cluster]))
        base_weighting_auto = {}
        base_weighting_history = {}
        base_weighting_auto_sum = 0.0
        base_weighting_history_sum = 0.0
        max_distance = 0.0
        min_distance = 999999999999999.0
        max_num_perturbation = 0
        min_num_perturbation = 999999999999999.0
        # prev_weighting: auto-correlation of perturbation events in space to MOST-RECENT previous ones
        # if +1 then 0% (in isolation from all_weighting) to choose a patch of minimal overall distance from most
        # recently-perturbed patches;
        # if -1 then 0% to choose a patch of maximal overall distance from most recently-perturbed patches
        #
        # all_weighting: probabilistic weighting to patches according to how many times previously chosen - do
        # they accumulate or spread out?
        # if +1 then 0% (in isolation from prev_weighting) to choose least-perturbed patch,
        # if -1 then 0% to choose a most-perturbed patch
        for patch_num in eligible_patch_nums:
            weighted_distance = 0.0
            for prev_patch_num in most_recently_perturbed_patch_nums:
                weighted_distance += np.linalg.norm(
                    system_state.patch_list[patch_num].position - system_state.patch_list[prev_patch_num].position)
            min_distance = min(min_distance, weighted_distance)
            max_distance = max(max_distance, weighted_distance)
            base_weighting_auto[patch_num] = weighted_distance
            max_num_perturbation = max(max_num_perturbation, system_state.patch_list[patch_num].num_times_perturbed)
            min_num_perturbation = min(min_num_perturbation, system_state.patch_list[patch_num].num_times_perturbed)
            base_weighting_history[patch_num] = system_state.patch_list[patch_num].num_times_perturbed

        for patch_num in eligible_patch_nums:
            # prev
            if prev_weighting is not None and len(prev_weighting) > 0:
                new_amount = weighting_amount_calc(max_val=max_distance, min_val=min_distance,
                                                   actual_val=base_weighting_auto[patch_num], weighting=prev_weighting)
            else:
                new_amount = 0.0
            base_weighting_auto[patch_num] = new_amount
            base_weighting_auto_sum += new_amount
            # all
            if all_weighting is not None and len(all_weighting) > 0:
                new_amount = weighting_amount_calc(max_val=max_num_perturbation, min_val=min_num_perturbation,
                                                   actual_val=base_weighting_history[patch_num],
                                                   weighting=all_weighting)
            else:
                new_amount = 0.0
            base_weighting_history[patch_num] = new_amount
            base_weighting_history_sum += new_amount
    else:
        base_weighting_auto = {x: 1 for x in eligible_patch_nums}
        base_weighting_auto_sum = max(1.0, len(eligible_patch_nums))
        base_weighting_history = {x: 1 for x in eligible_patch_nums}
        base_weighting_history_sum = max(1.0, len(eligible_patch_nums))

    # update combined weightings to obtain final weighting which gets passed to the 'draw'
    final_weighting = {}
    for patch_num in eligible_patch_nums:
        final_value_holder = 0.0
        if prev_weighting is not None and len(prev_weighting) > 0 and base_weighting_auto_sum > 0.0:
            # for the weighting lists, the first entry is their relative importance (to each other!)
            final_value_holder += prev_weighting[0] * base_weighting_auto[patch_num] / base_weighting_auto_sum
        if all_weighting is not None and len(all_weighting) > 0 and base_weighting_history_sum > 0.0:
            final_value_holder += all_weighting[0] * base_weighting_history[patch_num] / base_weighting_history_sum
        final_weighting[patch_num] = final_value_holder

    perturbation_list = cluster_builder(system_state=system_state,
                                        cluster_spec=patches_affected,
                                        eligible_patch_nums=eligible_patch_nums,
                                        clusters_must_be_separated=clusters_must_be_separated,
                                        probability_weighting=final_weighting, )
    return perturbation_list


def weighting_amount_calc(max_val, min_val, actual_val, weighting: list):
    # function called in perturbation_patch_selection to streamline the all_weighting and prev_weighting calculations
    #
    # the function is y = |alpha - x|^beta + gamma, where weighting = [relative_weight, alpha, beta, gamma]
    # thus for beta-order polynomial increase: alpha = 0, and beta-order polynomial decrease: alpha = 1
    # gamma is the vertical shift - can be positive to adjust the relative probability throughout so that
    # there is no longer a zero probability at the extreme end, but ordinarily we will keep this at zero.
    #
    # This is not checked for suitability until the probability distribution is passed to cluster_draw(), where errors
    # will halt the code if we have negative probabilities etc. Thus I recommend only passing weightings from the
    # following options (or None, of course!):
    # - Linear increase from (0,0) to (1,1): [0, 1, 0]
    # - Linear decrease from (0,1) to (1,0): [1, 1, 0]
    # - Constant value of 1: [?, 0, 0]
    # - Quadratic increase from (0,0) to (1,1): [0, 2, 0]
    # - Quadratic decrease from (1,1) to (0,0): [1, 2, 0]
    #
    if max_val - min_val > 0.0:
        normalised_val = (actual_val - min_val) / (max_val - min_val)
        new_amount = np.abs(weighting[1] - normalised_val) ** weighting[2] + weighting[3]
    else:
        new_amount = 1.0
    return new_amount


def cluster_builder(system_state, cluster_spec, eligible_patch_nums, clusters_must_be_separated,
                    probability_weighting=None):
    #
    # This is used both by perturbation_patch_selection() AND by reserve_construction()
    #
    cluster_list = []
    for cluster_num, cluster in enumerate(cluster_spec):
        current_cluster = []

        if cluster["num_patches"] >= 1:

            # remove any ineligible habitats
            actual_patch_nums = []
            if cluster["habitat_nums_permitted"] is None:
                actual_patch_nums = deepcopy(eligible_patch_nums)  # just to be sure we do not backtrace changes
            else:
                for patch_num in eligible_patch_nums:
                    if system_state.patch_list[patch_num].habitat_type_num in cluster["habitat_nums_permitted"]:
                        actual_patch_nums.append(patch_num)

            # remove patches adjacent to prior clusters if necessary (i.e. separation of 2 steps - one stepping stone!)
            if cluster_num > 0 and clusters_must_be_separated:
                actual_patch_nums_copy = deepcopy(actual_patch_nums)
                # need to use a copy as removing from the same list while iterating over it does not behave properly!
                current_patches = list(set([patch for cluster in cluster_list for patch in cluster]))
                for patch_num in actual_patch_nums_copy:
                    for prev_patch in current_patches:
                        if system_state.patch_adjacency_matrix[prev_patch, patch_num] \
                                > 0 or system_state.patch_adjacency_matrix[patch_num, prev_patch] > 0:
                            actual_patch_nums.remove(patch_num)
                            break
                del actual_patch_nums_copy  # clear this variable as no longer required

            # now choose the first reserve patch in this cluster
            type_patch_nums = initial_cluster_choice(system_state=system_state,
                                                     actual_patch_nums=actual_patch_nums,
                                                     cluster_initial=cluster["initial"])

            # draw one if possible
            draw_num = int(cluster_draw(type_patch_nums, actual_patch_nums, probability=probability_weighting))
            current_cluster.append(draw_num)
            actual_patch_nums.remove(draw_num)

            # Now generate the rest of the cluster
            if cluster["num_patches"] >= 2:
                for next_reserve in range(cluster["num_patches"] - 1):
                    type_patch_nums = cluster_next_element(system_state=system_state,
                                                           actual_patch_nums=actual_patch_nums,
                                                           current_cluster=current_cluster,
                                                           cluster_arch_type=cluster["arch_type"])

                    # draw one if possible
                    draw_num = int(cluster_draw(type_patch_nums, actual_patch_nums, probability=probability_weighting))
                    current_cluster.append(draw_num)
                    actual_patch_nums.remove(draw_num)

        cluster_list = cluster_list + [current_cluster]
    return cluster_list


def initial_cluster_choice(system_state, actual_patch_nums, cluster_initial: list):
    # Called only by cluster_builder().
    # Returns the list of patches eligible to be the first element of the cluster.
    x_max = system_state.dimensions[0]
    y_max = system_state.dimensions[1]
    type_patch_nums = []

    if cluster_initial[0] == "random":
        type_patch_nums = actual_patch_nums

    elif cluster_initial[0] == "patch_number":
        if type(cluster_initial[1]) == int and cluster_initial[1] in actual_patch_nums:
            type_patch_nums.append(cluster_initial[1])
        elif cluster_initial[1] == "high":
            type_patch_nums.append(max(actual_patch_nums))
        elif cluster_initial[1] == "low":
            type_patch_nums.append(min(actual_patch_nums))
        else:
            raise "Cluster initial value patch number not recognised."

    elif cluster_initial[0] == "position":

        if type(cluster_initial[1]) == np.ndarray:  # second argument of form (x,y)
            # iterate over eligible patches to find the one at this location
            for patch_num in actual_patch_nums:
                if (system_state.patch_list[patch_num].position == cluster_initial[1]).all():
                    type_patch_nums.append(patch_num)
                    break
        elif cluster_initial[1] == "edge":
            # for both this and below, at simulation_obj generation we should store in system_state the maximum
            # possible values of x and y as attributes that could be retrieved here to seek {0, max_x} for
            # edge or int(max_x/2) for center.
            for patch_num in actual_patch_nums:
                if system_state.patch_list[patch_num].position[0] in [0, x_max] or \
                        system_state.patch_list[patch_num].position[1] in [0, y_max]:
                    type_patch_nums.append(patch_num)

        elif cluster_initial[1] == "center":  # beware spelling
            for patch_num in actual_patch_nums:
                if system_state.patch_list[patch_num].position[0] in [np.floor(x_max / 2), np.ceil(x_max / 2)] \
                        and system_state.patch_list[patch_num].position[1] in [np.floor(y_max / 2),
                                                                               np.floor(y_max / 2)]:
                    type_patch_nums.append(patch_num)
        else:
            raise "Cluster initial value position not recognised."

    elif cluster_initial[0] in ['degree', 'centrality']:

        attr_name = cluster_initial[0]
        attr_value = cluster_initial[1]
        if type(attr_value) in [int, np.ndarray]:  # specific value desired
            # iterate over eligible patches to find all with this value
            for patch_num in actual_patch_nums:
                if getattr(system_state.patch_list[patch_num], attr_name) == attr_value:
                    type_patch_nums.append(patch_num)
        elif attr_value in ['high', 'low']:
            if attr_value == 'high':
                desired_value = max([getattr(x, attr_name) for x in system_state.patch_list])
            else:
                desired_value = min([getattr(x, attr_name) for x in system_state.patch_list])
            for patch_num in actual_patch_nums:
                if getattr(system_state.patch_list[patch_num], attr_name) == desired_value:
                    type_patch_nums.append(patch_num)
        else:
            raise f"Cluster initial value {cluster_initial[0]} not recognised."

    else:
        raise "Cluster initial value type not recognised."

    return type_patch_nums


def cluster_next_element(system_state, current_cluster: list, actual_patch_nums: list, cluster_arch_type: str):
    # Called only by cluster_builder().
    # Returns the list of patches eligible to be drawn as the next element of the cluster.
    type_patch_nums = []
    if cluster_arch_type == "random":
        type_patch_nums = actual_patch_nums

    elif cluster_arch_type == "star":
        # try to choose a neighbour of the initial node
        for patch_num in actual_patch_nums:
            if system_state.patch_adjacency_matrix[current_cluster[0], patch_num] > 0 \
                    or system_state.patch_adjacency_matrix[patch_num, current_cluster[0]] > 0:
                type_patch_nums.append(patch_num)

    elif cluster_arch_type == "chain":
        # try to choose a neighbour of the final node
        for patch_num in actual_patch_nums:
            if system_state.patch_adjacency_matrix[current_cluster[-1], patch_num] > 0 \
                    or system_state.patch_adjacency_matrix[patch_num, current_cluster[-1]] > 0:
                type_patch_nums.append(patch_num)

    elif cluster_arch_type == "disconnected":
        # try to choose a neighbour of None of the current nodes
        for patch_num in actual_patch_nums:
            for already_chosen_num in current_cluster:
                if system_state.patch_adjacency_matrix[already_chosen_num, patch_num] == 0 \
                        and system_state.patch_adjacency_matrix[patch_num, already_chosen_num] == 0:
                    type_patch_nums.append(patch_num)

    elif cluster_arch_type == "box":
        # try to choose a neighbour of multiple current nodes
        adjacency_counter = []
        for patch_num in actual_patch_nums:
            num_adjacent = 0
            for already_chosen_num in current_cluster:
                if system_state.patch_adjacency_matrix[already_chosen_num, patch_num] > 0 \
                        or system_state.patch_adjacency_matrix[patch_num, already_chosen_num] > 0:
                    num_adjacent += 1
            adjacency_counter.append([patch_num, num_adjacent])
        # what was the greatest adjacency?
        max_adjacency = max([x[1] for x in adjacency_counter])
        # choose from those with greatest adjacency only
        type_patch_nums = [x[0] for x in adjacency_counter if x[1] == max_adjacency]

    elif cluster_arch_type == "position_box":
        # choose a box consisting of positionally-close patches, regardless of actual connectivity
        distance_counter = []
        for patch_num in actual_patch_nums:
            summative_distance = 0.0
            for already_chosen_num in current_cluster:
                summative_distance += np.linalg.norm(system_state.patch_list[patch_num].position
                                                     - system_state.patch_list[already_chosen_num].position)
            distance_counter.append([patch_num, summative_distance])
        # what was minimum total positional distance of any acceptable patch to the existing members of the cluster?
        min_distance = min([x[1] for x in distance_counter])
        # choose from those with least distance only
        type_patch_nums = [x[0] for x in distance_counter if x[1] == min_distance]

    else:
        raise "Cluster arch-type not recognised."
    return type_patch_nums


def cluster_draw(type_patch_nums, actual_patch_nums, probability=None):
    # Called in two places of the loops in cluster_builder().
    # We pass in the desired and the possible patches to be chosen from, and this returns the next element for the
    # cluster if possible.
    if len(type_patch_nums) > 0:
        target = type_patch_nums
    elif len(actual_patch_nums) > 0:
        target = actual_patch_nums
    else:
        raise "Insufficient number of eligible cluster choices."
    if probability is not None:
        # now build the probability distribution and re-normalise as some patches may/not be included
        selected_probability = []
        probability_total = 0.0
        for patch_num in target:
            selected_probability.append(probability[patch_num])
            probability_total += probability[patch_num]
        if probability_total > 0.0:
            final_probability = [x / probability_total for x in selected_probability]
        elif len(target) > 0:
            # there are targets, but the total probability was zero?
            final_probability = [1 / len(target) for _ in target]
            print("Warning: cluster patch choices made but zero total probability. Due to extreme weighting {-1, +1} "
                  "and either on the first draw where baseline was e.g. zero previous perturbations everywhere, or "
                  "because cluster topology restricted available options which may overwrite weighting preferences.")
        else:
            # there are no targets?
            raise "Insufficient number of eligible cluster choices at the target probability stage."
    else:
        # uniform is no distribution is specified for the actual_patch_nums
        final_probability = [1 / len(target) for _ in target]
    draw_num = int(np.random.choice(target, 1, p=final_probability)[0])  # returns 1x1 array so use [0] to extract
    return draw_num
