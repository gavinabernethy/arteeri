import copy
import random
import numpy as np


def reset_temp_values(patch_list):
    # reset all movement and feeding values
    for patch in patch_list:
        patch.sum_competing_for_resources = 0.0
        for local_pop in patch.local_populations.values():
            local_pop.current_temp_change = 0.0
            local_pop.holding_population = local_pop.population
            local_pop.g_values = {}
            local_pop.kills = {
                "g0": {},
                "g1": {},
                "g2": {},
                "g3": {},
            }
            local_pop.killed = {
                "g0": {},
                "g1": {},
                "g2": {},
                "g3": {},
            }
    reset_dispersal_values(patch_list=patch_list)


def reset_dispersal_values(patch_list):
    # reset all movement values only
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            local_pop.population_leave = 0.0
            local_pop.leaving_array = np.zeros([len(patch_list), 1])
            local_pop.population_enter = 0.0


def temporal_function(parameter, time):
    # retrieves the current value of a potentially temporally-varying species-specific parameter
    if parameter["type"] is None:
        value = None
    elif parameter["type"] == "constant":
        value = parameter["constant_value"]
    elif parameter["type"] == "sine":
        value = parameter["amplitude"] * np.sin(time * (2.0 * np.pi / parameter["period"]) + parameter["phase_shift"]) \
                + parameter["vertical_shift"]
    elif parameter["type"] == "vector_exp":
        index = np.mod(time, np.floor(parameter["period"]))
        value = parameter["vector_exp"][index]
    elif parameter["type"] == "vector_imp":
        # find the current time in the cycle
        index = np.mod(time, np.floor(parameter["period"]))  # find largest key below or equal to the index.
        # Note that dictionaries are unordered, so filter the starting_values less than or equal to the index
        # and then choose the largest to get the key for our current time period
        key = max(filter(lambda x: x <= index, parameter["vector_imp"].keys()))
        # pass this back to the dictionary of values
        value = parameter["vector_imp"][key]
    else:
        raise Exception("Type not recognised.")  # Note that we accept 'None' as valid!
    return value


# ------------------------ DISPERSAL TYPES ------------------------ #

def dispersal_scheme_step_polynomial(species_from, movement_score, parameters):
    # Potentially two separate polynomials, whose order and coefficients are contained in lists in the species
    # parameters. These apply separately depending on whether the local population is above or below its
    # species-specific threshold fraction of its patch-specific carrying capacity.
    #
    # e.g. For simple diffusion that scales with local population density, could do: [0, M] [0, M] for some M.
    #
    # In the case of water vole and mink, we want: [0], [M0, M1] with large M0, M1 so that there is no emigration at all
    # under the carrying capacity (except random extra movement), but it is high (and rapidly linearly
    # increasing) when this IS exceeded.
    poly_para = species_from.species.current_coefficients_lists
    density = species_from.holding_population / species_from.carrying_capacity  # note that this depends on patch size
    if density <= poly_para["DENSITY_THRESHOLD"]:
        coefficients = poly_para["UNDER"]
    else:
        coefficients = poly_para["OVER"]
    leaver_proportion = 0.0
    for power, coefficient in enumerate(coefficients):
        leaver_proportion += coefficient * density ** power
    leaving_pop = parameters["pop_dyn_para"][
                      "MU_OVERALL"] * movement_score * leaver_proportion * species_from.carrying_capacity
    # note that the final amount is scaled by carrying capacity
    return leaving_pop


def dispersal_scheme_uniform_diffusion(species_from, movement_score, parameters):
    # simple diffusion weighted by the species-and-habitat-specific traversal score, their movement speed,
    # and the overall movement parameter
    # this is NOT dependent on patch size/carrying capacity in the source patch, only the size of the local population
    leaving_pop = parameters["pop_dyn_para"]["MU_OVERALL"] * \
                  movement_score * \
                  species_from.holding_population
    return leaving_pop


def dispersal_scheme_stochastic_quantity(species_from, movement_score, parameters):
    # stochastic diffusion weighted by the species-and-habitat-specific traversal score and their movement speed
    # and the overall movement parameter
    leaving_pop = np.random.rand() * \
                  parameters["pop_dyn_para"]["MU_OVERALL"] * \
                  movement_score * \
                  species_from.holding_population
    return leaving_pop


def dispersal_scheme_stochastic_binomial(species_from, movement_score, parameters):
    # stochastic binary migration, with probability of emigration weighted by the species mobility, and the amount
    # of movement weighted by the overall movement parameter and the species-and-habitat-specific traversal score.
    leaving_pop = np.random.binomial(1, max(0.0, min(1.0, species_from.species.current_dispersal_mobility))) * \
                  parameters["pop_dyn_para"]["MU_OVERALL"] * \
                  movement_score * \
                  species_from.holding_population
    return leaving_pop


def dispersal_scheme_adaptive(species_from, movement_score, parameters):
    # migration is scaled by the fractional decline in the local population as a direct result of
    # ecological_iterate (feeding and reproduction, after the last dispersal was completed), and scaled by
    # the species-and-habitat-specific traversal score and their movement speed and the overall movement parameter
    non_dispersal_change = species_from.local_growth + species_from.direct_impact_value + \
                            species_from.prey_gain - species_from.predation_loss
    current_pop = species_from.holding_population  # current population has been saved this step following ODE/growth.
    previous_pop = species_from.holding_population - non_dispersal_change
    # did foraging/reproduction/being predated upon result in a net decline for this local population?
    if non_dispersal_change < 0.0 and current_pop > 0.0 and previous_pop != 0.0:
        # instance of movement being enacted
        leaving_pop = parameters["pop_dyn_para"]["MU_OVERALL"] * movement_score * (-non_dispersal_change) / \
                      previous_pop * species_from.holding_population
    else:
        leaving_pop = 0.0
    return leaving_pop


# ----------------------------------------------------------------- #

def calculate_possible_movement(species_from, patch_to_num, parameters):
    leaving_pop = 0.0
    # the species-and-habitat-specific traversal score
    movement_score = species_from.actual_dispersal_targets[patch_to_num]
    # pass to the species-specific dispersal mechanism:
    dispersal_function = {
        "step_poly": dispersal_scheme_step_polynomial,
        "diffusion": dispersal_scheme_uniform_diffusion,
        "stochastic_quantity": dispersal_scheme_stochastic_quantity,
        "stochastic_binomial": dispersal_scheme_stochastic_binomial,
        "adaptive": dispersal_scheme_adaptive,
    }
    if species_from.species.current_dispersal_mechanism == "no_dispersal":
        pass
    else:
        leaving_pop = dispersal_function[species_from.species.current_dispersal_mechanism](species_from,
                                                                                           movement_score, parameters)
        # 1 if destination is a higher patch number, -1 otherwise
        directional_difference = 2 * int(patch_to_num - species_from.patch_num > 0) - 1
        if directional_difference * species_from.species.current_dispersal_direction > 0:
            # adding 1.0 here means that, at least for a patch with equal numbers of links in both directions, on
            # average there is no modification from this directional preference to the total amount of migration
            # although of course in principle that won't exactly work out.
            leaving_pop = (1.0 + abs(species_from.species.current_dispersal_direction)) * leaving_pop
        else:
            # this will just give leaving_pop unchanged for all destination patches if current_dispersal_direction = 0.
            leaving_pop = (1.0 - abs(species_from.species.current_dispersal_direction)) * leaving_pop
    species_min_amount_to_move = max(0.0, species_from.species.minimum_population_size)
    if leaving_pop < species_min_amount_to_move:
        # the amount that wants to move is smaller than the minimum population size - should it be allowed to move?
        if species_from.species.always_move_with_minimum:
            leaving_pop = species_min_amount_to_move
        else:
            leaving_pop = 0.0
    leaving_pop = min(leaving_pop, species_from.holding_population)
    species_from.population_leave += leaving_pop
    species_from.leaving_array[patch_to_num] = leaving_pop


def pre_dispersal_of_local_population(
        patch_list,
        parameters,
        local_pop,
        specified_fraction=None,
        is_dispersal_override=False,
        # if True, can call this function in a perturbation and force dispersal even if not normally permitted
):
    # Calculates movements FROM a given patch. But not actually enacted yet, as these all need to occur simultaneously.

    # Record current population as the potential dispersers (for source calculation)
    local_pop.potential_dispersal = local_pop.holding_population

    # First, check all conditions, including: is there somewhere that CAN actually be travelled to?
    if local_pop.holding_population > 0.0 and len(local_pop.actual_dispersal_targets) > 0 \
            and (is_dispersal_override or local_pop.species.is_dispersal):

        # possible
        for patch_to_num in local_pop.actual_dispersal_targets:
            calculate_possible_movement(species_from=local_pop,
                                        patch_to_num=patch_to_num,
                                        parameters=parameters)

        # is the amount to leave pre-determined?
        if specified_fraction is not None:
            # this is currently only for use in a population perturbation where we force dispersal (potentially of
            # specific species from specific habitats), but where dispersal mobility is NON-ZERO otherwise the
            # leaving array will be empty as the mobility score is used when determining the possible target list
            local_pop.leaving_array = local_pop.leaving_array * local_pop.holding_population * specified_fraction / \
                                      local_pop.population_leave
            local_pop.population_leave = sum(local_pop.leaving_array)
        else:
            # Otherwise (for normal dispersal):

            if local_pop.population_leave > local_pop.holding_population:
                # re-scale if necessary

                if local_pop.holding_population > len(local_pop.actual_dispersal_targets) \
                        * local_pop.species.minimum_population_size:
                    # in this case we can reduce all amounts uniformly and at least some will remain above minimum
                    local_pop.leaving_array = local_pop.leaving_array * local_pop.holding_population / \
                                              local_pop.population_leave
                else:
                    # otherwise, we incrementally reduce random amounts to random destinations so that some still have a
                    # significant amount of migrants (probably but not guaranteed in all circumstances, draw-dependent).
                    temp_pop_leave = local_pop.population_leave
                    while temp_pop_leave > local_pop.holding_population:  # while too many
                        # choose a current destination randomly
                        destination = random.choice(list(local_pop.actual_dispersal_targets.keys()))
                        current_amount = local_pop.leaving_array[destination]
                        # if the amount going there is at least as big as the minimum population, randomly reduce it!
                        if current_amount >= local_pop.species.minimum_population_size:
                            draw_reduction_to = np.random.uniform(0.0, current_amount)
                            if draw_reduction_to < local_pop.species.minimum_population_size:
                                # if we would be reducing below the minimum pop, just set that movement to zero
                                temp_pop_leave -= current_amount
                                local_pop.leaving_array[destination] = 0.0
                            else:
                                temp_pop_leave -= (current_amount - draw_reduction_to)  # actual reduction is difference
                                local_pop.leaving_array[destination] = draw_reduction_to  # set to new (reduced) amount
                local_pop.population_leave = float(sum(local_pop.leaving_array))
            else:
                # Then add optional stochastic modifiers to see if ONE more individual wanders out;
                # This must be done AFTER the normalisation so that it is one WHOLE individual (relative to species
                # minimum population size) and not subdivided;
                # Only do this if rescaling NOT necessary and >= minimum pop. still available after regular dispersal.
                binomial = np.random.binomial(1, parameters["species_para"][local_pop.species.name][
                    "DISPERSAL_PARA"]["BINOMIAL_EXTRA_INDIVIDUAL"])
                if local_pop.holding_population - local_pop.population_leave >= \
                        local_pop.species.minimum_population_size and binomial:
                    local_pop.population_leave += local_pop.species.minimum_population_size

                    # select one destination at random to receive the +1 member
                    destination = random.choice(list(local_pop.actual_dispersal_targets.keys()))
                    # we actually scale the amount added to leaving array so that the subsequently applied dispersal
                    # penalty (if non-zero) will NOT affect this individual (but will impact on ALL other dispersal)!
                    local_pop.leaving_array[destination] += \
                        local_pop.species.minimum_population_size / local_pop.species.dispersal_efficiency

        # then update the destinations with the final arrivals
        for patch_to_num in local_pop.actual_dispersal_targets:
            species_find = patch_list[patch_to_num].local_populations[local_pop.name]

            # the final destination populations are updated with the confirmed leavers, scaled by an efficiency
            # penalty if this is specified for the species or for all in the simulation
            species_find.population_enter += species_find.species.dispersal_efficiency * float(
                local_pop.leaving_array[patch_to_num])


def find_best_actual_scores(local_pop, target, query_attr, max_path_attr, mobility_scaling_attr,
                            heaviside_threshold_attr, is_heaviside_manual, heaviside_manual_value,
                            home_patch_traversal_score, is_scaled_target_size):
    # used by both build_actual_dispersal_targets() and build_interacting_populations_list() to determine the
    # local_population's access to distant patches given their CURRENT mobility scores and path length restrictions
    #
    # We already have the potential travel COSTs, so we choose the best available one, then add search cost (for prey
    # or nest sites) in the target patch, take reciprocal and scale by mobility (and size if required).
    path_length = 0
    if home_patch_traversal_score == 0.0:
        # this is not strictly needed as a separate consideration, but it speeds up the calculation
        score = 0.0
    else:
        if getattr(local_pop.species, query_attr):
            # path length restricted, but is the best one within the allowed range anyway?
            if target["routes"]["best"][0] <= getattr(local_pop.species, max_path_attr):
                cost = target["routes"]["best"][1]
                path_length = target["routes"]["best"][0]
            else:
                cost = float('inf')
                # otherwise look for best reachable
                for try_path_length in range(1, getattr(local_pop.species, max_path_attr) + 1):
                    if try_path_length in target:
                        cost = min(cost, target[try_path_length][0])
                        path_length = try_path_length
        else:
            # path unrestricted so just return the best overall cost
            cost = target["routes"]["best"][1]
            path_length = target["routes"]["best"][0]

        # now pass into Heaviside step function
        if is_heaviside_manual:
            heaviside_threshold = heaviside_manual_value
        else:
            heaviside_threshold = getattr(local_pop.species, heaviside_threshold_attr)
        if cost <= heaviside_threshold:
            filtered_cost = 0.0
        else:
            filtered_cost = cost - heaviside_threshold

        # add search cost in target patch
        filtered_cost += target["target_patch_size"] / target["target_patch_traversal"]

        # scale by mobility attribute (and target_size if appropriate) and return
        score = getattr(local_pop.species, mobility_scaling_attr) * filtered_cost ** (-1.0)

        if is_scaled_target_size:
            score = score * target["target_patch_size"]

    return score, path_length


def build_actual_dispersal_targets(patch_list, species_list, is_dispersal, time):
    # Determine a dictionary of which locations can ACTUALLY be reached and with what movement score, given
    # the current dispersal mobility score, minimum link strength, and path length restriction
    if is_dispersal:

        # initialise if running for the first time
        for species in species_list:
            if None in [species.current_dispersal_mechanism, species.current_dispersal_mobility,
                        species.current_max_dispersal_path_length, species.current_minimum_link_strength_dispersal,
                        species.current_dispersal_direction, species.current_coefficients_lists]:
                species.current_dispersal_mechanism = temporal_function(
                    species.dispersal_para["DISPERSAL_MECHANISM"], time)
                species.current_dispersal_mobility = temporal_function(
                    species.dispersal_para["DISPERSAL_MOBILITY"], time)
                species.current_max_dispersal_path_length = temporal_function(
                    species.dispersal_para["MAX_DISPERSAL_PATH_LENGTH"], time)
                species.current_minimum_link_strength_dispersal = temporal_function(
                    species.dispersal_para["MINIMUM_LINK_STRENGTH_DISPERSAL"], time)
                species.current_dispersal_direction = temporal_function(
                    species.dispersal_para['DISPERSAL_DIRECTION'], time)
                species.current_coefficients_lists = temporal_function(
                    species.dispersal_para['COEFFICIENTS_LISTS'], time)
        for patch in patch_list:
            for local_pop in patch.local_populations.values():
                # reset them
                local_pop.actual_dispersal_targets = {}
                temp_dict = {}

                if local_pop.species.is_dispersal:
                    # Iterate through the patches that can IN PRINCIPLE be reached, as determined during
                    # system_state.build_paths_and_adjacency_lists() based on more fundamental properties of the
                    # spatial network topology/adjacency, habitat types, and species-habitat traversal scores.
                    this_patch_species_traversal = patch.this_habitat_species_traversal[local_pop.species.name]
                    for reachable_patch_num in patch.adjacency_lists[local_pop.name]:

                        # don't include same patch
                        if reachable_patch_num != patch.number:

                            z = patch.species_movement_scores[local_pop.name][reachable_patch_num]
                            target_score, unused_path_length = find_best_actual_scores(
                                local_pop=local_pop,
                                target=z,
                                query_attr="is_dispersal_path_restricted",
                                max_path_attr="current_max_dispersal_path_length",
                                mobility_scaling_attr="current_dispersal_mobility",
                                is_heaviside_manual=True,
                                heaviside_manual_value=0.0,
                                heaviside_threshold_attr=None,
                                home_patch_traversal_score=this_patch_species_traversal,
                                is_scaled_target_size=True,
                            )
                            if target_score >= local_pop.species.current_minimum_link_strength_dispersal:
                                temp_dict[reachable_patch_num] = target_score
                local_pop.actual_dispersal_targets = temp_dict


def build_interacting_populations_list(patch_list, species_list, is_nonlocal_foraging, is_local_foraging_ensured, time):
    # this function should NOT be only looking at non-zero population sizes, as the list will not be rebuilt
    # if they are populated at a later time. It also involves the within-patch predator-prey interactions, so
    # do NOT skip this method if is_nonlocal_foraging is false

    # initialise if running for the first time
    for species in species_list:
        if None in [species.current_prey_dict, species.current_foraging_mobility, species.current_foraging_kappa,
                    species.current_max_foraging_path_length, species.current_minimum_link_strength_foraging,
                    species.current_predation_rate, species.current_predation_efficiency,
                    species.current_predation_focus]:
            species.current_prey_dict = temporal_function(species.predation_para["PREY_DICT"], time)
            species.current_foraging_mobility = temporal_function(species.predation_para['FORAGING_MOBILITY'], time)
            species.current_foraging_kappa = temporal_function(species.predation_para["FORAGING_KAPPA"], time)
            species.current_max_foraging_path_length = temporal_function(
                species.predation_para["MAX_FORAGING_PATH_LENGTH"], time)
            species.current_minimum_link_strength_foraging = temporal_function(
                species.predation_para["MINIMUM_LINK_STRENGTH_FORAGING"], time)
            species.current_predation_efficiency = temporal_function(
                species.predation_para["PREDATION_EFFICIENCY"], time)
            species.current_predation_focus = temporal_function(
                species.predation_para["PREDATION_FOCUS"], time)
            species.current_predation_rate = temporal_function(species.predation_para["PREDATION_RATE"], time)

    # reset interacting population lists
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            local_pop.interacting_populations = []
    # Build list of dictionaries of all the populations that each population can interact with (in either way)
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            this_patch_species_traversal = patch.this_habitat_species_traversal[local_pop.species.name]
            this_patch_species_feeding = patch.this_habitat_species_feeding[local_pop.species.name]

            for patch_to_num in patch.adjacency_lists[local_pop.name]:
                patch_to = patch_list[patch_to_num]

                # only permit including local_populations of other patches if stated
                if is_nonlocal_foraging or patch_to_num == patch.number:

                    for local_pop_to in patch_to.local_populations.values():
                        # we need to include both, as there may be a local population who can reach but cannot be
                        # reached, and we need to note them as interacting with the other the adjacency lists may not
                        # be symmetric, so it is insufficient to just have one of these

                        # now account for species-specific limitations and foraging path length limits
                        local_pop_score = 0.0
                        local_pop_to_score = 0.0
                        patch_to_species_traversal = patch_to.this_habitat_species_traversal[local_pop_to.species.name]
                        patch_to_species_feeding = patch_to.this_habitat_species_feeding[local_pop_to.species.name]

                        path_to_length = 0
                        path_from_length = 0
                        if patch_to_num == patch.number:
                            # for WITHIN-PATCH FEEDING:
                            if is_local_foraging_ensured:
                                # with this global option set to true, within-patch feeding is always set to
                                # precisely 1.0 for any species, as in earlier versions of Artemis
                                local_pop_score = 1.0
                                local_pop_to_score = 1.0
                            else:
                                # normally, within-patch score is now master_mu * species_mu * this_habitat_traversal
                                if local_pop.species.current_foraging_mobility is None:
                                    local_pop_score = None
                                else:
                                    local_pop_score = \
                                        local_pop.species.current_foraging_mobility * this_patch_species_traversal

                                if local_pop_to.species.current_foraging_mobility is None:
                                    local_pop_to_score = None
                                else:
                                    local_pop_to_score = \
                                        local_pop_to.species.current_foraging_mobility * patch_to_species_traversal
                        else:
                            if local_pop.species.is_nonlocal_foraging:
                                # score dictionary for THIS species' local population to THAT patch
                                z = patch.species_movement_scores[local_pop.name][patch_to_num]
                                local_pop_score, path_to_length = find_best_actual_scores(
                                    local_pop=local_pop, target=z,
                                    query_attr="is_foraging_path_restricted",
                                    max_path_attr="current_max_foraging_path_length",
                                    mobility_scaling_attr="current_foraging_mobility",
                                    is_heaviside_manual=False,
                                    heaviside_manual_value=None,
                                    heaviside_threshold_attr="current_foraging_kappa",
                                    home_patch_traversal_score=this_patch_species_traversal,
                                    is_scaled_target_size=False,
                                )
                                if local_pop_score < local_pop.species.current_minimum_link_strength_foraging:
                                    local_pop_score = 0.0

                            if local_pop_to.species.is_nonlocal_foraging:
                                # score dictionary for THAT species' local population to THIS patch
                                z = patch_to.species_movement_scores[local_pop_to.name][patch.number]
                                local_pop_to_score, path_from_length = find_best_actual_scores(
                                    local_pop=local_pop_to, target=z,
                                    query_attr="is_foraging_path_restricted",
                                    max_path_attr="current_max_foraging_path_length",
                                    mobility_scaling_attr="current_foraging_mobility",
                                    is_heaviside_manual=False,
                                    heaviside_manual_value=None,
                                    heaviside_threshold_attr="current_foraging_kappa",
                                    home_patch_traversal_score=patch_to_species_traversal,
                                    is_scaled_target_size=False,
                                )
                                if local_pop_to_score < local_pop_to.species.current_minimum_link_strength_foraging:
                                    local_pop_to_score = 0.0

                        # Only actual interactions (at least one-way) get added to the list, and only if the in-patch
                        # feeding score of that species was non-zero
                        if (local_pop_score != 0.0 and this_patch_species_feeding > 0.0) or \
                                (local_pop_to_score != 0.0 and patch_to_species_feeding > 0.0):
                            local_pop.interacting_populations.append({"object": local_pop_to,
                                                                      "score_to": local_pop_score,
                                                                      "score_from": local_pop_to_score,
                                                                      "is_same_patch": patch_to_num == patch.number,
                                                                      "path_to_length": path_to_length,
                                                                      "path_from_length": path_from_length,
                                                                      })
                            local_pop_to.interacting_populations.append({"object": local_pop,
                                                                         "score_to": local_pop_to_score,
                                                                         "score_from": local_pop_score,
                                                                         "is_same_patch": patch_to_num == patch.number,
                                                                         "path_to_length": path_from_length,
                                                                         "path_from_length": path_to_length,
                                                                         })

    # check that there are no duplicates - this will be caused by two populations who CAN both reach each other,
    # and this leads to inconsistencies in the predation calculations
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            temp_list = []
            for x in local_pop.interacting_populations:
                if x not in temp_list:
                    temp_list.append(x)
            local_pop.interacting_populations = temp_list


def checker(output_attribute, input_temporal, is_change):
    # used by change_checker() to compare new and previous values for changes
    if output_attribute != input_temporal:
        is_change = True
    return input_temporal, is_change


def change_checker(species_list, patch_list, time, step, is_dispersal, is_nonlocal_foraging, is_local_foraging_ensured):
    # this function deals with the temporal variation of species parameters.
    #
    # update temporary values of species properties and check if anything has changed
    is_change = False
    for species in species_list:
        # check if any species-dependent properties have changed due to temporal variation.
        # if so, update the current values and mark a change so that the lists can be rebuilt
        # for each entry: [to update, base attribute, nested attribute (if needed)]
        update_and_check = [
            ['current_r_value', 'growth_para', 'R'],
            ['current_prey_dict', 'predation_para', 'PREY_DICT'],
            ['current_predation_efficiency', 'predation_para', 'PREDATION_EFFICIENCY'],
            ['current_predation_focus', 'predation_para', 'PREDATION_FOCUS'],
            ['current_predation_rate', 'predation_para', 'PREDATION_RATE'],
            ['current_foraging_mobility', 'predation_para', 'FORAGING_MOBILITY'],
            ['current_foraging_kappa', 'predation_para', 'FORAGING_KAPPA'],
            ['current_minimum_link_strength_foraging', 'predation_para', 'MINIMUM_LINK_STRENGTH_FORAGING'],
            ['current_max_foraging_path_length', 'predation_para', 'MAX_FORAGING_PATH_LENGTH'],
        ]
        # don't bother checking dispersal parameters if dispersal is not enabled (remember that is_dispersal is NOT
        # allowed to vary with time!
        if species.is_dispersal:
            update_and_check = update_and_check + [
                ['current_dispersal_mobility', 'dispersal_para', 'DISPERSAL_MOBILITY'],
                ['current_dispersal_direction', 'dispersal_para', 'DISPERSAL_DIRECTION'],
                ['current_dispersal_mechanism', 'dispersal_para', 'DISPERSAL_MECHANISM'],
                ['current_coefficients_lists', 'dispersal_para', 'COEFFICIENTS_LISTS'],
                ['current_minimum_link_strength_dispersal', 'dispersal_para', 'MINIMUM_LINK_STRENGTH_DISPERSAL'],
                ['current_max_dispersal_path_length', 'dispersal_para', 'MAX_DISPERSAL_PATH_LENGTH']
            ]
        for _ in update_and_check:
            if _[2] is not None:
                # i.e. if we need to access a dictionary-nested attribute
                result, is_change = checker(getattr(species, _[0]),
                                            temporal_function(getattr(species, _[1])[_[2]], time), is_change)
            else:
                # or if we can access the attribute for comparison directly
                result, is_change = checker(getattr(species, _[0]),
                                            temporal_function(getattr(species, _[1]), time), is_change)
            setattr(species, _[0], result)

        # check predation efficiency and focus are suitable only when set (instead of for every predation loop)
        if species.current_predation_efficiency is not None and (
                species.current_predation_efficiency < 0.0 or species.current_predation_efficiency > 1.0):
            raise Exception("ERROR: species predation_efficiency should be in the interval [0, 1] or None")

        if species.current_predation_focus is not None and species.current_predation_focus < 0.0:
            raise Exception("ERROR: species predation_focus should be non-negative or None")

    if is_change or step == 0:
        # something changed (or we are at the start of the simulation) - need to rebuild the lists!
        # IMPORTANT: this is okay, only so long as we do not change any of:
        # - the fundamental ability of a given species to pass through a given habitat type,
        # - the spatial topology and adjacency of the patch network,
        # - the size or habitat type of any patch.
        # Such alterations would instead be handled by "perturbation.py" and require us to re-determine the possible
        # pathing structure in a more fundamental manner by executing the (very intensive) function:
        # system_state.build_all_patches_species_paths_and_adjacency()
        #
        # You must also be aware that this occurs by TIME (i.e. day) rather than step - thus, for example, changes that
        # are specified to change daily by sine function will only update every 10 steps if main_para:steps_to_days=10.
        # This gives us some modulo control over how often to expend computational time updating.
        print(f" ...Step {step}: species behaviour change identified.")
        build_actual_dispersal_targets(
            patch_list=patch_list, species_list=species_list,
            is_dispersal=is_dispersal, time=time)
        build_interacting_populations_list(
            patch_list=patch_list, species_list=species_list,
            is_nonlocal_foraging=is_nonlocal_foraging, is_local_foraging_ensured=is_local_foraging_ensured,
            time=time)


def foraging_calculator(patch_list, time, current_patch_list):
    # this function is responsible for how we decide what predator populations will eat using the g0-g3 system

    # calculate idealised feeding (g0 -> g1)
    for patch_num in current_patch_list:
        habitat_reproduction_scores = patch_list[patch_num].this_habitat_species_feeding
        for local_pop in patch_list[patch_num].local_populations.values():
            if local_pop.holding_population > 0.0 and habitat_reproduction_scores[local_pop.species.name] > 0.0:
                # no predation gain if habitat-feeding score is 0.0, and distance-foraging metrics set to -1 for visual
                if local_pop.species.current_prey_dict is not None and len(local_pop.species.current_prey_dict) != 0:
                    local_pop.calculate_predation(time=time)
                # g0 is the total prey available for that local population
                # g1 is the actual number of kills it would make if there are no other competitors to share with
            else:
                local_pop.weighted_foraging_distance = -1.0  # default to -1 rather than 0 to distinguish absence or
                local_pop.maximum_foraging_distance = -1.0  # lack of predation by this species

    # competition and prey allocation (g1 -> g2)
    for patch_num in current_patch_list:
        for local_pop in patch_list[patch_num].local_populations.values():
            if local_pop.holding_population > 0.0 and len(local_pop.species.predator_list) != 0:
                local_pop.predator_allocation()
                # g2 is the total prey available after sharing allocation with other predators

    # calculate top-up feeding (g2 -> g3)
    for patch_num in current_patch_list:
        habitat_reproduction_scores = patch_list[patch_num].this_habitat_species_feeding
        for local_pop in patch_list[patch_num].local_populations.values():
            if local_pop.holding_population > 0.0 and habitat_reproduction_scores[local_pop.species.name] > 0.0:
                if local_pop.species.current_prey_dict is not None and len(local_pop.species.current_prey_dict) != 0:
                    local_pop.predator_shortfall_distribution()
                    # g3 is the actual number of kills that will be implemented


def growth_caller(parameters, patch_list, time, alpha, is_dispersal, current_patch_list):
    # this implements local growth (i.e. reproduction and mortality) across the entire system, looking at
    # the .holding_population's and adding the resulting changes to the .current_temp_change's
    for patch_num in current_patch_list:
        patch_list[patch_num].sum_competing_for_resources = 0.0
        for local_pop in patch_list[patch_num].local_populations.values():
            if local_pop.holding_population > 0.0:
                # sum up the total resource competition first
                patch_list[patch_num].sum_competing_for_resources += local_pop.resource_usage_conversion \
                                                                     * local_pop.holding_population
        for local_pop in patch_list[patch_num].local_populations.values():
            if local_pop.holding_population > 0.0:
                local_pop.growth(parameters=parameters, time=time, alpha=alpha,
                                 patch_competitors=patch_list[patch_num].sum_competing_for_resources)
            else:
                local_pop.local_growth = 0.0


def foraging_caller(parameters, patch_list, time, alpha, is_dispersal, current_patch_list):
    # this implements predation across the entire system, calling the functional responses twice and looking at
    # the .holding_population's for predator and prey population values to calculate from.
    #
    # first we calculate the desired feeding for all local populations
    foraging_calculator(patch_list=patch_list, time=time, current_patch_list=current_patch_list)
    # then enact the feeding result - adding the resulting changes to the .current_temp_change's
    for patch in patch_list:
        for local_population in patch.local_populations.values():
            this_patch_species_feeding = patch.this_habitat_species_feeding[local_population.species.name]
            local_population.foraging(this_patch_species_feeding=this_patch_species_feeding)


def direct_impact_caller(parameters, patch_list, time, alpha, is_dispersal, current_patch_list):
    # this implements direct impact across the entire system, looking at the .holding_population's and adding the
    # resulting changes to the .current_temp_change's
    for patch in patch_list:
        for local_population in patch.local_populations.values():
            local_population.direct_impact(time=time)


def dispersal_caller(parameters, patch_list, time, alpha, is_dispersal, current_patch_list):
    # this implements dispersal across the entire system, looking at the .holding_population's and adding the resulting
    # changes to the .current_temp_change's
    if is_dispersal:
        # build the temporary movement first
        for patch_from_num in current_patch_list:
            for local_pop in patch_list[patch_from_num].local_populations.values():
                # this nested for loop is calculating the leaving values from this patch, which are all zero by default
                # so if there is no population then it can be safely skipped for efficiency (happens within function)
                pre_dispersal_of_local_population(patch_list=patch_list, parameters=parameters, local_pop=local_pop)
        # finally enact all movement
        for patch in patch_list:
            for local_pop in patch.local_populations.values():
                local_pop.current_temp_change += local_pop.population_enter - local_pop.population_leave


def update_populations(patch_list, species_list, time, step, parameters, current_patch_list, is_ode_recordings):
    # this is the full function for a single standard iteration of the ecological model - including growth
    # (reproduction and mortality), predation and being predated upon, any special direct impacts or pure direct
    # impacts, and dispersal.
    # The order in which these sub-stages of a single step are enacted is specified according to their priority in
    # the main_para.
    alpha = parameters["pop_dyn_para"]["COMPETITION_ALPHA_SCALING"]
    is_nonlocal_foraging = parameters["pop_dyn_para"]["IS_NONLOCAL_FORAGING_PERMITTED"]
    is_local_foraging_ensured = parameters["pop_dyn_para"]["IS_LOCAL_FORAGING_ENSURED"]
    is_dispersal = parameters["pop_dyn_para"]["IS_DISPERSAL_PERMITTED"]

    function_priority_dictionary = parameters["main_para"]["ECO_PRIORITIES"]
    function_name_to_actual = {
        "foraging": foraging_caller,
        "direct_impact": direct_impact_caller,
        "growth": growth_caller,
        "dispersal": dispersal_caller,
    }

    # reset temporary and previous values
    reset_temp_values(patch_list=patch_list)

    # did any species parameters change?
    change_checker(species_list=species_list, patch_list=patch_list, time=time, step=step,
                   is_dispersal=is_dispersal, is_nonlocal_foraging=is_nonlocal_foraging,
                   is_local_foraging_ensured=is_local_foraging_ensured)

    # iterate over the 4 possible priorities (if all four sub-stages have their own separate timing within a step)
    for priority in range(4):
        if len(function_priority_dictionary[priority]) > 0:
            # are there any ecological functions at this level of priority?
            for function in function_priority_dictionary[priority]:
                # each of the four functions makes changes to the local_population.current_temp_change, WITHOUT
                # checking if the population would have gone negative at that point.
                # thus all functions that share the same priority are enacted concurrently
                function_name_to_actual[function](
                    parameters=parameters, patch_list=patch_list, time=time, alpha=alpha,
                    is_dispersal=is_dispersal, current_patch_list=current_patch_list,
                )
            # now we will update all local populations with the total result of all functions that were applied at this
            # priority (i.e. this explicit sub-step within the step)
            for patch in patch_list:
                for local_pop in patch.local_populations.values():

                    if parameters["main_para"]["MODEL_TIME_TYPE"] == "discrete":
                        new_population = local_pop.holding_population + local_pop.current_temp_change
                    elif parameters["main_para"]["MODEL_TIME_TYPE"] == "continuous":
                        new_population = local_pop.holding_population + parameters[
                            "main_para"]["EULER_STEP"] * local_pop.current_temp_change
                    else:
                        raise Exception("Model type not recognised in 'main_para[MODEL_TIME_TYPE]' -"
                                        " discrete or continuous?")

                    # only at the end of a priority sub-step do we NOW check that the net result has not made this
                    # population go negative
                    if 0.0 != new_population < local_pop.species.minimum_population_size:
                        # # CHECK: when is the overwriting occurring?
                        # print(f"Step {step}: patch {patch.number} species {local_pop.species.name } "
                        #       f"overwriting {new_population}")
                        new_population = 0.0
                    # finally, reset the current_temp_change running total for the next priority sub-step, and
                    # update the holding_population which is what the next stage's functions will look at
                    local_pop.holding_population = copy.deepcopy(new_population)
                    local_pop.current_temp_change = 0.0

    # finally, for all populations we need to check a special case, calculate the net change that has occurred due
    # to internal (i.e. non-dispersal) effects, update the .population to conclude the step, and record all histories
    for patch in patch_list:
        for local_pop in patch.local_populations.values():
            if local_pop.species.is_predation_only_prevents_death:
                # cannot gain new members due to predation (only when permitted by the growth function) however,
                # they can keep themselves alive due to predation counting the natural death rate or carrying cap.
                # thus we must always overwrite these in the end for such species! However this does NOT mean that we
                # have treated as though all processes are concurrent - the actual values of each component have
                # been computed differently based on the sub-stage sequencing.

                #
                # check if they grew
                if local_pop.local_growth < 0:
                    # allow predation to only save the population by undoing the negative impact of growth

                    local_pop.holding_population = min(0.0, local_pop.local_growth + local_pop.prey_gain) + \
                                                   local_pop.direct_impact_value - local_pop.predation_loss + \
                                                   local_pop.population_enter - local_pop.population_leave
                else:
                    # discount the beneficial effect of predation on this predator altogether
                    local_pop.holding_population = local_pop.local_growth + local_pop.direct_impact_value - \
                                                   local_pop.predation_loss + local_pop.population_enter - \
                                                   local_pop.population_leave
                if local_pop.holding_population < local_pop.species.minimum_population_size:
                    local_pop.holding_population = 0.0

            # Record new populations, set new population from holding_population, and call ODE recording for this step
            # Note that "local_pop.holding_population - local_pop.population" would give the TOTAL change in any case.
            local_pop.internal_change = local_pop.holding_population - local_pop.population - \
                                        local_pop.population_enter + local_pop.population_leave

            # # CHECK: do the discrepancies here match when excess deaths are overwritten by the min function?
            # expected_change = local_pop.local_growth + local_pop.direct_impact_value + local_pop.population_enter + \
            #     local_pop.prey_gain - local_pop.predation_loss - local_pop.population_leave
            # discrepancy = np.abs(local_pop.internal_change - expected_change)
            # if discrepancy > 0.0000001:
            #     print(f"Step {step}: patch {patch.number} species {local_pop.species.name } "
            #           f"discrepancy is {discrepancy}")

            local_pop.population = copy.deepcopy(local_pop.holding_population)
            if is_ode_recordings:
                local_pop.ode_recordings(time=time, step=step)
            local_pop.record_population_history()
