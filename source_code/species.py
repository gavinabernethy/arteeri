class Species:

    def __init__(self,
                 name=None,
                 species_num=None,
                 core_para=None,
                 initial_population_para=None,
                 growth_para=None,
                 predation_para=None,
                 dispersal_para=None,
                 pure_direct_impact_para=None,
                 direct_impact_on_me=None,
                 perturbation_para=None,
                 population_dynamics_para=None,  # not stored, only used in construction
                 ):

        # fundamentals:
        self.name = name
        self.species_num = species_num

        # core:
        self.minimum_population_size = set_default_value(core_para, "MINIMUM_POPULATION_SIZE", 0.0)
        self.lifespan = set_default_value(core_para, "LIFESPAN", 0.0)
        self.seasonal_period = set_default_value(core_para, "SEASONAL_PERIOD", 0)
        self.predator_list = set_default_value(core_para, "PREDATOR_LIST", None)

        # initial:
        self.initial_population_para = initial_population_para
        self.initial_population_mechanism = set_default_value(
            initial_population_para, "INITIAL_POPULATION_MECHANISM", None)

        # growth:
        self.growth_para = growth_para
        self.resource_usage_conversion = set_default_value(
            growth_para, "RESOURCE_USAGE_CONVERSION", None)
        self.growth_function = set_default_value(growth_para, "GROWTH_FUNCTION", None)

        # growth - offset:
        self.is_growth_offset = set_default_value(
            growth_para, ["ANNUAL_OFFSET", "IS_GROWTH_OFFSET"], False)
        self.growth_annual_duration = set_default_value(
            growth_para, ["ANNUAL_OFFSET", "ANNUAL_DURATION"], False)
        self.growth_vector_offset_species = set_default_value(
            growth_para, ["ANNUAL_OFFSET", "GROWTH_OFFSET_SPECIES"], False)
        self.is_growth_offset_local = set_default_value(
            growth_para, ["ANNUAL_OFFSET", "IS_GROWTH_OFFSET_LOCAL"], False)
        self.growth_vector_offset_local = set_default_value(
            growth_para, ["ANNUAL_OFFSET", "GROWTH_OFFSET_LOCAL"], False)

        # predation:
        self.predation_para = predation_para
        self.is_predation_only_prevents_death = set_default_value(
            predation_para, "IS_PREDATION_ONLY_PREVENTS_DEATH", False)
        self.is_nonlocal_foraging = set_default_value(
            predation_para, "IS_NONLOCAL_FORAGING", False)
        self.is_foraging_path_restricted = set_default_value(
            predation_para, "IS_NONLOCAL_FORAGING_PATH_RESTRICTED", False)
        self.predation_focus_type = set_default_value(
            predation_para, "PREDATION_FOCUS_TYPE", False)


        # dispersal:
        self.dispersal_para = dispersal_para
        self.is_dispersal = set_default_value(dispersal_para, "IS_DISPERSAL", None)
        self.is_dispersal_path_restricted = set_default_value(
            dispersal_para, "IS_DISPERSAL_PATH_RESTRICTED", False)
        self.always_move_with_minimum = set_default_value(
            dispersal_para, "ALWAYS_MOVE_WITH_MINIMUM", False)
        dispersal_penalty = max(set_default_value(dispersal_para, "SS_DISPERSAL_PENALTY", 0.0),
                                population_dynamics_para["GENERAL_DISPERSAL_PENALTY"])
        self.dispersal_efficiency = 1.0 - min(0.999999999, dispersal_penalty)  # convert penalty to efficiency in [0, 1]


        # direct impact:
        self.pure_direct_impact_para = pure_direct_impact_para
        self.is_pure_direct_impact =  set_default_value(pure_direct_impact_para, "IS_PURE_DIRECT_IMPACT", False)
        self.is_direct_offset = set_default_value(pure_direct_impact_para, ["ANNUAL_OFFSET", "IS_DIRECT_OFFSET"], False)
        self.direct_annual_duration = set_default_value(
            pure_direct_impact_para, ["ANNUAL_OFFSET", "ANNUAL_DURATION"], 0)
        self.direct_vector_offset_species = set_default_value(
            pure_direct_impact_para, ["ANNUAL_OFFSET", "DIRECT_OFFSET_SPECIES"], None)
        self.is_direct_offset_local = set_default_value(
            pure_direct_impact_para, ["ANNUAL_OFFSET", "IS_DIRECT_OFFSET_LOCAL"], False)
        self.direct_vector_offset_local = set_default_value(
            pure_direct_impact_para, ["ANNUAL_OFFSET", "DIRECT_OFFSET_LOCAL"], None)
        self.direct_impact_on_me = direct_impact_on_me

        # perturbation:
        self.perturbation_para = perturbation_para
        self.is_perturbs_environment = set_default_value(perturbation_para, "IS_PERTURBS_ENVIRONMENT", False)

        # CURRENT holding values - growth:
        self.current_r_value = None
        self.current_cml_para = None  # parameters for sine, tent, shift maps

        # CURRENT holding values - foraging
        self.current_prey_dict = None
        self.current_predation_pragmatism = None
        self.current_predation_focus = None
        self.current_predation_rate = None
        self.current_foraging_mobility = None
        self.current_foraging_kappa = None
        self.current_max_foraging_path_length = None
        self.current_minimum_link_strength_foraging = None

        # CURRENT holding values - dispersal
        self.current_dispersal_mobility = None
        self.current_dispersal_direction = None
        self.current_dispersal_mechanism = None
        self.current_coefficients_lists = None
        self.current_max_dispersal_path_length = None
        self.current_minimum_link_strength_dispersal = None

def set_default_value(parameter_dict, target_name, default_value):
    # This is used throughout the Species() initialisation, allowing us to simply pass None for entire _para sections,
    #  such as dispersal_para, perturbation_para etc. when specifying the species parameter dictionary, and the
    #  defaults (mostly no effect) will be used. Not recommended for key parameters such as core_para and growth_para.
    #
    # It checks the key (which can be one or two levels) and sets the default value if it does not exist.
    if parameter_dict is None:
        return default_value
    else:
        if type(target_name) == str:
            try:
                output = parameter_dict[target_name]
            except KeyError:
                output = default_value
        elif type(target_name) == list and len(target_name) == 2:
            try:
                output = parameter_dict[target_name[0]][target_name[1]]
            except KeyError:
                output = default_value
        else:
            raise Exception(f"Expecting {target_name} to be either a string or a two-element list.")
        return output
