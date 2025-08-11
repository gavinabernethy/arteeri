# Store library of species parameter descriptions that can then saved in the parameters script if needed
import numpy as np

# ----------------------------------------- DEFAULT ----------------------------------------- #

default_species = {
    "CORE_PARA": {
        "MINIMUM_POPULATION_SIZE": 1.0,
        "LIFESPAN": 1,
        "PREDATOR_LIST": [],  # What species are predators of this species?
        "SEASONAL_PERIOD": 0,  # used in the growth and direct impact offsets
    },

    "INITIAL_POPULATION_PARA": {
        "INITIAL_POPULATION_MECHANISM": "random_binomial",
        "IS_ENSURE_MINIMUM_POPULATION": False,  # note that this overwrites even an explicit patch_vector
        "CONSTANT_VALUE": None,
        "GAUSSIAN_MEAN": None,
        "GAUSSIAN_ST_DEV": None,
        "BINOMIAL_MAXIMUM_MULTIPLIER": 1.5,
        "BINOMIAL_PROBABILITY": 1.0,  # probability of occurrence in a given patch under any '_binomial' scheme
        "HABITAT_TYPE_NUM_BINOMIAL_DICT": {},  # probability of occurrence in patch of given habitat_type_num
        "PATCH_VECTOR": [],
    },

    "GROWTH_PARA":
        {
            "GROWTH_FUNCTION": "logistic",
            "R": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 5.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "RESOURCE_USAGE_CONVERSION": 1.0,  # in [0, 1] - how much resource do you use relative to other species?
            # Essentially this allows a conversion scale for different carrying capacities of species.
            "CARRYING_CAPACITY": 10.0,  # relative to a standard size and quality of habitat of a given type
            "ANNUAL_OFFSET": {
                "IS_GROWTH_OFFSET": False,  # is there an annual offset to early/late seasonal
                # behaviour, for example to delay mating season or late spring etc.
                "ANNUAL_DURATION": 0.0,  # should usually match the period of the R value. This is how long each year is
                "GROWTH_OFFSET_SPECIES": [],  # list - each entry is the annual offset. Can be stochastic!
                "IS_GROWTH_OFFSET_LOCAL": False,  # is there an annual offset that varies by patch?
                "GROWTH_OFFSET_LOCAL": [],  # list of lists - each entry is list of annual offsets per patch
            },
            "CML_PARA": [
                    # This should be an ordered list of the 1 (shift), 3 (tent) or 4 (sine) parameters used in these maps:
                    {"type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
            ],
        },
    "PREDATION_PARA":
        {
            "PREDATION_FUNCTION": "lotka_volterra",
            "PREY_DICT": {
                "type": 'constant',  # {'constant', 'vector_exp', 'vector_imp'}
                "constant_value": {},  # dictionary with name: [z-score, preference boost]
                "period": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
            },
            "IS_NONLOCAL_FORAGING": False,
            "MINIMUM_LINK_STRENGTH_FORAGING": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "IS_NONLOCAL_FORAGING_PATH_RESTRICTED": False,
            "MAX_FORAGING_PATH_LENGTH": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "FORAGING_MOBILITY": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "FORAGING_KAPPA": {  # should typically be zero, unless you want to EFFORTLESSLY forage over a range
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "PREDATION_RATE": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "PREDATION_PRAGMATISM": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0.5,  # in the range [0,+INF). Determines when preferences considered.
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "PREDATION_FOCUS_TYPE": "best_score",  # must be either "best_score" or "best_yield" - determines the
                # precise form of the effort function, so that rho -> +infty leads to either local (if non-zero prey)
                # foraging or focus on only the best return single prey population (from both score and population).
            "PREDATION_FOCUS": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 1.0,  # should be non-negative
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "ECOLOGICAL_EFFICIENCY": 0.0,  # fraction of killed/consumed prey biomass
            # that converts directly to new biomass of this species
            "IS_PREDATION_ONLY_PREVENTS_DEATH": False,  # cant gain new members due to pred. (only to r)
            # If false, this represents a long-lived species who need to eat to live, but they only reproduce and grow
            # the population when explicitly allowed to breed using R and the growth function.
        },
    "DISPERSAL_PARA":
        {
            "IS_DISPERSAL": True,  # This MUST NOT be temporally-varied.
            "DISPERSAL_MECHANISM": {
                # CHOOSE FROM: 'step_poly', 'diffusion', 'stochastic_quantity', 'stochastic_binomial', 'adaptive'
                "type": 'constant',  # {'constant', 'vector_exp', 'vector_imp'}
                "constant_value": "step_poly",
                "period": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
            },
            "ALWAYS_MOVE_WITH_MINIMUM": False,  # this should certainly be false if using stochastic_binomial
            "SS_DISPERSAL_PENALTY": 0.0,  # this is a fraction of movement that dies/never arrives. It can overwrite
            # the system-wide general value in pop_dyn_para but ONLY IF IT IS LARGER!
            #
            # Below is the cut-off, whilst dispersal mobility (or mechanism-specific equivalents) raise
            # the amounts going everywhere accessible:
            # large mobility / small threshold = much dispersal to far (and close) places
            # large mobility / high threshold = much dispersal only to close places
            # small mobility / small threshold = low level of dispersal to far (and close) places
            # small mobility / high threshold = low level of dispersal only to close places
            #
            # Be clear that this minimum link strength is about the perceived "closeness" or accessibility of the
            # location - it is NOT the minimum amount of the population required to register movement (this is just the
            # minimum population size).
            "MINIMUM_LINK_STRENGTH_DISPERSAL": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            # THIS IS REDUNDANT FOR STEP_POLY DISPERSAL IF CF_LISTS SCALED
            # Also, note that a typical movement_score is ~0.1 - 0.5 for adjacent patches.
            # Thus, for about 10% emigration, we want dispersal_mobility * mu_overall = 0.2
            "DISPERSAL_MOBILITY": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 1.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            # Is there a preference for the permitted direction of movement? 1.0 = ONLY move up, -1.0 = ONLY move down.
            "DISPERSAL_DIRECTION": {  # This should be a numerical value from [-1.0, 1.0]
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 0.0,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            # Do we insist that species may only migrate to immediate neighbouring patches in a given step?
            # This should always be TRUE and MUST NOT BE GREATER than the main_para["ASSUMED_MAX_PATH_LENGTH"]!!!
            # Otherwise, a removed/perturbed patch may yet be remembered in the species movement scores for distant
            # local populations who can continue to access it - their paths and scores will NOT have been updated as
            # the lower assumed_max_path_length will mean the removed patch wasn't recorded in this distant patch's
            # stepping_stone_lists and thus did not seem to need updated following the perturbation.
            "IS_DISPERSAL_PATH_RESTRICTED": True,  # Keep this TRUE, as explained above!
            "MAX_DISPERSAL_PATH_LENGTH": {
                "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                "constant_value": 1,
                "period": None,
                "amplitude": None,
                "phase_shift": None,
                "vertical_shift": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                "logistic_r": None,  # r-value of the logistic map to generate the time-series
                "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
            },
            "BINOMIAL_EXTRA_INDIVIDUAL": 0.5,
            # these coefficients controls WHEN the population leaves, and by how much, but not WHERE they go (near/far)!
            "COEFFICIENTS_LISTS": {
                "type": 'constant',  # {'constant', 'vector_exp', 'vector_imp'}
                "constant_value": {
                    "DENSITY_THRESHOLD": 1.0,  # uses UNDER if x/K <= this value, OVER otherwise
                    "UNDER": [0],
                    "OVER": [0, 1.0],
                },
                "period": None,
                "vector_exp": None,  # [value_0, value_1, ..., value_period]
                "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
            },
        },
    "PURE_DIRECT_IMPACT_PARA":
        {
            "IS_PURE_DIRECT_IMPACT": False,  # direct impact but not from any species. e.g. culling
            "TYPE": "vector",
            "IMPACT": 0.0,
            "PROBABILITY": 0.0,
            "DIRECT_VECTOR": [],
            "ANNUAL_OFFSET": {
                "IS_DIRECT_OFFSET": False,  # is there an annual offset to early/late seasonal
                # behaviour, for example to delay mating season or late spring etc.
                "ANNUAL_DURATION": 0.0,
                # should usually match the period of the R value (seasonal_duration). This is how long each year is
                "DIRECT_OFFSET_SPECIES": [],  # list - each entry is the annual offset. Can be stochastic!
                "IS_DIRECT_OFFSET_LOCAL": False,  # is there an annual offset that varies by patch?
                "DIRECT_OFFSET_LOCAL": [],
                # list of lists - each entry is list of offsets per each patch for that season
            },
        },
    "DIRECT_IMPACT_ON_ME": {},  # dictionary of species names (including self) and linear impact scores
    "PERTURBATION_PARA": {
        "IS_PERTURBS_ENVIRONMENT": False,  # does this species induce perturbations in the physical environment?
        "TO_IMPACT": [],  # list containing some of 'same', 'adjacent', 'xy-adjacent'
        "IMPLEMENTATION_PROBABILITY_COEFFICIENTS": {
            # for each potentially-impacted patch, occurs with probability = X_0*chi(x) + X_1*x + X_2*x^2 + X_3*x^3
            # dependent upon the density, that is x = local population / carrying_capacity
            "SAME": [0, 0, 0, 0],
            "ADJACENT": [0, 0, 0, 0],
            "XY_ADJACENT": [0, 0, 0, 0],
        },
        "PERTURBATION": {
            "IS_REMOVAL": False,
            "IS_HABITAT_TYPE_CHANGE": False,
            "HABITAT_TYPE_NUM_TO_CHANGE_TO": None,  # integer habitat type number
            "IS_QUALITY_CHANGE": False,
            "RELATIVE_QUALITY_CHANGE": None,  # ± float amount?
            "IS_ADJACENCY_CHANGE": False,
            "ABSOLUTE_ADJACENCY_CHANGE": None,  # 1 or 0
        },
    },

    # # template for a parameter that may vary over time
    # "temporally_varying_parameter":
    #     {
    #         "type": None,  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
    #         "constant_value": None,
    #         "period": None,
    #         "amplitude": None,
    #         "phase_shift": None,
    #         "vertical_shift": None,
    #         "vector_exp": None,  # [value_0, value_1, ..., value_period]
    #         "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
    #               NOTE DICT. UNORDERED! The keys indicate when a given behaviour begins in the cycle [0, period)
    #          "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
    #          "logistic_r": None,  # r-value of the logistic map to generate the time-series
    #          "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
    #     }
}

# --------------------------- ARTEMIS SAMPLE (PREY AND PREDATOR) ---------------------------- #

ARTEMIS_SAMPLE_MASTER = {
    "prey": {
        "CORE_PARA":{
            "MINIMUM_POPULATION_SIZE": 0.000001,
            "LIFESPAN": 100,
            "PREDATOR_LIST": ['predator'],
            "SEASONAL_PERIOD": 0,
        },
        "INITIAL_POPULATION_PARA": {
            "INITIAL_POPULATION_MECHANISM": "gaussian",
            "IS_ENSURE_MINIMUM_POPULATION": True,
            "CONSTANT_VALUE": None,
            "GAUSSIAN_MEAN": 0.1,
            "GAUSSIAN_ST_DEV": 0.01,
            "BINOMIAL_MAXIMUM_MULTIPLIER": None,
            "BINOMIAL_PROBABILITY": None,
            "HABITAT_TYPE_NUM_BINOMIAL_DICT": {},
            "PATCH_VECTOR": None,
        },
        "GROWTH_PARA":
            {
                "GROWTH_FUNCTION": "logistic",
                "R": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 4.0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "RESOURCE_USAGE_CONVERSION": 0.5,  # in [0, 1] - how much resource do you use relative to other species?
                # Essentially this allows a conversion scale for different carrying capacities of species.
                "CARRYING_CAPACITY": 1.0,
                "ANNUAL_OFFSET": {
                    "IS_GROWTH_OFFSET": False,  # is there an annual offset to early/late seasonal
                    # behaviour, for example to delay mating season or late spring etc.
                    "ANNUAL_DURATION": 0.0,  # This is how long each year is
                    "GROWTH_OFFSET_SPECIES": [],  # list - each entry is the annual offset. Can be stochastic!
                    "IS_GROWTH_OFFSET_LOCAL": False,  # is there an annual offset that varies by patch?
                    "GROWTH_OFFSET_LOCAL": [],  # list of lists - each entry is list of annual offsets per patch
                },
                "CML_PARA": [
                    # This should be an ordered list of the 1 (shift), 3 (tent) or 4 (sine) parameters used in these maps:
                    {"type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                    },
                ],
            },
        "PREDATION_PARA": None,
        "DISPERSAL_PARA":
            {
                "IS_DISPERSAL": True,
                "DISPERSAL_MECHANISM": {
                    "type": 'constant',  # {'constant', 'vector_exp', 'vector_imp'}
                    "constant_value": "diffusion",
                    "period": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                },
                "ALWAYS_MOVE_WITH_MINIMUM": False,  # this should certainly be false if using stochastic_binomial
                "SS_DISPERSAL_PENALTY": 0.0,  # this is a fraction of movement that dies/never arrives. It can
                # overwrite the system-wide general value in pop_dyn_para but ONLY IF IT IS LARGER!
                "MINIMUM_LINK_STRENGTH_DISPERSAL": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.01,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "DISPERSAL_MOBILITY": {
                    # THIS IS REDUNDANT FOR STEP_POLY DISPERSAL IF CF_LISTS SCALED
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.025,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                # Is there a preference for permitted movement direction? 1.0 = ONLY move up, -1.0 = ONLY move down.
                "DISPERSAL_DIRECTION": {  # This should be a numerical value from [-1.0, 1.0]
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "IS_DISPERSAL_PATH_RESTRICTED": True,
                "MAX_DISPERSAL_PATH_LENGTH": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 1,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "BINOMIAL_EXTRA_INDIVIDUAL": 0.0,
                "COEFFICIENTS_LISTS": {
                    "type": None,  # {'constant', 'vector_exp', 'vector_imp'}
                    "period": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                },
            },
        "PURE_DIRECT_IMPACT_PARA": None,
        "DIRECT_IMPACT_ON_ME": {},  # dictionary of species names (including self) and linear impact scores
        "PERTURBATION_PARA": None,
    },

    "predator": {
        "CORE_PARA":{
            "MINIMUM_POPULATION_SIZE": 0.0001,
            "LIFESPAN": 100,
            "PREDATOR_LIST": [],
            "SEASONAL_PERIOD": 0,
        },
        "INITIAL_POPULATION_PARA": {
            "INITIAL_POPULATION_MECHANISM": "gaussian",
            "IS_ENSURE_MINIMUM_POPULATION": True,  # note that this applies even for an explicit patch_vector
            "CONSTANT_VALUE": None,
            "GAUSSIAN_MEAN": 0.01,
            "GAUSSIAN_ST_DEV": 0.001,
            "BINOMIAL_MAXIMUM_MULTIPLIER": None,
            "BINOMIAL_PROBABILITY": None,
            "HABITAT_TYPE_NUM_BINOMIAL_DICT": {},
            "PATCH_VECTOR": None,
        },
        "GROWTH_PARA":
            {
                "GROWTH_FUNCTION": "logistic",
                "R": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.5,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "RESOURCE_USAGE_CONVERSION": 1.0,  # in [0, 1] - how much resource do you use relative to other species?
                # Essentially this allows a conversion scale for different carrying capacities of species.
                "CARRYING_CAPACITY": 1.0,
                "ANNUAL_OFFSET": {
                    "IS_GROWTH_OFFSET": False,  # is there an annual offset to early/late seasonal
                    # behaviour, for example to delay mating season or late spring etc.
                    "ANNUAL_DURATION": None,  # This is how long each year is
                    "GROWTH_OFFSET_SPECIES": [],  # list - each entry is the annual offset. Can be stochastic!
                    "IS_GROWTH_OFFSET_LOCAL": False,  # is there an annual offset that varies by patch?
                    "GROWTH_OFFSET_LOCAL": [],  # list of lists - each entry is list of annual offsets per patch
                },
                "CML_PARA": [
                    # This should be an ordered list of the 1 (shift), 3 (tent) or 4 (sine) parameters used in these maps:
                    {"type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                    },
                ],
            },
        "PREDATION_PARA":
            {
                "PREDATION_FUNCTION": "lotka_volterra",
                "PREY_DICT": {
                    "type": 'constant',  # {'constant', 'vector_exp', 'vector_imp'}
                    "constant_value": {'prey': [1.0, 0.0]},  # dictionary with name: [z-score, preference boost]
                    "period": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": {},  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                },
                "IS_NONLOCAL_FORAGING": True,
                "MINIMUM_LINK_STRENGTH_FORAGING": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.0001,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "IS_NONLOCAL_FORAGING_PATH_RESTRICTED": True,
                "MAX_FORAGING_PATH_LENGTH": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 2,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "FORAGING_MOBILITY": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 1.0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "FORAGING_KAPPA": {  # should typically be zero, unless you want to EFFORTLESSLY forage over a range
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 1.0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "PREDATION_RATE": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 10.0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "PREDATION_PRAGMATISM": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 1.0,  # in the range [0,+INF). Determines when preferences considered.
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "PREDATION_FOCUS_TYPE": "best_score",  # must be either "best_score" or "best_yield" - determines the
                  # precise form of the effort function, so that rho -> +infty leads to either local (if non-zero prey)
                  # foraging or focus on only the best return single prey population (from both score and population).
                "PREDATION_FOCUS": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 1.0,  # should be non-negative
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "ECOLOGICAL_EFFICIENCY": 0.3,
                "IS_PREDATION_ONLY_PREVENTS_DEATH": False,  # cant gain new members due to pred. (only to r)
            },
        "DISPERSAL_PARA":
            {
                "IS_DISPERSAL": True,
                "DISPERSAL_MECHANISM": {
                    "type": 'constant',  # {'constant', 'vector_exp', 'vector_imp'}
                    "constant_value": "diffusion",
                    "period": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                },
                "ALWAYS_MOVE_WITH_MINIMUM": False,  # this should certainly be false if using stochastic_binomial
                "SS_DISPERSAL_PENALTY": 0.0,  # this is a fraction of movement that dies/never arrives.
                # It can overwrite the system-wide general value in pop_dyn_para but ONLY IF IT IS LARGER!
                "MINIMUM_LINK_STRENGTH_DISPERSAL": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.01,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "DISPERSAL_MOBILITY": {
                    # THIS IS REDUNDANT FOR STEP_POLY DISPERSAL IF CF_LISTS SCALED
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.025,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                # Is there a preference for permitted movement direction? 1.0 = ONLY move up, -1.0 = ONLY move down.
                "DISPERSAL_DIRECTION": {  # This should be a numerical value from [-1.0, 1.0]
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 0.0,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "IS_DISPERSAL_PATH_RESTRICTED": True,
                "MAX_DISPERSAL_PATH_LENGTH": {
                    "type": 'constant',  # {'constant', 'sine', 'vector_exp', 'vector_imp', 'logistic_map'}
                    "constant_value": 1,
                    "period": None,
                    "amplitude": None,
                    "phase_shift": None,
                    "vertical_shift": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                    "logistic_initial": None,  # initial value of this parameter (WILL BE RE-SCALED BY THE MAXIMUM)
                    "logistic_r": None,  # r-value of the logistic map to generate the time-series
                    "logistic_max": None,  # theoretical maximum value of this parameter, to re-scale all to [0, 1]
                },
                "BINOMIAL_EXTRA_INDIVIDUAL": 0.0,
                "COEFFICIENTS_LISTS": {
                    "type": None,  # {'constant', 'vector_exp', 'vector_imp'}
                    "constant_value": None,
                    "period": None,
                    "vector_exp": None,  # [value_0, value_1, ..., value_period]
                    "vector_imp": None,  # { 0 : value_0, ... , lower_time_limit_N : value_N }
                },
            },
        "IS_PURE_DIRECT_IMPACT": False,  # direct impact but not from any species
        "PURE_DIRECT_IMPACT_PARA":
            {
                "TYPE": None,
                "IMPACT": None,
                "PROBABILITY": None,
                "DIRECT_VECTOR": [],
                "ANNUAL_OFFSET": {
                    "IS_DIRECT_OFFSET": False,  # is there an annual offset to early/late seasonal
                    # behaviour, for example to delay mating season or late spring etc.
                    "ANNUAL_DURATION": None,  # This is how long each year is
                    "DIRECT_OFFSET_SPECIES": [],  # list - each entry is the annual offset. Can be stochastic!
                    "IS_DIRECT_OFFSET_LOCAL": False,  # is there an annual offset that varies by patch?
                    "DIRECT_OFFSET_LOCAL": [],  # list of lists - each entry is list of annual offsets per patch
                },
            },
        "DIRECT_IMPACT_ON_ME": {},  # dictionary of species names (including self) and linear impact scores
        "IS_PERTURBS_ENVIRONMENT": False,  # does this species induce perturbations in the physical environment?
        "PERTURBATION_PARA": {
            "TO_IMPACT": [],  # list containing some of 'same', 'adjacent', 'xy-adjacent'
            "IMPLEMENTATION_PROBABILITY_COEFFICIENTS": {
                # for each potentially-impacted patch, occurs with probability = X_0*chi(x) + X_1*x + X_2*x^2 + X_3*x^3
                # dependent upon the density, that is x = local population / carrying_capacity
                "SAME": [0, 0, 0, 0],
                "ADJACENT": [0, 0, 0, 0],
                "XY_ADJACENT": [0, 0, 0, 0],
            },
            "PERTURBATION": {
                "IS_REMOVAL": False,
                "IS_HABITAT_TYPE_CHANGE": False,
                "HABITAT_TYPE_NUM_TO_CHANGE_TO": None,  # integer habitat type number
                "IS_QUALITY_CHANGE": False,
                "RELATIVE_QUALITY_CHANGE": None,  # ± float amount?
                "IS_ADJACENCY_CHANGE": False,
                "ABSOLUTE_ADJACENCY_CHANGE": None,  # 1 or 0
            },
        },
    },
}
