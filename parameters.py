from parameters_species_repository import *

meta_para = {
    "IS_NEW_PROGRAM": True,
    "REPEAT_PROGRAM_CODE": None,  # what is the simulation number to be repeated?
    "REPEAT_PROGRAM_PATH": None,  # what is the output path of the simulation to be repeated? We need to know where to
    # find their meta_data and parameter files in the results/sub_folder/folder structure.
    "NUM_REPEATS": 1,  # how many simulations should be executed with the current parameter set?
    "IS_RUN_SAMPLE_SPATIAL_DATA_FIRST": True,  # should we execute sample_spatial_data() before running the batch set?
    # if false then we will try to load the SPATIAL_TEST_SET below. So if you want to do several batches with the same
    # spatial set then generate it separately by executing sample_spatial_data.py then run the batches with this FALSE.
    # Note that if TRUE then any existing data will be overwritten and new adjacency sets will need to be generated
    # regardless of the save/load variables parameters.
}

master_para = {
"main_para":
        {
            # ---------- CONTROL PARAMETERS: ---------- #
            "IS_SIMULATION": True,  # if False then we init Simulation_obj but do not execute .full_simulation()
            "NUM_TRANSIENT_STEPS": 1000,
            "NUM_RECORD_STEPS": 1000,
            "NUM_PATCHES": 20,
            # ----------------------------------------- #

            "MODEL_TIME_TYPE": "discrete",  # continuous ODEs ('continuous') or discrete maps ('discrete')?
            "EULER_STEP": 0.1,  # ONLY used if continuous - solve ODEs by Euler method
            "STEPS_TO_DAYS": 1,  # be aware that this affects how often temporal functions are updated!

            "ECO_PRIORITIES": {0: {'foraging', 'direct_impact', 'growth', 'dispersal'}, 1: {}, 2: {}, 3: {}},
            # What is the order or concurrence in which foraging (predation), growth (reproduction and mortality),
            # direct impact, and dispersal should be resolved?

            "MAX_CENTRALITY_MEASURE": 10,  # Max amount of NxN matrix multiplication when determining patch.centrality
            "ASSUMED_MAX_PATH_LENGTH": 3,  # used for shortcuts in rebuilding paths AND multiplying adjacency matrix!
            # THIS VALUE NEEDS TO BE AT LEAST EQUAL TO THE MAXIMUM MAX_DISPERSAL_PATH_LENGTH ACROSS ALL SPECIES!!!
            # and IF YOU CHANGE THIS then "IS_LOAD_ADJ_VARIABLES" BELOW MUST BE "FALSE" AS WE NEED TO REBUILD THEM!!!
            #
            # Note that saving the adjacency variables does seem to be extremely slow in DEBUG mode.
            "IS_SAVE_ADJ_VARIABLES": False,  # Save patch.stepping_stone_list,.species_movement_scores,.adjacency_lists?
            "IS_LOAD_ADJ_VARIABLES": False,  # Load patch.stepping_stone_list,.species_movement_scores,.adjacency_lists?

            # ------------- Generation data - needs to be set before spatial habitat generation ------------- #
            "SPECIES_TYPES": {
                # This is the broader dictionary of what types of species are "known" to the system, even if they do
                # not actually feature in the simulation at the beginning (or at all!).
                # Key (index) must be non-negative integers without gaps, and the value must be a recognised species
                # name loaded (below) from the parameters_species_repository.py:
                0: "prey",
                1: "predator",
            },  # key numbering must remain consistent with column ordering of the loaded arrays

            "HABITAT_TYPES": {
                # Key (indexing) must be non-negative integers without gaps. Value can be any given name.
                0: 'habitat_type_0',
                1: 'habitat_type_1',
            },
            "GENERATED_SPEC": {
                #
                # NOTE THAT IF YOU WISH TO CHANGE THESE YOU MUST RE-RUN SPATIAL HABITAT GENERATION BEFORE MAIN.PY
                #
                "FEEDING": {
                    "IS_SPECIES_SCORES_SPECIFIED": True,  # if false, then randomly generated from the uniform
                    # distribution over [MIN_SCORE, MAX_SCORE]
                    "MIN_SCORE": None,
                    "MAX_SCORE": None,
                    # specify habitat scores for generation
                    # If used, this needs to have keys from 0, ...,  total_possible_habitats, indexing lists with
                    # length equal to the total possible number of scores (i.e. the number of species)
                    "HABITAT_SCORES": {0: [1.0, 1.0],
                                       1: [0.0, 1.0],
                                       },
                },
                "TRAVERSAL": {
                    "IS_SPECIES_SCORES_SPECIFIED": True,  # if false, then randomly generated from the uniform
                    # distribution over [MIN_SCORE, MAX_SCORE]
                    "MIN_SCORE": None,
                    "MAX_SCORE": None,
                    "HABITAT_SCORES": {0: [1.0, 1.0],
                                       1: [1.0, 1.0],
                                       },
                },
            },

            # --- Initial choice of the above species and habitats that should actually be present at the start --- #

            # Now indicate which species from the "SPECIES_TYPES" keys that you want to be present at the beginning.
            "INITIAL_SPECIES_SET": {0, 1},  # each must exist as a key in the types dictionary, ordering not needed

            # each must be present in the types dictionary, ordering not needed
            # THIS ALSO NEEDS TO BE SET BEFORE SPATIAL HABITAT GENERATION!
            "INITIAL_HABITAT_SET": {0, 1},
            # if the following is None then probabilities are treated as uniform when combined with auto-correlation
            "INITIAL_HABITAT_BASE_PROBABILITIES": None,

            # do we load the hurst module and attempt to calculate Hurst exponents?
            "IS_CALCULATE_HURST": False,

            # When conducting distance metric, network and complexity analyses that include linear regressions, do we
            # record the vectors of values, to reconstruct the raw data scatter plots against the fitted models later?
            "IS_RECORD_METRICS_LM_VECTORS": False,

            # Parameters for conducting the complexity analysis - increasing these can be computationally-expensive.
            "COMPLEXITY_ANALYSIS": {
                "NUM_CLUSTER_DRAWS": 100,  # how many samples do we try to draw for each delta?
                "NUM_CLUSTER_DRAW_ATTEMPTS": 20,  # how many attempts to draw each sample (may fail if disconnected)?
                "MAX_DELTA": 64,  # maximum delta is min(num_patches_in_subnetwork/2 , this_value)
                "IS_PARTITION_ANALYSIS": True,  # toggle partition analysis to determine resolution for max complexity
                "NUM_EVO_PARTITIONS": 10,
                "NUM_REG_PARTITIONS": 10,
                "PARTITION_SUCCESS_THRESHOLD": 0.8,  # what fraction (of the max possible) patches must be in clusters
                    # of precisely the required size, for the partition to be considered a success?
            },
        },
    "graph_para":
        {
            # This set of parameters is accessed during the execution of sample_spatial_data.py for generating
            # the spatial network files that are then used in the simulation.
            "SPATIAL_TEST_SET": 1,
            "SPATIAL_DESCRIPTION": "artemis_01",
            #
            # Network topology:
            "GRAPH_TYPE": "ring", # choices are: "manual", "lattice", "line", "ring",
            #  "star", "random", "small_world", "scale_free", "cluster", "balanced_tree", "power_law_tree"
            #
            # Drawing the visual layout of the network - the only simulation impact that take any account of the
            # patch.position are "IS_LATTICE_WRAPPED" and the cluster type "position_box".
            "GRAPH_LAYOUT": "spiral",  # choices are "grid", "tree", "space_filling_curve, "spiral"

            "ADJACENCY_MANUAL_SPEC": None,  # should be None if we want to generate the patch adjacency matrix by
            # other means, and a list (length = num_patches) of lists (length = num_patches) if we want to use it
            "LATTICE_GRAPH_CONNECTIVITY": 1.0,
            "IS_LATTICE_INCLUDE_DIAGONALS": False,
            "IS_LATTICE_WRAPPED": False,  # should only be used for GRAPH_TYPE 'lattice' and GRAPH_LAYOUT 'grid'.
            "RANDOM_GRAPH_CONNECTIVITY": None,
            "SMALL_WORLD_NUM_NEIGHBOURS": None,  # if num_patches > 2, this need to be at least 2 or graph fails
            "SMALL_WORLD_SHORTCUT_PROBABILITY": None,
            "CLUSTER_NUM_NEIGHBOURS": None,
            "CLUSTER_PROBABILITY": None,
            "TREE_BRANCHING": 3,  # for balanced_tree
            "TREE_POWER": None,  # for power_law_tree
            #
            # Cliquey Network options:
            "WITHIN_CLIQUE_PROBABILITY": 0.8,  # for cliquey networks
            "BETWEEN_CLIQUE_PROBABILITY": 0.01,  # for cliquey networks
            "NUMBER_OF_CLIQUES": 5,  # for cliquey networks
            "IS_HABITAT_CLIQUE": True,  # if true then cliques will have uniform habitat type
            #
            # Determine habitat type:
            "IS_HABITAT_PROBABILITY_REBALANCED": True,  # are habitat probabilities sequentially biased to recover?
            "HABITAT_TYPE_MANUAL_ALL_SPEC": None, # should be None if we want to generate habitats by probability,
                # otherwise a list of the habitat nums
            "HABITAT_SPATIAL_AUTO_CORRELATION": 0.8,  # in range [-1, 1]
            "HABITAT_TYPE_MANUAL_OVERWRITE": None,  # set this to None or empty dict, unless you want
            # to manually specify the habitat types of only certain patches in an otherwise randomly-generated system.
            # If you want to specify ALL patches then use the MANUAL_ALL_SPEC option instead.
            "IS_HABITAT_CLUSTERS": False,  # if True then all else ignored except the following two options...
            "HABITAT_CLUSTER_SIZE": [],  # if IS_HABITAT_CLUSTERS, then what size of clusters should be drawn?
                # This should be either a single value for all clusters to be the same, or a list of the desired cluster
                # sizes to be sequentially alternated through. If IS_BIND_HABITAT_TO_CLUSTER_SIZE then this needs to be
                # the same length as the number of habitats at generation. Otherwise, this is not bound to specific
                # habitat types, so simply list the range of same-habitat sizes desired and each cluster is assigned
                # a habitat type to minimise overlapping boundaries.
            "HABITAT_CLUSTER_TYPE_STR": None,  # cluster topology (box, star, chain, random, disconnected).
            # But also select "chess_box" to use "box" but with constraints to resemble a chess board. Note that this
            # has limited effectiveness if used in conjunction with different cluster sizes, so most effective if
            # content with a single cluster size across all habitat types.
            "IS_CHESS_BIND_HABITAT_TO_SIZE": False,  # if TRUE then HABITAT_CLUSTER_TYPE_STR needs to be "chess_box" and
            # HABITAT_CLUSTER_SIZE needs to be a list of length equal to the size of "main_para"["INITIAL_HABITAT_SET"].
            # Set to FALSE if not a chess-board design or just one cluster size chosen, or is a chessboard with multiple
            # cluster sizes but you do not want to swap habitat types to minimise overlap.
            #
            # Patch size (scales the carrying capacity for all local populations):
            "PATCH_SIZE_MANUAL_SPEC": None,  # should be None if we want to generate size by probability
            "MIN_SIZE": 1.0,
            "MAX_SIZE": 1.0,
            # Patch quality (scales the reproductive rate (r-parameter) for all local populations):
            "QUALITY_TYPE": "gradient",  # quality types are: 'manual', 'random', 'auto_correlation', 'gradient'
            "QUALITY_MANUAL_SPEC": None,  # should be None if we want to generate quality by other means
            "QUALITY_SPATIAL_AUTO_CORRELATION": 0.5,  # in range [-1, 1]
            "MIN_QUALITY": 1.0,
            "MAX_QUALITY": 1.0,
            "QUALITY_FLUCTUATION": 0.0,
            "QUALITY_AXIS": "x+y",  # 'x' or 'y' or 'x+y'
            "IS_ENVIRONMENT_NATURAL_RESTORATION": False,  # is there a tendency to return to certain characteristics?
            "RESTORATION_PARA": {
                "IS_QUALITY_CHANGE": False,
                "QUALITY_CHANGE_PROBABILITY": 0.0,
                "QUALITY_CHANGE_SCALE": 0.0,  # should be in range (0, 1] - how much do we move towards the target?
                "QUALITY_DESIRED": 0.0,
                "IS_HABITAT_CHANGE": False,
                "HABITAT_CHANGE_PROBABILITY": None,
                "HABITAT_TYPE_NUM_DESIRED": None,
            },
        },
    "plot_save_para":
        {
            "IS_ALLOW_FILE_CREATION": True,  # prevents creation of any files for running on remote clusters
            "IS_SUB_FOLDERS": False,  # do we nest the results folders in sub-folders?
            "SUB_FOLDER_CAPACITY": 3,  # if so, what should the maximum capacity of these be?
            "IS_PRINT_KEY_OUTPUTS_TO_CONSOLE": True,  # prints final and average local populations to console
            "IS_PRINT_DISTANCE_METRICS_TO_CONSOLE": True,  # JSON.dumps() of species and community distribution analysis
            "IS_SAVE": True,  # do you save ANY data files?
            "IS_PLOT": True,  # do you plot ANY final graphs? Must be enabled to save any subsets controlled below.
            "IS_PLOT_DISTANCE_METRICS_LM": False,  # Do you plot the distance-metric linear models (with
            # scatter plots of the base data, IF COLLECTED - see option in main_para)?
            "IS_PLOT_COMPLEXITY": True,  # Do you plot the complexity power laws?
            "IS_RECORD_AND_PLOT_LESSER_LM": False,  # if True, then also attempt to fit, store, and then plot _base,
            # _nz, and shifted (less reliable) fitted linear models for inter_species_predictions.
            "PLOT_INIT_NETWORK": True,  # do we plot the initial network before the simulation (before species pathing)?
            "MANUAL_SPATIAL_NETWORK_SAVE_STEPS": [],  # LIST of integer steps during which to plot the spatial network:
            # - include 0 to plot early state of the network (AFTER first step 0 iterates) before patch perturbations;
            # - include -1 to plot the initialised system before ANY steps or perturbations executed whatsoever,
            # including the initial distributions of the species populations.
            #
            # Data control options (requires IS_SAVE to be true):
            "IS_SAVE_LOCAL_POP_HISTORY_CSV": False,  # produce individual .csv file with only the core time series
            # (population size, internal change, dispersal in and out) for each local_pop object.
            "IS_SAVE_SYSTEM_STATE_DATA": True,  # produce JSON of system state, including, for example, histories of
            # perturbation, biodiversity, and the mean and s.d. of quality and size of patches present at that time in
            # the network, and the full set of network distance-metric analysis.
            "IS_SAVE_PATCH_DATA": False,  # produce JSON file of every patch object (including, for example,
            # the full history per-patch of different local clustering values).
            "IS_SAVE_PATCH_LOCAL_POP_DATA": False,  # produce a JSON file of every local_population object - however
            # note that this requires IS_SAVE_PATCH_DATA to be true first.
            "IS_ODE_RECORDINGS": False,  # do we save the history of each iteration of the ODE details as an attribute
            # of each local population object (it would then be printed as part of IS_SAVE_PATCH_LOCAL_POP_DATA)?
            # This is memory-intensive and mainly intended for debugging.
            "IS_SAVE_DISTANCE_METRICS": True,  # produce JSON of species and community distribution analysis.
            "IS_PICKLE_SAVE": False,  # save the Python objects.
            "IS_SAVE_CURRENT_MOVE_SCORES": False,  # writes the final movement scores to the simulation-specific folder.
            #
            # Plot control options (requires IS_PLOT to be true):
            "IS_PLOT_LOCAL_TIME_SERIES": False,  # produce plots of all local time-series of species properties, imposed
            # on the same figure.
            "LOCAL_PLOTS": False,  # produce separate plot files for each patch showing the time-series.
            "IS_PLOT_ACCESSIBLE_SUB_GRAPHS": False,  # patch plots showing the fully-connected network sub-graphs from
            # the POV of each species (can be memory intensive and cause crashes if repeated too often).
            "IS_PLOT_ADJACENCY_SUB_GRAPHS": False,  # plots the undirected adjacency-based sub-graphs, regardless of any
            # species actual ability to traverse them.
            "IS_PLOT_INTERACTIONS": False,  # this draws at least num_species * num_patches plots
            "IS_PLOT_DEGREE_DISTRIBUTION": False,
            "IS_PLOT_UNRESTRICTED_PATHS": False,  # # plots all reachable foraging/direct dispersal paths per species
            # WITHOUT threshold and path-length restriction (crowded plot if species able to traverse most habitats).
            "IS_BIODIVERSITY_ANALYSIS": False,  # produce a species-area curve
            "MIN_BIODIVERSITY_ATTEMPTS": 100,  # used in the estimation of the SAR for very high numbers of patches
        },
    "pop_dyn_para":
        {
            "MU_OVERALL": 1.0,  # scales dispersal for all species and movements
            "GENERAL_DISPERSAL_PENALTY": 0.1,  # in [0, 1], this sets a baseline fractional LOSS for all movement,
            # and is overwritten ONLY if the species-specific cost is HIGHER.
            "COMPETITION_ALPHA_SCALING": 0.25,  # Sets relative strength [0 - 1] of inter-specific competition
            # to intra-specific competition for natural abiotic resources.
            "IS_LOCAL_FORAGING_ENSURED": False,  # if true, all feeding scores within-patch set to 1 (NOT m_i * h_i,k)
            #
            # These act like safety valves - overriding species-specific options if set to False to turn off such acts.
            "IS_NONLOCAL_FORAGING_PERMITTED": True,
            "IS_DISPERSAL_PERMITTED": True,
            "IS_PURE_DIRECT_IMPACT": False,  # allows custom ±c*x impact on a species for harvesting or culling etc.
            # Can be specified to occur on a periodic cycle, and the cycle can be offset on certain years or habitats.
            "IS_DIRECT_IMPACT": False,  # allows custom ±c*x*y impacts between species to be specified, inc. mutualism
            "IS_DIRECT_IMPACT_NONLOCAL": False,  # occurs between ALL interacting_populations, not just locally.
            "IS_SPECIES_PERTURBS_ENVIRONMENT": False,  # determines whether species permitted to induce perturbations.
        },
    "perturbation_para":
        {
            # a list of integers from [0 : num_patches - 1] of patches designated reserves and cannot be perturbed:
            "IS_RESERVE_LOAD": False,  # Want to load a previously-drawn reserve set for this spatial network?
            "IS_RESERVE_SAVE": False,  # Want to save this reserve set for the spatial network (overwrites previous)?
            "RESERVE_PATCH_CLUSTERS": [[]],  # default is empty list of lists (per cluster)
            "IS_RESERVE_GENERATION": False,
            "CLUSTERS_MUST_BE_SEPARATED": False,
            "RESERVE_CLUSTERS": None,
            "IS_PLOTS": False,
            "IS_OUTPUT_DATAFILES": False,
            # NO MORE THAN ONE PERTURBATION ARCHETYPE SHOULD BE ENACTED IN A GIVEN TIME-STEP - AS THEY WILL NOT BE
            # SIMULTANEOUS AND AS SYSTEM_STATE.PERTURBATION_HISTORY IS A DICTIONARY AND ONLY HOLDS ONE VALUE PER STEP!
            "PERT_STEP_DICTIONARY": {13: "a"},  # {100 * x + 50: 'b' for x in range(100)}
            # e.g. {x: 'a' for x in range(2000)} so that pert type 'a' occurs at step x
            "PERT_ARCHETYPE_DICTIONARY": {
                "a": {
                    "perturbation_type": "patch_perturbation",
                    "perturbation_subtype": "change_parameter",
                    # "change_habitat" or "change_parameter" or "remove_patch" or ""change_adjacency"
                    "patch_list_overwrite": None,  # list if desired
                    "patches_affected": [{"num_patches": 1,
                                          "habitat_nums_permitted": None,  # if None then all by default
                                          "initial": ["random"],
                                          "arch_type": "box",
                                          }],
                    "is_pairs": False,
                    "habitat_nums_to_change_to": None,
                    "parameter_change": 0.5,  # numerical value
                    "parameter_change_type": "relative_multiply",  # 'absolute', 'relative_add', or 'relative_multiply'
                    "parameter_change_attr": "quality",  # 'quality' or 'size'
                    "adjacency_change": None,
                    "is_reserves_overwrite": False,
                    "clusters_must_be_separated": True,
                    "proximity_to_previous": 0,  # 0, 1, 2, 3 - will not apply to removal perturbations
                    # [relative_weight, alpha, beta, gamma] or None - preference of total distance from last pert.
                    "prev_weighting": [1, 1, 5, 0],
                    # [relative_weight, alpha, beta, gamma] or None - preference function of previously pert. patches
                    "all_weighting": None,
                    "rebuild_all_patches": False,
                    "contagion_probability": 0.95,  # applies to patch perturbation only - either float or
                    # list of length equal to the number of habitats
                    "contagion_delay": 10,  # applies to patch perturbation only. Must be >0 for pert to take effect.
                    "contagion_cooldown": 20,  # applies to patch perturbation only. Prevents immediate reinfection if >=2*delay.
                },
            },  # see comment examples at base of script
        },
    # ------------------------------------- COPY FROM SPECIES REPOSITORY ------------------------------------- #
    # convert these to attributes of the species class, rather than storing them in the parameters
    # species should then not be named directly, but just stored in a dictionary with the names given here
    #
    # editing should be done in the repository - create new copies there if needed with unique names for an experiment,
    # so that everything is saved and not overwritten
    "species_para":
        {
            x: ARTEMIS_SAMPLE_MASTER[x] for x in ARTEMIS_SAMPLE_MASTER
        },
}

# "PERT_ARCHETYPE_DICTIONARY": {
#     "a": {
#         "perturbation_type": "patch_perturbation",
#         "perturbation_subtype": "change_parameter",
#         # "change_habitat" or "change_parameter" or "remove_patch" or ""change_adjacency"
#         "patch_list_overwrite": None,  # list if desired
#         "patches_affected": [{"num_patches": 2,
#                               "habitat_nums_permitted": None,  # if None then all by default
#                               "initial": ["random"],
#                               "arch_type": "box",
#                               }],
#         "is_pairs": False,
#         "habitat_nums_to_change_to": None,
#         "parameter_change": 0.5,  # numerical value
#         "parameter_change_type": "relative_multiply",  # 'absolute', 'relative_add', or 'relative_multiply'
#         "parameter_change_attr": "quality",  # 'quality' or 'size'
#         "adjacency_change": None,
#         "is_reserves_overwrite": False,
#         "clusters_must_be_separated": True,
#         "proximity_to_previous": 0,  # 0, 1, 2, 3 - will not apply to removal perturbations
#         # [relative_weight, alpha, beta, gamma] or None - preference of total distance from last pert.
#         "prev_weighting": [1, 1, 5, 0],
#         # [relative_weight, alpha, beta, gamma] or None - preference function of previously pert. patches
#         "all_weighting": None,
#         "is_restoration": False,
#         "rebuild_all_patches": False,
#         "contagion_probability": 0.0,  # applies to patch perturbation only - either float or list of length
#               # equal to num habitats
#         "contagion_delay": 0,  # applies to patch perturbation only. Must be >0 for pert to take effect.
#         "contagion_cooldown": 0,  # applies to patch perturbation only. Prevents immediate reinfection if >=2*delay.
#     },
# },
#
#
#     "b": {
#         "perturbation_type": "population_perturbation",
#         "perturbation_subtype": "extinction",  # 'dispersal' or 'extinction' or 'displacement'
#         "probability": 1.0,
#         "fraction_of_population": 0.2,
#         "patch_list_overwrite": None,
#         "is_reserves_overwrite": False,
#         "species_affected": None,  # if None (list of names of species) then default is all
#         "all_patches_habitats_affected": None,  # if None (list of numbers of habitat types) then default
#         # is all - but note that this option is only in the case that all patches are being selected from.
#         "patches_affected": None,  # if None (list of patch numbers) default is all except reserves
#         "clusters_must_be_separated": False,
#         "proximity_to_previous": 0,
#         # [relative_weight, alpha, beta, gamma] or None - preference of total distance from last pert.
#         "prev_weighting": None,
#         # [relative_weight, alpha, beta, gamma] or None - preference function of previously pert. patches
#         "all_weighting": [1, 1, 3, 0],
#         "is_restoration": False,
#     },
