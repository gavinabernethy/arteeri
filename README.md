# A.R.T.E.E.R.I.
## Arteeri spatial population dynamics tool (formerly 'Artemis').

Dr Gavin M Abernethy,
University of Stirling.

DOI: 10.5281/zenodo.16794869

## What is this software for?

A.R.T.E.E.R.I. (A Ranged Terrestrial, Ecological and Evolutionary model with Reserves Iterator) is designed to be a
general-purpose Python tool enabling the simulation and analysis of highly-generalised 
mathematical ecological models.

### The models may potentially involve:
- Multiple species (community modelling).
- Multiple patches (explicit spatial modelling).
- Arbitrary specification of many features, including:
  - Temporal variation of species properties.
  - Predator-prey relationships, and preferred prey.
  - The graph of the spatial network.
  - The distribution of habitat types on the spatial network.
  - The ability of each species to traverse and reside in the different habitat types.
  - Perturbations induced upon the space or upon particular species.
  - Numerical simulation of either continuous or discrete methods.
  - Ordering of the ecological sub-steps (predation, reproduction, dispersal) with an iteration.
- Future releases will incorporate eco-evolutionary capabilities.

The tools are sufficiently general that non-ecological dynamical systems, such as simple coupled maps, can be 
specified and analysed.

## Getting started:
To use Artemis, only the files parameters.py and parameters_species_repository.py should be edited.

1) Within parameters_species_repository.py create a named dictionary corresponding to your experiment.
   - In this dictionary, each key should be the name of a species.
   - The value should be a dictionary copied from "default_species". Edit the properties of each species as desired:
     - CORE_PARA - Basic properties of lifespan, and the list of predators of this species.
     - INITIAL_POPULATION_PARA - Choice of mechanisms and parameters governing the initial distribution of the species.
        Can be dependent upon habitat or patch.
     - GROWTH_PARA - Set the growth/reproductive function (e.g. logistic map) for the species, and associate properties
        such as carrying capacity and growth rate R. These can be varied according to arbitrary temporal patterns.
     - PREDATION_PARA - If this species predates others, set: the prey species and preferences for them, the predation 
        functional response and associated parameters, including those governing the spatial foraging strategy.
     - DISPERSAL_PARA - Controls for the species dispersal, including the basic mechanism.
     - PURE_DIRECT_IMPACT_PARA - This provides a venue for any bespoke factors you wish to impact the population, such
        as culls (by % or raw amount) or increase of the population by external human intervention.
     - DIRECT_IMPACT_ON_ME - Similarly, any additional impacts that need to be factored, which are caused by the 
        presence of particular other species not modelling by predation or competition.
     - PERTURBATION_PARA - Enables specifying whether the species' presence induces change to the habitat of the same
        or proximate patches, such as reducing the quality or changing the type of habitat with a given probability.
2) Within parameters.py we set all the options for the spatial network, the simulation, and import the species:
   - meta_para - Whether a new simulation should be instigated, and how many repeated simulation to execute.
   - master_para - Holds the parameters governing a single simulation.
     - main_para - Properties including the number of patches, number of time-steps, number of habitat types and their 
        frequencies and associated reproduction and traversal scores for each species, which species and habitat types
        should be present at habitat generation (rather than being introduced subsequently via a perturbation).
     - graph_para - Further specifying the spatial network (which may have already been generated separately),
       including the graph type and associated connectivity parameters, the size and quality of patches, the habitat 
        distribution (manual, clusters, auto-correlation), and restoration characteristics against perturbations.
     - plot_save_para - Options controlling the post-simulation analysis and the amount and detail of files and figures 
        to be generated and saved.
     - pop_dyn_para - Parameters governing the population dynamics of all species, including 'valves' to prevent or 
        allow overall behaviours such as dispersal or ranged foraging.
     - perturbation_para - Enables specifying individual or sequences of population- or patch-level perturbations, 
        which may recur in waves and impact clusters of patches. Perturbation event specifications are stored in the
        PERT_ARCHETYPE_DICTIONARY and then referenced to occur at specified time-steps in PERT_STEP_DICTIONARY. Also 
        allows prescribing of reserves protected against selection.
     - species_para - Imports the desired species dictionary from species_repository_para.py, by replacing 'temp' in
        the following with the name of the dictionary: "x: temp[x] for x in temp".
   
## Overview of scripts and folders:
### Executables:
- main.py - execute to run the main ecological simulation program. You do not need to edit it unless you
   wish to exactly repeat a previous simulation (in this case, comment out the penultimate line new_program() and 
   uncomment repeat_program(sim=?), passing in the number of the prior simulation to repeat).
- sample_spatial_data.py - execute this before the main program if we need to generate a spatial network that the 
   simulation will then take place on. By default, the main program will call this to generate a new network.
- re_analysis.py - used to reload a previously-completed simulation (without re-running it) for examination.

### Advanced use only:
#### Core:
- local_population.py - contains the local_population class and methods for calculating the growth and feeding within 
   the patch. 
- patch.py - contains the patch class. 
- perturbation.py - contains all the functions for conducting patch and habitat perturbations at a given point during 
   the simulation. 
- population_dynamics.py - handles interactions between populations in different patches, and dispersal between 
   patches. 
- simulation_obj.py - handles the initialisation of the program, loading in required data files and putting 
   everything in the necessary format to conduct the simulation. Also controls what happens during the main loop, 
   such as when to call the perturbations.
- species.py - contains the species class. 
- system_state.py - holds the current species and patches, and information on the current spatial network and how it 
    is connected, along with associated file system_state_functions.py.

#### Data handling:
- data_manager.py - handles saving data output files and producing graphs. Functions are stored in the related files:
    - data_core_functions.py
    - data_plot_functions.py
    - data_save_functions.py

#### Additional functions for specialist usage:
- cluster_functions.py - functions used to partition the network for complexity analysis, or to generate connected 
   clusters of patches for perturbation waves or as reserves if required.
- degree_distribution.py - functions for fitting power laws to frequency distributions.

### Generated folders:
- spatial_data_files - stores the files generated by sample_spatial_data.py for use in the main program.
- results - stores results of a simulation in the following directory:
  - [unique simulation number] - each execution of the main program will generate the first available integer
     starting at 100 and assign to this simulation. 
    - parameters.json - contains parameters of the simulation, which can be consulted or reloaded for a re-analysis. 
      This is a copy of the parameters as it was at the beginning of the simulation. 
    - metadata.json - contains information such as the random seeds to repeat the whole simulation, and runtimes.
    - Copies of the parameters.py and parameter_species_repository.py files are stored with the results.
    - [initial_network] - records analytics of the initialised network, after distribution of initial populations but
     prior to the execution of the first time-step.
    - [integer time-step] - at the point when outputs are produced (usually the final time-step), a folder is created 
     with name given by the current step. 
      - It contains further copies of the metadata and parameters, as they were at the time of the output step.
      - Subfolder [data] - contains any output .txt or .JSON files.
      - Subfolder [figures] - contains any output .PNG figures.