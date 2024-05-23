import numpy as np
from population_dynamics import temporal_function


# ------------------------ FUNCTIONAL RESPONSES ------------------------ #

def functional_response(predation_func, predation_para, predator_population, prey_population, attack_rate, time):
    predation_function = {
        "lotka_volterra": functional_response_lotka_volterra,
        "holling_II": functional_response_holling_II,
        "beddington_deangelis": functional_response_beddington_deangelis,
        "ratio_dependent": functional_response_ratio_dependent,
    }
    # include the min function to prevent an individual predator species feeding on more than the entire prey population
    func_res = min(prey_population,
                   predator_population * predation_function[predation_func](attack_rate=attack_rate,
                                                                            predation_para=predation_para,
                                                                            predator_population=predator_population,
                                                                            prey_population=prey_population)
                   )
    return func_res


def functional_response_lotka_volterra(attack_rate, predation_para, predator_population, prey_population):
    func_res = attack_rate * prey_population
    return func_res


def functional_response_holling_II(attack_rate, predation_para, predator_population, prey_population):
    func_res = attack_rate * predation_para["B"] * prey_population / \
               (1.0 + predation_para["B"] * prey_population)
    return func_res


def functional_response_beddington_deangelis(attack_rate, predation_para, predator_population, prey_population):
    func_res = attack_rate * predation_para["B"] * prey_population / \
               (1.0 + predation_para["B"] * prey_population + predation_para["C"] * predator_population)
    return func_res


def functional_response_ratio_dependent(attack_rate, predation_para, predator_population, prey_population):
    func_res = attack_rate * predation_para["B"] * prey_population / \
               (predation_para["B"] * prey_population + predation_para["C"] * predator_population)
    return func_res


# ------------------------ CLASS: LOCAL POPULATION ------------------------ #

class Local_population:

    def __init__(self, species, patch, parameters, current_patch_list=None):
        self.species = species
        self.patch_num = patch.number
        self.parameters = parameters
        self.name = species.name
        self.ode_recording = {}
        self.occupancy = 0
        self.population = 0.0
        self.set_initial_population(patch, current_patch_list=current_patch_list)
        self.temp_population = 0.0
        self.actual_dispersal_targets = {}
        self.g_values = {}
        self.kills = {
            "g1": {},
            "g2": {},
            "g3": {},
        }
        self.killed = {
            "g1": {},
            "g2": {},
            "g3": {},
        }
        self.leaving_array = None
        self.interacting_populations = []
        self.population_history = []
        self.population_leave = 0.0
        self.population_enter = 0.0
        self.net_enter = 0.0
        self.internal_change = 0.0
        self.net_internal = 0.0
        self.sink = 0.0
        self.source = 0.0
        self.internal_change_history = []
        self.population_enter_history = []
        self.population_leave_history = []
        self.record_population_history()  # need this so that the initial population is recorded
        self.population_history_hurst_exponent = 0.0
        self.average_population = 0.0
        self.population_period = 0.0
        self.st_dev_population = 0.0
        self.max_abs_population = 0.0
        self.average_population_leave = 0.0
        self.average_population_enter = 0.0
        self.average_net_enter = 0.0
        self.average_internal_change = 0.0
        self.average_net_internal = 0.0
        self.average_sink = 0.0
        self.average_source = 0.0
        self.recent_occupancy_change_frequency = 0  # how many times in averaging period/num_steps did occupancy change?
        self.stored_pop_values = {}  # used for comparisons during perturbation studies
        self.stored_change = {}  # used for comparisons during perturbation studies
        self.growth_function = {
            "malthusian": self.growth_malthusian,
            "logistic": self.growth_logistic,
        }
        self.current_growth_vector_offset = 0
        self.current_direct_vector_offset = 0
        # Patch-dependent and stored as local population attributes (but may need updated if patch altered)
        self.r_mod = 0.0
        self.carrying_capacity = 0.0
        self.resource_usage_conversion = 0.0
        self.rebuild_patch_dependent_properties(patch=patch)
        self.current_temp_change = 0.0  # for the ODE sub-stages
        self.holding_population = 0.0  # for the ODE sub-stages
        self.r_value = 0.0
        self.r_final = 0.0
        self.l_final = 0.0
        self.k_final = 0.0
        self.competitors_final = 0.0
        self.local_growth = 0.0
        self.direct_impact_value = 0.0
        self.prey_gain = 0.0
        self.predation_loss = 0.0
        self.survivors = 0.0
        self.prey_shortfall = 0.0
        self.predator_shortfall = 0.0

    def set_initial_population(self, patch, current_patch_list):
        # how to specify (or reset) the initial population in the patch based on species-specific scheme
        if self.species.initial_population_mechanism == "constant":
            self.population = self.species.initial_population_para["VALUE"]
        elif self.species.initial_population_mechanism == "random_binomial":
            probability = self.species.initial_population_para["BINOMIAL_PROBABILITY"]
            self.population = self.species.initial_population_para[
                                  "MAXIMUM_MULTIPLIER"] * np.random.rand() * np.random.binomial(n=1, p=probability)
        elif self.species.initial_population_mechanism == "constant_binomial":
            probability = self.species.initial_population_para["BINOMIAL_PROBABILITY"]
            self.population = self.species.initial_population_para[
                                  "MAXIMUM_MULTIPLIER"] * np.random.binomial(n=1, p=probability)
        elif self.species.initial_population_mechanism == "habitat_binomial":
            # A habitat-specific probability of spawning (with fixed size)
            probability = self.species.initial_population_para["HABITAT_TYPE_NUM_BINOMIAL_DICT"][patch.habitat_type_num]
            self.population = self.species.initial_population_para[
                                  "MAXIMUM_MULTIPLIER"] * np.random.binomial(n=1, p=probability)
        elif self.species.initial_population_mechanism == "patch_vector":
            patch_vector = self.species.initial_population_para["PATCH_VECTOR"]
            # check for errors
            if len(patch_vector) != len(current_patch_list) or patch.number not in range(len(patch_vector)):
                raise Exception(f'Error with the length of the initial population patch vector '
                                f'of species {self.species} in patch {patch.number}.')
            else:
                self.population = patch_vector[patch.number]
        else:
            raise Exception(f'Unrecognised (or not given) mechanism for initial population '
                            f'of species {self.species} in patch {patch.number}.')

    def rebuild_patch_dependent_properties(self, patch):
        # necessary as the patch properties may be changed, and this would not change the attributes of a patch object
        # copy stored as an attribute of this local population. MUST BE CAREFUL WITH THIS!
        self.r_mod = patch.quality * patch.this_habitat_species_feeding[self.name]
        self.resource_usage_conversion = self.species.resource_usage_conversion * patch.this_habitat_species_feeding[
            self.name]
        self.carrying_capacity = self.species.growth_para["CARRYING_CAPACITY"] * patch.size

    def record_population_history(self):
        self.population_history.append(self.population)
        self.internal_change_history.append(self.internal_change)
        self.population_leave_history.append(self.population_leave)
        self.population_enter_history.append(self.population_enter)

    def growth_malthusian(self, r_value, patch_competitors, alpha):
        r_ = r_value * self.r_mod
        l_ = self.species.lifespan
        k_ = k_ = self.carrying_capacity
        competitors = 0.0
        gain = r_ * self.holding_population
        mortality = (1 / l_) * self.holding_population
        growth = gain - mortality
        return growth, r_, l_, k_, competitors

    def growth_logistic(self, r_value, patch_competitors, alpha):
        r_ = r_value * self.r_mod
        l_ = self.species.lifespan
        k_ = self.carrying_capacity
        # This is based on RP's model, noting that the maximum function is required to prevent a low r-value
        # allowing one of the "death" terms to flip sign and cause infinite growth.
        if self.resource_usage_conversion == 0.0:
            competitors = self.holding_population
        else:
            # account for all local populations in patch who may need natural resources and/or nesting sites
            # alpha accounts for the absolute scaling of inter-specific competition
            # whilst the .resource_usage_conversion is only for relative scaling of the effect between species.
            competitors = (alpha * patch_competitors + (1.0 - alpha) *
                           self.resource_usage_conversion * self.holding_population) / self.resource_usage_conversion
        growth = self.holding_population * (r_ - 1.0 / l_ - np.max([1.0, (r_ - 1.0 / l_)]) * competitors / k_)
        return growth, r_, l_, k_, competitors

    def set_current_vector_offset(self, time, vector_statement):
        # For both growth and direct impact, we have the additional options of an offset to the annual periodic vectors
        # which can be specified for both individual years AND individual patches!
        #
        # Setup
        if vector_statement == "GROWTH":
            species_vector = self.species.growth_vector_offset_species
            local_vector = self.species.growth_vector_offset_local
            local_is = self.species.is_growth_offset_local
        elif vector_statement == "DIRECT":
            species_vector = self.species.direct_vector_offset_species
            local_vector = self.species.direct_vector_offset_local
            local_is = self.species.is_direct_offset_local
        else:
            raise Exception("Unknown.")
        # main
        max_seasons = len(species_vector)  # if we only specify a few seasons' offsets, repeat the seasons afterwards
        offset_species = species_vector[np.mod(int(time / self.species.seasonal_period), max_seasons)]
        if local_is:
            max_seasons = len(local_vector)
            offset_local = local_vector[np.mod(int(time / self.species.seasonal_period), max_seasons)][self.patch_num]
        else:
            offset_local = 0
        # Update target
        if vector_statement == "GROWTH":
            self.current_growth_vector_offset = offset_species + offset_local
        elif vector_statement == "DIRECT":
            self.current_direct_vector_offset = offset_species + offset_local
        else:
            raise Exception("Unknown.")

    def foraging(self):
        # -------------------------- PREDATION -------------------------- #
        # Predation by this population
        if self.species.current_prey_dict is not None and len(self.species.current_prey_dict) == 0:
            prey_gain = 0.0
        else:
            prey_gain = self.species.predation_para["ECOLOGICAL_EFFICIENCY"] * \
                        sum([x for x in self.kills["g3"].values()])

        # Predation upon this population
        if len(self.species.predator_list) == 0:
            predation_loss = 0.0
        else:
            predation_loss = sum(self.killed["g3"].values())
        self.current_temp_change += prey_gain - predation_loss
        self.predation_loss = predation_loss
        self.prey_gain = prey_gain

    def direct_impact(self, time):
        # ------------------------ DIRECT IMPACT ------------------------ #
        direct_impact = 0.0
        if self.parameters["pop_dyn_para"]["IS_DIRECT_IMPACT"]:
            direct_impact += self.calculate_direct_impact()
        if self.parameters["pop_dyn_para"]["IS_PURE_DIRECT_IMPACT"]:
            if self.species.is_pure_direct_impact:

                if self.species.pure_direct_impact_para["TYPE"] == "binomial":
                    # binomial application to this local population with fixed probability and impact if drawn
                    direct_impact += self.species.pure_direct_impact_para["IMPACT"] * \
                                     np.random.binomial(n=1, p=self.species.pure_direct_impact_para["PROBABILITY"])

                elif self.species.pure_direct_impact_para["TYPE"] == "vector":
                    # if applicable, we check for (species-specific) offset at the beginning of this new season
                    if self.species.is_direct_offset and np.mod(time, self.species.seasonal_period) == 0:
                        self.set_current_vector_offset(time=time, vector_statement="DIRECT")
                    # then the annual duration tells us the length of the year (may be different to the season) and
                    # where we are in it - subject to the current season's (possibly patch-specific) offset if necessary
                    year_time = np.mod(time + self.current_direct_vector_offset, self.species.direct_annual_duration)
                    direct_impact += self.species.pure_direct_impact_para["DIRECT_VECTOR"][year_time]

                else:
                    raise Exception("Error - (pure) direct impact type is not yet coded.")
        self.direct_impact_value = direct_impact
        self.current_temp_change += direct_impact

    def growth(self, parameters, time, patch_competitors, alpha):
        # ------------------- REPRODUCTION / MORTALITY ------------------- #
        # Retrieve current R-value
        r_value = self.species.current_r_value

        # IF APPLICABLE retrieve the species-specific annual offset for this year and location
        # Note that this requires running the temporal function again locally for the patch population.
        if self.species.growth_para["R"]["type"] in ["vector_exp", "vector_imp"] and self.species.is_growth_offset:
            self.set_current_vector_offset(time=time, vector_statement="GROWTH")
            year_time = np.mod(time + self.current_growth_vector_offset, self.species.growth_annual_duration)
            r_value = temporal_function(self.species.growth_para["R"], year_time)

        # Final growth function (accounting for reduced r due to predation, and including mortality which is impacted
        # by carrying capacity if relevant).
        # Note that this can be negative due to overcrowding (intra-specific competition, only for this sub-population),
        # and/or if the growth value is not large enough to compensate for the natural mortality.
        # Within this function, the final effective value of r will also take into account the patch quality
        # and species-specific habitat feeding score.

        local_growth, r_final, l_final, k_final, competitors_final = self.growth_function[
            self.species.growth_function](
            r_value=r_value,
            patch_competitors=patch_competitors,
            alpha=alpha,
        )
        if parameters["main_para"]["MODEL_TIME_TYPE"] == "discrete":
            local_growth_change = local_growth - self.holding_population
        elif parameters["main_para"]["MODEL_TIME_TYPE"] == "continuous":
            local_growth_change = local_growth
        else:
            raise Exception("Model type not recognised in 'main_para[MODEL_TIME_TYPE]' - discrete or continuous?")
        self.local_growth = local_growth_change
        self.current_temp_change += local_growth_change
        self.r_value = r_value
        self.r_final = r_final
        self.l_final = l_final
        self.k_final = k_final
        self.competitors_final = competitors_final

    def ode_recordings(self, time, step):
        # permanently record all aspects for error-checking
        self.ode_recording[step] = {
            "time": time,
            "new_population": self.population,
            "r_value": self.r_value,
            "r_mod": self.r_mod,
            "r_final": self.r_final,
            "l_value": self.l_final,
            "k_value": self.k_final,
            "competitors": self.competitors_final,
            "local_growth": self.local_growth,
            "direct_impact": self.direct_impact_value,
            "prey_gain": self.prey_gain,
            "predation_loss": self.predation_loss,
            "g_values": self.g_values,
            "kills": {y: [x for x in self.kills[y].values()] for y in ["g0", "g1", "g2", "g3"]},
            "killed": {y: [x for x in self.killed[y].values()] for y in ["g0", "g1", "g2", "g3"]},
        }

    def calculate_direct_impact(self):
        # Is there a direct linear impact (if with self, growth function should be malthusian so as not to account for
        # it twice)? For interference competition, mutualistic interactions, or simple CML-style multi-species models.
        # also for general impact from non-species special factors (e.g. culls)
        #
        # include scores (may be positive or negative) for self (0 if already accounted for in logistic growth function)
        # in the parameters["species_para"][species_name]["DIRECT_IMPACT_ON_ME"]
        direct_impact = 0.0
        for population in self.interacting_populations:
            if (self.parameters["pop_dyn_para"]["IS_DIRECT_IMPACT_NONLOCAL"] or self.patch_num == population[
                "object"].patch_num) and population["object"].species.name in self.species.direct_impact_on_me:
                direct_impact += self.species.direct_impact_on_me[population["object"].species.name] * \
                                 self.holding_population * population["object"].holding_population
        return direct_impact

    def calculate_predation(self, time):
        # this function:
        # - sums up the total prey available to this local predator population
        # - calculates its functional response
        # - then distributes that functional response amongst the hunted prey populations in accordance with their size
        # It should be run first to get the predation by the predator if it had sole access to all prey
        prey_sum_shortfall = 0.0
        self.g_values["g2"] = 0.0
        if self.holding_population > 0.0:

            total_prey_hunted = 0.0
            total_effort = 0.0
            predation_efficiency = self.species.current_predation_efficiency
            predation_focus = self.species.current_predation_focus

            for population in self.interacting_populations:
                # this accounts already if nonlocal feeding is allowed in general, for this species, at this range
                if population["object"].holding_population > 0.0:
                    if population["object"].name in self.species.current_prey_dict:
                        # What is predator's preference for this species?
                        this_preference = self.species.current_prey_dict[population["object"].name]
                        # Calculate the effort function. This is designed so that:
                        # - Low efficiency => prioritise prey based on preferred prey species and population size.
                        # - High efficiency => prioritise prey based on accessibility score and population size.
                        #

                        this_effort = population["object"].holding_population * (
                                predation_efficiency * population["score_to"] + (1.0 - predation_efficiency) *
                                this_preference) ** predation_focus

                        # Allocate relative prey hunted:
                        this_prey_hunted = this_effort * population["score_to"] * population[
                            "object"].holding_population
                        # Update running totals:
                        total_effort += this_effort
                        total_prey_hunted += this_prey_hunted
                        # Store info
                        self.kills["g0"][population["object"]] = (this_prey_hunted, this_effort)  # g0

            # Now rescale efforts to sum to 1.0 (may be above or below this):
            if total_effort > 0.0:
                effort_rescale = 1.0 / total_effort
                rescaled_total_prey_hunted = total_prey_hunted * effort_rescale
                # Outcome of g0 just for recording purposes:
                self.g_values["g0"] = rescaled_total_prey_hunted
                # Now pass this actual amount of hunted prey to the functional response, to determine how much of the
                # hunted prey can actually be eaten:
                self.g_values["g1"] = functional_response(
                    predation_func=self.species.predation_para["PREDATION_FUNCTION"],
                    predation_para=self.species.predation_para,
                    predator_population=self.holding_population,
                    prey_population=rescaled_total_prey_hunted,
                    attack_rate=self.species.current_predation_rate,
                    time=time,
                )
                # Now we know the maximum amount that can be eaten, the idea is to proportionally distribute this
                # desired feeding amongst the different available prey.
                for population in self.interacting_populations:
                    if population["object"].holding_population > 0.0:
                        if population["object"].name in self.species.current_prey_dict:
                            if rescaled_total_prey_hunted > 0.0:
                                # effect the rescaling for the predation on this particular prey
                                final_effort = effort_rescale * self.kills["g0"][population["object"]][1]
                                final_prey_eaten = effort_rescale * self.kills["g0"][population[
                                    "object"]][0] * self.g_values["g1"] / rescaled_total_prey_hunted
                            else:
                                final_effort = 0.0
                                final_prey_eaten = 0.0
                            self.kills["g1"][population["object"]] = (final_prey_eaten, final_effort)
                            # Now register the effects on the prey's object directly
                            population["object"].killed["g1"][self] = final_prey_eaten  # g1
                            # record shortfall in prey hunted vs. killed
                            prey_sum_shortfall += max(0.0, self.kills["g0"][population["object"]][0] - final_prey_eaten)
            else:
                self.g_values["g0"] = 0.0
                self.g_values["g1"] = 0.0
        self.prey_shortfall = prey_sum_shortfall

    def predator_allocation(self):
        # this function re-scales the maximum kills/killed if this local prey population is over-hunted
        # by multiple local predator populations
        predator_sum_shortfall = 0.0
        if self.holding_population > 0.0:
            total_predated = sum(self.killed["g1"].values())
            if total_predated > self.holding_population:
                # rescaling the allocations of this prey is necessary
                for population in self.killed["g1"]:
                    # iterate over the potential predators of self
                    scaled_predation = population.kills["g1"][self][0] * self.holding_population / total_predated
                    self.killed["g2"][population] = scaled_predation
                    population.kills["g2"][self] = (scaled_predation, population.kills["g1"][self][1])
                    # build g2
                    population.g_values["g2"] += scaled_predation
                    # running total of predator shortfalls
                    predator_sum_shortfall += max(0.0,
                                                  population.kills["g0"][self][0] - population.kills["g2"][self][0])
                # record remaining prey
                self.survivors = 0.0
            else:
                # no change
                self.killed["g2"] = self.killed["g1"]
                for population in self.killed["g1"]:
                    population.kills["g2"][self] = population.kills["g1"][self]
                    # build g2
                    population.g_values["g2"] += population.kills["g1"][self][0]
                    # running total of predator shortfalls
                    predator_sum_shortfall += max(0.0,
                                                  population.kills["g0"][self][0] - population.kills["g2"][self][0])
                # record remaining prey
                self.survivors = self.holding_population - total_predated
        self.predator_shortfall = predator_sum_shortfall

    def predator_shortfall_distribution(self):
        if self.holding_population > 0.0 and self.g_values['g0'] > 0:  # make sure you are hunting something
            g3_running_total = 0.0
            for population in self.interacting_populations:
                prey = population["object"]
                if prey.holding_population > 0.0:
                    if prey.name in self.species.current_prey_dict:

                        # disparity is d_{i,j,k,l}
                        disparity = max(0.0, self.kills["g0"][prey][0] - self.kills["g2"][prey][0])
                        if disparity > 0:
                            denominator_1 = self.prey_shortfall
                            denominator_2 = prey.predator_shortfall

                            if denominator_1 * denominator_2 == 0.0:
                                top_up = 0.0
                            else:
                                top_up = max(0.0,
                                             disparity * (self.g_values["g2"] - self.g_values["g1"]) / denominator_1 *
                                             min(1.0, prey.survivors / denominator_2))
                                g3_running_total += top_up
                            total_of_this_prey_killed = self.kills["g2"][prey][0] + top_up
                            self.kills["g3"][prey] = total_of_this_prey_killed
                            # Now register the effects on the prey's object directly
                            prey.killed["g3"][self] = total_of_this_prey_killed
                        else:
                            # no disparity, so we skip the top-up stage
                            self.kills["g3"][prey] = self.kills["g2"][prey][0]
                            prey.killed["g3"][self] = self.kills["g2"][prey][0]
            self.g_values["g3"] = self.g_values["g2"] + g3_running_total

    def build_recent_time_averages(self, current_step, back_steps):
        # Can be called at any step to calculate the recent averages of population changes that are being stored in
        # full arrays of the history (but it would be needlessly inefficient to calculate them every time-step, so we
        # only call this in anticipation of upcoming plots - i.e. mainly at the end of the simulation)
        #
        # Mean and standard deviation of recent population history
        self.average_population = np.sum(self.population_history[current_step - back_steps: current_step]) / back_steps
        self.st_dev_population = np.std(self.population_history[current_step - back_steps: current_step])

        # Recent variations in the local population - periodicity, maximum absolute variation, occupancy change:
        period = 0
        period_epsilon = max(0.00000000001, self.st_dev_population * 0.0000001)
        max_abs_var = np.abs(self.population_history[current_step] - self.average_population)
        num_occupancy_changes = 0
        min_pop = self.species.minimum_population_size
        for n in range(1, back_steps):
            # did the occupancy change?
            new_pop = self.population_history[current_step - n]
            old_pop = self.population_history[current_step - n + 1]
            if (new_pop < min_pop <= old_pop) or (new_pop >= min_pop > old_pop):
                num_occupancy_changes += 1
            # update greatest absolute deviation from the mean
            max_abs_var = max(max_abs_var, np.abs(self.population_history[current_step - n] - self.average_population))
            # periodicity check
            if period == 0:
                # only check if not already determined
                if abs(self.population_history[current_step - n] -
                       self.population_history[current_step]) < period_epsilon:
                    # if we think we have a period M, only record it if we can confirm against X_{n-2M} ~ X_{n} also
                    if len(self.population_history) >= 2 * n:
                        if abs(self.population_history[current_step - 2 * n] - self.population) < period_epsilon:
                            period = n
        self.population_period = period
        self.max_abs_population = max_abs_var
        if back_steps > 1:
            self.recent_occupancy_change_frequency = num_occupancy_changes / (back_steps - 1.0)
        else:
            self.recent_occupancy_change_frequency = 0.0

        # Average population change due to the internal ODE/Difference Equation (including possibly distant foraging by
        # this species and distant predation upon this species) AND direct impact (i.e. everything except dispersal):
        self.average_internal_change = np.sum(
            self.internal_change_history[current_step - back_steps: current_step]) / back_steps
        # Average population emigrated during dispersal
        self.average_population_leave = np.sum(
            self.population_leave_history[current_step - back_steps: current_step]) / back_steps
        # Average population immigrated during dispersal
        self.average_population_enter = np.sum(
            self.population_enter_history[current_step - back_steps: current_step]) / back_steps
        # Average net immigration during dispersal
        self.average_net_enter = np.sum(
            np.array(self.population_enter_history[current_step - back_steps: current_step]) -
            np.array(self.population_leave_history[current_step - back_steps: current_step])) / back_steps
        # Average net internal (see description below)
        total_change = np.sum(
            np.abs(self.internal_change_history[current_step - back_steps: current_step])) + np.sum(
            np.abs(np.array(self.population_enter_history[current_step - back_steps: current_step]) -
                   np.array(self.population_leave_history[current_step - back_steps: current_step])))
        if total_change == 0.0:
            self.average_net_internal = 0.0
        else:
            self.average_net_internal = np.sum(np.abs(self.internal_change_history[
                                                      current_step - back_steps: current_step])) / total_change
        # Average sink detection (i.e. fraction of current population that entered by immigration last step)
        sum_sink_change = 0.0
        for n in range(back_steps):
            if self.population_history[current_step - n] != 0.0:
                sum_sink_change += self.population_enter_history[current_step - n] / \
                                   self.population_history[current_step - n]
        self.average_sink = (1.0 / back_steps) * sum_sink_change
        # Average source detection (i.e. fraction of previous population that left by emigration)
        sum_source_change = 0.0
        for n in range(back_steps):
            if self.population_history[current_step - n - 1] != 0.0:
                sum_source_change += self.population_leave_history[current_step - n] / \
                                     self.population_history[current_step - n - 1]
        self.average_source = (1.0 / back_steps) * sum_source_change

    def update_local_nets(self):
        # occupancy of patch
        if self.population >= self.species.minimum_population_size:
            self.occupancy = 1
        else:
            self.occupancy = 0
        # net values:
        #
        # Net immigration (vs. emigration)
        self.net_enter = self.population_enter - self.population_leave
        #
        # Net internal:
        # A measure of the fraction of the MAGNITUDE of the population's change that is due to the internal ODE vs.
        # the total combination of the internal ODE plus net migration
        total_change = np.abs(self.internal_change) + np.abs(self.net_enter)
        if total_change == 0.0:
            self.net_internal = 0.0
        else:
            self.net_internal = np.abs(self.internal_change) / total_change
        # Sink detection - (non-native entered compared to net CURRENT population):
        # Fraction of the current population that entered via immigration at the most recently completed step
        # Note that this can be greater than 1.0, if sufficient population died, left, or was predated simultaneously.
        if self.population == 0.0:
            self.sink = 0.0
        else:
            self.sink = self.population_enter_history[-1] / self.population_history[-1]
        # Source detection - (emigration compared to PREVIOUS population):
        # Fraction of the previous population that emigrated in most recently completed step
        if self.population_history[-2] == 0.0:
            self.source = 0.0
        else:
            self.source = self.population_leave_history[-1] / self.population_history[-2]
        # Source notes:
        # This is called after (not 'within') a time-step. At the end of the step histories are updated so
        # self.population_leave is amount that left in the step, and self.population & self.population_history[-1]
        # will BOTH be the current population at the end of this same step. Thus use [-2] to access the previous
        # population recorded before this most recent time step. Because of the different possible sub-stage
        # configurations, the most consistent across them is just to record the proportion leaving in this
        # step (whether before, during, or after reproduction) compared to population present at the end of the
        # previous full time-step.
        # Timeline:
        # - record population_history[-2] etc.
        # - complete one entire time-step of reproduction, foraging, direct impact, dispersal
        # - update internal_change, population_leave, population_enter and FROM ALL OF THESE population and
        # population_history[-1], internal_change_history[-1], population_leave[-1], population_enter[-1].
        # So we do not actually record the population state between ecological sub-stages.
