##########################################################################
## A genetic algorithm for the graph coloring problem
## Author: Atilio Gomes Luiz
## Date: May 28, 2024
##
## This is an implementation of the algorithm described in the technical report:
## "An ordering-based genetic approach to graph coloring"
## By Cornelius Croitoru, Ovidiu Gheorgies and Adriana Gheorghies
## Link: https://typeset.io/pdf/an-ordering-based-genetic-approach-to-graph-coloring-1ddn24whgb.pdf
#########################################################################
import random
import heapq

# ----------------------
# Class chromosome.
# A chromosome is an object that encapsulates 
# two pieces of information: 
# a permutation of the graph's vertices and 
# the fitness value of this permutation.
# ----------------------
class chromosome:
    def __init__(self, permutation_list, fitness_value):
        self.permutation_list = permutation_list.copy()
        self.fitness_value = fitness_value
    
    def get_fitness(self):
        return self.fitness_value
    
    def set_fitness(self, fitness_value):
        self.fitness_value = fitness_value
    
    def get_permutation(self):
        return self.permutation_list
    
    def __lt__(self, other):
        return self.fitness_value <= other.fitness_value
    
    def __str__(self):
        #return f"chromosome(permutation_list={self.permutation_list}, fitness={self.fitness_value})"
        return f"number of colors = {self.fitness_value}"
    

# --------------------------------
# Class ga_graph_coloring
# --------------------------------
class ga_graph_coloring:
    ##
    # Constructor
    ##
    def __init__(self, graph, size_population, crossover_rate, 
                 mutation_rate_block, mutation_rate_neighbors, elitism_rate):
        # graph to be colored
        self.G = graph
        # size of the population (number of chromosomes)
        self.size_population = size_population 
        # probability of the crossover operation being carried out
        self.crossover_rate = crossover_rate
        # probability of the block_move_mutation operation being carried out
        self.mutation_rate_block = mutation_rate_block
        # probability of the swap_neighbors_mutation operation being carried out
        self.mutation_rate_neighbors = mutation_rate_neighbors
        # elitism rate: determines the quantity of best chromossomes that survive from on
        # generation to another
        self.elitism_rate = elitism_rate
        # list that saves all the chromosomes
        self.population = list()
        # maximum degree of the graph G
        self.G_maximum_degree = self.__max_degree()
        

    ##
    # calculates the maximum degree of the graph
    ##
    def __max_degree(self):
        maximum_degree = -1
        for v in self.G:
            if len(self.G[v]) > maximum_degree:
                maximum_degree = len(self.G[v])
        return maximum_degree

    ##
    # This function is responsible for running the genetic algorithm.
    ##
    def run(self, number_of_generations):
        # Step 1: Create initial population of randomly generated chromosomes
        self.__create_initial_population()
        best_so_far_cost = self.__get_best_chromosome().get_fitness()

        # Step 2: Run the generations until stopping criteria is reached.
        # The loop will run until either the maximum number of generations is reached,
        # or if for 30 generations no better coloring is found.
        count_unchangeness = 0
        counter = 0
        while counter < number_of_generations and count_unchangeness <= 30:
            self.__update_population()
            iteration_best_cost = self.__get_best_chromosome().get_fitness()
            if(iteration_best_cost < best_so_far_cost):
                best_so_far_cost = iteration_best_cost
                count_unchangeness = 0
            else:
                count_unchangeness = count_unchangeness + 1

            counter = counter + 1
        
        return self.__get_best_chromosome()


    ##
    # Create an initial population of chromosomes
    # Each chromosome contains a randomly generated 
    # permutation of the sequence [0..n-1]
    ##
    def __create_initial_population(self):
        for i in range(0, self.size_population):
            permutation = [x for x in range(0,len(self.G))]
            random.shuffle(permutation)
            individual = chromosome(permutation, self.__fitness(permutation))
            self.population.append(individual)


    ##
    # Returns the best chromosome in the population
    ##
    def __get_best_chromosome(self):
        best_chromosome = self.population[0]
        for individual in self.population[1:]:
            if individual.get_fitness() < best_chromosome.get_fitness():
                best_chromosome = individual 
        return best_chromosome


    ##
    # Function that calculates the fitness value of a vertex ordering of the graph G.
    # Let n be the number of vertices o the graph G.
    # The vertices of G are integers in the range [0..n-1].
    # This function receives as input a list called 'permutation_list' 
    # containing a permutation of the numbers [0..n-1]
    # and returns the maximum color used by a Greedy Coloring of G, according to 
    # the given permutation of the vertices of the graph.
    ##
    def __fitness(self, permutation_list):
        vertex_colors = [0] * len(self.G)
        maximum_color_used = 0
        for v in permutation_list:
            for color in range(1, self.G_maximum_degree + 2):
                color_is_forbidden = False
                for neighbor in self.G[v]:
                    if vertex_colors[neighbor] == color:
                        color_is_forbidden = True
                        break
                if color_is_forbidden:
                    continue
                else:
                    vertex_colors[v] = color
                    if color > maximum_color_used:
                        maximum_color_used = color
                    break

        return maximum_color_used
    

    ##
    # Crossover operator: Permutation One Point Crossover (POP)
    # This procedure receives two parent chromossomes as input and 
    # returns two offsprings.
    # First, it randomly selects a cut-point.
    # The first offspring is created as follows:
    # The first portion of parent 1 up to the cut-point becomes the
    # first part of offspring 1, and the remainder of offspring 1 is 
    # obtained by copying the vertices absent from the first portion
    # of the offspring 1 in the same ordering as they occur in parent 2.
    # The second offspring is generated in a symmetrical way.
    ##
    def __crossover_pop(self, chromosome1, chromosome2):
        # creates the first offspring
        cut_point = random.randint(0,len(chromosome1.get_permutation())-2)
        permutation_list = list()
        for i in range(0, cut_point+1):
            permutation_list.append(chromosome1.get_permutation()[i])
        sublist = chromosome1.get_permutation()[cut_point+1:]
        for v in chromosome2.get_permutation():
            if v in sublist:
                permutation_list.append(v)
        offspring1 = chromosome(permutation_list, self.__fitness(permutation_list))

        # creates the second offspring
        cut_point = random.randint(0,len(chromosome2.get_permutation())-2)
        permutation_list = list()
        for i in range(0, cut_point+1):
            permutation_list.append(chromosome2.get_permutation()[i])
        sublist = chromosome2.get_permutation()[cut_point+1:]
        for v in chromosome1.get_permutation():
            if v in sublist:
                permutation_list.append(v)
        offspring2 = chromosome(permutation_list, self.__fitness(permutation_list))

        return [offspring1, offspring2]
    

    ##
    # First mutation operator (it is not used)
    ##
    def __random_swap_mutation(self, chromosome, number_of_interchanges):
        counter = number_of_interchanges
        while counter > 0:
            p1 = random.randint(0, len(chromosome.get_permutation())-1)
            p2 = p1
            while p2 == p1:
                p2 = random.randint(0, len(chromosome.get_permutation())-1)
            aux = chromosome.get_permutation()[p1]
            chromosome.get_permutation()[p1] = chromosome.get_permutation()[p2]
            chromosome.get_permutation()[p2] = aux
            counter = counter - 1

        chromosome.set_fitness(self.__fitness(chromosome.get_permutation()))
        

    ##
    # Second mutation operator
    ##
    def __block_move_mutation(self, chromosome, block_size):
        n = len(chromosome.get_permutation())
        i = random.randint(0, n-block_size)
        list_of_values = [x for x in range(0,i)]
        for x in range(i+block_size, n-1):
            list_of_values.append(x)
        j = random.choice(list_of_values)
        if j < i:
            sublist = chromosome.get_permutation()[j:i]
            for l in range(0,block_size):
                chromosome.get_permutation()[j+l] = chromosome.get_permutation()[i+l]
            for l in range(0, i-j):
                chromosome.get_permutation()[j+block_size+l] = sublist[l]
        if j >= i + block_size:
            sublist1 = chromosome.get_permutation()[i:i+block_size]
            sublist2 = chromosome.get_permutation()[i+block_size:j]
            sublist = sublist2 + sublist1
            counter = 0
            while counter < len(sublist):
                chromosome.get_permutation()[i+counter] = sublist[counter]
                counter = counter + 1

        chromosome.set_fitness(self.__fitness(chromosome.get_permutation()))

    ##
    # Third mutation operator
    ##
    def __neighbors_swap_mutation(self, chromosome):
        # Step 1: randomly choose a vertex v of self.G
        v = random.choice([x for x in range(0,len(self.G))])
        # Step 2: randomly choose pairs of neighbors of v and swap them
        neighbors = [w for w in self.G[v]]
        random.shuffle(neighbors)
        for i in range(0, len(neighbors), 2):
            if i < len(neighbors)-1:
                aux = chromosome.get_permutation()[neighbors[i]]
                chromosome.get_permutation()[neighbors[i]] = chromosome.get_permutation()[neighbors[i+1]]
                chromosome.get_permutation()[neighbors[i+1]] = aux

   
    ##
    # Returns a list containing the k chromosomes with best 
    # fitness value from the current population.
    # The number k is given by:
    # k = self.size_population * self.elitism_rate
    #
    # In order to select the k best chromossomes in a fast way,
    # this function uses a data structure called priority queue.
    # Here, I am using the priority queue implementation called 'heapq' 
    ##
    def __elitism(self):
        priority_queue = []

        for individual in self.population:
            number = int(individual.get_fitness())
            heapq.heappush(priority_queue, (number, individual))
            
        number_of_individuals = int(self.size_population * self.elitism_rate)

        new_population = list()
        for i in range(0, number_of_individuals):
            fitness, individual = heapq.heappop(priority_queue)
            new_population.append(individual)

        return new_population
    
    ##
    # Tournament Selection
    # This function is responsible for selecting t chromosomes
    # from the current population (where t = num_parents).
    # The basic idea behind tournament selection is to randomly 
    # choose a subset of k individuals from the population (where k = tournament_size), 
    # compare their fitness values, and select the best individual from 
    # this subset. This process is repeated in order to select the required number of parents.
    ##
    def __tournament_selection(self, population, tournament_size, num_parents):
        parents = list()
        while len(parents) < num_parents:
            # Step 1: randomly select 'tournament_size' different individuals
            tournament = random.sample(population, tournament_size)
            # Step 2: evaluate fitness of each individual in the tournament
            best_individual = tournament[0]
            for individual in tournament[1:]:
                if individual.get_fitness() > best_individual.get_fitness():
                    best_individual = individual
            # Step 3: Select the best individual as a parent
            parents.append(best_individual)
        
        return parents


    ##
    # This function generates a new population based on the 
    # application of the crossover, mutation and selection operators.
    ## 
    def __update_population(self):
        # Step 1: save the k best individuals from the current generation
        new_population = self.__elitism()

        # Step 2: generate n-k new child chromosomes
        while len(new_population) < self.size_population:
            # selects two parents from a tournament of 3 individuals
            parents = self.__tournament_selection(self.population, 3, 2)
            # undertake crossover and mutation operations according to the specified rates
            if random.random() < self.crossover_rate:
                children = self.__crossover_pop(parents[0],parents[1])
            else:
                children = parents
            if(random.random() < self.mutation_rate_block):
                self.__block_move_mutation(children[0], 3)
            if(random.random() < self.mutation_rate_block):
                self.__block_move_mutation(children[1], 3)
            if(random.random() < self.mutation_rate_neighbors):
                self.__neighbors_swap_mutation(children[0])
            if(random.random() < self.mutation_rate_neighbors):
                self.__neighbors_swap_mutation(children[1])
            
            # add the children to the new population
            new_population.append(children[0])
            if len(new_population) < len(self.population):
                new_population.append(children[1])

        self.population = new_population


