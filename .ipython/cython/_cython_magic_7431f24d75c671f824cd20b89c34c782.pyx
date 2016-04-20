
##
import pandas as pd
import numpy as np
import igraph as ig
import random
import os
import sys



class SIModel(object):
    """SI epidemic model for networks."""
    
    def __init__(self, graph, nodes_to_watch, p=0.1, tests=False, 
                 starting_infections=1, seed=None):
        """Constructs an SI model on the given `graph` with
        infection rate `p`"""
        self.graph = graph
        self.p= float(p)
        self.infected = set()
        self.uninfected = set(range(graph.vcount()-1))
        self.watched_infected = set()
        self.watched_uninfected = nodes_to_watch.copy()
        self.nodes_to_watch = nodes_to_watch
        self.tests = tests
        self.share_infected = 0.

        if seed is not None:
            np.random.seed(seed)
            random.seed(seed)

        
        self.set_initial_infections(starting_infections)
        
    def step(self, iteration_counter=None):
        """Runs a single step of the SI model simulation."""
        # Infected vertices spread the infeciton
        newly_infected = set()
        
        binomial = np.random.binomial
        choice = np.random.choice
        get_neighbors = self.graph.neighbors
        
        for vertex in self.infected:
            neighbors = get_neighbors(vertex)
            sample_size = binomial(len(neighbors), self.p)    
            newly_infected.update(choice(neighbors,sample_size, replace=False))

        self.move_to_infected(newly_infected)
        
        if self.tests:
            print('iteration {}'.format(iteration_counter))
            print("infected, then uninfected")
            print(self.watched_infected)
            print(self.watched_uninfected)
            print("neighbors of infected")
            ns = set()
            for i in self.watched_infected:
                ns.update(get_neighbors(i))
            print(ns)
        
        self.share_infected = self.share_watched_infected()

    def move_to_infected(self, newly_infected):
        # Move into infected sets, remove if in uninfected
        self.infected.update(newly_infected)
        self.uninfected.difference_update(newly_infected)
        
        
        # Deal with subset of nodes we watch. 
        newly_infected.intersection_update(self.nodes_to_watch)

        self.watched_infected.update(newly_infected)
        self.watched_uninfected.difference_update(newly_infected)

        if self.tests:
            assert len(self.infected.intersection(self.uninfected)) == 0        
            assert newly_infected.issubset(self.nodes_to_watch)
            assert len(self.watched_infected.intersection(self.watched_uninfected)) == 0


    def share_watched_infected(self):
        if self.tests:
            assert len(self.watched_infected)+len(self.watched_uninfected) == len(self.nodes_to_watch)
        
        return len(self.watched_infected) / len(self.nodes_to_watch)

    def set_initial_infections(self, starting_infections):
        initials = set(random.sample(self.watched_uninfected, starting_infections))
        self.move_to_infected(initials)
        self.share_infected = self.share_watched_infected()



def diffusion_run(nodes_to_watch, graph, p, district, tests=False, max_iter=1000,
                  seed=None):
                      
    si_model = SIModel(graph=graph, nodes_to_watch=nodes_to_watch, 
                    tests=tests, p=p, starting_infections=1, seed=seed)
    

    values_to_record = [0.1, 0.25, 0.5, 0.75]
    result = pd.Series([np.nan, np.nan, np.nan, np.nan], index=values_to_record)
    result.name = district

    for i in range(1,max_iter):
        si_model.step(i)
        for threshold in values_to_record:
            if pd.isnull(result.loc[threshold]) and si_model.share_infected >= threshold:                
                result.loc[threshold] = i
        if pd.notnull(result.loc[values_to_record[-1]]): break

    return result

