# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 12:10:37 2016

@author: Nick
"""
from __future__ import division


##
import numpy as np
import igraph as ig
import random
import pandas as pd


cdef list sample_from_neighbors(list neighbors, float p):
    to_return = []    
    for i in neighbors:
        if random.random() < p:
            to_return.append(i)
    return to_return


cdef set stand_alone_step(object graph, set infected, float p):
        # Infected vertices spread the infeciton
        newly_infected = set()
                
        for vertex in infected:
            neighbors = graph.neighbors(vertex)
            choices = sample_from_neighbors(neighbors, p)
            newly_infected.update(choices)

        return newly_infected
        


class SIModel(object):
    """SI epidemic model for networks."""
    
    def __init__(self, graph,  nodes_to_watch,  p=0.1,  tests=False, 
                  starting_infections=1,  seed=None,  debug=False, 
                  initial_infection=None):
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
        self.debug = debug

        if seed is not None:
            np.random.seed(seed)
            random.seed(seed)

        if initial_infection is not None:
            if len(initial_infection) == starting_infections:
                self.move_to_infected(initial_infection)
                
            else:
                raise ValueError("If pass initial infection set, starting_infection must be same size!")
        else:
            self.set_initial_infections(starting_infections)
        
    def step(self, iteration_counter=None):
        """Runs a single step of the SI model simulation."""
        newly_infected = stand_alone_step(self.graph, self.infected, self.p)
        self.move_to_infected(newly_infected)
        
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

    @property
    def share_infected(self):
        if self.tests:
            assert len(self.watched_infected)+len(self.watched_uninfected) == len(self.nodes_to_watch)
        
        return len(self.watched_infected) / len(self.nodes_to_watch)

    def set_initial_infections(self, starting_infections):
        initials = set(random.sample(self.watched_uninfected, starting_infections))
        self.move_to_infected(initials)



def diffusion_run(nodes_to_watch, graph, p, district, tests=False, max_iter=1000,
                  seed=None, thresholds=[0.1, 0.25, 0.5, 0.75], starting_infections=1,
                    debug=False, initial_infection=None):
    
    if tests:
        test_suite()
    
    # Sanity check
    assert max(nodes_to_watch) < graph.vcount()
                  
    si_model = SIModel(graph=graph, nodes_to_watch=nodes_to_watch, 
                    tests=tests, p=p, starting_infections=starting_infections, 
                    seed=seed, debug=debug, initial_infection=initial_infection)
    

    result = pd.Series(np.nan, index=thresholds)
    result.name = district

    for i in range(0,max_iter):
        for threshold in thresholds:
            if pd.isnull(result.loc[threshold]) and si_model.share_infected >= threshold:                
                result.loc[threshold] = i
        if pd.notnull(result.loc[thresholds[-1]]): break
        si_model.step(i)


    return result
    

def test_suite():
    from pandas.util.testing import assert_series_equal
    line = ig.Graph()
    line.add_vertices(range(2))
    line.add_edges([(0,1)])    
    
    # Start with 1 of 2, should spread from half (starting seed) to full from 0 to 1. 
    test_result = diffusion_run(graph=line, nodes_to_watch=set(range(2)), p=1, tests=False,
                      district=0)
    assert_series_equal(test_result, pd.Series([0.0,0,0,1],index = [0.1, 0.25, 0.5, 0.75], name=0))

    # If one watched and one seeded, should all trigger at 0.
    test_result = diffusion_run(graph=line, nodes_to_watch=set(range(1)), p=1, tests=False,
                      district=0)
    assert_series_equal(test_result, pd.Series([0.0,0,0,0],index = [0.1, 0.25, 0.5, 0.75],name=0))                      

    # If watch 2 and start with 2, trigger right away
    test_result = diffusion_run(graph=line, nodes_to_watch=set(range(2)), p=1, tests=False,
                      district=0, starting_infections=2)
    assert_series_equal(test_result , pd.Series([0.0,0,0,0],index = [0.1, 0.25, 0.5, 0.75],name=0))                    


    # Immediate spread in full graph
    full = ig.Graph.Full(10)
    test_result = diffusion_run(graph=full, nodes_to_watch=set(range(10)), p=1, tests=False,
                      district=0)
    assert_series_equal(test_result , pd.Series([0.0,1,1,1],index = [0.1, 0.25, 0.5, 0.75],name=0))

    # No spread in disconnected
    sparse = ig.Graph()
    sparse.add_vertices(4)
    sparse.add_edges([(0,1), (2,3)])
    test_result = diffusion_run(graph=sparse, nodes_to_watch=set(range(4)), p=1, tests=False,
                      district=0)
    assert_series_equal(test_result , pd.Series([0.0,0,1,np.nan],index = [0.1, 0.25, 0.5, 0.75],name=0))                      
    
    
    # No spread if p=0
    sparse = ig.Graph.Full(2)
    test_result = diffusion_run(graph=sparse, nodes_to_watch=set(range(2)), p=0, tests=False,
                      district=0)
    assert_series_equal(test_result , pd.Series([0.0,0,0,np.nan],index = [0.1, 0.25, 0.5, 0.75],name=0))                      
    
    

    # If make line and start at tail, should diffuse one per step
    line = ig.Graph()
    line.add_vertices(range(10))
    line.add_edges([(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9)])    
    
    test_result = diffusion_run(graph=line, nodes_to_watch=set(range(10)), p=1, initial_infection={0}, district=0)
    assert_series_equal(test_result, pd.Series([0.0,2,4,7],index = [0.1, 0.25, 0.5, 0.75], name=0))

    # Modify thresholds
    test_result = diffusion_run(graph=line, nodes_to_watch=set(range(10)), p=1, initial_infection={0}, district=0, thresholds=[0,0.1,1])
    assert_series_equal(test_result, pd.Series([0.0,0,9],index = [0,0.1,1], name=0))

    ######
    # Test probabilities
    ####

    # negative binomial -- num failures before finishing, so subtract 1 from step count
    # For p=0.2
    # Mean 4, sd 4.5 -- 3se is 0.42. 
    run = 1000
    output = pd.Series(np.nan*run)
    for i in range(run):
        r = diffusion_run(graph=ig.Graph.Full(2), nodes_to_watch=set(range(2)), p=0.2, tests=False,
                      district=0, thresholds=[1], initial_infection={0})
        output.loc[i] = r.loc[1]

    output.mean()-1
    try:       
        assert output.mean()-1 > 3.5 and output.mean()-1 <4.5
    except:
        raise ValueError("Test on 1-link diffusion should have taken avg 5 steps (se 0.14), took {}".format(output.mean()))


    # More complex. Failures before 9 successes, so steps minus 9.
    # mean 3.86, sd2.3, 3 standard errors thus ~.7
    # for p=0.7
    line = ig.Graph()
    line.add_vertices(range(10))
    line.add_edges([(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9)])    
    run = 100
    output = pd.Series(np.nan*run)
    for i in range(run):
        r = diffusion_run(graph=line, nodes_to_watch=set(range(10)), p=0.7, tests=False,
                      district=0, thresholds=[1], initial_infection={0})
        output.loc[i] = r.loc[1]

    output.mean() - 9
    try:       
        assert output.mean() - 9 < 4.5 and output.mean() - 9 > 3.1
    except:
        raise ValueError("Test on 9-line diffusion should have taken avg 12.86 steps (se 0.22), took {}".format(output.mean()))

    print('test graphs ok')
