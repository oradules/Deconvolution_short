#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 09:24:58 2021

@author: mdouaihy

Code written by Ovidiu Radulescu, University of Montpellier, June 2019
this program implements the local optimisation algorithm for polymerase positions
"""

import numpy as np
from sumSignalDroso import sumSignal1_par

def optimize_local_nocol_par(target,guess,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq,
                        TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym):

   ###pattern perturbation to search local minimum

    def GD_y_fitness(x):
            return sum((sumSignal1_par(x,FreqEchSimu, FreqEchImg, TaillePreMarq,
                                        TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)-target)**2)  ## positions of poly
             
    positions = np.where(guess==1)[0]
    Nbr_poly_estimate = len(positions)
    shift_window = round(num_possible_poly/Nbr_poly_estimate)+20 ## one position can move [-s_w,s_w]
    
    Min_fit = GD_y_fitness(positions)
       

        
    for posi_i in range(len(positions)):
        new_pos = positions.copy()
        to_test=positions[posi_i]+np.arange(-shift_window,shift_window+1,3) ### positions to test
        
        ### perform the collision test before the loop;
        to_remove=np.where((to_test<0) | (to_test>num_possible_poly-1) | np.isin(to_test,positions))[0]
        to_test=np.delete(to_test,to_remove)

        ####compute the fitness for each admissible position 
        for j in range(len(to_test)):   
            new_pos[posi_i]=to_test[j] ## admissible positions
            fitness=GD_y_fitness(new_pos)
            if fitness<Min_fit:
                positions[posi_i] = to_test[j]
                Min_fit = fitness
    return positions, Min_fit