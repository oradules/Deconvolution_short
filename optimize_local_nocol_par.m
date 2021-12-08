%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% this program implements the local optimisation algorithm for polymerase positions %%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% new version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [positions,Min_fit]=optimize_local(target,guess,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)
%%%%% pattern perturbation to search local minimum


        GD_y_fitness = @(x) sum((sumSignalv_par(x,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)-target).^2); % x: positions of poly
         
        positions = find(guess==1);
        Nbr_poly_estimate = length(positions);
        shift_window = round(num_possible_poly/Nbr_poly_estimate)+20;%one position can move [-s_w,s_w]
       
        Min_fit = GD_y_fitness(positions);
        
        for posi_i = 1: length(positions)
            new_pos = positions;
            to_test=positions(posi_i)+(-shift_window:3:shift_window); %%% positions to test
            
            %%%%% perform the collision test before the loop;
            to_test(ismember(to_test,positions) | (to_test <=0) | (to_test > num_possible_poly) )=[];
            %%%%% compute the fitness for each admissible position 
            
            
            for j = 1:length(to_test)   
                new_pos(posi_i)=to_test(j);%%%% admissible positions
                fitness=GD_y_fitness(new_pos);
                if fitness<Min_fit
                    positions(posi_i) = to_test(j);
                    Min_fit = fitness;
                end    
            end
        end
%Min_fit