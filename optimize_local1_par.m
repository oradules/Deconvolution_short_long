%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% this program implements the local optimisation algorithm for polymerase positions %%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [positions,Min_fit]=optimize_local(target,guess,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)
%%%%% pattern perturbation to search local minimum


        GD_y_fitness = @(x) sum((sumSignal1_par(x,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)-target).^2); % x: positions of poly
         
        positions = find(guess==1);
        Nbr_poly_estimate = length(positions);
        shift_window = round(num_possible_poly/Nbr_poly_estimate)+20;%one position can move [-s_w,s_w]
       
        Min_fit = GD_y_fitness(positions);
        
        for posi_i = 1: length(positions)
            new_pos = positions;
            for j = -shift_window:shift_window
                new_pos(posi_i) = positions(posi_i) + j; %%% displaced position
                if new_pos(posi_i) <= 0 || new_pos(posi_i) > num_possible_poly
                    continue
                end
                fitness=GD_y_fitness(new_pos);
                if fitness<Min_fit
                    positions = new_pos;
                    Min_fit = fitness;
                end
            end
        end
%Min_fit