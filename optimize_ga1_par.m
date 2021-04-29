%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% this program implements the genetic algorithm for polymerase positions %%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pattern_fit=optimize_ga(target,Nbr_poly_estimate,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym, generations)

GD_fitness = @(x) sum((sumSignal1_par(find(x==1),FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)-target).^2); % x: positions of poly
               
        
Nbr_simu_DNA = 500; % number of "chromosome" in GA 
Pattern_polys = zeros(Nbr_simu_DNA,num_possible_poly);
for i = 1:Nbr_simu_DNA
    %[num_possible_poly,Nbr_poly_estimate]
    Pattern_polys(i,randperm(num_possible_poly,Nbr_poly_estimate)) = 1; % randomly choose poly position 
end
       
options = gaoptimset;
options = gaoptimset(options,'PopulationType', 'bitstring', 'CreationFcn',@gacreationuniform,'MutationFcn',@mutationuniform,...
            'TolFun',1e-12,'Paretofraction',0.35,'Generations',generations);
options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
options = gaoptimset(options,'InitialPopulation', Pattern_polys);
options = gaoptimset(options,'Display','none'); %%% 'iter' or 'none' 

	
%options = gaoptimset(options,'ga','UseParallel', false);
     
pattern_fit = ga(GD_fitness,num_possible_poly,[],[],[],[],[],[],[],[],options);
       
end
     
 