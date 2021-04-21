%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% genetic algorithm optimization for polymerase positions
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% uses optimize_ga1_par.m, sumSignal1_par.m, optimize_local1_par.m  and 
%%%% several short movie datasets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% computation of polymerase positions is performed in parallel, one cell per processor 
%%%% the result is written in the file result_name_cPosPred.mat   %%%%%%%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TaillePreMarq = 700; % 700 bases (ref publi)
TailleSeqMarq = 5800; % 2900 bases (ref publi)
TailleSeqMarq = 2900; % 2900 bases (ref publi)
TaillePostMarq = 1600 + 67*100; % 1600 bases (ref publi)
TaillePostMarq = 5300 + 67*100; % 1600 bases (ref publi)
EspaceInterPolyMin = 10; % 
Polym_speed = 67; % average speed bases par seconde (Ref publi)
FrameLen = 3; %%% frame length in seconds
Intensity_for_1_Polym = 1; %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% file names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
names = {'tat_latency'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




nthr='05'; %%%% the low intensity value considered to be zero

for iii=1:length(names)


name = names{iii};


fname = ['data_',name,'_short_',nthr,'.mat'];
load(fname);   %%%% load datafile





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sd=size(DataExp); %%%%% size matrix data
nloops = sd(2);
frame_num=sd(1); %%%%%% number of frames
FreqEchImg = (1/FrameLen); %%%%% 1/3 image per seconde data time sampling
DureeSimu = frame_num*FrameLen; %%% film duration in s
DureeSignal = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed; % (s)
DureeAnalysee = DureeSignal + DureeSimu ; % (s)
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); % how many interval(possible poly start position) in 1s
num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DataPred=zeros(sd);
PosPred=zeros(num_possible_poly,sd(2));
cPosPred=cell(1,sd(2));


Fit=zeros(sd(2),1);
%%%%%% parallel computing %%%%%%%%%%%%%%%%%%%
mycluster = parcluster; %%% create cluster
set(mycluster,'NumWorkers',50) ; %%% change max number of workers
parpool(mycluster,min([nloops,50])); %%% open a pool of workers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parpool




tic
%%%% iterate on cells


parfor iexp=1:nloops
    
    [num2str(iexp),'/',num2str(nloops)]
    
    DataExpSmooth = (DataExp(:,iexp))';
    
  
%%%%%%%%%%%%%%% GA optimisation 
%%%%%%%%%%%%%%%%%%%% initial estimate number of polymerases
area= FreqEchImg/Polym_speed*Intensity_for_1_Polym*(TaillePostMarq+TailleSeqMarq/2);
Nbr_poly_estimate = min([round( sum(DataExpSmooth) / area ),num_possible_poly]);

generations=400; %%% number of generations in the genetic algorithm

x_GA_art = optimize_ga1_par(DataExpSmooth,Nbr_poly_estimate,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym, generations);
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
[positions_fit,Min_Fit]=optimize_local1_par(DataExpSmooth,x_GA_art,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym); %%%% set of indices
       
prediction=sumSignal1_par(positions_fit,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym);


%%%%% store predictions
DataPred(:,iexp)=prediction';
cPosPred{iexp}=sort(positions_fit);



Fit(iexp)=Min_Fit;
end %%% iexp
toc



fname = ['result_',name,'_cPosPred',num2str(EspaceInterPolyMin),'.mat'];
save(fname,'DataExp','DataPred','cPosPred','Fit');

delete(gcp('nocreate'))


end %%%% for iii


