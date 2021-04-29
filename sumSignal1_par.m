%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%% computes signal given polymerase positions
% input: transcription start positions(Trans_position), Parameters
% output: sum of signals
% call function: getSignal()
% Parameters = {FreqEchSimu, FreqEchImg,DureeSimu,NbrSondeFluo,...
%            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed};

function [Sum_signals] = sumSignal(Trans_positions,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)
%%%%%%% compute signal from positions 
        
    Sum_signals = zeros(1,frame_num);
    
    ximage = (1:frame_num)/FreqEchImg*Polym_speed; %%%% frame positions in bp
    
    for i = 1:length(Trans_positions)
        
        xpos = Trans_positions(i)/FreqEchSimu*Polym_speed-(TaillePreMarq+TailleSeqMarq+TaillePostMarq); %%%% polymerase start in bp 
        
        ypos = ximage - xpos - TaillePreMarq; %%% frame position relative to signal raise
        
        ind = find(ypos > 0 & ypos < TailleSeqMarq + TaillePostMarq); 
        
        Sum_signals(ind) = Sum_signals(ind) +Signal_par(ypos(ind),Intensity_for_1_Polym,TailleSeqMarq);
    
    end        
        
end 