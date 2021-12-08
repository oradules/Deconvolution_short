%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%% computes signal given polymerase positions
% input: transcription start positions(Trans_position) in min spacings, Parameters
% output: sum of signals
% call function: getSignal()
% Parameters = {FreqEchSimu, FreqEchImg,DureeSimu,NbrSondeFluo,...
%            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed};

function [Sum_signals] = sumSignal(Trans_positions,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)
%%%%%%% compute signal from positions 
    Taille = (TaillePreMarq+TailleSeqMarq+TaillePostMarq);
    Sum_signals_matrix = zeros(frame_num,length(Trans_positions));
    ximage=(repelem([1:frame_num],length(Trans_positions),1)')/FreqEchImg*Polym_speed; %%%% frame positions in bp
    xpos=(Trans_positions/FreqEchSimu)*Polym_speed-Taille;
      t1=repelem((xpos+TaillePreMarq)',frame_num,1);
      
      ypos=ximage-t1;
      ind=(ypos > 0)&(ypos < (TailleSeqMarq + TaillePostMarq));
      
    Sum_signals_matrix(ind) = Sum_signals_matrix(ind) + Signal_par(ypos(ind),Intensity_for_1_Polym,TailleSeqMarq);
    Sum_signals=sum(Sum_signals_matrix');
    
end
