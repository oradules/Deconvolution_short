%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
function S=Signal(ypos,Intensity_for_1_Polym,TailleSeqMarq) 
%%%%%% signal from one polymerase
    S = ones(size(ypos))*Intensity_for_1_Polym;
    ind2=find(ypos < TailleSeqMarq);
    S(ind2) = ypos(ind2)/TailleSeqMarq*Intensity_for_1_Polym;        
end 