function [Pol_speed,SpaceInterPolyMin,Pre,Seq,Post,EspaceInter,FrameLen,Intensity_for_1_Polym]=parameters_droso();
%%%%%%%%%%%%%%%%%%%%%%%%%%% -------Parameters DROSO
Pol_speed = 45; % average speed bases per second (Ref publi)
SpaceInterPolyMin = 30; % in bp 
Pre = 41; % 700 bases (ref publi)
Seq = 1292; % 2900 bases (ref publi)
Post = 4526 + Pol_speed*0; % 1600 bases + 50s polya signal
EspaceInter = 30; % en base (40 bases)
FrameLen = 3.9; %%% frame length in seconds
Intensity_for_1_Polym = 1; 
end

