%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% Copyright : This is published under 3-clause BSD
%%%%% Last change March 2021
%%%%% computes polymerase positions by deconvolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%% -------Parameters DROSO
Polym_speed = 45; % average speed bases par seconde (Ref publi)
TaillePreMarq = 41; % 700 bases (ref publi)
TailleSeqMarq = 1292; % 2900 bases (ref publi)
TaillePostMarq = 4526 + Polym_speed*0; % 1600 bases + 0s polya signal
EspaceInterPolyMin = 30; % minimum inter POLII space in bp
FrameLen = 3.86; %%% frame length in seconds
Intensity_for_1_Polym = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% hardware parameters
nworkers = 50;  %%% max number of workers (less than the number of cores on your computer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% all the generated mat files will be placed here
MatFilePath = '../matresultfiles/';
mkdir(MatFilePath);

%%%% load the .mat files generated by Read_data
%%% where are the files
DataFilePath = '../matdatafiles/';

data_list=struct2cell(dir(fullfile(DataFilePath)));
file_name_list = data_list(1,3:end);

n_exp=length(file_name_list);
    
for iii=1:n_exp
    
    
    
fnameinput = file_name_list{iii};


['data',fnameinput,' ',num2str(iii),'/',num2str(n_exp)] 

%%%% select only those mat files whose name are prefixed by "data_" 
if strfind(fnameinput,'data_') == 1

%%% full path input file name
ffname = [DataFilePath,fnameinput];

load(ffname);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sd=size(DataExp); %%%%% size matrix data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_num=sd(1); %%%%%% number of frames
FreqEchImg = (1/FrameLen); %%%%% data time sampling
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


nloops = sd(2);

%%%%%% parallel computing %%%%%%%%%%%%%%%%%%%
mycluster = parcluster; %%% create cluster
set(mycluster,'NumWorkers',nworkers) ;
parpool(mycluster,min([nloops,nworkers])); %%% open a pool of workers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
%%%% iterate on cells
parfor iexp=1:nloops

    
    [num2str(iexp),'/',num2str(nloops)]
    
    DataExpSmooth = (DataExp(:,iexp))';
    
  
%%%%%%%%%%%%%%% GA optimisation 
%%%%%%%%%%%%%%%%%%%% estimate number of poly
area= FreqEchImg/Polym_speed*Intensity_for_1_Polym*(TaillePostMarq+TailleSeqMarq/2);
Nbr_poly_estimate = min([round( sum(DataExpSmooth) / area ),num_possible_poly]);

generations=400; %%%% number of generation in the GA step

x_GA_art = optimize_ga1_par(DataExpSmooth,Nbr_poly_estimate,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym, generations);
     
     
%%%%%%%%%%%%%%%% local optimisation, version without overlap    
[positions_fit,Min_Fit]=optimize_local1_par(DataExpSmooth,x_GA_art,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym); %%%% set of indices
%Min_Fit        
prediction=sumSignal1_par(positions_fit,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym);

%%%%% store predictions
DataPred(:,iexp)=prediction';
%PosPred(:,iexp)=column;
cPosPred{iexp}=positions_fit;



Fit(iexp)=Min_Fit;  
end %%% iexp
toc


for i=1:sd(2)
    PosPred(cPosPred{i},i)=1;
end

%%%%% extract short name from file name
iend=strfind(fnameinput,'Calibrated');
name = fnameinput(6:iend-1);


fnameoutput = [MatFilePath,'result_',name,'_00.mat'];
save(fnameoutput,'DataExp','DataPred','PosPred','Fit');

delete(gcp('nocreate'))

end %%%% if strfind
end %%%% for iii


