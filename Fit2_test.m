%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier
%%%%% Copyright : This is published under 3-clause BSD
%%%%% last change november 2021
%%%%% performs 2 exp fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


global FreqEchSimu FreqEchImg DureeAnalysee TaillePreMarq ...
            TailleSeqMarq TaillePostMarq  Polym_speed frame_num num_possible_poly EspaceInterPolyMin ...
            DureeSimu Intensity_for_1_Polym;
global dt;
tic

%%%% load parameters
[ Polym_speed, EspaceInterPolyMin,TaillePreMarq,TailleSeqMarq,TaillePostMarq,EspaceInterPolyMin,FrameLen,Intensity_for_1_Polym] = parameters_droso();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FreqEchImg = (1/FrameLen); % image per second data time sampling   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); % how many interval(possible poly start position) in 1s
DureeSignal = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed; % (s)

fsz=16;lw=2;
sz=10;

%%% list of movies to process
names={{'krINRmut_2states'},...
    {'Sna+INR_M2'},...
    {'SnaTATALight+INR_M2'},...
    {'SnaTATAlight_2states'},...
    {'Sna_2states'},...
    {'Kr-INRBrk_2states'},...
    {'Kr_M2'},...
    {'SnaTATAMut_2states'},...
    {'eef1a1wt_2states'},...
    {'eef1a1G4mut_M2'},...
    {'eef1a1TATAmut_M2'},...
    {'HIVwt_M2'},...
    {'HIVlo_M2'},...
    {'HIVno_M2'}};

names={{'krINRmut_2states'},...
    {'Sna+INR_M2'},...
    {'SnaTATALight+INR_M2'},...
    {'SnaTATAlight_2states'},...={{'Kr_M2'},{'Sna_2states'}};

    {'Sna_2states'},...
    {'Kr-INRBrk_2states'},...
    {'Kr_M2'},...
    {'SnaTATAMut_2states'}};

names
likelihood=0; %%%%% if this is true use max likelihood for parametric survival function fit

if likelihood
    DataFilePath0 = '../modelfit_results_likelihood'; %%%%% where to write fit results
else
    DataFilePath0 = '../modelfit_results_test'; %%%%% where to write fit results    
end
DataFilePath='../matresultfiles/'; %%%%% where are the deconvolution results

mkdir(DataFilePath0);
xlsfilename = [DataFilePath0,'/results.xlsx'];
%xlswrite(xlsfilename,{'Data','k1p','k1m','k2','lambda1','lambda2','A1','A2','p1','p2','Obj','KS test','Nuclei','Frames'},1,'A1');

writecell({'Data','lambda1','lambda2','A1','A2','k1p','k1m','k2','Obj','KS test'},xlsfilename,'Sheet',1,'Range',['A1']);
res=[];
par_noise=[];
for ifile = 1:length(names)
    
    close all

lnames = names{ifile};
 


if length(lnames) > 1
%%%%% lump files, Procustes method 
dataExp=[];
dataPred=[];
posPred=[];
tmax=[];
%%%% first compute max dimension
nmax=1;nmaxpos=1;
for iii=1:length(lnames)
    fname = [DataFilePath,'result_',lnames{iii},'_artificial_00.mat'];
    load(fname);
    n2=size(DataExp);
    n3=size(PosPred);
    if n2(1) > nmax
        nmax=n2(1);
    end
    if n3(1) > nmaxpos
        nmaxpos=n3(1);
    end
end



for iii=1:length(lnames)
    fname = [DataFilePath,'result_',lnames{iii},'_artificial_00.mat'];
    load(fname);
    
    n2=size(DataExp);
    n3=size(PosPred);    
    
%%%%%%%%%%%%%% padding 
DataExp=[DataExp;zeros(nmax-n2(1),n2(2))];
DataPred=[DataPred;zeros(nmax-n2(1),n2(2))]; 
PosPred=[PosPred;zeros(nmaxpos-n3(1),n3(2))];
    

%%%%%%%%%%%%% concatenating
     dataExp=[dataExp,DataExp]; 
     dataPred=[dataPred,DataPred];
     posPred=[posPred,PosPred];
        
        
 tmax=[tmax;(n2(1)-1)/FreqEchImg*ones(n2(2),1)]; %%%% movie length       
        
end

%%Maria: should we also extract the file name like that:
%%%%% extract short name from result file name
%iend = strfind(fname,'_00');
%name = fname(8:iend-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ovidiu: No, because names contains also the information on how to 
%% lump the files


%%%%% padded and concatenated data%%%%
DataExp=dataExp;
DataPred=dataPred;
PosPred=posPred;
else
 fname = [DataFilePath,'result_',lnames{1},'_artificial_00.mat'];
    load(fname);   
 n2=size(DataExp);
 n3=size(PosPred);
 tmax= (n2(1)-1)/FreqEchImg*ones(n2(2),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% if only one file per phenotype this is just the containt of the file
%%%% with no padding and concatenation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end











%%%% name of the phenotype
name=lnames{1}; %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


dirwrite=[DataFilePath0,'/',name,'_result'];
mkdir(dirwrite);

n=size(DataExp);
nexp=n(2);
% -------Parameters----------
DureeSimu = n(1)*FrameLen; %%% in s
frame_num = n(1);

DureeAnalysee = DureeSignal + DureeSimu ; % (s)
num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)); 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=[];  %%%%% will contain the start of the analyzed region
%%%%%%%%%%%% 
for data_i=1:nexp

max_intensity=max(DataPred(:,data_i));
%%%% find first hit
ihit=min(find(DataPred(:,data_i) > max_intensity/5 ));
if isempty(ihit)
    ihit = n(1);
end

ihit=0; %%%%% for artificial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=(ihit-1)/FreqEchImg; %%%% t0 
T0=[T0,t0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end





%%%% show data, check that positions predict data
isel_vis=[4,13,115];
isel_vis=1:10;
n=size(DataExp); 
stime=(0:n(1)-1)*FrameLen;
frame_num=n(1);
for i=1:length(isel_vis)
    subplot(length(isel_vis),1,i)
    hold off
 plot(stime,DataExp(:,isel_vis(i)),'r')
 hold on
 plot(stime,DataPred(:,isel_vis(i)),'k')
positions   =  find(PosPred (:,isel_vis(i)) == 1);  
prediction=sumSignal1_par(positions,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym);
plot(stime,prediction,'xk')

bppositions =  positions/FreqEchSimu*Polym_speed-(TaillePreMarq+TailleSeqMarq+TaillePostMarq); %%% positions in bp
fpositions  =  bppositions /Polym_speed; %%%% positions in minutes
m= max(DataPred(:,isel_vis(i)));
for j=1:length(positions)
  plot([fpositions(j),fpositions(j)],[0,m/10],'c')
  hold on
end
times = positions / FreqEchSimu -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed ; %%% starting times of polymerases in seconds
%    eliminate events that don't have effect between T0 and tmax
times = times ( times + (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed > T0(isel_vis(i)) & ...
  times + (TaillePreMarq)/Polym_speed < tmax(isel_vis(i)) );  
for j=1:length(times)
  plot([times(j),times(j)],[m/10,m/10],'co')
  hold on
end
%%%% add T0 and tmax
plot([T0(isel_vis(i)),T0(isel_vis(i))],[0,m],'k','linewidth',2)
plot([tmax(isel_vis(i)),tmax(isel_vis(i))],[0,m],'k','linewidth',2)
%%%% add signal from individual polymerases
for ipol=1:length(positions)
    position=positions(ipol);
    onepol=onepolSignal1_par(position,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym);
    plot(stime,onepol,'linewidth',1,'color',[82, 190, 128]/256)  
end
% axis([0, max(stime), 0, m]);
 if i==3
     xlabel('Time [min]','fontsize',fsz)
 end
 if i==2
     ylabel('Florescence intensity','fontsize',fsz)
 end
 set(gca,'fontsize',fsz)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% compute distribution of spacings %%%%%%%%%%%%%
nn=size(PosPred);
dt=[]; dtc=[]; 
for i=1:nn(2) %%% for all cells  
  times = find (PosPred (:,i) == 1) / FreqEchSimu -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed ; %%% starting times of polymerases in seconds
%    eliminate events that don't have effect between T0 and tmax
  times = times ( times + (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed > T0(i) & ...
  times + (TaillePreMarq)/Polym_speed < tmax(i) );  
lt = length(times);  
    if lt > 1
        dtimes = diff(times);
        dt=[dt;dtimes];
    if tmax(i) > times(end) 
        dtc=[dtc;tmax(i)-times(end)]; %%%% the very last time
    end
    else
       if lt == 1
            dtc=[dtc;tmax(i)-times(end)]; %%%% the very last time
       end
       if lt==0
            dtc=[dtc;tmax(i)-T0(i) + (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed]; %%%% the very last time
       end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%% define matrices that store parameters and objective functions
store=[];

if ~isempty(dt)

for cens=0:1

if cens
    [fs,xs,flo,fup]=ecdf([dt;dtc],'censoring',[zeros(size(dt));ones(size(dtc))]);
else
    [fs,xs,flo,fup]=ecdf(dt);
end
%%%% fit distribution of spacings using combination of two exponentials
xs=xs(1:end-1);
fs=fs(1:end-1);
flo=flo(1:end-1);
fup=fup(1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%

sN = sqrt(length(xs));

%%%% try several cost functions
if ~likelihood
%%%% 1 : log scale only
exp_fitness = @(k) (abs(log(k(3)*exp(k(1)*xs)+(1-k(3))*exp(k(2)*xs))-log(1-fs)))/sN; % k: parameters
opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, ...
        'MaxFunEvals',1e6,'TolX', 1e-10);
%%%% 2 : mixture linear/log with parameter alpha = 0.6
%alpha=0.6;
%NS = length(xs);fact1=sqrt((1-alpha)/NS);fact2=sqrt(alpha/NS);eshort=(1-fs);
%exp_fitness = @(k)( [ log( ( k(3)*exp(k(1)*xs)+(1-k(3))*exp(k(2)*xs) )./eshort )* fact1;...
%                    ( k(3)*exp(k(1)*xs)+(1-k(3))*exp(k(2)*xs)- eshort )*fact2]); %%%%% mixed
%opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, ...
%        'MaxFunEvals',1e6,'TolX', 1e-10);                
else
%%%% 3 : log likelihood, directly on spacings
dt = [dt;dtc];
%exp_fitness = @(k) (sum( - log(  -k(3)*k(1)*exp(k(1)*dt) - (1-k(3))*k(2)*exp(k(2)*dt) )));
options = optimoptions('fmincon','GradObj','on','TolFun', 1e-8,'TolX', 1e-10);
%%% programs in C:\Users\ovidiu\Dropbox\artificial_data_short\test1_artificial_maxlikelihood_2states.m
end

k00=[-0.01,-0.001,1];    
amp = [log(100),log(100)]; 
NbIterationinFit=100;    

if likelihood
    lb=[-inf;-inf;0];ub=[0;0;1]; %%%% combined with fmincon
end


    for mc = 1:NbIterationinFit
        
        %%%% Change k00
        factor=exp(amp.*(2*rand(1,2)-1)); 
        k0 = k00;
        k0(1:2)=k0(1,2).*factor;
        k0(3) =  2*rand(1,1)-1; %%% A1
        %%%% sort k0(1:2)
        k0(1:2)=sort(k0(1:2),'ascend');
        
        %%% impose constraints
        if ~( k0(1)*k0(3)+k0(2)*(1-k0(3)) < 0 )
            while ~( k0(1)*k0(3)+k0(2)*(1-k0(3)) < 0 )
                k0(3)=2*rand(1,1)-1;  %%% A1 
            end
        end
         
if ~likelihood        
        %%%% Use the fcn lsqnonlin
            [k, obj] = lsqnonlin(exp_fitness,k0,[],[],opts);
        %%%% of use fmin con
else
        [k, obj] = fmincon(@maxlikelihood2withgradient,k0,[],[],[],[],lb,ub,[],options);
end        
        %%%% sort k
        A=[k(3),1-k(3)]; %%% A values before sorting 
        [kk,IX]=sort(k(1:2),'ascend'); %%% sort lambdas
        k(1:2)=kk;
        A=A(IX);
        k(3)=A(1);
        
        store = [store;[k, obj, cens]];
        
        disp(mc)
        
    
    end    
 


end % cens 


%%%% select optimal and suboptimal obj < 2 objmin
ind = find  (max(abs(imag(store(:,1:3))),[],2) < 1e-10);
[objmin,indmin]=min(store(ind,4));
imin=ind(indmin);
overflow=1;
ind=find(store(:,4) < (1+overflow)*objmin & max(abs(imag(store(:,1:3))),[],2) < 1e-10 );
ksel=real(store(ind,1:3));
kmin=real(store(imin,1:3) );
censmin=store(imin,5);

if censmin
    [fs,xs,flo,fup]=ecdf([dt;dtc],'censoring',[zeros(size(dt));ones(size(dtc))]);
else
    [fs,xs,flo,fup]=ecdf([dt;dtc]);% no censoring
end

%%%%% plot survival function 
h=figure(70)
hold off
semilogy(xs,1-fs,'or') %%% empirical function
hold on
semilogy(xs,1-flo,'--r') %%%% lower confidence 
semilogy(xs,1-fup,'--r') %%%% upper confidence
pred=kmin(3)*exp(kmin(1)*xs)+(1-kmin(3))*exp(kmin(2)*xs);
semilogy(xs,pred,'k','linewidth',2) %%% predicted 2 exp
axis([0, 250, 1e-6, 1])
xlabel('Time [s]','fontsize',fsz)
ylabel('Survival function','fontsize',fsz)


%%%% compute 3 rates k1p,m k2 from the 3 parameters %%%
S1 = kmin(3)*kmin(1)+(1-kmin(3))*kmin(2);
k2 = -S1; 
S2 = kmin(3)*(kmin(1))^2+(1-kmin(3))*(kmin(2))^2;
S3 = kmin(3)*(kmin(1))^3+(1-kmin(3))*(kmin(2))^3;
k1m = S1-S2/S1; 
k1p = (S3*S1-S2^2)/S1/(S1^2-S2); 
title(['k_1^-=',num2str(k1m,1),'k_1^+=',num2str(k1p,1),'k_2=',num2str(k2,1)])




figfile=[dirwrite,'/Fit_',lnames{1},'.pdf'];
print(h,'-dpdf',figfile)



%%%%%%% KS test

thr= 20;
ind = find(xs>thr); %%% perform the test on times larger than 20 s
dist = max(abs( fs(ind)-1+pred(ind)));
nsample=length(dt(dt>thr));
N=10;
c=sqrt(nsample)*dist;
r=[1:N];
aa=2*sum(exp(-2*c.^2*r.^2).*((-1).^(r-1)));

h=figure(80)
hold off
semilogy(xs,fs,'kx')
hold on
semilogy(xs,1-pred,'ro')
xlabel('Time [s]','fontsize',fsz)
ylabel('CDF','fontsize',fsz)
figfile=[dirwrite,'/Fit_CDF_',lnames{1},'.pdf'];
print(h,'-dpdf',figfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%% compute intervals
S1 = ksel(:,3).*ksel(:,1)+(1-ksel(:,3)).*ksel(:,2);
K2 = -S1; 
S2 = ksel(:,3).*(ksel(:,1)).^2+(1-ksel(:,3)).*(ksel(:,2)).^2;
S3 = ksel(:,3).*(ksel(:,1)).^3+(1-ksel(:,3)).*(ksel(:,2)).^3;
K1m = S1-S2./S1; 
K1p = (S3.*S1-S2.^2)./S1./(S1.^2-S2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nn=size(DataExp);
%%%% optimal
%res= [k1p,k1m,k2,kmin(1),kmin(2),kmin(3),1-kmin(3),k1m/(k1m+k1p),k1p/(k1m+k1p),objmin,aa,nn(2),nn(1)];

res= [kmin(1),kmin(2),kmin(3),1-kmin(3),k1p,k1m,k2,objmin,aa];

P1=K1m./(K1m+K1p);
P2=K1p./(K1m+K1p);

%resl= [min(K1p),min(K1m),min(K2),min(ksel(:,1)),min(ksel(:,2)),min(ksel(:,3)),min(1-ksel(:,3)),min(P1),min(P2)];

%resl = max([resl;zeros(1,length(resl))]); 

resl= [min(ksel(:,1)),min(ksel(:,2)),min(ksel(:,3)),min(1-ksel(:,3)),min(K1p),min(K1m),min(K2)];

%resh= [max(K1p),max(K1m),max(K2),max(ksel(:,1)),max(ksel(:,2)),max(ksel(:,3)),max(1-ksel(:,3)),max(P1),max(P2)];

resh= [max(ksel(:,1)),max(ksel(:,2)),max(ksel(:,3)),max(1-ksel(:,3)),max(K1p),max(K1m),max(K2)];

% xlswrite(xlsfilename,{strrep(name,'result_','')},1,['A',num2str(4*ifile-2)]); %%% filename
% 
% xlswrite(xlsfilename,res,1,['B',num2str(4*ifile-2)]); %%% best result
% xlswrite(xlsfilename,resl,1,['B',num2str(4*ifile-1)]); %%% low
% xlswrite(xlsfilename,resh,1,['B',num2str(4*ifile)]); %%%  high

writecell({strrep(name,'result_','')},xlsfilename,'Sheet',1,'Range',['A',num2str(4*ifile-2)]); %%% filename

writematrix(res,xlsfilename,'Sheet',1,'Range',['B',num2str(4*ifile-2)]); %%% best result
%xlswrite(xlsfilename,res2,1,['B',num2str(4*ifile-1)]); %%% second best result
writematrix(resl,xlsfilename,'Sheet',1,'Range',['B',num2str(4*ifile-1)]); %%% low
writematrix(resh,xlsfilename,'Sheet',1,'Range',['B',num2str(4*ifile)]); %%%  high


else
    
    
xlswrite(xlsfilename,{strrep(name,'result_','')},1,['A',num2str(4*ifile-2)]); %%% filename
nn=size(DataExp);
%%%% 
xlswrite(xlsfilename,res,1,['B',num2str(4*ifile-2)]); %%% best result
    
end

end % ifile
toc
