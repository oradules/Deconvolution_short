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
    {'SnaTATAlight_2states'},...
    {'Sna_2states'},...
    {'Kr-INRBrk_2states'},...
    {'Kr_M2'},...
    {'SnaTATAMut_2states'}};

names={{'Kr_M2'},{'Sna_2states'}};


likelihood=0; %%%%% if this is true use max likelihood for parametric survival function fit

if likelihood
    DataFilePath0 = '../modelfit_results_3States_likelihood'; %%%%% where to write fit results
else
    DataFilePath0 = '../modelfit_results_3States_testV'; %%%%% where to write fit results    
end
DataFilePath='../matresultfiles/'; %%%%% where are the deconvolution results

mkdir(DataFilePath0);
xlsfilename = [DataFilePath0,'/results.xlsx'];
%xlswrite(xlsfilename,{'Data','k1p','k1m','k2','lambda1','lambda2','A1','A2','p1','p2','Obj','KS test','Nuclei','Frames'},1,'A1');

writecell({'Data','lambda1','lambda2','lambda3','A1','A2','A3',...
    'M1:k1p','M1:k1m','M1:k2p','M1:k2m','M1:k3','M1:p1','M1:p2','M1:p3',...
    'M2:k1p','M2:k1m','M2:k2p','M2:k2m','M2:k3','M2:p1','M2:p2','M2:p3',...
    'Obj','KS test'},xlsfilename,'Sheet',1,'Range',['A1']);
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
prediction=sumSignalv_par(positions,FreqEchSimu, FreqEchImg, TaillePreMarq, ...
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
exp_fitness = @(k) (abs(log(k(4)*exp(k(1)*xs)+k(5)*exp(k(2)*xs)+(1-k(4)-k(5))*exp(k(3)*xs))-log(1-fs)))/sN; % k: parameters

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

%%%Maria: here two states and three states don't deffer
dt = [dt;dtc];
%exp_fitness = @(k) (sum( - log(  -k(3)*k(1)*exp(k(1)*dt) - (1-k(3))*k(2)*exp(k(2)*dt) )));
options = optimoptions('fmincon','GradObj','on','TolFun', 1e-8,'TolX', 1e-10);
%%% programs in C:\Users\ovidiu\Dropbox\artificial_data_short\test1_artificial_maxlikelihood_2states.m
end

k00=[-0.1,-0.01,-0.001,0.25,0.25];    
amp = [log(100),log(100),log(100),log(1),log(1)]; 
NbIterationinFit=100;    

if likelihood
    lb=[-inf;-inf;0];ub=[0;0;1]; %%%% combined with fmincon
end


    for mc = 1:NbIterationinFit
        
         %%%% Change k00
        factor=exp(amp.*(2*rand(1,5)-1)); 
        
        k0 = k00.*factor;
        k0(4:5)=2*rand(1,2)-1; 
        
        %%%% sort k0(1:3)
        k0(1:3)=sort(k0(1:3),'ascend');
        
        if ~( sum(k0(4:5)) < 1 && k0(1)*k0(4)+k0(2)*k0(5)+k0(3)*(1-sum(k0(4:5))) < 0 )
            while ~(sum(k0(4:5)) < 1 && k0(1)*k0(4)+k0(2)*k0(5)+k0(3)*(1-sum(k0(4:5))) < 0)
                k0(4:5)=2*rand(1,2)-1;  %%% A1,A2 values
            end
        end    
         
if ~likelihood        
        %%%% Use the fcn lsqnonlin
            [k, obj] = lsqnonlin(exp_fitness,k0,[],[],opts);
        %%%% of use fmin con
else
    %Maria: what is maxlikelihood2withgradient
    [k, obj] = fmincon(@maxlikelihood2withgradient,k0,[],[],[],[],lb,ub,[],options);
end        
        %%%% sort k
        A=[k(4),k(5),1-k(4)-k(5)]; %%% A values before sorting 
        [kk,IX]=sort(k(1:3),'ascend'); %%% sort lambdas
        k(1:3)=kk;
        A=A(IX);
        k(4:5)=A(1:2);
        
        store = [store;[k, obj, cens]];
        
        disp(mc)
        
    
    end    
 


end % cens 


%%%% select optimal and suboptimal obj < 2 objmin
ind = find  (max(abs(imag(store(:,1:3))),[],2) < 1e-10); % takes the part where we have the imaginary part of lambda i almost 0
[objmin,indmin]=min(store(ind,end-1)); %minimum of the objectives for all the previous lambdas
imin=ind(indmin); % finding the positions of the lambda where we have the lowest objective
overflow=1;
ind=find(store(:,end-1) < (1+overflow)*objmin & max(abs(imag(store(:,1:3))),[],2) < 1e-10 ); % taking suboptimal objectives such that they are <2*optimal objective and they satisfy the fact that the imaginary part of lambda is almost 0
ksel=real(store(ind,1:5)); % taking an array of lambda i's and Ai's where we have the the objective function less than <2*times the minimum objective function
kmin=real(store(imin,1:5) ); %taking the lambda i's and A i's that are the fittest
censmin=store(imin,end);


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
pred=kmin(4)*exp(kmin(1)*xs)+kmin(5)*exp(kmin(2)*xs)+(1-kmin(4)-kmin(5))*exp(kmin(3)*xs);
semilogy(xs,pred,'k','linewidth',2) %%% predicted 2 exp
axis([0, 250, 1e-6, 1])
xlabel('Time [s]','fontsize',fsz)
ylabel('Survival function','fontsize',fsz)


%%%% compute 5 rates k1p,m k2p,m k3 from the 5 parameters %%%
l1=kmin(1);
l2=kmin(2);
l3=kmin(3);
A1=kmin(4);
A2=kmin(5);
A3=1-A1-A2;
L1=l1+l2+l3;
L2=l1.*l2+l1.*l3+l2.*l3;
L3=l1.*l2.*l3;
S1=A1.*l1+A2.*l2+A3.*l3;
S2=A1.*l1.^2+A2.*l2.^2+A3.*l3.^2;
S3=A1.*l1.^3+A2.*l2.^3+A3.*l3.^3;
%%%%% Model M2
kk3= -S1;
kk1p = 1/2 * ( -L1+S2/S1 + sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
kk2p = 1/2 * ( -L1+S2/S1 - sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
kk1m = 1/2 * (S1-S2/S1 - (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));
kk2m = 1/2 * (S1-S2/S1 + (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));
pp1=kk1m*kk2p/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m);
pp2=kk1p*kk2m/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m);
pp3=kk1p*kk2p/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m);
%%%% model M1
k1p=-L3.*(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1p
k2p=-(S2.^2-S1.*S3)./S1./(S1.^2-S2); %%% k2pk
k3=-S1; %%%% k3
k2m=(S1.^2-S2)./S1; %%% k2m
k1m=-A1.*A2.*A3.*(l1-l2).^2.*(l1-l3).^2.*(l2-l3).^2.*S1./(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1m
p1=k1m*k2m/(k1p*k2m+k1m*k2m+k1p*k2p);
p2=k1p*k2m/(k1p*k2m+k1m*k2m+k1p*k2p);
p3=k1p*k2p/(k1p*k2m+k1m*k2m+k1p*k2p);

title(['k_1^-=',num2str(k1m,1),'k_1^+=',num2str(k1p,1),'k_2^-=',num2str(k2m,1),...
    'k_2^+=',num2str(k2p,1),'k_3=',num2str(k3,1)])



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

nn=size(DataExp);
%%%% optimal
res= [kmin(1),kmin(2),kmin(3),kmin(4),kmin(5),1-kmin(5),k1p,k1m,k2p,k2m,k3, p1,p2,p3,...
    kk1p,kk1m,kk2p,kk2m,kk3,pp1,pp2,pp3,objmin,aa];




%%%%% compute intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1=ksel(:,1);l1sel=l1;
l2=ksel(:,2);l2sel=l2;
l3=ksel(:,3);l3sel=l3;
A1=ksel(:,4);A1sel=A1;
A2=ksel(:,5);A2sel=A2;
A3=1-A1-A2;A3sel=A3;
L1=l1+l2+l3;
L2=l1.*l2+l1.*l3+l2.*l3;
L3=l1.*l2.*l3;
S1=A1.*l1+A2.*l2+A3.*l3;S1sel=S1;
S2=A1.*l1.^2+A2.*l2.^2+A3.*l3.^2;
S3=A1.*l1.^3+A2.*l2.^3+A3.*l3.^3;
    
%%%% model M2
KK1p = 1/2 * ( -L1+S2./S1 + sqrt((S1.*L1-S2).^2-4*L3.*S1)./S1 );
KK2p = 1/2 * ( -L1+S2./S1 - sqrt((S1.*L1-S2).^2-4*L3.*S1)./S1 ); 
KK1m = 1/2 * (S1-S2./S1 - (-S1.^2.*L1+S1.*S2+S1.*L2-L3+S2.^2./S1-S3)./sqrt((S1.*L1-S2).^2-4*L3.*S1)); 
KK2m = 1/2 * (S1-S2./S1 + (-S1.^2.*L1+S1.*S2+S1.*L2-L3+S2.^2./S1-S3)./sqrt((S1.*L1-S2).^2-4*L3.*S1)); 
KK3=-S1; %%%% k3
PP1=KK1m.*KK2p./(KK1p.*KK2p+KK1m.*KK2p+KK1p.*KK2m);
PP2=KK1p.*KK2m./(KK1p.*KK2p+KK1m.*KK2p+KK1p.*KK2m);
PP3=KK1p.*KK2p./(KK1p.*KK2p+KK1m.*KK2p+KK1p.*KK2m);
%%%% model M1
K1p=-L3.*(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1p
K2p=-(S2.^2-S1.*S3)./S1./(S1.^2-S2); %%% k2p
K3=-S1; %%%% k3
K2m=(S1.^2-S2)./S1; %%% k2m
K1m=-A1.*A2.*A3.*(l1-l2).^2.*(l1-l3).^2.*(l2-l3).^2.*S1./(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1m
P1=K1m.*K2m./(K1p.*K2m+K1m.*K2m+K1p.*K2p);
P2=K1p.*K2m./(K1p.*K2m+K1m.*K2m+K1p.*K2p);
P3=K1p.*K2p./(K1p.*K2m+K1m.*K2m+K1p.*K2p);


resl= [min(ksel(:,1)),min(ksel(:,2)),min(ksel(:,3)),min(ksel(:,4)),min(ksel(:,5)) ,min(1-ksel(:,5)),...
    min(K1p),min(K1m),min(K2p),min(K2m),min(K3),min(P1),min(P2),min(P3),...
    min(KK1p),min(KK1m),min(KK2p),min(KK2m),min(KK3),min(PP1),min(PP2),min(PP3)];


resh= [max(ksel(:,1)),max(ksel(:,2)),max(ksel(:,3)),max(ksel(:,4)),max(ksel(:,5)) ,max(1-ksel(:,5)),...
    max(K1p),max(K1m),max(K2p),max(K2m),max(K3),max(P1),max(P2),max(P3),...
    max(KK1p),max(KK1m),max(KK2p),max(KK2m),max(KK3),max(PP1),max(PP2),max(PP3)];
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
