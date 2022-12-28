% Foreset-Bottom set Transition Model
% Input data (comparable to Get_Data & )

clear all
% Some constants

g = 9.81; 
timeday = 86400;    %sec
timeyr = 31557600;  %sec
ro = 1; % do not yet know what this is  near bed mud concentration. // eqn19 // eqn16 purely depos
Rip = 1/(0.6)^2; %2.78
drforeset =10;  %10
intM=0.35;
pi = 3.1415926; nu = 1e-6;
% alp = 11.2; np = 4.5; taustarsGc = 0.03;	%Constants in Parker approximation of Einstein bedload relation
% aleh = 0.05; neh = 2.5; taustarsSc = 0; 	%Constants in Engelund-Hansen load relation

% Main model parameters
Cza = 20;			% O: Dimensionless Chezy resistance coefficient subaerial
Czs = 20;           % O: Dimensionless Chezy resistnace coefficient subaqueous 30
Ds = 0.1;			% mm: Grain size of sand
Rs = 1.65;          % O: submerged specific gravity of sand  1.65
lps = 0.6;          % Bed porosity, sand-bed reach
Dm = 0.015;			% mm: Grain size of mud 0.04
Rm = 1.65;          % O: submerged specific gravity of mud   1.65
lpm = 0.6;          % Bed porosity, mud-bed reach
Lamms=0.0;
% 


al =  0.05*(Cza)^2;              % coeff in sand load relation = 0.05*(Cza)^2 for E-H
nl = 2.5;               % exponent in sand load relation = 2.5 for E-H

tsc = 0;                % critical Shields stress in load relation = 0 for E-H



%qs 2 xil 7.4 

etabI = 0;              % m: initial elevation of bottom of the foreset
Sfbi = 2E-05;			% O: Slope of fluvial bed 0.00073
Stbi = 0;             % O: Slope of subaqueous basement  2.0E-04
Sa = 0.002;               % O: slope of foreset face  0.2

rsi = 1000;             % m: initial length of fluvial zone
rmax = 50000;       % m: maximum length  12000

dt = 0.000001;   	
% dt = 0.001;   			% yr: Time step 0.0002
Ntoprint =500000;     % N: Number of steps until a printout of results is made  500000
%dt*Ntoprint = ??,   means one printout is for ?? year
 
N =20;			% N: Number of intervals in each zone (excluding ghost node) 20
Ns =100; % N of settling region

Nprint = 20;		% N: Number of printouts after the initial one 15 
xiddot = 5;			% mm/yr: Rate of base level rise 1000

xiddot = xiddot / 1000 / timeyr;
% dt = timeyr * dt/100;
dt = timeyr * dt;

dsbar = 1/N;
dstbar = 1/Ns;

Cfa = 1 / Cza ^ 2;
Cfs = 1 / Czs ^ 2; 

Ds = Ds / 1000;
Dm = Dm / 1000;

M = Ns+N+1;

%Defining parameters
etafluv = zeros(N+1, 1);
sfbar = zeros(N+1, 1);
rfluv = zeros(N+1, 1);

etaturb = zeros(Ns+1, 1);
stbar = zeros(Ns+1, 1);
rturb = zeros(Ns+1, 1);

Hfluv = zeros(N+1, 1);
qs = zeros(N+1, 1);
qm = zeros(Ns+1, 1);
dqm=zeros(Ns+1,1);
Sturb = zeros(Ns+1, 1);
cibar=zeros(Ns+1,1);
cibardev=zeros(Ns+1,1);
deta=zeros(Ns+1,1);
qmdev=zeros(Ns+1,1);
qms=zeros(N+1,1);
topsetlength=zeros(Nprint,1);
foresetlength=zeros(Nprint,1);
eta = zeros(M+1, 1);
r = zeros(M+1, 1);

ystrata = zeros(M+1, Nprint+1);
xstrata = zeros(M+1, Nprint+1);
shoreline = zeros(Nprint+1,1);
fbtransition = zeros(Nprint+1,1);
timeplot = zeros(Nprint+1,1);


etasdot=zeros(Nprint,1);
rbdot=zeros(Nprint,1);
etabdot=zeros(Nprint,1);
rsdot=zeros(Nprint,1);






%%
qw = 5; 			% m2/s: water discharge at flood over width
% qso = 5E-05;	% m2/s: Feed rate of sand %default 0.001
% qmo = qso*0.2; 	% m2/s: Feed rate of mud 0.00050
for etatI=5:5:15
    for sr=80:10:100
       
      
% sr=70;
qt=0.0005; % qso = 8.0124e-04 qmo = 0.0016  // qt=0.0024
qso=qt*sr*0.01;
qmo=qt-qso;

%% Compute Normal Flow
% This subroutine computes the Shields stress, slope and depth of the normal flow
% of the fluvial region associated with the volume sand feed rate per unit wifth qso,
% the water discharge per unit width qw and the sand grain size Ds.

tsn = (qso / ((Rs * g * Ds) ^ 0.5 * Ds * al)) ^ (1 / nl) + tsc;
Sn = Cza * (Rs * Ds * tsn) ^ 1.5 * g ^ 0.5 / qw;
Hn = Rs * Ds * tsn / Sn;


% etatI = 3;        % m: initial elevation of top of the foreset  70
xil =etatI+Hn;              % m: water surface elevation o fht lake 73
Hturb = etatI*ones(Ns+1, 1);
intM=0.35;







poster_0='DimentsionLess_depth';
poster_1='Dimension_depth';
poster1='Sratio';
figurename1=sprintf('%s %d %s %d',poster_0,etatI,poster1,sr); % dimensionless
figurename2=sprintf('%s %d %s %d',poster_1,etatI,poster1,sr); % dimension 




washload=0;
capturemud=0;
fbt_delta_settling_Loop
% close all
%% results of main code are saved here 
results.shoreline=shoreline;
results.etashoreline=etashoreline;
results.fbb=fbb;
results.etafbb=etafbb;
results.rsdot=rsdot;
results.etasdot=etasdot;
results.rbdot=rbdot;
results.etabdot=etabdot;
results.finalshoreline=shoreline(Nprint+1);
% results.washload=washload;
results.topsetlength=topsetlength;
results.foresetlength=foresetlength;
% results.topsetlength=topsetlength;

R.r(etatI,sr) = results;


    end
end
save('datat.mat','R')