
% Foreset-Bottom set Transition Model
% Input data (comparable to Get_Data & )

clear all
% Some constants

% 2D parameters
intM=0.35; % intermettency
Bc = 300; % river width = 100;
thetaF = 120; % delta fan angle, opening angle [degree]
thetaF = thetaF / 360 * (2 * pi() );
thetaP = 60; % advection settling plume angle [degree]
thetaP = thetaP / 360 * (2 * pi() );
Qbfw=4100; % m3/s
Qbfs=38.40968; % Mt/yr
Fsand=0.183175;
qw=Qbfw/Bc; % m2/s
annualBedloadSand = Fsand * Qbfs;
annualBedloadSand = annualBedloadSand/(1.65+1)*1000000/(60*60*24*365.25);
FloodBedload = annualBedloadSand/intM;
Ndiversion = 1; Fdiversion=1;
Qso=FloodBedload/Ndiversion/Fdiversion;
qso=Qso/Bc; % m2/sec
% qso=Qbfs*(10^12)/31557600/2.65/1000000/Bc*Fsand;
g = 9.81; 
timeday = 86400;    %sec ( of day )
timeyr = 31557600;  %sec (of year)
ro = 1; % do not yet know what this is  near bed mud concentration. // eqn19 // eqn16 purely depos
Rip = 1/(0.6)^2; %2.78
drforeset =10;  %10
 
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
Lamms=0.49;            % mud captured per sand 

al =  0.05*(Cza)^2;              % coeff in sand load relation = 0.05*(Cza)^2 for E-H
nl = 2.5;               % exponent in sand load relation = 2.5 for E-H

tsc = 0;                % critical Shields stress in load relation = 0 for E-H



%qs 2 xil 7.4 

Sfbi = 2E-05;			% O: Slope of fluvial bed 0.00073
% Stbi = 1.8E-04;             % O: Slope of subaqueous basement  2.0E-04
Stbi=0.00018;
Sa = 0.002;               % O: slope of foreset face  0.2

rsi = 4300;             % m: initial length of fluvial zone
rmax = 25000;       % m: maximum length  12000

% dt = 0.00000001;   	% yr
dt = 0.000002*365.25; % day % 4 instead of 2
% dt = 0.001;   			% yr: Time step 0.0002
Ntoprint =2500000;     % N: Number of steps until a printout of results is made  1000000
%dt*Ntoprint = ??,   means one printout is for ?? year
 
N =20;			% N: Number of intervals in each zone (excluding ghost node) 20
Ns =20; % N of settling region
Nprint =20;		% N: Number of printouts after the initial one 15 
xiddot = 7;			% mm/yr: Rate of base level rise 1000

dt*Nprint*Ntoprint/365.25;

xiddot = xiddot / 1000 / timeyr; % m / sec
% dt = timeyr * dt/100;
dt = timeday * dt; % dt in seconds

dsbar = 1/N;
dstbar = 1/Ns;

Cfa = 1 / Cza ^ 2;
Cfs = 1 / Czs ^ 2; 

Ds = Ds / 1000;
Dm = Dm / 1000;

M = Ns+N+1;

%Defining parameters
etafluv = zeros(N, 1);
sfbar = zeros(N, 1);
rfluv = zeros(N, 1);

etaturb = zeros(Ns+1, 1);
stbar = zeros(Ns+1, 1);
rturb = zeros(Ns+1, 1);
% rturbratio = zeros(Ns+2,1);
Hfluv = zeros(N, 1);
% qso = zeros(N+1, 1);
qm = zeros(Ns+1, 1);
dqm=zeros(Ns+1,1);
Sturb = zeros(Ns+1, 1);
cibar=zeros(Ns+1,1);
cibardev=zeros(Ns+1,1);
deta=zeros(Ns+1,1);
qmdev=zeros(Ns+1,1);
qms=zeros(N,1);
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






%% Compute Normal Flow
% This subroutine computes the Shields stress, slope and depth of the normal flow
% of the fluvial region associated with the volume sand feed rate per unit wifth qso,
% the water discharge per unit width qw and the sand grain size Ds.

tsn = (qso / ((Rs * g * Ds) ^ 0.5 * Ds * al)) ^ (1 / nl) + tsc;
% tsn = 1.86;
Sn = Cza * (Rs * Ds * tsn) ^ 1.5 * g ^ 0.5 / qw;
Hn = Rs * Ds * tsn / Sn;

aleh = 0.1;
tausform = 1.86;
loadcoef = aleh * tausform *Cza / Rs  ;

%%
etabI = -1.5;              % m: initial elevation of bottom of the foreset  1.5

etatI=0; % m
% qmo=3e-06;
retentionrate= 0.0; % wax lake delta application; 0, 0.1, 0.2, 0.3, 0.4, 0.5 
qmo = ((qso * (1-Fsand)/Fsand) - Lamms*qso)*retentionrate ; % mud dirscharge to bottomset

xil =etatI+Hn*1;              % m: water surface elevation
Hturb = etatI*ones(Ns+1, 1);

% interM=0.35;
% 
% % figure saving
% poster_0='DimentsionLess_depth';
% poster_1='Dimension_depth';
% poster1='Sratio';
% figurename1=sprintf('%s %d %s %d',poster_0,etatI,poster1,sr); % dimensionless
% figurename2=sprintf('%s %d %s %d',poster_1,etatI,poster1,sr); % dimension 




washload=0;
capturemud=0;
fbt_delta_settling_fieldapp

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
% results.Lamms=capturemud/(qso*dt*Ntoprint*Nprint);
% results.mudratio=((qmo*dt*Ntoprint*Nprint)-capturemud) / (qso*dt*Ntoprint*Nprint);
results.foresetlength=foresetlength;
results.topsetlength=topsetlength;
results.topsetlength=topsetlength;
% 
 R.r(1,1) = results;

% save('New_RR_20_floc.mat','R') % change to - retention rate - and save

