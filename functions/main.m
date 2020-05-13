% clear all
% is =4 
addpath('./functions')

%% DEFINITIONS
include100=1; % 0 to exclude events between 100-250 km
savefig=0; %1 to save figure, 0 not to.

% station params
% PO stations: CHGQ, MATQ, WEMQ, YOSQ
stas=['SOLO';'TBO ';'GTOO';'SILO';'TIMO';'PNPO';'SUNO';'MATQ';'PKLO';...
    'KILO';'MALO';'WEMQ';'VLDQ';'VIMO';'EPLO';'CHGQ';'ULM ';'KAPO';'MUMO';'KASO'];
nets=[repmat('CN',7,1);'PO';'CN';'CN';'CN';'PO';repmat('CN',3,1);'PO';'CN';'CN';'CN';'CN'];

% is=1;
sta=strtrim(stas(is,:));
net=nets(is,:);
% sta='YOSQ'
% net='PO'
fdir=['/data2/home/genevieve/research/metalearth/MGBmethod/data/',sta]
% channelcode='H'
ndt0=0.01; %Station sampling currently all stations selected by Dave are sampled at 0.01s (100 Hz)

% Lengths and decimation parameters
rlength=50;  % desired length of receiver function.
%  ndecimate=10;
ndecimate=2;
hlength=200;
ndt=ndt0*ndecimate;
dtap=3.0; %taper? % default: 3.0
% Parameter rlength is desired length of receiver function.
% This is required in order to determine lengths of time series
% to be used in the deconvolution in order to avoid wrap around
% effects.
ndec=ndecimate;
npts=fix(hlength/(ndt0));
irpts=fix(rlength/(ndt0));
    
% filtering params
% Moho: .02 - .8
% 410-670: .02 .2
lf=0.02;
hf=0.8;
nyq=0.5/(ndecimate*ndt0);
wn=[lf/nyq,hf/nyq];
[b,a]=butter(2,wn);

% other
tshift=0.0; % used for plotting time to depth in makeplot

% Near surface crustal parameters for free surface correction.
% Vp within 5.9-6.1, Vs within 3.7-3.8
a0=6;
b0=3.75;
b02=b0*b0; % Convenient definitions.
a02=a0*a0;


%% . . Read in data.
readdata;

%% . . Define the p/BAZ bins
disp('Making the p/baz grid')
makegrid;
% makegrid_kmeans;

%% . . Normalize traces, assemble in bins.
disp('Normalizing')
normalize;

%% . . Deconvolution.
disp('Deconvolution & plot')
% makeplot2
lp=0.75;
depthmap
% lp=2.00;
% depthmap;

