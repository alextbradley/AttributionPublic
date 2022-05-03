% 03/05/22
%
% Generate mitgcm driver/data files for a given simulation at a given timeslice
%
% This code produces a folder within the current directory which includes:
% - bathymetry
% - lev_s (initial salinity guess, set to boundary restoring profile everywhere)
% - lev_t (initial temperature guess, set to boundary restoring profile everywhere)
% - pload (pressure loading file)
% - shelfice_topo (ice topography file)
% - RBCt, RBCs (restoring files)
% - RBCt_mask, RBCs_mask (restoring masks)

% Preliminaries 
clear;
addpath('..')
iceout_folder = '/data/icesheet_output/shelf/aleey/wavi/';

%
% Input info
% 
timeslices = [25]; %array of timeslices to generate
run_nos    = ["08101"]; %array of simulations to generate timeslices for

%
% Grid info
%
deltaZ = 10; %vertical grid spacing
nz     = 110; %number of grid cells in z
dz = deltaZ*ones(1,nz);
zz = [0,cumsum(dz)];

%
% Parameters
%
tl = 1.2; %lower layer temperature
tu = -1;  %upper layer temperature
sl = 34.6; %lower layer salinity 
su = 34.0; %upper layer salinity


%
% Find the correct forcing from the run name
% force_no as follows:
% 0 <-> constant forcing
% 1 <-> anthropogenic forcing
% 2 <-> natural forcing
%

iR = 1; %will loop over iR (run_nos)
split_run_no = split(run_nos(iR), ''); 
force_no = str2num(split_run_no(4));
exp_no   = str2num(split_run_no(6));

%
% get the anomaly
%
if force_no == 1 %anthropogenic forcing
	anomaly = load(strcat("../ocean_forcing/anthropogenic_random_forcing/random_forcing_" , num2str(exp_no), "/anomaly.mat"));
elseif force_no == 2 %natural forcing
	anomaly = load(strcat("../ocean_forcing/natural_random_forcing/random_forcing_" , num2str(exp_no), "/anomaly.mat"));
elseif force_no == 1 && exp_no == 1 %cold forcing, pycnocline centre at -600m
	anomaly.t = 0;
	anomaly.pc = -600; 
elseif force_no == 1 && exp_no == 2 %intermediate forcing, pycnocline centre at -500m
	anomaly.t = 0; 
	anomaly.pc = -500;
elseif force_no == 1 && exp_no == 3 %warm forcing, pycnocline centre at -400m
	anomaly.t = 0;
	anomaly.pc = -400; 
else
	error("Check file number. No forcing of this type found")
end

 
% 
% Loop over timeslices. For each timeslice, find the nearest pycnocline depth from the anomaly file and use this to generate lev_s, lev_t, RBC_s, RBC_t. Use the nearest solution file to return the geometry.
%
for iT = 1:length(timeslices)
	[~,idx] = min(abs(timeslices(iT) - anomaly.t)); %anomalies have quite high resolution (0.1 yrs) so OK to just take the nearest entry
	pc = anomaly.pc(idx); %pycnocline centre at this timepoint
        ta = ambient_profile(pc, pw, -zz, tl, tu);
        sa = ambient_profile(pc, pw, -zz, sl, su);





end %end loop over timeslice

