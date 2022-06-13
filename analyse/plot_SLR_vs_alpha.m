% 31/05/22
% Produce a plot of SLR as a function of alpha for different realizations of forcing. 

%
% Preliminaries
%
clear
iceout_folder = '/data/icesheet_output/aleey/wavi/';
density_ice   = 918.0;
dx = 1000;
dy = 1000;


%
% Ensemble
%
big_ensemble = 0; 	       %set to one for whole ensemble
gamma_T    = [8, 9, 10, 11,12];   %gamma_T values, all *1e-3
ensembles  = [1,2]; 	       %1: anthropogenic, 2: natural. New plot for each
SLR_time   = 0; 		%time at which to compute SLR
timeslices = [SLR_time];

%ensemble member numbers
if big_ensemble
	members = 1:20;
else
	members = 1:5;
end

%
% genndata
%
SLR = zeros(length(ensembles),length(members),length(gamma_T));
numplot = 1;

%get the data
ss = get_wavi_results(gamma_T, ensembles, members, timeslices,iceout_folder,dx,dy);

for ie = 1:length(ensembles)
figure(numplot); clf; suptitle(strcat("ensemble number ", num2str(ensembles(ie))));
count = 1;
for im = 1:length(members)
%plot gamma_T vs SLR for this ensemble member
subplot(ceil(length(members)/5),5,count); hold on; box on
plot(gamma_T,[ss(:,ie,im).slr], 'ro--', 'markersize', 5)
xlabel('$\gamma_T$', 'interpreter', 'latex');
ylabel('SLR (mm)', 'interpreter', 'latex');
count = count + 1;
end %end loop over members
numplot = numplot + 1;
end %end loop over ensembles
