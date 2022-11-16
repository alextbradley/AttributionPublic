% Generate the data file for figure 5, which shows the probability enhancement 
%
% 
% Preliminaries
%
addpath('functions')


% 
% Parameters
%
x       = -1000:10:4000; %range of plausible SLR values
gammas  = 8:12; 	%gamma_T values (x 10^3)
mu      = 10; 		%prior hyperparameter
sigma_m = 0.1;          %error covariance in melting 
sigma_a = 1;            %error covariance in alpha

 
% 
% Run info
%
ensembles  = [3;2];       %ensemble 3: 100m trend per century, ensemble 2: 50m trend per century
members    = [1:20;1:20]; %ensemble member numbers
timeslices = 0:25:100;    %calibration timeslices (must be subset of [0,25,50,75,100]) (only relevant for uniformD = false)
SLR_time   = 1:1:100;     %times at which to output the pdfs
uniformD   = true;        %set to true to assume all gamma§ equally likely (no mitgcm info)

%
% Generate the data
%
total_runs = length(ensembles)*length(members)*length(SLR_time); count = 1;
cdf = zeros(length(x), length(ensembles), length(members), length(SLR_time));
for ie = 1:length(ensembles)  %loop over ensembles
ensemble = ensembles(ie);
for im = 1:length(members)    %loop over ensemble members
for it = 1:length(SLR_time)   %loop over time output points
        [cc,pp]  = get_cdf(x,gammas,ensemble,members(ie,im),timeslices,mu,sigma_m,sigma_a, SLR_time(it),uniformD);
        cdf(:,ie,im,it) = cc;

        fprintf('completed %.3d of %.3d \n', count, total_runs); count = count + 1;

end %end loop over SLR_time
end
end

%
% Save the data
% 
save('data-for-figures/figure5-data.mat', 'cdf', 'ensembles', 'members', 'timeslices', 'SLR_time', 'x', 'uniformD', 'sigma_a', 'sigma_m', 'mu')

