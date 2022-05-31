% 10/05/22
%
% Basic analysis of an mitgcm run.
%
% Plots:
% (1) plot mean melt rate and total flux as a function of time
% (2) contour plot of final melt rate pattern

%
% Plot flags
%
melt_evolution = 1; %plot mean melt rate and total flux as a function of time
contour_melt = 1; %contour plot of final melt rate pattern


caseno = "11101000";
main_dir = '/data/oceans_output/shelf/aleey/mitgcm/ATTRmit_'; %not in git repo
result_dir = strcat(main_dir, caseno, '/run/');
input_dir = strcat(main_dir, caseno, '/input/');


%
% Parameters
%
secs_per_year = 365.25*24*60*60;
secs_per_month = secs_per_year/12;
density_ice = 918.0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%melt rates and co-ordinates
state2D_fname = strcat(result_dir, 'state2D.nc');
melt = ncread(state2D_fname, 'SHIfwFlx', [1, 1, 1], [Inf, Inf, Inf]);
melt = -melt * secs_per_year / density_ice;
x    = ncread(state2D_fname, 'LONGITUDE');
y    = ncread(state2D_fname, 'LATITUDE');
nx   = length(x);
ny   = length(y);
dx   = diff(x); dx = dx(1);
dy   = diff(y); dy = dy(1);
time = ncread(state2D_fname, 'TIME');

%ice topo
fid = fopen(strcat(input_dir, 'shelfice_topo.bin'));
icetopo = fread(fid, 'real*8', 'b');
icetopo = reshape(icetopo, nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numplot = 1;
if melt_evolution
    total_flux = zeros(1,length(time));
    mean_melt  = zeros(1,length(time));
    for iT = 1:length(time)
       m = squeeze(melt(:,:,iT));
       total_flux(iT) = sum(sum(m))*dx*dy;
       mean_melt(iT)  = mean(mean(m(icetopo ~= 0)));
    end
    
    figure(numplot); clf; 
    subplot 121; box on
    plot(time/secs_per_month, total_flux/1e9, 'ro-', 'markerfacecolor', 'r');
    xlabel('time (mo)');
    ylabel('total flux (km^3)');
    
    subplot 122; box on
    plot(time/secs_per_month, mean_melt, 'ro-', 'markerfacecolor', 'r');
    xlabel('time (mo)');
    ylabel('mean shelf melt rate (m/yr)');
    
    numplot = numplot+1;
end

if contour_melt
    figure(numplot); clf; 
    m = squeeze(melt(:,:,end));
    m(icetopo == 0) = nan;
    contourf(x,y,m', 20, 'linestyle', 'none') 
    c = colorbar;
    c.Label.String = 'melt rate (m/yr)';
    axis equal
end
