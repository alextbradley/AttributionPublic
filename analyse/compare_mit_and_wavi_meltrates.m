% 09/06/22
%
% Produce maps of WAVI (left) and MITgcm (right) melt rates for a given forcing and a given melt parameter
%
%

%
% Preliminaries
%
clear
addpath('../plot_tools')
iceout_folder = '/data/icesheet_output/aleey/wavi/';
oceanout_folder = '/data/oceans_output/shelf/aleey/mitgcm/';
density_ice   = 918.0;
dx = 1000;
dy = 1000;

% 
% ensemble values (new plot for each ensemble member and gamma_T value)
%
gamma_T    = [12];   %gamma_T values, all *1e-3
ensembles  = [1];         %1: anthropogenic, 2: natural. 
members    = [2,3,4,5];         %ensemble member number. New plot for each different forcing
timeslices = [0,25,50,75,100]; lt = length(timeslices);

%
% Get data
%
ss_wavi = get_wavi_results(gamma_T, ensembles, members, timeslices,iceout_folder,dx,dy);
ss_mit = get_mitgcm_results(gamma_T, ensembles, members, timeslices,oceanout_folder);


%
% make plots
%
numplot = 1;
for ig = 1:length(gamma_T)
for ie = 1:length(ensembles)
for im = 1:length(members)
figure(numplot); clf;
count = 1;
for it = 1:length(timeslices)
	%get mit and wavi melt rate
        m_wavi = ss_wavi(ig,ie,im,it).m;
        m_mit  = ss_mit(ig,ie,im,it).m;
	idx    = (m_wavi ~= 0) & (m_mit ~= 0);
	m_wavi_mean = mean(m_wavi(idx));
	m_mit_mean  = mean(m_mit(idx));

	m_mit(m_mit == 0) = nan;
	m_wavi(m_wavi == 0) = nan;
	
	%for working out x lims
	sz = size(m_mit);
	[xidx,~] = find(~isnan(m_mit)); %indices of non-nan entries in mitgcm
	xidx = min(xidx);
	%saturate values
	sv = 50;
	m_wavi_sat = saturate(m_wavi,sv,0);
	m_mit_sat  = saturate(m_mit, sv,0);

	%plot them
	subplot(lt,2,count);
	contourf(m_wavi_sat', 20, 'linestyle','none')
	colorbar;
	title(strcat("WAVI melt rate at $t = ", num2str(timeslices(it)), "$, $\bar{m} = ", num2str(m_wavi_mean), "$ m/yr" ), 'interpreter', 'latex' ); 
	xlim([xidx - 10,sz(1)])
	count = count + 1;

	subplot(lt,2,count);
	contourf(m_mit_sat', 20, 'linestyle','none');
	colorbar;
	title(strcat("mitgcm melt rate at $t = ", num2str(timeslices(it)), "$, $\bar{m} = ", num2str(m_mit_mean), "$ m/yr" ), 'interpreter', 'latex' ); 
	xlim([xidx - 10,sz(1)])
	count = count + 1;
	
	%set x scale
end %end loop over timeslices
numplot = numplot + 1;
end %end loop over members
end %end loop over ensembles
end %end loop over gamma_T
