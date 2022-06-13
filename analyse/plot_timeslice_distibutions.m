% 02/06/22
% Plot distributions for one particular realization of forcing.
%
% Produce a plot with (4 x N) subplots. N is the number of timeslices. Each column is a different timeslice 
% First row:  mean melt rate as a function of alpha for wavi and mitgcm
% Second row: P(mitgcm|alpha)
% Third row:  P(alpha|mu) [nb same for each timeslice)
% Fourth row: P(alpha|mitgcm) [normalized product of rows three and four]


%
% Preliminaries
%
clear
iceout_folder = '/data/icesheet_output/aleey/wavi/';
oceanout_folder = '/data/oceans_output/shelf/aleey/mitgcm/';
density_ice   = 918.0;
dx = 1000;
dy = 1000;
sigma_a = 1e-3*[1,2,3];   %sigma_a values
mu = 10e-3;	   %prior parameter
sigma_m = [1,2,3];

%
% Ensemble
%
gamma_T    = [8, 9, 10, 11,12];   %gamma_T values, all *1e-3
ensembles  = 1;         %1: anthropogenic, 2: natural. 
members    = 1; 	%ensemble member number. New plot for each different forcing
timeslices = [0,25,50,75,100];

%
% Get data
%
ss_wavi = get_wavi_results(gamma_T, ensembles, members, timeslices,iceout_folder,dx,dy);
ss_mit = get_mitgcm_results(gamma_T, ensembles, members, timeslices,oceanout_folder);

%
% make plots
%
numplot = 1;
for ie = 1:length(ensembles)
for im = 1:length(members)
figure(numplot); clf;

for it = 1:length(timeslices)
    for ig = 1:length(gamma_T)

	%get mit and wavi melt rate
	m_wavi = ss_wavi(ig,ie,im,it).m;
	m_mit  = ss_mit(ig,ie,im,it).m;

	%restrict to those values in shared domain
	idx = (m_mit ~= 0) & (m_wavi ~= 0); %where we have non zero entries in both (wavi has non-zero gl entries which mitgcm doesn't) 	
	
	wavi_mean_melt(ig) = mean(m_wavi(idx));
	mit_mean_melt(ig) = mean(m_mit(idx));
    end %end loop over gamma
    
    %
    % first row: mean melt rate as a function of gamma for mit and wavi
    %
    subplot(5,length(timeslices),it); hold on; box on
    plot(gamma_T, wavi_mean_melt, 'ko-', gamma_T, mit_mean_melt, 'ro-');
    xlabel('$\gamma_T \times 10^3$', 'interpreter', 'latex');
    ylabel('mean melt rate (m/yr)')
    title(strcat('t = ', num2str(timeslices(it))))
    if it == 1
	legend('WAVI', 'MITgcm')
    end 

    %
    % second row: difference between mean melt rates
    %
    
 
    %
    % third row: P(mitgcm|alpha) = exp(-D^2/(2sigma_m^2))
    %
    colmap = parula(length(sigma_m)+1);
    subplot(5,length(timeslices), it + length(timeslices)); hold on; box on;
    for is = 1:length(sigma_m)
	pp = exp(-(wavi_mean_melt - mit_mean_melt).^2/2/sigma_m(is)^2);

        plot(gamma_T, pp, 'o-', 'color', colmap(is,:));
	legendinfo{is} = ['$\sigma_m = ', num2str(sigma_m(is)), '$'];
    end %end loop over sigma_m values
    if it == 1
	legend(legendinfo, 'interpreter', 'latex')
    end
    xlabel('$\gamma_T \times 10^3$', 'interpreter', 'latex');
    ylabel('$P(\mathrm{mitgcm}|\gamma_T)$', 'interpreter', 'latex')

    %
    % fourth row: P(alpha|mu) = exp(-(alpha - mu)^2/(2sigma_a^2))
    %
    colmap_a = lines(length(sigma_a)+1);
    subplot(5,length(timeslices), it + 2*length(timeslices)); hold on; box on;
    for is = 1:length(sigma_a)
	pp = exp(-(gamma_T*1e-3 - mu).^2/2/sigma_a(is)^2);

        plot(gamma_T, pp, 'o-', 'color', colmap_a(is,:));
	legendinfo{is} = ['$\sigma_a = ', num2str(sigma_a(is)), '$'];
    end %end loop over sigma_m values
    if it == 1
	legend(legendinfo, 'interpreter', 'latex')
    end
    xlabel('$\gamma_T \times 10^3$', 'interpreter', 'latex');
    ylabel('$P(\alpha|\mu)$', 'interpreter', 'latex')

    %
    % fifth row: P(alpha|mit,mu)
    %
    subplot(5,length(timeslices), it + 3*length(timeslices)); hold on; box on;
    clear legendinfo
    count = 1;
    for isa = 1:length(sigma_a)
    for ism = 1:length(sigma_m)
	p_alpha_given_mu = exp(-(gamma_T*1e-3 - mu).^2/2/sigma_a(isa)^2);
	p_mit_given_alpha = exp(-(wavi_mean_melt - mit_mean_melt).^2/2/sigma_m(ism)^2);
	pp = p_alpha_given_mu .* p_mit_given_alpha;
        plot(gamma_T, pp, 'o-', 'color', colmap_a(isa,:), 'markeredgecolor', colmap(ism,:));
	legendinfo{count} = ['$\sigma_a = ', num2str(sigma_m(ism)), ', \sigma_a = ' num2str(sigma_a(isa)), '$'];
	count = count + 1;
    end %end loop over sigma_m values
    end %end loop oover sigma_a values
    if it == 1
	ll = 	legend(legendinfo, 'interpreter', 'latex');
    end
    xlabel('$\gamma_T \times 10^3$', 'interpreter', 'latex');
    ylabel('$P(\alpha|\mathrm{mitgcm})$', 'interpreter', 'latex')
    

end %end loop over timeslices 
numplot = numplot + 1;
end %end loop over ensembles
end %end 


