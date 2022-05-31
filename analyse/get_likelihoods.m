% 12/05/2022
%
% Generate an array of likelihoods.
%
% Note: it is hard coded into this script that the ice front is at x = 300
% with the ocean restoring boundary at x = 340km, and that the wavi domain
% has resolution 1km (i.e. the grid is 300 x 50)

%
% Preliminaries
%
clear
iceout_folder = '/data/icesheet_output/aleey/wavi/';
oceanout_folder = '/data/oceans_output/shelf/aleey/mitgcm/';

plot_stuff       = 0; %plot melt rate in each case as we run through
plot_heatmaps    = 0; %plot heatmap of correlations in a specific case
plot_single_entry= 0; %contour plot melt rate for a single case
SLR_vs_alpha     = 0; %plot SLR as a function of alpha for each enseble member
pdf_plots        = 1; %plot the pdfs of slr for each ensemble members
addpath('../plot_tools');

%
% Ensemble Values
%
gamma_T    = [8, 9, 10, 11];   %gamma_T values, all *1e-3
%gamma_T    = [8];
ensembles  = [1,2];                %one for anthropogenic, two for natural
%ensembles  = [1];
members    = [1, 2, 3, 4, 5];      %ensemble member numbers
timeslices = [0, 25, 50, 75,100]; %timeslice

SLR_time   = 100; %which time to take sea level rise by



%
% Parameters
%
secs_per_year = 365*24*60^2;
density_ice   = 918.0;
dx = 1000; 
dy = 1000; %grid sizes

correlations = zeros(length(gamma_T), length(ensembles), length(members), length(timeslices));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAVI runs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss_wavi = struct;
wavi_mean_melt = zeros(length(gamma_T), length(ensembles), length(members), length(timeslices));

numplot = 1;

for ig = 1:length(gamma_T)
    for ie = 1:length(ensembles)
        
        if plot_stuff
        figure(numplot); clf; suptitle(['WAVI, \gamma_T = ' , num2str(gamma_T(ig), '%02.f'), '*1e-3, ensemble number ', num2str(ensembles(ie))]);
        end
        
        count = 1; %initize run counter 
        for im = 1:length(members)
            %generate run number
            run_no = strcat(num2str(gamma_T(ig), '%02.f'), num2str(ensembles(ie)), num2str(members(im), '%02.f'));
            
            %if run exists, load it and loop over timeslices
            fname = strcat(iceout_folder, 'ATTR_', run_no, '/run/outfile.nc');
       
            if exist(fname, 'file') %if file exists, load info
                h = ncread(fname, 'h');
                t = ncread(fname, 'TIME');
                gr = ncread(fname ,'grfrac');
                m = ncread(fname, 'm');
                x = ncread(fname, 'x');
                y = ncread(fname, 'y');
                b = ncread(fname, 'b');
                
                %get the SLR associated 
                haf = h - (1028.0/918.0)*(-b);
                haf0 = haf(:,:,1); %haf at first time
                vaf0 = sum(sum(haf0(haf0 > 0))) * dx *dy;
                [~,idxT] = min(abs(t - SLR_time));
                hafT = haf(:,:,idxT);
                vafT = sum(sum(hafT(hafT > 0)))*dx *dy;
                dmaf = (vaf0 - vafT)*density_ice;   %change in mass above floatation
                slr = dmaf/1e9/361.8;               %associated sea level rise
                
             
                for it = 1:length(timeslices)
                    ss_wavi(ig, ie, im, it).run_no = run_no;
                    ss_wavi(ig, ie, im, it).slr = slr;   %slr same for all timeslices
                    
                    [~, idx] = min(abs(timeslices(it) - t)); %get index of wavi outfile closest to specified timeslice

                    ss_wavi(ig, ie, im, it).h = squeeze(h(:,:,idx)); %ice thickness at this timeslice 
                    ss_wavi(ig, ie, im, it).m = squeeze(m(:,:,idx)); %melt at this timeslice 
                    ss_wavi(ig, ie, im, it).gr = squeeze(gr(:,:,idx)); %grounded fraction at this timeslice
                
                    mm =  squeeze(m(:,:,idx));
                    wavi_mean_melt(ig,ie,im,it) = mean(mm(mm~=0));
                    
                    
                    %make plot
                    if plot_stuff
                       subplot(length(members),length(timeslices), count); box on
                       mm = ss_wavi(ig, ie, im, it).m;
                       
                       mm = saturate(mm,50,0);
                       contourf(x/1e3,y/1e3,mm', 20, 'linestyle', 'none'); 
                       title(['member no ' , num2str(members(im), '%02.f') , ', timeslice ',  num2str(timeslices(it), '%03.f')])

                    end
                    count = count + 1;
                end %end loop over timeslice  
    
            end %end if statement
        end %end loop over ensemble members
        
        if plot_stuff
            numplot = numplot + 1;
        end
        
    end %end loop over ensembles
end %end loop over gamma_T


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MITgcm runs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss_mit = struct;
mit_mean_melt = zeros(length(gamma_T), length(ensembles), length(members), length(timeslices));

for ig = 1:length(gamma_T)
    for ie = 1:length(ensembles)
        
        if plot_stuff
        figure(numplot); clf; suptitle(['MITgcm, \gamma_T = ' , num2str(gamma_T(ig), '%02.f'), '*1e-3, ensemble number ', num2str(ensembles(ie))]);
        end
        
        count = 1; %initize run counter 
        for im = 1:length(members) %loop over ensemble members
            
            for it = 1:length(timeslices)
                %filename 
                run_no = strcat(num2str(gamma_T(ig), '%02.f'), num2str(ensembles(ie)), num2str(members(im), '%02.f'), num2str(timeslices(it), '%03.f'));
                fname = strcat(oceanout_folder, 'ATTRmit_', run_no, '/run/state2D.nc');
                
                %if simulation exists, get the melt rate at the final time
                if exist(fname, 'file')
                    melt = ncread(fname, 'SHIfwFlx', [1, 1, 1], [Inf, Inf, Inf]);
                    melt = -melt * secs_per_year / density_ice;
                    melt = squeeze(melt(:,:,end));
                    szm  = size(melt);
                    %x    = ncread(state2D_fname, 'LONGITUDE');
                    %y    = ncread(state2D_fname, 'LATITUDE');
                    
                    %put the melt on the same grid as wavi. Final 40 grid
                    %points are ocean only (restoring boundary at 340km,
                    %ice front at 300km)
                    m = zeros(300,50);
                    nf = szm - 40; %40 ocean ocean only grid points 
                    m(end - nf + 1:end,:) = melt(1:nf, :);                   
                    ss_mit(ig,ie, im, it).m = m;
                    
                    mit_mean_melt(ig,ie,im,it) = mean(m(m~=0));
                    
                    %compute the correlation coefficient
                    m_mit = ss_mit(ig,ie, im,it).m;
                    m_wavi = ss_wavi(ig, ie, im, it).m;
                    idx = (m_mit ~= 0) & (m_wavi ~= 0); %where we have non zero entries in both (wavi has non-zero gl entries which mitgcm doesn't)
                    mv_mit = m_mit(idx); mv_wavi = m_wavi(idx); %melt entries in the non-zero domain for both
                    c     = corr([mv_mit, mv_wavi]);     %correlation of these entries
                    correlations(ig, ie, im, it) = c(1,2); %take the off diagonal entry
                    %pause
                    
                end
                
                %make plot
                if plot_stuff
                    subplot(length(members),length(timeslices), count); box on
                    mm = ss_mit(ig, ie, im, it).m;
                    mm = saturate(mm,50,0);
                    contourf(x/1e3,y/1e3,mm', 20, 'linestyle', 'none');
                    title(['member no ' , num2str(members(im), '%02.f') , ', timeslice ',  num2str(timeslices(it), '%03.f')])
                    
                end
                count = count + 1;
                
                
            end %end loop over timelice

        end %end loop over ensemble members
        
        if plot_stuff
        numplot = numplot + 1;
        end
        
    end %end loop over ensembles
end %end loop over gamma_T

rsquared = correlations.^2;

%% Plots
%
%Make a heatmap comparing mean melt rate
%
if plot_heatmaps
ig = 1; 
ie = 1; 

figure(numplot); clf; 

subplot(1,3,1); heatmap(timeslices, members, squeeze(wavi_mean_melt(1,1,:,:)));
xlabel('timeslices'); ylabel('ensemble members')
title('WAVI mean melt rate');

subplot(1,3,2); heatmap(timeslices, members, squeeze(mit_mean_melt(1,1,:,:)));
xlabel('timeslices'); ylabel('ensemble members')
title('MITgcm mean melt rate');

subplot(1,3,3); heatmap(timeslices, members, squeeze(rsquared(1,1,:,:)));
xlabel('timeslices'); ylabel('ensemble members')
title('melt rate correlations');

suptitle(['mean melt rate for \gamma_T = ' , num2str(gamma_T(ig), '%02.f'), '*1e-3 ensemble no ' , num2str(ensembles(ie)) ])
numplot = numplot + 1;
end
%%
%
% Plot a single melt rate and look where there are entries in one and not the other
% 
if plot_single_entry
%should we check each entry?
ig = 4; 
ie = 1; 
im = 3;
it = 5;

figure(numplot); clf; 
subplot(3,1,1); 
mm_wavi = ss_wavi(ig,ie, im, it).m;
mm_wavi = saturate(mm_wavi, 50, 0);
contourf(x,y,mm_wavi', 20, 'linestyle', 'none');
colorbar;
%spy(mm_wavi);
title('WAVI melt rate');

subplot(3,1,2); 
mm_mit = ss_mit(ig,ie, im, it).m;
mm_mit = saturate(mm_mit, 50, 0);
contourf(x,y,mm_mit', 20, 'linestyle', 'none');
colorbar;
%spy(mm_mit);
title('MITgcm melt rate');

subplot(3,1,3);
swav = mm_wavi ~= 0;
smit = mm_mit ~= 0;
surf(x,y, (smit - swav)'); %-1 where melt entries in wavi and not mitgcm, +1 where melt entries in mitgcm and not wavi
view([0 90])
colorbar
title('-1 where melt entries in wavi and not mitgcm, +1 where melt entries in mitgcm and not wavi')

tt = suptitle(['MITgcm, \gamma_T = ' , num2str(gamma_T(ig), '%02.f'), '*1e-3, ensemble number ', num2str(ensembles(ie)), ...
    ', member no ' , num2str(members(im)) , ', timeslice ' ,  num2str(timeslices(it)) , '. rsquared: ', num2str(rsquared(ig, ie, im, it))]);
tt.Position(2) = tt.Position(2) + 0.02;
numplot = numplot + 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Probability Stuffs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmam = 1; %error covariance in MITgcm
sigmaa = 1e-3;
mu     = 10e-3;
a      = gamma_T*1e-3;

P_mmit_given_alpha = exp(-(mit_mean_melt - wavi_mean_melt).^2 / 2 /sigmam^2);
P_alpha_given_mu = zeros(length(gamma_T), length(ensembles), length(members), length(timeslices)); 
%fill in P(alpha|mu), which doesn't know about naything expcet for gamma_T
%( = alpha)
for ie = 1:length(ensembles)
    for im = 1:length(members)
        for it = 1:length(timeslices)       
            P_alpha_given_mu(:, ie, im,it)   = exp(-(a - mu).^2 / 2 / sigmaa^2);
        end
    end
end

P_alpha_given_mit_mu = P_mmit_given_alpha .* P_alpha_given_mu; %no normalization?
P_alpha_given_mit_mu = mean(P_alpha_given_mit_mu, 4); %average over timeslices
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot SLR as a funnction of alpha. Different ensemble on each row of subplot
if SLR_vs_alpha
    numplot = 1;
figure(numplot); clf; 
count = 1;
slr_out = linspace(0,3000);
alpha_test = linspace(min(gamma_T)/1e3, max(gamma_T)/1e3);

melt_diff = mit_mean_melt - wavi_mean_melt;
mean_melt_diff = mean(melt_diff, 4); %timeslice average mean melt


iT = 1; %SLR calculations always use same timeslice
for ie = 1:length(ensembles)
    for im = 1:length(members)
        subplot(length(ensembles), length(members), count)
        slrdata = [ss_wavi(:,ie, im,it).slr];
        
        plot(gamma_T/1e3, slrdata, 'ro')
        xx = linspace(min(gamma_T), max(gamma_T))*1e-3;
        yy = interp1(gamma_T/1e3, slrdata, xx);
        hold on; 
        plot(xx, yy, 'k--');
        
        % add an example showing that the intersections code works
        X1 = gamma_T/1e3; Y1 = slrdata;
        X2 = gamma_T/1e3; Y2 = 500*ones(size(gamma_T));
        [x0,y0,iout,jout] = intersections(X1,Y1,X2,Y2);
        if ~isempty(x0)
            plot(x0, y0, 'ko', 'markerfacecolor', 'k')
        end
        
        
        count = count +1;
    end
end
figure(numplot); suptitle('slr as a function of alpha')
numplot = numplot + 1;

%%
figure(numplot);  clf;
count = 1;
for ie = 1:length(ensembles)
    for im = 1:length(members)
  
        subplot(length(ensembles), length(members), count); hold on
        plot(gamma_T/1e3, mean_melt_diff(:, ie, im), 'ro');
        mm = mean_melt_diff(:, ie, im);
        mmi = interp1(gamma_T/1e3, mm, xx);
        plot(xx, mmi, 'k--');
        
        
        count = count +1;
    end
end

 suptitle('mean melt diff')
numplot = numplot +1;
end

%%%%%%%%%%%%%%%%%%%%%%% generate pdfs %%%%%%%%%%%%%%%%%%%%


figure(3); clf; 
count = 1;
slr_out = linspace(0,4000,1e4);%needs to be regularly spaced
ds = diff(slr_out); ds = ds(1);
alpha_test = linspace(min(gamma_T)/1e3, max(gamma_T)/1e3);

melt_diff = mit_mean_melt - wavi_mean_melt;
mean_melt_diff = mean(melt_diff, 4); %timeslice average mean melt difference

sigma_a = 1e-3;
sigma_m = 2; 
mu = 1e-2;

iT = 1; %SLR calculations always use same timeslice

pslr_all = zeros(length(ensembles), length(members), length(slr_out));
for ie = 1:length(ensembles)
    for im = 1:length(members)
	im
        slrdata = [ss_wavi(:,ie, im,it).slr]; 
        mm = mean_melt_diff(:, ie, im);
        
        %generate the pdfs
        pslr = zeros(1,length(slr_out));
        P_mmit_given_alpha = zeros(1,length(slr_out));
        P_alpha_given_mu = zeros(1,length(slr_out));
        for is = 1:length(slr_out)
            %for every value of slr, find the associated value(s) of alpha
            Y2 = slr_out(is)*ones(size(gamma_T));
            [alpha,slr,iout,jout] = intersections(gamma_T/1e3,slrdata,gamma_T/1e3,Y2);
            
	    for ia = 1:length(alpha)
            %for every value of alpha, compute P(alpha|mu) = exp((-alpha - mu).^2 / 2 / sigma_m^2 )and
            %P(m_mit|alpha) = exp((m_mit(alpha) - m_wavi(alpha))
                P_alpha_given_mu = exp(-(alpha(ia) - mu).^2 / 2 / sigma_a^2 );
		
        	mean_melt_diff_at_this_alpha = interp1(gamma_T/1e3, mm, alpha(ia)); %get the difference in the mean melt rates at this value of alpha
		P_mmit_given_alpha = exp(-(mean_melt_diff_at_this_alpha)^2 /2 / sigma_m^2);
		p_this_alpha = P_alpha_given_mu * P_mmit_given_alpha;
		pslr(is) = pslr(is) +  p_this_alpha; %add the probability associated with this alpha
                
		%pslr(is) =  p_this_alpha; %add the probability associated with this alpha	
            end
        end
	%normalize the distribution        
	pslr = pslr / sum(pslr*ds);	
 
        subplot(length(ensembles), length(members), count)
	plot(slr_out,pslr, 'k');
	xlabel('slr (mm)'); ylabel('prob') 
        count = count +1;
	pslr_all(ie,im,:)=pslr;
    end
end

numplot = numplot + 1;
%plot anthrop vs natural
figure(4); clf; hold on
colmap = lines(5);
for ie = 1:length(ensembles)
	pslr = squeeze(pslr_all(ie,:,:));
	pslr = mean(pslr,1); %take the mean over the runs
	meandist = sum(slr_out .* pslr); %mean of the distibution
	vardist = sum((slr_out-meandist).^2 .* pslr)/2; %and variance
	normal_distn = 1 / sqrt(2*pi*vardist) * exp(-(slr_out - meandist).^2 / 2 /vardist);

%	plot(slr_out, pslr, 'linewidth', 2, 'color', colmap(ie,:));
	plot(slr_out, normal_distn,'linewidth', 2, 'color', colmap(ie,:));
end
title('Distributions modelled as Gaussians')
legend('anthropogenic', 'natural'); 


