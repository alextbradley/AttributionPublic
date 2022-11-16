% Make figure 5 of the ms, showing the (a) raw and (b) percentage enhancement in likelihood as a function of time and slr.

%
% Load in the data
% 
f       = load('data-for-figures/figure5-data.mat');
cdf     = squeeze(mean(f.cdf, 3)); %ensemble mean cdf
genfits = 0; %set to 1 to pass through the fitting loop 

%
% Other stuff
%
xl = [-100, 3000]; %xlimits in the plot
yl = [10,100];      %y limits
fs = 11; 	   %font size

%
% Perform fits
%
ft = fittype('1/2 + 1/2*erf(b*sign(x-a)*sqrt(abs(x-a)) + c)', 'independent', 'x', 'dependent', 'y'); %fit type to be applied

%initialize space
if genfits
cdf_mean_fit = nan(size(cdf));
pdf_mean_fit = nan(length(f.x)-1, length(f.ensembles), length(f.SLR_time));

%count number of sims
count = 1;
total_fits = length(f.ensembles)*length(f.SLR_time);

for ie = 1:length(f.ensembles)
for it = 2:length(f.SLR_time);
mean_profile = squeeze(cdf(:,ie,it)); %ensemble mean cdf at this time
mean_profile = reshape(mean_profile, [length(mean_profile),1]); % make a col vector
mean_profile_fit = fit(f.x', mean_profile, ft);
cdf_mean_fit(:,ie,it) = mean_profile_fit(f.x);

%differentiate to get the pdf
mm = mean_profile_fit(f.x);
dx = diff(f.x); dx = dx(1);
pdf_mean_fit(:,ie,it) = diff(mm)/ dx;

fprintf('completed %.3d of %.3d fits \n', count, total_fits); count = count + 1;

end %end loop over slr time
end %end loop over ensembles
end %end genfits flag

%
% Make the plots
%
fig = figure(1); clf; fig.Position(3:4) = [1200,400];
addpath('plottools');

% (a) density enhancement
density_enhancement = squeeze(pdf_mean_fit(:,1,:) - pdf_mean_fit(:,2,:));
xm = (f.x(2:end) + f.x(1:end-1))/2;; 
ax(1) = subplot(1,2,1); %contourf(xm,f.SLR_time, density_enhancement', 50, 'linecolor', 'none'); 
imagesc(xm,f.SLR_time, density_enhancement')
caxis(5e-4*[-1,1])
colormap(ax(1), cmocean('balance'));
c(1) = colorbar;
c(1).Label.String = 'anthropogenic density enhancement';
c(1).Label.FontSize = fs;
xlim(xl);
ylim(yl);
xlabel('sea level rise (mm)');
ylabel('time (years)');
ax(1).FontSize = fs;
set(gca, 'YDir', 'normal')

% (b) percentage enhancement
percentage_enhancement = squeeze((pdf_mean_fit(:,1,:) - pdf_mean_fit(:,2,:))./pdf_mean_fit(:,2,:))*100;
percentage_enhancement(abs(density_enhancement) < 1e-7) = nan;
ax(2) = subplot(1,2,2);
m = log10(percentage_enhancement').*sign(percentage_enhancement'); %what to plot

pos_enhancement = nan(size(percentage_enhancement));
pos_enhancement(percentage_enhancement > 0) =percentage_enhancement(percentage_enhancement > 0);

neg_enhancement = nan(size(percentage_enhancement));
neg_enhancement(percentage_enhancement < 0) =percentage_enhancement(percentage_enhancement <0);
h = imagesc(xm,f.SLR_time, log10(pos_enhancement)');

set(h, 'AlphaData', ~isnan(pos_enhancement')); %make the nan entries transparent
caxis([-1,4])
xlim(xl)
ylim(yl);
colormap(ax(2), cmocean('tempo'));
%c(2) = colorbar;
set(gca, 'YDir', 'normal')
xlabel('sea level rise (mm)');
ylabel('time (years)');
ax(2).FontSize = fs;

axnew = axes;
axnew.Position = ax(2).Position;
h = imagesc(xm,f.SLR_time, log10(abs(neg_enhancement))');
set(h, 'AlphaData', ~isnan(neg_enhancement')); %make the nan entries transparent
caxis([-1,4])
colormap(axnew, (cmocean('amp')));
%c(3) = colorbar;
set(gca, 'YDir', 'normal')
xlim(xl)
axnew.Visible = 'off';
%drawnow; pause; aa

% make a new axis for the colorbar
axn = axes();

nanmat = 1e3*percentage_enhancement';
cmap=  [flipud(cmocean('amp')); cmocean('tempo')];
axn.Position = ax(2).Position;
hh = imagesc(xm,f.SLR_time, nanmat);
caxis([-1,1])
xlim(xl)
colormap(axn, cmap);
set(gca, 'YDir', 'normal')
cc = colorbar;
set(hh, 'AlphaData', zeros(size((nanmat))));
axn.Visible = 'off';
cc.Position(1) = 0.92;
cc.Ticks = linspace(-1,1,13);
cc.TickLabels = {"-10^{4}", "-10^{3}", "-10^{2}", "-10^{1}", "-10^0", "-10^{-1}", "0", "10^{-1}", "10^0", "10^1", "10^2", "10^3", "10^4"};
cc.Label.String = 'anthropogenic enhancement ratio';
cc.Label.FontSize = fs;
