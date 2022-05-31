%Basic case analysis script

%% Preliminaries
clear
addpath('../plot_tools');
%% Case info
%specify case NB: check grid resolution
folder = "/data/icesheet_output/aleey/wavi/ATTR_10002/run/";
fname = strcat(folder,'outfile.nc');
dx = 1000; dy = 1000;

%% plot flags
grv_plot   = 1; %grounded fraction, mean melt and total melt flux as a function of time
slice_plot = 1; %slice taken along the centreline
gl_pos_plot= 1; %plot bathymetry and contour of grounding line
melt_vid   = 0; %video of melting profile
vel_plot   = 1; %contourf of velocity at final time
xgl_plot   = 1; %line plot of centreline gl position as a function of time

numplot = 1;
%% get solution
time = ncread(fname, 'TIME');
h_stack = ncread(fname, 'h');
x = ncread(fname, 'x');
y = ncread(fname, 'y');
b = ncread(fname, 'b');
u = ncread(fname, 'u');
v = ncread(fname, 'v');

%grfrac_stack = ncread(fname, 'grfrac');
melt = ncread(fname, 'm');

nx = length(x);
ny = length(y);
bathy = squeeze(b(:,:,1));

grounded_volume = zeros(1,length(time));
mean_melt = zeros(1,length(time));
total_melt_flux = zeros(1,length(time));


for i = 1:length(time)
    h = squeeze(h_stack(:,:,i));
    grfrac_stack(:,:,i) = (h > -1028/918 * bathy);
    grfrac = squeeze(grfrac_stack(:,:,i));
    grounded_volume(i)  = sum(sum(h .* grfrac .* dx * dy));
    
    idx = grfrac == 0;
    mm = squeeze(melt(:,:,i));
    mean_melt(i) = mean(mm(idx));
    total_melt_flux(i) = sum(sum(mm))*dx*dy;
end

%% Plots
% Plot 1: grounded fraction, mean melt and total melt flux as a function of
% time
if grv_plot
    figure(numplot); clf;
    subplot 131; hold on 
    plot(time, grounded_volume/1e9, 'o-')
    xlabel('time (years)');
    ylabel('grounded volume ($\mathrm{km}^3$)')
    
    subplot 132; hold on 
    plot(time, mean_melt, 'o-')
    xlabel('time (years)');
    ylabel('mean melt rate (m/yr)')
    
    subplot 133; hold on 
    plot(time, total_melt_flux/1e9, 'o-')
    xlabel('time (years)');
    ylabel('total melt flux (cubic km/yr)')

    numplot = numplot +1;
end
%
%% Plot 2: evolution of the gl
if gl_pos_plot
    figure(numplot); clf; hold on; box on
    %  add base layer of bathymetry
    contourf(x/1e3,y/1e3, bathy', 30, 'linestyle', 'none');
    ax = gca;
    colormap(ax, autumn(100));
    xlabel('x (km)');
    ylabel('y (km)');
    title('bathymetry and grounding line position')
    
    
    % add gronding lines
    axnew = axes;
    axnew.Position = ax.Position;
    hold on
    
    colmap = parula(length(time));
    for i = 1:length(time)
        contour(x/1e3,y/1e3,squeeze(grfrac_stack(:,:,i))', [0.5, 0.5], 'linecolor', colmap(i,:));
        % drawnow
        % pause
    end
    
    c = colorbar(axnew);
    c.TickLabels = {'0', num2str(time(end))};
    c.Ticks = [min(c.Ticks), max(c.Ticks)];
    c.Label.String = 'time (years)';
    axnew.Position = ax.Position;
    axnew.Visible = 'off';
    numplot = numplot +1;
end

%% Plot 3: slice along the centre
if slice_plot
    figure(numplot); clf; hold on; box on
    idx = floor(ny/2);
    %idx = 5
    b = bathy(:,idx);
    plot(x, b, 'm', 'linewidth', 2);
    colmap = parula(length(time));
    for i = 1:length(time)
        h = squeeze(h_stack(:,idx,i));
        grfrac = squeeze(grfrac_stack(:,idx,i));
        isfloat = (grfrac == 0);
        
        base = zeros(size(x));
        base(~isfloat) = b(~isfloat);
        base(isfloat) = -918/1028 * h(isfloat);
        s = h + base;
        plot(x,s, 'color', colmap(i,:));
        plot(x,base, '--', 'color', colmap(i,:), 'linewidth', 1.5);
        %pause
        %time(i)
    end
    c = colorbar;
    c.TickLabels = {'0', num2str(time(end)/2), num2str(time(end))};
    c.Ticks = [min(c.Ticks), (min(c.Ticks) + max(c.Ticks))/2,  max(c.Ticks)];
    c.Label.String = 'time (years)';
    xlim([min(x), max(x)]);
    numplot = numplot +1;
end
%% Make a video of melting
if melt_vid
    figure(numplot); clf; hold on
    colmap = parula(length(time));
    for i = 1:length(time)
        clf; hold on
        mm = squeeze(melt(:,:,i));
        
        idx = squeeze(grfrac_stack(:,:,i) == 1);
        mm(idx) = nan;
        %mm = saturate(mm, 20, 0);
        contourf(x/1e3,y/1e3, mm', 30, 'linestyle', 'none');
        ax = gca;
        colormap(ax, autumn(100));
        c = colorbar;
        xlabel('x (km)');
        ylabel('y (km)');
        title(['t = ' num2str(time(i))])
        contour(x/1e3,y/1e3,squeeze(grfrac_stack(:,:,i))', [0.5, 0.5], 'linecolor', 'k');
        %xlim([400, 640])
        xlim([200, 300])
        drawnow
        pause
        
    end
    numplot = numplot +1;
end

%% Plot the velocities at final timestep
if vel_plot
    figure(numplot); clf; hold on;
    idx = floor(length(time)/2);
    vv = squeeze(v(:,:,idx));
    uu = squeeze(u(:,:,idx));
    vel = sqrt(uu.^2 + vv.^2);
    contourf(x/1e3, y/1e3, vel', 20, 'linestyle', 'none');
    xlabel('x (km)');
    ylabel('y (km)');
    contour(x/1e3,y/1e3,squeeze(grfrac_stack(:,:,end))', [0.5, 0.5], 'linecolor', 'k');
    colorbar
    %xlim([250, 300])
    numplot = numplot + 1;
end

%% Plot the centreline grounding line position as a function of time
if xgl_plot
    figure(numplot); clf; hold on;
    xgl = zeros(1,length(time));
    for i = 1:length(time)
        centreline_grfrac = squeeze(grfrac_stack(:,floor(ny/2),i));
        idx = find(centreline_grfrac < 1, 1, 'first');
        xgl(i) = x(idx);
    end
    
    plot(time, xgl/1e3, 'ro-');
    xlabel('t (yrs)');
    ylabel('$x_{gl}$ (km)', 'interpreter', 'latex');
    ylim([180, 280])
    
    numplot = numplot + 1;

end
