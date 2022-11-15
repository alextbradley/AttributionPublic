% Make figure 3 of the attribution manuscript
%
% Stochastic retreat: grounded volume as a function of time for gamma_T = 10e-3 with natural and anthro and corresponding forcing profiles
%
%
% Preliminaries
%
addpath('plottools')
outfile_path = '/data/icesheet_output/aleey/wavi/'; %change to full path of result location
gendata = 1; %gendata loop flag

%
% Constants
%
rhoi = 918;  %ice density
rhow = 1028; %water density
dx   = 1000; 
dy   = 1000; %grid resolution

%
% Run info
%
gammas        = [10; 10];
ensembles     = [3; 2]; %anthro (with trend) first, then natural (no trend)
members       = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20"];
members       = repmat(members, [2,1]);
forcing_paths = ["../gendata/ocean/ocean_forcing/A100_T100";
		 "../gendata/ocean/ocean_forcing/A100_T000"];

sz      = size(members);

%
% Initialize storage
%
if gendata
vaf  = cell(sz);
t    = cell(sz); %times corresponding to vaf (different to those for forcing)
pc   = cell(sz); %pycnocline centre (i.e. forcing)
tpc  = cell(sz); %times of those pycnocline centre position

%
% Get the results for VAF and associate time and forcing and associated time
%
for i = 1:sz(1) %forcing type
for j = 1:sz(2) %member number
	%forcing
	force_path = strcat(forcing_paths(i), "/random_forcing_", num2str(j),"/anomaly.mat");
	f = load(force_path);
	pc{i,j} = f.pc;
	tpc{i,j} = f.t;

	%vaf evolution
	fnum = strcat(num2str(gammas(i)),num2str(ensembles(i)),members(i,j));
        fname =  strcat(outfile_path, 'ATTR_', fnum, '/run/outfile.nc');

	hh = ncread(fname, 'h', [1, 1, 1], [Inf,Inf,Inf]); %ice thickness
        gg = ncread(fname, 'grfrac', [1, 1, 1], [Inf,Inf,Inf]); %ice grounded fraction
	tt = ncread(fname, 'TIME');
	t{i,j} = tt;
 	vv = hh.*gg; 
	vv = sum(sum(vv,2),1)*dx*dy;
	vv = squeeze(vv);	
	vaf{i,j} = vv;
	
end %end loop over member number
end %end loop over ensembles
end %end gendata flag

%
% Make plot
%
figure(1); clf;
colmap = [70,154,158; 187,102,168]/255;

% 
% forcing profiles
%
subplot(1,2,1); hold on; box on
pc_mean = zeros(2001,2);
for i = 1:sz(1)
for j = 1:sz(2)
	tt = cell2mat(tpc(i,j));
	pp = cell2mat(pc(i,j));
	p = plot(tt,pp,'linewidth', 1);
	p.Color = [colmap(i,:),0.15];
	pc_mean(:,i) = pc_mean(:,i) + pp;
end
pc_mean = pc_mean/sz(2);
plot(tt,pc_mean(:,i), 'linewidth', 1.5, 'color', colmap(i,:));
end
plot(tt,-500 + tt, '--','linewidth', 1.25, 'color', colmap(1,:)); %add the anthro trend
plot(tt,-500*ones(size(tt)), '--','linewidth', 1.25, 'color', colmap(2,:)); %add the natural trend
 
xlim([0, 100]);
xlabel('time (yrs)');
ylabel('pycnocline depth (m)');
ylim([-750, -250]);
ax(1) = gca; ax(1).FontSize = 11;
ax(1).YTick = -700:100:-300;

%
% VAF evolution 
% 
subplot(1,2,2); hold on; box on
for i = 1:sz(1)
for j = 1:sz(2)
	tt = cell2mat(t(i,j));
	vv = cell2mat(vaf(i,j));
	volume_change = vv(1) - vv; 
	SLR = volume_change / 395 / 1e9; %SLR in mm
	p = plot(tt,SLR,'linewidth', 1.25);
	p.Color = [colmap(i,:),0.5];
end
end



xlim([0, 100]);
xlabel('time (yrs)');
ylabel('sea level rise (mm)');
ax(2) = gca;
ax(2).FontSize = ax(1).FontSize;
ax(2).YLim = [-1,25];
grid on


fig = gcf;
fig.Position(3:4) = [950, 400];
