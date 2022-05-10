% 28/04/22
% 
% Plot the grounding line position as a function of time for both anthropogenic and natural forcing
%
%
% Preliminaries
%
clear
set(0,'DefaultTextInterpreter','latex','DefaultAxesFontSize',12,'DefaultTextFontSize',12);
gendata = 1; %flag to specify dats
outdir ="/data/icesheet_output/aleey/wavi/"; %outputs main dir

%
% Grid info
%
ny = 50;
nx = 300;

%
% Run Spec
%
run_nos = ["11101", "11102", "11103", "11104", "11105";  %anthropogenic
           "11201", "11202", "11203", "11204", "11205"]; %natural
forcing = ["1/", "2/", "3/","4/", "5/"];     
forcing_files_folder = "../gendata/ocean/ocean_forcing/";
forcing_files = [strcat("anthropogenic_random_forcing/random_forcing_", forcing);strcat("natural_random_forcing/random_forcing_", forcing)];
forcing_files = strcat(forcing_files_folder, forcing_files);
sz = size(run_nos);       
   
%
% Gendata loop
ss = struct;
for iF = 1:2 %forcing type
    for iN = 1:sz(2) %ensemble member no
        fname = strcat(outdir, "ATTR_", run_nos(iF, iN), "/run/outfile.nc");
       
        %load raw quantities
        ss(iF, iN).t = ncread(fname, 'TIME');
        ss(iF, iN).h = ncread(fname, 'h');
        ss(iF, iN).x = ncread(fname, 'x');
        ss(iF, iN).y = ncread(fname, 'y');
        ss(iF, iN).b = ncread(fname, 'b');
        ss(iF, iN).u = ncread(fname, 'u');
        ss(iF, iN).v = ncread(fname, 'v');
        ss(iF, iN).m = ncread(fname, 'm');
        ss(iF, iN).grfrac = ncread(fname, 'grfrac');
        
        %compute grounding line position
        ss(iF, iN).xgl = zeros(1,length(ss(iF, iN).t));
        for it = 1:length(ss(iF,iN).t)
            centreline_grfrac = squeeze(ss(iF,iN).grfrac(:,floor(ny/2),it));
            idx = find(centreline_grfrac < 1, 1, 'first');
            ss(iF,iN).xgl(it) = ss(iF,iN).x(idx);
        end
        
        %load the pycnocline depth
        pc = load(strcat(forcing_files(iF, iN), 'anomaly.mat'));
        ss(iF,iN).pc = pc.pc;
        ss(iF,iN).pct = pc.t;
    end
end

%
% Make plot
%
figure(1); clf;
subplot 211 ; hold on; box on
colmap = [255,165,0; 0,90,255]/255;

for iF = 1:2
    mean_profile = zeros(size(ss(1,1).t));
    for iN = 1:sz(2)
        lh = plot(ss(iF, iN).t, ss(iF,iN).xgl/1e3);
        lh.Color=[colmap(iF,:),0.2]; %fourth entry sets the alpha
        lh.HandleVisibility = 'off';

        
        mean_profile = mean_profile + ss(iF, iN).xgl'/sz(2);
    end
    %add the mean retreat profile
    plot(ss(iF, iN).t, mean_profile/1e3, 'color', colmap(iF, :));
    
end


ylabel('$x_{gl}$ (km)')
legend({'Anthropogenic', 'Natural'}, 'Interpreter', 'latex', 'location', 'southwest')

% add the forcing files
subplot 212; hold on; box on
for iF = 1:2
    mean_profile = zeros(size(ss(1,1).pct));
    for iN = 1:sz(2)
        lh = plot(ss(iF, iN).pct, ss(iF,iN).pc);
        lh.Color=[colmap(iF,:),0.2]; %fourth entry sets the alpha
        lh.HandleVisibility = 'off';

        
        mean_profile = mean_profile + ss(iF, iN).pc/sz(2);
    end
    %add the mean forcing profile
    plot(ss(iF, iN).pct, mean_profile, 'color', colmap(iF, :));
    
end
ylabel('pycnocline center depth (m)')

for i = 1:2
    subplot(2,1,i);
    xlabel('time (yrs)')
    %xlim([0, 100])
end

