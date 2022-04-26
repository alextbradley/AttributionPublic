%%% 04/04/22 
%%% Generate "random" ambient forcing anomaly files.
%%% 
%%% Generate ambient forcing files with pycnocline centres oscillated
%%% according to AR1 red noise. This function creates N unseeded forcing
%%% profiles. For each forcing profile, we generate the associated ref
%%% temp, salinity, depth and time vectors used in the quadratic forcing
%%% melt rate parameterization.
%%% 

clear
%% Input parameters
N  = 20;             %Number of forcing files
dt = 0.1;           %timestep of output
r = 0.8;            %autocorrelation of red noise
pyc_amp = 100;      %amplitude of the oscillation of the pycnocline
pyc_mean = -500;    %mean of the untrended pycnocline centre 
t_end = 200;        %end time
t = 0:dt:t_end;     %reference time points
depth = 0:10:1000;  %reference depth points

pc_inc_per_century =  50;      %specify how much the pycnocline mean depth increases per century if anthropogenic applied 

for iN = 1:N
    %make natural files 
    folder = strcat("./natural_random_forcing/random_forcing_", num2str(iN), "/");
     if ~exist(folder, 'dir'); mkdir(folder); end
    
    % Parameters of profile
    pc = pyc_mean + pyc_amp*generate_random_forcing_anomaly(t,r);  %pycnocline center
        
    pw = 400*ones(length(t));                 %pycnocline width  (same size as time)  
    % constant upper and lower temperatures and salinities
    tl = 1.2;
    tu = -1;
    sl = 34.6;
    su = 34.0;
    
    % Generate profiles
    ta = ambient_profile(pc, pw, -depth, tl, tu);
    sa = ambient_profile(pc, pw, -depth, sl, su);
    
%    % Make video
%     figure(1); clf;
%     for i = 1:length(t)
%         subplot 121
%         plot(ta(i,:), -depth, 'r');
%         
%         subplot 122
%         plot(sa(i,:), -depth, 'b');
%         
%         sgtitle(['t = ' num2str(t(i))])
%         drawnow
%         shg
%         %pause
%     end
    
    
    % Save files
    acc = 'real*8';
    
    %reference time
    fname_out = strcat(folder, 'ref_time.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,t,acc);fclose(fid);
    
    %reference depth
    fname_out = strcat(folder, 'ref_depth.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,depth,acc);fclose(fid);
    
    %reference ambient temp
    fname_out = strcat(folder, 'forcing_ambient_temperature.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,ta,acc);fclose(fid);
    
    %reference ambient salt
    fname_out = strcat(folder, 'forcing_ambient_salinity.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,sa,acc);fclose(fid);
    
    %save mat files of the pycnocline centre
    save(strcat(folder, 'anomaly.mat'), 't', 'pc');
    
    %pause
    
    %make anthropognic files 
    folder = strcat("./anthropogenic_random_forcing/random_forcing_", num2str(iN), "/");
    if ~exist(folder, 'dir'); mkdir(folder); end

      
    % Parameters of profile
    pc = pc +  pc_inc_per_century*(t'/100); %add the anthropogenic contribution  
    pw = 400*ones(length(t));                 %pycnocline width  (same size as time)  
    
    % constant upper and lower temperatures and salinities
    tl = 1.2;
    tu = -1;
    sl = 34.6;
    su = 34.0;
    
    % Generate profiles
    ta = ambient_profile(pc, pw, -depth, tl, tu);
    sa = ambient_profile(pc, pw, -depth, sl, su);
   
    
    % Save files
    acc = 'real*8';
    
    %reference time
    fname_out = strcat(folder, 'ref_time.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,t,acc);fclose(fid);
    
    %reference depth
    fname_out = strcat(folder, 'ref_depth.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,depth,acc);fclose(fid);
    
    %reference ambient temp
    fname_out = strcat(folder, 'forcing_ambient_temperature.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,ta,acc);fclose(fid);
    
    %reference ambient salt
    fname_out = strcat(folder, 'forcing_ambient_salinity.bin');
    fid=fopen(fname_out,'w','b');
    fwrite(fid,sa,acc);fclose(fid);
    
    %save mat files of the pycnocline centre
    save(strcat(folder, 'anomaly.mat'), 't', 'pc');
    
end
%% Functions
function va = ambient_profile(pc, pw, z, vl, vu)
%loop over
va = zeros(length(pc), length(z));
for i = 1:length(pc) %loop over time points
    for j = 1:length(z) %loop over grid points
        
        va(i,j) = pointwise_ambient_profile(pc(i), pw(i), z(j), vl, vu);
        
    end
end
end


function va = pointwise_ambient_profile(pc, pw, z, vl, vu)
%return the ambient value at a depth z with the pycnocline centre at pc,
%width pw
d_low = pc - pw/2;
d_hi  = pc + pw/2;
if z< d_low
    va = vl;
    
elseif z > d_hi
    va = vu;
    
else
    va = vl + (vu - vl)/(d_hi - d_low) * (z - d_low);
    
end
end


