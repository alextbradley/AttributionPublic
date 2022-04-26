% 25/04/22
% Create a plot of the pycnocline centre position in the anthropogenic and
% natural forcing cases.

clear
set(0,'DefaultTextInterpreter','latex','DefaultAxesFontSize',14,'DefaultTextFontSize',14);
set(0, 'defaultaxeslinewidth', 1.25);
set(0, 'defaultlinelinewidth', .9);


anth_dir = "./anthropogenic_random_forcing/"; 
nat_dir = "./natural_random_forcing/";

N = length(dir(nat_dir)) - 2; %minus two for ".." and "." directories
%N = 5;
%N = 1;
nt = 2001;  %number of time points

%for each directory, load the pycnocline centre and plot
anth_data = zeros(N,nt);
nat_data  = zeros(N,nt);
for iN= 1:N
    %anthropogenic
    anth = load(strcat(anth_dir, "random_forcing_", num2str(iN), "/anomaly.mat"));
    anth_data(iN,:) = anth.pc;
    
    %natural
    nat = load(strcat(nat_dir, "random_forcing_", num2str(iN), "/anomaly.mat"));
    nat_data(iN,:) = nat.pc;
    
end

t = nat.t; %same time output for all forcing files

%make plot
figure(1); clf; 

%anthropogenic forcing
subplot 211 ; hold on; box on
for i = 1:1
lh = plot(t, anth_data(i,:));
lh.Color=[0,0,1,0.2]; %fourth entry sets the alpha
end
plot(t, mean(anth_data, 1), 'b') %add the mean taken at each time point
plot(t, -500 + 50*(t/100), 'b--'); %add the mean trend (hard coded, rather than using fit)
plot(t, -500*ones(size(t)), 'r--') %zero trend (i.e. natural forcing) (again hard coded)
title('Anthropogenic forcing');

subplot 212 ; hold on; box on
for i = 1:N
lh = plot(t, nat_data(i,:));
lh.Color=[1,0,0,0.2]; %fourth entry sets the alpha
end
plot(t, mean(nat_data, 1), 'r') %add the mean taken at each time point
plot(t, -500 + 50*(t/100), 'b--'); %anthropogenic trend
plot(t, -500*ones(size(t)), 'r--') %natural trend 
title('Natural forcing');


for i = 1:2
    subplot(2,1,i);
    xlabel("$t$");
    ylabel("$p_c$ (m)");
end
















