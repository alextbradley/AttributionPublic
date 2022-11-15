% 25/04/22
% Create a plot of the pycnocline centre position in the anthropogenic and
% natural forcing cases.

clear
set(0,'DefaultTextInterpreter','latex','DefaultAxesFontSize',14,'DefaultTextFontSize',14);
set(0, 'defaultaxeslinewidth', 1.25);
set(0, 'defaultlinelinewidth', .9);


anth_dir = "./anthropogenic_random_forcing/"; 
anth_dir = "./A100_T100/"; 
nat_dir = "./natural_random_forcing/";

%N = length(dir(nat_dir)) - 2; %minus two for ".." and "." directories
%N = 5;
N = 20; %first 20 entries
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

natcol = [1,165/255,0,.35];
anthcol  = [0,90/255,1,.35];
%anthropogenic forcing
subplot 211 ; hold on; box on
for i = 1:10
lh = plot(t, anth_data(i,:));
lh.Color=anthcol; %fourth entry sets the alpha
end
plot(t, mean(anth_data, 1), 'Color', anthcol(1:3), 'linewidth', 1.5) %add the mean taken at each time point
%plot(t, -500 + 50*(t/100), '--', 'Color', anthcol(1:3)); %add the mean trend (hard coded, rather than using fit)
%plot(t, -500*ones(size(t)), '--', 'Color', natcol(1:3)); %zero trend (i.e. natural forcing) (again hard coded)
title('Anthropogenic forcing');

subplot 212 ; hold on; box on
for i = 1:10
lh = plot(t, nat_data(i,:));
lh.Color = natcol; %fourth entry sets the alpha
end
plot(t, mean(nat_data, 1), 'Color', natcol(1:3), 'linewidth', 1.5) %add the mean taken at each time point
%plot(t, -500 + 50*(t/100), '--', 'Color', anthcol(1:3)); %add the mean trend (hard coded, rather than using fit)
%plot(t, -500*ones(size(t)), '--', 'Color', natcol(1:3)); %zero trend (i.e. natural forcing) (again hard coded)

title('Natural forcing');

subplot 211
plot(t, mean(nat_data, 1), 'Color', natcol(1:3), 'linewidth', 1.5) %add the mean taken at each time point

subplot 212
plot(t, mean(anth_data, 1), 'Color', anthcol(1:3), 'linewidth', 1.5) %add the mean taken at each time point


for i = 1:2
    subplot(2,1,i);
    xlabel("$t$");
    ylabel("$p_c$ (m)");
end
















