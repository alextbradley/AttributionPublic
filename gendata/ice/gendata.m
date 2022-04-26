% 12/4/21
% Generate data for a case similar to De Rydt and Gudmundsson 2016, JGR
% Earth Surface. 
% The valley is 50 km-wide, 300km long and 1100 m-deep. It features a ridge
% towards the downstream end of the valley. 

clear
numplot = 1;
%
% Grid info
%
dx = 1000; 
dy = 1000;
nx = floor(300*1e3/dx);
ny = floor(50*1e3/dy);
x  = 0:dx:(nx-1)*dx;
y  = (-25*1e3 + dx/2):dx:(25*1e3 - dx/2);
[xx,yy] = meshgrid(x,y);

%
% Bathymetry
%
bathy_valley = -500-600*sin(pi/2 +  pi*yy / 5e4);
sigma_b      = 1.1e4;
bathy_bump   =  -1100 + 400*exp(-(xx - 2.65e5).^2 ./2 /sigma_b^2);
bathy = max(bathy_valley,bathy_bump);
bathy = bathy'; %nx x ny

%save
fid=fopen('bathy.bin','w','b');
fwrite(fid,bathy,'real*8');fclose(fid);

%
% Rate factor
%
% the rate factor takes the value A = 2.94 × 10−9 a−1 kPa−3 for -20km < y <
% 20km (corresponding to -20C ice) and 5.04 × 10−9 a−1 kPa−3 for -25km < y
% < -20km and 20km < y < 25km (corresponding to -10C ice)
%
rate_factor = 2.94*1e-9*1e-9*ones(nx,ny);   %extra factor 1e-9 for the kPa conversion
idx_margin  = (yy < -20*1e3) | (yy > 20*1e3);
rate_factor(idx_margin') = 5.04 * 1e-9 * 1e-9;

%save
fid=fopen('glen_a_ref.bin','w','b');
fwrite(fid,rate_factor,'real*8');fclose(fid);


%
% Weertman C
%
% Paper says "Constant Weertman C of 20 m a−1 kPa−3". This corresponds to C
% = (20*1e-9)^(-1/3) = 368.4 oin WAVI. But it's constant so this can just
% be entered explicitly in the driver code
%
% n = 3;
% weertman_c = 20^n*ones(nx,ny); %is this right?
% fid=fopen('weertman_c.bin','w','b');
% fwrite(fid,weertman_c,'real*8');fclose(fid);

%
% Accumulation
%
% Accumulation varies from 15ma−1 at the ice divide (x=0km) to 1ma−1 at
% x=150km and is set to a constant value of 1ma−1 between x = 150km and the
% ice front at x = 300km
accumulation = ones(nx,ny);
[~,idx] = min(abs(x - 150*1e3));
accumulation(1:idx,:) = 15 - (15 - 1)/150e3 *(xx(:,1:idx)');
total_accumulation = sum(sum(accumulation))*dx*dy;
fprintf('total accumulation is %.1f Gt/yr \n', total_accumulation/1e9);

fid=fopen('accumulation.bin','w','b');
fwrite(fid,accumulation,'real*8');fclose(fid);

%
% initial condition
%
% take an initial condition corresponding to grounding at the ridge crest
% with surface Sm above sea level and shelf thickness equal to grounding
% thickness
[~,idx_ridge] = min(abs(x - 265*1e3));
S = 125; 
h_init = zeros(nx,ny);
h_init(1:idx_ridge+1,:) = -bathy(1:idx_ridge+1,:) + S;
for j = (idx_ridge+2:nx)
    h_init(j,:) = -bathy(idx_ridge+2,:) + S ;
end

fid=fopen('h_init.bin','w','b');
fwrite(fid,h_init,'real*8');fclose(fid);

%% %make a picture
fig =figure(1); clf; 
fig.Position(3:4) = [1400, 650];
subplot 231
contourf(x,y,bathy', 30, 'linestyle', 'none' ); 
c =  colorbar; c.Label.String = 'bathymetry';
xlabel('x');
ylabel('y');

subplot 232
idx = floor(ny/2);
plot(x,bathy(:,idx));
xlabel('x');
ylabel('centreline bed');
hold on
plot(x,S*ones(1,length(x)), 'k');
plot(x,S*ones(1,length(x)) - h_init(:,idx), 'r');
legend({'bathy', 'surface', 'ice base'}, 'location', 'northwest')

subplot 233
contourf(x,y,rate_factor',20, 'linestyle', 'none'); c =  colorbar; c.Label.String = 'rate factor';
xlabel('x');
ylabel('y');

subplot 234
contourf(x,y,accumulation', 20, 'linestyle', 'none'); c =  colorbar; c.Label.String = 'accumulation (m/yr)';
xlabel('x');
ylabel('y');


