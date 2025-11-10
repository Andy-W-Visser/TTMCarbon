load('CTL.mat') % lower resolution - less than a minute
firebrick = '#B22222';
pale_golden_rod = '#EEE8AA';
%% Preparation of the transport matrix
grid = output.grid;
TR = output.TR; % yr^-1
msk = output.msk;
M3d = output.M3d; % land = 0. ocean = 1
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
V = VOL(msk.pkeep);

[lp1,fp1,lp2,fp2] = eqage(TR,grid,M3d); % age in years

ibottom = sum(M3d,3)
lp1(isnan(lp1)) = 0;
fp1(isnan(fp1)) = 0;
DM3dp = M3d(:,:,1:end-1) - M3d(:,:,2:end);
DM3d = cat(3,DM3dp,M3d(:,:,end));
LP1bottom = sum(lp1.*DM3d,3);
FP1bottom = sum(fp1.*DM3d,3);


%% 

LAT = grid.YT; LAT = [LAT(:,1) LAT LAT(:,end)];
LON = grid.XT; LON = [LON(:,1)-2 LON LON(:,end)+2];
LP1bottom = [LP1bottom(:,end) LP1bottom LP1bottom(:,1)];

figure(2);
clf;
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%worldmap world
surfm(LAT, LON, FP1bottom);
colormap(turbo)
geoshow('landareas.shp','FaceColor','#EEE8AA')
colorbar
title('First passage time from sea floor (year)')

figure(3);
clf;
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%worldmap world
surfm(LAT, LON, LP1bottom);
colormap(turbo)
geoshow('landareas.shp','FaceColor','#EEE8AA')
colorbar
title('Last passage time from sea floor (year)')