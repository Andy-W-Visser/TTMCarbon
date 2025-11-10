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

m = size(TR,1);
sink = zeros(m,1);
sink(1:length(msk.hkeep)) = 1e10; % instantaneous Surface SINK
SSINK = spdiags(sink,0,m,m);
A = TR-SSINK; % transport + sink in surface

y = -90:2:90; ny = length(y);
x = 0:2:360; nx = length(x);
z = grid.zt; nz = length(z);
[latg,long] = meshgrid(y, x);

%% Input flux

bottomindex = sum(M3d,3); % bottom at zwstar(bottomindex + 1), bottomindex = 0 => land 

latx = 60;
lonx = -12;
lati = (latx + 90)/2;
loni = (360 + lonx)/2;

TFluxcopepod = 80*1E3*1E6; %g/year
area = output.grid.Areat(lati,loni);
qcopepod = TFluxcopepod/area/100;


boti = bottomindex(lati,loni);

Q = zeros(size(M3d));
Q(lati,loni,boti) = qcopepod;

q = sum(Q,3);

figure(1);
clf;
imagesc(q); axis xy;

figure(5); clf;
imagesc(bottomindex); axis xy;



%%

q_ocim = Q(msk.pkeep); 
q_ocim(isnan(q_ocim)) = 0;
export = V'*q_ocim / 1e15 % [PgC / yr]
tic
cseq = -A \q_ocim;
toc
TotCseq = V'*cseq / 1e15 % grams to Pg
seqtime = TotCseq / export

% Re interpolate to Q
qocim = 0*M3d+NaN; % make a 3-d array of NaN
qocim(msk.pkeep) = q_ocim;
qocim(isnan(qocim)) = 0;

C_eq = 0*M3d+NaN;
C_eq(msk.pkeep) = cseq;
C_eq(isnan(C_eq)) = 0;

%% A simple plot
cc = sum(C_eq.*grid.DZT3d,3); cc = [cc, cc(:,end)]; % gC/m2



figure(2);
clf;
cc(cc<0) = 0;

axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
title('Sequestered DIC (gC/m^{-2})')
surfm(latg,long,cc');
colormap(cool);
geoshow('landareas.shp','FaceColor','#EEE8AA')
colorbar
