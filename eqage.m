function [lp1,fp1,lp2,fp2] = eqage(TR,grid,M3d)
% [lp1,fp1,lp2,fp2]=eqage(TR,grid,M3d);
%
% Compute the first and centred second moments of the last passage time
% distribution (LP1 & P2) and of the first passage time distribution (FP1
% & FP2)
%
% d(age)/dt -TR*age = 1
% age =0 at sea surface
%
% d  [age(isrf)]   [ TR(isrf,isrf)  TR(isrf,iint)][age(isrf)]   [1(isrf)]
% -  [         ] - [                             ][         ] = [       ]
% dt [age(iint)]   [ TR(iint,isrf)  TR(iint,iint)][age(iint)]   [1(iint)]
% subject to   age(isrf) = 0
%
% At steady state d/dt -->0
%
% TR(iint,iint)*age(iint) = 1(iint)  ==> age(iint) = TR(iint,iint)\1(iint)
%
% NOTE: The units of the transport operator is yr^-1
%
[ny,nx,nz] = size(M3d);
% land sea mask
iocn = find(M3d(:)==1);
ilnd = find(M3d(:)==0);
% find the surface points
Q = 0*M3d; Q(:,:,1) = 1;
iint = find(Q(iocn)==0);
%
A = TR(iint,:); A = A(:,iint);
%
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
W = d0(VOL(iocn(iint)));
%
FA = mfactor(A);
%
n = length(iint);
rhs = ones(n,1);
%
% mean last passage time (ideal age)
%
lp1 = zeros(ny,nx,nz);
lp1(iocn(iint)) = -mfactor(FA,rhs);
lp1(ilnd) = NaN;
% centered second moment of the last passage time distribution
lp2 = zeros(ny,nx,nz);
lp2(iocn(iint)) = -2*mfactor(FA,lp1(iocn(iint)));
lp2 = sqrt(lp2-lp1.^2);
lp2(ilnd) = NaN;
%
% mean first passage time
%
fp1 = zeros(ny,nx,nz);
fp1(iocn(iint)) = -W\mfactor(FA,W*rhs,'transpose');
fp1(ilnd) = NaN;
% centered second moment of the first passage time distribution
fp2 = zeros(ny,nx,nz);
fp2(iocn(iint)) = -2*mfactor(FA,fp1(iocn(iint)));
fp2(ilnd) = NaN;
fp2 = sqrt(fp2-fp1.^2);

