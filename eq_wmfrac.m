function F = eq_wmfrac(TR,REG,M3d,grid)
% Want to solve dF/dt = A*F subject to  F prescribed at the surface
%
% Partition F and A into surface and interior grid points such that
%
% d [F_s] = [A_ss A_si][F_s]
% dt[F_i]   [A_is A_ii][F_i]
%
% If F_s is prescribed then  d/dt F_i = A_ii*F_i +A_is*F_s
% Note that A_is*F_s is the effective source from having prescribed
% the surface sst.
% 
% To compute A_is*F_s we use the sbc=0 transport operator and compute
% A*F where
% F = [F_s] so that A*F = [A_ss A_si][F_s] = [A_ss*F_s]
%     [ 0 ]               [A_is A_ii][ 0 ]   [A_is*F_s]
% thus the interior grid points of A*F give the needed quantity.

% partition the transport operator
n = size(TR,1);
[ny,nx,nz] = size(M3d);
iocn = find(M3d(:)==1);
isurf = find(M3d(:,:,1)==1);
iint = setdiff(iocn,isurf);
ns = length(isurf);
Aii = TR(ns+1:end,ns+1:end);
Ais = TR(ns+1:end,1:ns);

% solve for the interior distribution given surface b.c.
f = -Aii\(Ais*REG(isurf));

% make a 3-d array
F = 0*M3d+NaN;
F(isurf) = REG(isurf);
F(iint) = f;
