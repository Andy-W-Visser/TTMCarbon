function A = d0(r)
  
  m = length(r(:));
  
  K = speye(m);
  
  [i,j] = find(K);
  
  A = sparse(i,j,r(:),m,m);