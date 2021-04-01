function [A,b,c,info] = presolve(A,b,c,lbounds,ubounds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Presolve Input Data.                          
% .
% Usage:    [A,b,c,info] = presolve(A,b,c,lbounds,ubounds)
% Problem:   min c'x, s.t. Ax=b, lbounds<=x<=ubounds. 
% Input:    Data:            (A, b, c) 
%           LowerBound:      lbounds 
%           UpperBound:      ubounds
% 
% Author: Tianyi Lin, UC Berkeley. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% store original data.
[m, ~] = size(A);
if ~issparse(A)
    A = sparse(A); 
end;
b = sparse(b); 
c = sparse(c);
lbounds = sparse(lbounds);
info.b = b; 
info.c = c;

%% initialize statue.
info.feasible      = 1;
info.data_changed  = 0;
info.fixed_exist   = 0;  
info.zrcols_exist  = 0;
info.sgtons_exist  = 0;
info.free_exist    = 0;
info.lbounds_non0  = 0; 
info.nub           = 0;
info.ubounds_exist = 0; 

%% detect infeasible variable.
if any(lbounds > ubounds)
   info.feasible = 0; 
   return;
end;

%% delete fixed variables.
fixed = (lbounds == ubounds);
info.fixed_exist = any(fixed); 
if (info.fixed_exist)
   ifix = find(fixed); 
   infx = find(1 - fixed);
   xfix = lbounds(ifix);
   c    = c(infx);
   b    = b - A(:,ifix)*sparse(xfix);
   A    = A(:,infx);
   lbounds = lbounds(infx);
   ubounds = ubounds(infx);
   
   % record the information.
   info.data_changed = 1;
   info.ifix = ifix; 
   info.infx = infx; 
   info.xfix = xfix;
end 

%% calculate row statistics.
rnnzct = sum(spones(A'));

%% delete zero rows.
if any(rnnzct == 0)
   izrows = rnnzct == 0;
   if any(b(izrows) ~= 0) 
      info.feasible = 0; 
      return;
   end
   
   inzrows = find(rnnzct > 0);
   A = A(inzrows,:); 
   b = b(inzrows); 
   rnnzct = rnnzct(inzrows);
   info.data_changed = 1; 
   [m, ~] = size(A);
end

%% remove all linearly dependent rows.
sprk = sprank(A');
if (sprk < m)
   [dmp, ~] = dmperm(A);
   irow = dmp(1:sprk);
   A = A(irow,:); 
   b = b(irow); 
   rnnzct = rnnzct(irow);
   info.data_changed = 1;
end

%% delete zero columns.
zrcol = (max(spones(A)) == 0)';                    % the maximum absolute value of the element in each column
if any(zrcol == 1)
   info.zrcols_exist = 1;
   izrcol = find(zrcol);
   temp = ubounds(izrcol);
   if any(c(izrcol) < 0 & temp == Inf)
      info.feasible = 0; 
      return;
   end
   
   % variable x associated with zero columns.
   temp(isinf(temp)) = 0; 
   xzrcol = zeros(size(izrcol)) + (c(izrcol) < 0).*temp + (c(izrcol) > 0).*lbounds(izrcol);
   inzcol = find(1 - zrcol);
   A = A(:,inzcol);
   c = c(inzcol);
   lbounds = lbounds(inzcol);
   ubounds = ubounds(inzcol);
   
   % record the information. 
   info.izrcol = izrcol; 
   info.xzrcol = xzrcol; 
   info.inzcol = inzcol; 
   info.data_changed = 1;
end

%% delete singleton rows.
singleton = (rnnzct == 1);
nsgrows = nnz(singleton);                       % the number of singleton rows. 
if nsgrows >= max(1, 0.01*size(A,1))
   info.sgtons_exist = 1;
   isgrows = find(singleton);                   % the index of singleton rows. 
   iothers = find(1 - singleton);               % the index of other rows. 

   Atmp  = A(isgrows,:);                        % the singleton rows.  
   btmp  = b(isgrows);                          % the right-hand number. 
   Atmp1 = spones(Atmp);               
   
   if nsgrows == 1 
      isolved  = find(Atmp1);                   % the index of the singleton. 
      insolved = find(Atmp1 == 0);              % all the other components. 
      xsolved  = btmp/Atmp(isolved);
   else
      colnnzct = sum(Atmp1);                    % summing all rows to see if dependent rows exist. 
      isolved  = find(colnnzct);
      insolved = find(colnnzct == 0);
      [ii, jj] = find(Atmp); 
      Atmp = Atmp(ii,jj);                       % the diagonal elements are non-zero singletons in Atmp. 
      btmp = btmp(ii);                          % just the value of btmp. 
      xsolved  = btmp./diag(Atmp);              % solve it. 
      if any(colnnzct >  1)
         repeat = diff([0; jj]) == 0;           % diff() is the difference and approximation derivatives. 
         for i = 1:length(xsolved) - 1
            if repeat(i+1) && xsolved(i+1) ~= xsolved(i)
               info.feasible = 0;               % if repeat==1, then this means that two elements should 
               return;                          % be the same. If different, then infeasible. 
            end; 
         end;
         ii = find(~repeat); 
         jj = ii;
         Atmp = Atmp(ii, jj); 
         btmp = btmp(ii);
         xsolved  = btmp./diag(Atmp);           % this is the solution we want. 
      end;
   end;

   if any(xsolved < lbounds(isolved)) || any(xsolved > ubounds(isolved))
      info.feasible = 0; 
      return;
   end
   
   b = b(iothers) - A(iothers,isolved)*xsolved;
   A = A(iothers, insolved);
   c = c(insolved);
   lbounds = lbounds(insolved);
   ubounds = ubounds(insolved);
   
   % record the information.
   info.isolved = isolved; 
   info.insolved = insolved; 
   info.xsolved = xsolved; 
   info.data_changed = 1;
end;

%% delete free variables.
free = (lbounds == -Inf & ubounds == Inf);
info.free_exist = any(free); 
if (info.free_exist)
   ifree = find(free);
   nfree = length(ifree);
   infree = find(1 - free);
   c = [c(infree); c(ifree); -c(ifree)]; 
   A = [A(:,infree) A(:,ifree) -A(:,ifree)]; 
   lbounds = [lbounds(infree); zeros(2*nfree, 1)];
   ubounds = [ubounds(infree); Inf(2*nfree, 1)];
   
   % record the information.
   info.data_changed = 1;
   info.ifree = ifree; 
   info.infree = infree; 
   info.nfree = nfree;
end

%% shift nonzero lower bounds.
info.lbounds_non0 = any(lbounds ~= 0);
if (info.lbounds_non0)
   b = b - A*lbounds;
   info.data_changed = 1;
end
info.lbounds = lbounds; 

%% find upper bounds.
iubounds = (ubounds ~= Inf);
info.ubounds_exist = full(any(iubounds));
if (info.ubounds_exist)
   b_ubounds = ubounds(iubounds)-lbounds(iubounds); 
   nub = length(b_ubounds);
   info.data_changed = 1;
end

[m, n] = size(A); 
info.m = m; 
info.n = n;

%% reformulate LP.
if (info.ubounds_exist)
    c = [c; sparse(nub, 1)]; 
    b = [b; b_ubounds];  
    E = sparse(nub, n);
    E(sub2ind(size(E),1:nub,find(iubounds>0)')) = 1; 
    A = [A sparse(m, nub);E speye(nub)];
    info.data_changed = 1; 
end
    