function [data, w] = normalize_data(data, scale, w)
%--------------------------------------------------------------------------
% data normalization
%
% Author: Tianyi Lin, UC Berkeley, 2018.  
%--------------------------------------------------------------------------
A     = data.A;
b     = data.b; 
c     = data.c; 
[m,n] = size(A);
 
min_scale   = 1e-3;
max_scale   = 1e3;
minRowScale = min_scale * sqrt(n);
maxRowScale = max_scale * sqrt(n);
minColScale = min_scale * sqrt(m);
maxColScale = max_scale * sqrt(m);

%% E scale
E = sqrt(sum(A.^2, 1))';
E(E < minColScale) = 1;
E(E > maxColScale) = maxColScale;
A = A*sparse(diag(1./E));

%% D scale:
D = sqrt(sum(A.^2, 2));
D(D < minRowScale) = 1;
D(D > maxRowScale) = maxRowScale;
A = sparse(diag(1./D))*A;
  
%% b and c scale
nmrowA = mean(sqrt(sum(A.^2, 2)));  
nmcolA = mean(sqrt(sum(A.^2, 1)));  

A = A*scale;

c = c./E;
cnorm = norm(c); 
% cnorm = sum(abs(c));
sc_c = nmrowA/max(cnorm, min_scale);
c = c * sc_c * scale;

b = b./D;
bnorm = norm(b); 
% bnorm = sum(abs(b)); 
sc_b = nmcolA/ max(bnorm, min_scale);
b = b * sc_b * scale;



%% normalized (A,b,c) record
data.A = A; 
data.b = b;
data.c = c;

w.D = D;
w.E = E;
w.sc_b = sc_b;
w.sc_c = sc_c;

end