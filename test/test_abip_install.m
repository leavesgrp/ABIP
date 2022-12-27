function [] = test_abip_install(test_qcp)
% Test ABIP installation
if nargin == 0
    test_qcp = false;
end % End if

m = 50;
n = 2000;

% Test LP
rng(24);
A = [sparse(sprand(m, n, 0.3)), speye(m)];
b = rand(m, 1);
c = rand(m + n, 1);

data.A = A;
data.b = b;
data.c = c;

K.l = m + n;
[xlp, ylp, slp, infolp] = abip(data, K); %#ok


if test_qcp
    
    param.solver = 1;
    [xlp2, ylp2, slp2, infolp2] = abip(data, K, param); %#ok
    
    clear data;
    
    % Test QCP
    data.A = sparse([1, 2, 3, 4, 5, 6, 7, 8;
        0, 1, 2, 1, 2, 3, 1, 2]);
    data.b = [4, 3]';
    data.c = [1, 0, 2, 1, 4, 2, 3, 0]';
    data.Q = speye(8);
    
    K.q = 3;
    K.rq = 3;
    K.f = 1;
    K.l = 1;
    
    [xqcp, yqcp, sqcp, infoqcp] = abip(data, K, param); %#ok
    
end % End if