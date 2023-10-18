function [A,b,lambda] = get_lasso_simu_data(m_lasso,n_lasso,seed)

    rng(seed);
    A = sprandn(m_lasso ,n_lasso, 0.15);
    v = zeros(n_lasso,1);
    for i = 1:n_lasso
        if rand() > 0.5
            v(i) = normrnd(0,1/n_lasso);
        end
    end
    epsilon = randn(m_lasso,1);
    b = A*v + epsilon;
    lambda_max = norm(A'*b,Inf);
    lambda = lambda_max / 5;
    
end
