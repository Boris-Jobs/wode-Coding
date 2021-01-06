function [L, S] = RobustPCA_laplacian(X, lambda, mu ,tol ,max_iter, beta)
    % - X is a data matrix (of the size N x M) to be decomposed
    %   X can also contain NaN's for unobserved values
    % - lambda - regularization parameter, default = 1/sqrt(max(N,M))
    % - mu - the augmented lagrangian parameter, default = 10*lambda
    % - tol - reconstruction error tolerance, default = 1e-6
    % - max_iter - maximum number of iterations, default = 1000

    [M, N] = size(X);
    unobserved = isnan(X);
    X(unobserved) = 0; 
    normX = norm(X, 'fro');
     
    % default arguments
    if nargin < 2 
        lambda = 1 / sqrt(max(M,N));
    end
    if nargin < 3
        mu = 10*lambda;
    end
    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5
        max_iter = 1000;
    end
    if nargin < 6
        beta=1.1;
    end
    
    % initial solution
    L = zeros(M, N);  %X=L+S
    S = zeros(M, N);
    Y1 = zeros(M, N);
    Y2 = zeros(M, N);
    w = zeros(M, M);
    MF =zeros(M, M);
    p=1.1;
    for i=1:M
        for j=1:M
            delta=norm(X(i,:)-X(j,:),'fro')/2;
            w(i,j)=exp(-delta);
        end
    end
    for i=1:M  
        for j=1:M
            if i~=j
                MF(i,j)=-w(i,j);
            end
            if i==j
                for k=1:M
                    MF(i,j)=MF(i,j)+w(i,k);
                end
                MF(i,j)=MF(i,j)-w(i,j);
            end
        end
    end
    for iter = (1:max_iter)
        % ADMM step: update L and S
        L = Do(1/mu, X - S +1/mu*Y1); 
        H =(2*beta*MF+mu*eye(M))^(-1)*(mu*S+Y2);
        S = So(lambda/mu, X - L + (1/mu)*(Y1-Y2)+H); 
        % and augmented lagrangian multiplier
        Z = X - L - S;  
        Z(unobserved) = 0; % skip missing values
        Z1=S-H;
        Z1(unobserved) = 0;
        Y1 = Y1 + mu*Z; 
        Y2 = Y2 + mu*Z1;
        mu=p*mu;
        err = norm(Z, 'fro') / normX; 
        if (iter == 1) || (mod(iter, 10) == 0) || (err < tol)
            fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
                    iter, err, rank(L), nnz(S(~unobserved)));
        end
        if (err < tol)
            break; 
        end
    end
end

function r = So(tau, X)  
    % shrinkage operator
    r = sign(X) .* max(abs(X) - tau, 0);
end

function r = Do(tau, X)
    % shrinkage operator for singular values
    [U, S, V] = svd(X, 'econ');
    r = U*So(tau, S)*V'; 
end
