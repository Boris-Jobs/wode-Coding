function [L, S] = RobustPCA(X, lambda, mu, tol, max_iter)%5个参数
    % - X 是将要被分解的N×M的数据矩阵
    %   X can also contain NaN's for unobserved values
    % - lambda - regularization parameter, default = 1/sqrt(max(N,M))
    % - mu - the augmented lagrangian parameter, default = 10*lambda
    % - tol - reconstruction error tolerance, default = 1e-6
    % - max_iter - maximum number of iterations, default = 1000
    %X=L+S, lambda, Y, mu, ρ=1.5

    [M, N] = size(X);
    unobserved = isnan(X);%返回一个与A相同维数的数组
    %若A的元素为NaN（非数值），在对应位置上返回逻辑1（真），否则返回逻辑0（假）
    X(unobserved) = 0;
    normX = norm(X, 'fro');%返回矩阵 X 的 Frobenius 范数

    % default arguments，nargin为输出参数的个数，一共5个参数
    %下列的if与end这一段是在设定默认值
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
    
    % initial solution
    L = zeros(M, N);
    S = zeros(M, N);
    Y = zeros(M, N);%相当于RPCA的λ
    
    for iter = (1:max_iter)
        % ADMM step: update L and S
        L = Do(1/mu, X - S + (1/mu)*Y);
        S = So(lambda/mu, X - L + (1/mu)*Y);
        % and augmented lagrangian multiplier
        Z = X - L - S;
        Z(unobserved) = 0; % skip missing values跳过不存在的值
        Y = Y + mu*Z;
        mu=1*mu;%【原程序增加项，加了这项以后速度加了个数量级，但sparse就不那么明显】
        err = norm(Z, 'fro') / normX;%用误差值的二范数来衡量分解程度的
        if (iter == 1) || (mod(iter, 10) == 0) || (err < tol)
            fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
                    iter, err, rank(L), nnz(S(~unobserved)));
                %Number of nonzero matrix elements
        end%第一次||每10次||误差满足时，fprintf一次
        
        if (err < tol)
            break; 
        end
        
    end
end

function r = So(tau, X)
    % shrinkage operator
    r = sign(X) .* max(abs(X) - tau, 0);%返回X的阈值函数取值
    %S单变量时就是取这个
end

function r = Do(tau, X)
    % shrinkage operator for singular values
    [U, S, V] = svd(X, 'econ');%economy size,
    %For m < n, only the first m columns of V arecomputed and S is m-by-m.
    r = U*So(tau, S)*V';
    %L单变量时就是取这个
end
