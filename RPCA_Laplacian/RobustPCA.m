function [L, S] = RobustPCA(X, lambda, mu, tol, max_iter)%5������
    % - X �ǽ�Ҫ���ֽ��N��M�����ݾ���
    %   X can also contain NaN's for unobserved values
    % - lambda - regularization parameter, default = 1/sqrt(max(N,M))
    % - mu - the augmented lagrangian parameter, default = 10*lambda
    % - tol - reconstruction error tolerance, default = 1e-6
    % - max_iter - maximum number of iterations, default = 1000
    %X=L+S, lambda, Y, mu, ��=1.5

    [M, N] = size(X);
    unobserved = isnan(X);%����һ����A��ͬά��������
    %��A��Ԫ��ΪNaN������ֵ�����ڶ�Ӧλ���Ϸ����߼�1���棩�����򷵻��߼�0���٣�
    X(unobserved) = 0;
    normX = norm(X, 'fro');%���ؾ��� X �� Frobenius ����

    % default arguments��narginΪ��������ĸ�����һ��5������
    %���е�if��end��һ�������趨Ĭ��ֵ
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
    Y = zeros(M, N);%�൱��RPCA�Ħ�
    
    for iter = (1:max_iter)
        % ADMM step: update L and S
        L = Do(1/mu, X - S + (1/mu)*Y);
        S = So(lambda/mu, X - L + (1/mu)*Y);
        % and augmented lagrangian multiplier
        Z = X - L - S;
        Z(unobserved) = 0; % skip missing values���������ڵ�ֵ
        Y = Y + mu*Z;
        mu=1*mu;%��ԭ������������������Ժ��ٶȼ��˸�����������sparse�Ͳ���ô���ԡ�
        err = norm(Z, 'fro') / normX;%�����ֵ�Ķ������������ֽ�̶ȵ�
        if (iter == 1) || (mod(iter, 10) == 0) || (err < tol)
            fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
                    iter, err, rank(L), nnz(S(~unobserved)));
                %Number of nonzero matrix elements
        end%��һ��||ÿ10��||�������ʱ��fprintfһ��
        
        if (err < tol)
            break; 
        end
        
    end
end

function r = So(tau, X)
    % shrinkage operator
    r = sign(X) .* max(abs(X) - tau, 0);%����X����ֵ����ȡֵ
    %S������ʱ����ȡ���
end

function r = Do(tau, X)
    % shrinkage operator for singular values
    [U, S, V] = svd(X, 'econ');%economy size,
    %For m < n, only the first m columns of V arecomputed and S is m-by-m.
    r = U*So(tau, S)*V';
    %L������ʱ����ȡ���
end
