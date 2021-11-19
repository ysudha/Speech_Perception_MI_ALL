function [L, S] = RobustPCA_time(X, lambda, mu, tol, max_iter)
    % - X is a data matrix (of the size N x M x ntime) to be decomposed
    %   X can also contain NaN's for unobserved values
    % - lambda - regularization parameter, default = 1/sqrt(max(N,M))
    % - mu - the augmented lagrangian parameter, default = 10*lambda
    % - tol - reconstruction error tolerance, default = 1e-6
    % - max_iter - maximum number of iterations, default = 1000

    [M, N, nt] = size(X); %nt = num of time windows 
    unobserved = isnan(X);
    X(unobserved) = 0;
    normX=0;
    for i = 1:nt
        normX = normX + norm(squeeze(X(:,:,i)), 'fro');
    end
    normX = normX/nt; %average fro norm of X

    % default arguments
    if nargin < 2
        lambda =  1 / sqrt(max(M,N));
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
    S = zeros(M, N, nt);    %looking for multiple sparse matrices, single L
    Y = zeros(M, N, nt);    % lagrange multipliers
    
    for iter = (1:max_iter)
        % ADMM step: update L and S
        %updating L using each S one at a time, average them all.
        for i = 1:nt
            L = L + Do(1/mu, squeeze(X(:,:,i)) - squeeze(S(:,:,i)) + (1/mu)*squeeze(Y(:,:,i)));
        end
        L = L./nt;
        
        %update each S using the single L
        for i = 1:nt
            S(:,:,i) = So(lambda/mu, squeeze(X(:,:,i)) - L + (1/mu)*squeeze(Y(:,:,i)));
        end
        % and augmented lagrangian multiplier
        Z=zeros(M,N,nt);
        for i = 1:nt
            Z(:,:,i) = squeeze(Z(:,:,i)) + (squeeze(X(:,:,i) - L - squeeze(S(:,:,i))));
        end
        Z = Z./nt;
        
        Z(unobserved) = 0; % skip missing values
        Y = Y + mu*Z;
        
        normZ=0;
        for i = 1:nt
            normZ = normZ + norm(squeeze(Z(:,:,i)), 'fro');
        end
        normZ = normZ/nt; %average fro norm of Z
    
        err = normZ/ normX;
  %      if (iter == 1) || (mod(iter, 10) == 0) || (err < tol)
 %           fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
  %                  iter, err, rank(L), nnz(S(~unobserved)));
  %      end
        if (err < tol) break; end
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