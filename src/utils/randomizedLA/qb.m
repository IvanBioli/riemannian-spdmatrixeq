function [Q, B]=qb(A, num_queries, args)
% trace_est = qb(A, num_queries, dimension, args)
% 
% Computes the a QB decomposition.
% 
% Required Inputs:
% - A: Matrix
% 
% - num_queries: Total number of matrix-vector products to compute
% 
% Name-Value Pair Inputs:
% 
% - sketch_dist: Function_handle that, given two inputs (m,n), returns a m by n
% matrix that will be used in the sketching matrix. (default value: random sign matrix)
% 
    % Set default values
    arguments
        A;
        num_queries;
        args.sketch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
    end

    % Sketch the matrix, and take the QR
    S = args.sketch_dist(size(A,2), num_queries);
    [Q, ~] = qr(A*S, 0); % 0 here means 'use economic qr'
    B = Q'*A;

end
