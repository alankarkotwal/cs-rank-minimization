function z = solveSCAe(A, y, K, t, eps, maxIter) %#ok<*INUSL>

%    oldZ = Inf*double(ones(size(A, 2), 1));
    z = double(zeros(size(A, 2), 1));
    
    iter = 0;
    
    while norm(y - A*z) > eps && iter < maxIter
        
        obj = @(zv) norm(A*zv-y);
        cons = @(zv) deal((norm(z, 1) - (gFunc(z, t) - xi(z, t)*(zv - z)))/t - K, []);
        
        %options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter', 'MaxIter', 10000);
        opts = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000, 'MaxFunEvals', 10000);

        [zv, fval, exitflag, output] = fmincon(obj, z, [], [], [], [], [], [], cons, opts);        
        
        z = zv;
        iter = iter + 1;
        
    end
    
    %disp(norm(y - A*z));
    %disp(l0norm(z, 0.1));
    
end

% CVX Optimization framework:

%         oldZ = z;
        
%         cvx_begin quiet
%         
%             variable zv(size(A, 2))
%             minimize norm(y - A*zv)
%             subject to
%                 (norm(z, 1) - (gFunc(z, t) - xi(z, t)*(zv - z)))/t <= K %#ok<*NOPRT>
%         
%         cvx_end
