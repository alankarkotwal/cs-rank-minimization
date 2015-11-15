function z = solveSCA(A, y, K, t, eps, maxIter) %#ok<*INUSL>

    z = double(zeros(size(A, 2), 1));
    
    iter = 0;
    
    while norm(y - A*z) > eps && iter < maxIter
        
        obj = @(zv) norm(A*zv-y);
        cons = @(zv) deal((norm(z, 1) - (gFunc(z, t) - xi(z, t)*(zv - z)))/t - K, []);
        
        opts = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000, 'MaxFunEvals', 10000);

        [zv, ~, ~, ~] = fmincon(obj, z, [], [], [], [], [], [], cons, opts);        
        
        z = zv;
        iter = iter + 1;
        
        disp(norm(y - A*z));
        
    end
    
end

% Old CVX stuff
%         oldZ = z;
%         
%         cvx_begin quiet
%         
%             variable zv(size(A, 2))
%             minimize norm(y - A*zv)
%             subject to
%                 (norm(zv, 1) - (gFunc(z, t) - xi(z, t)*(zv - z)))/t <= K %#ok<*NOPRT>
%         
%         cvx_end
%         z = zv;
%         iter = iter + 1;