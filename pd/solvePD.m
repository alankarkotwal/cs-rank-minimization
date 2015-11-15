function z = solvePD(A,y)
% This function solves Az = y s.t. l0norm(z) <= r, where
% A is the data matrix (mx64).
% y is the measurement(observation) vector (mx1).
% z (64x1) is a "r"-sparse reconstruction (1 <= r <=64).

% Constants

r = ceil(64*0.3); % We are penalizing the reconstructed signal (z) to be no more than 30% sparse.
sigma = 1.1; % greater then 1.
rho = 1; % positive.
eps_outer = 15;
eps_inner = 0.05;
eps_perturb = 0.5;
z_feas = zeros(64,1); % one feasible solution.
w = zeros(64,1);
w(randi(64,1,r)) = 1; % Initial choice of w (one r-sparse vector with ones being the only non-0 values).
gamma = 1.1*max(0.5*(norm(y,2))^2, 0.5*(norm(A*((A'*A + rho*eye(64))\(A'*y + rho*w))-y,2))^2 + (rho/2)*(norm(((A'*A + rho*eye(64))\(A'*y + rho*w))-w,2))^2);
z = z_feas; % First z initialized to z_feas.

while(1) 
    % Block Coordinate Descent :-
    while(1) 
    
        z_prev = z;
        w_prev = w;
        z = (A'*A + rho*eye(64))\(A'*y + rho*w); % Using Gradient Descent
    
        w = zeros(64,1);
        [~,indices] = sort(z,'descend');
        w(indices(1:r)) = z(indices(1:r));
	
        % Inner stopping criterion
        if(max(max(z-z_prev)/max(max(z),1), max(w-w_prev)/max(max(w),1)) <= eps_inner)
            rtemp = sum(logical(w));
            % Perturbation to check if the stationary point is a minimum or a
            % saddle.
            if(rtemp>1)
                while(1)
                    
                    wtemp = zeros(64,1);
                    [~,idx] = sort(w,'descend');
                    wtemp(idx(1:(rtemp-1))) = w(indices(1:(rtemp-1)));
           
                    ztemp = (A'*A + rho*eye(64))\(A'*y + rho*wtemp); % Using Gradient Descent

                    wtemp = zeros(64,1);
                    [~,indices] = sort(ztemp,'descend');
                    wtemp(indices(1:r)) = ztemp(indices(1:r));
        
                    qtemp = (0.5*(norm(A*ztemp-y,2))^2 + (rho/2)*(norm(ztemp-wtemp,2))^2);
                    q = (0.5*(norm(A*z-y,2))^2 + (rho/2)*(norm(z-w,2))^2);
            
                    if((qtemp/q) <= eps_perturb)
                        z = ztemp;
                        w = wtemp;
                        rtemp = sum(logical(w));
                        if(rtemp>1)
                            continue;
                        else
                            break;
                        end
                    else
                        break;
                    end
                end
                break;
            else
                break;
            end
        end
    end
    
	rho = sigma*rho;
	if((0.5*(norm(A*((A'*A + rho*eye(64))\(A'*y + rho*w))-y,2))^2 + (rho/2)*(norm(((A'*A + rho*eye(64))\(A'*y + rho*w))-w,2))^2) > gamma)
		w = z_feas;
	end
	
	% Outer stopping criterion
	if(max(abs(z-w)) <= eps_outer)
		break;
	end
end

z = w;

end