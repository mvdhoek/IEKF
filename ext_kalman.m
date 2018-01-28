%% ITERATIVE EXTENDED KALMAN FILTER
% 
% Implementation of the iterative extended kalman filter
% system dynamics and jacobian matrices of f(x,u) and h(x,u) must be
% defined in the main script and are given as input to this function.
%
% The Jacobian matrix is calculated numerically by numjacobian.m. The IEKF
% is enabled by default (IEKF = 1) with a maximum of 100 iterations
% (max_iter).
%
% M.A. van den Hoek
% Delft University of Technology
% Faculty of Aerospace Engineering
% Department of Control & Simulation
%

function [ x_k1,P_k1,K_k1 ] = ext_kalman( sys,initial_param,Q_k,R_k,inputs,z_k1,dt )

    % iterative extended kalman
    IEKF = 1;
    
    % simulation parameters
    max_iter = 100;

    % calculate size of time vector and state vector
    N = size(z_k1,2);
    M = size(initial_param.x0,2);
    
    % preallocate matrices
    x_k1 = zeros(M,N); x_k1(:,1) = initial_param.xk0;
    P_k1 = zeros(M,M,N); P_k1(:,:,1) = diag(initial_param.P0); % can append with timeseries, squeeze to plot
    K_k1 = zeros(M,size(z_k1,1),N);

    for i=1:N-1
        % one-step ahead state estimation
        x_k1(:,i+1) = x_k1(:,i) + rk4(sys.f,x_k1(:,i),inputs.u(:,i),dt);
        
        % Calculate Jacobian for one-step ahead estimation
        %F = J.f(x_k1(i+1));
        %H = J.h(x_k1(i+1)); -- deprecated, done numerically by numjacobian
        
        % Calculate Jacobian numerically using complex step
        F = numjacobian(sys.f,x_k1(:,i+1)');
        H = numjacobian(sys.h,x_k1(:,i+1)');
        
        % Discretize state matrices
        [phi,gamma] = c2d(F,sys.G,dt);
        
        % calculate covariance matrix estimation
        P_k1(:,:,i+1) = phi*P_k1(:,:,i)*phi' + gamma*diag(Q_k)*gamma';
    
        if (IEKF)
            % ITERATIVE EXTENSION OF KALMAN (Iterated EKF)
            % define tolerance
            e = 1e-10;

            % initialize iterator using state estimation
            eta = x_k1(:,i+1);
            eta_l = zeros(size(eta));
            iter = 1;
            
            while (norm(eta_l - eta,Inf)/norm(eta_l,Inf) >= e) && (iter <= max_iter)

                % save old eta at begin of loop, except at first iretation
                if iter ~= 1
                    eta = eta_l;
                end

                % Calculate measurement equation Jacobian Hx
                H = numjacobian(sys.h,eta');

                % calculate Kalman gain
                K_k1(:,:,i+1) = P_k1(:,:,i+1)*H' / (H * P_k1(:,:,i+1) * H' + diag(R_k));

                % measurement update using estimated state, kalman gain and
                % measurement
                eta_l = x_k1(:,i+1) + K_k1(:,:,i+1)*(z_k1(:,i+1)-sys.h(eta)-H*(x_k1(:,i+1)-eta));

                %display(norm(eta_l - eta,Inf)/norm(eta_l,Inf));
                
                % increase iteration
                %display(iter)
                iter = iter + 1;
            end

            % update optimal state estimate
            x_k1(:,i+1) = eta_l;
           
        else   
            % calculate kalman gain -- NON-ITERATIVE
            K_k1(:,:,i+1) = P_k1(:,:,i+1)*H.'/(H * P_k1(:,:,i+1) * H.' + diag(R_k));

            % measurement update using estimated state, kalman gain and measurement
            x_k1(:,i+1) = x_k1(:,i+1) + K_k1(:,:,i+1)*(z_k1(:,i+1)-sys.h(x_k1(:,i+1)));
                                   %-- NON-ITERATIVE
        end
        
        % update covariance estimation
        P_k1(:,:,i+1) = (eye(M)-K_k1(:,:,i+1)*H)*P_k1(:,:,i+1)*(eye(M)-K_k1(:,:,i+1)*H)' + K_k1(:,:,i+1)*diag(R_k)*K_k1(:,:,i+1)';
       
        % print progress
        %fprintf('iteration: %d \n',i)
    
    end

end

