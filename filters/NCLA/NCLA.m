% paper: A novel cooperative localization algorithm using enhanced particle filter technique in maritime search and rescue wireless sensor network.
% author: Yeong Jun Kim(colson)
% email: colson@korea.ac.kr || dud3722000@naver.com
% date: 2020-07-30
% Utilezed code: https://www.github.com/YeongJunKim/localization_sim
% This class for the linear system.

classdef NCLA < handle
    properties
        
        % x^t_n = A * [x^{t-1}_n + (v_a + v_d) delta(t)] + B * Q [u_{1,t} u_{2,t}]'
        
        % A : transition matrix
        % [ 1 1 0 0;
        %   0 1 0 0;
        %   0 0 1 1;
        %   0 0 0 1]
        % B : noise matrix
        % [ 0.5 0;
        %   1   0;
        %   0   0.5;
        %   0   1];
        % Q : variance of the noise.
        % u_{1,t}, u_{2,t} : transition noise.
        
        A; B; Q;            % mobility model
        H1; H2;             % combine H_ns & H_nm for simultaneous model.
        H_after;            % H_after(H1 - H2)
        
        H_ns; H_nm;         % measurement update
        
        P;                  % particle random point
        
        ns;                 % number of particles
        
        w;                  % weights
        
        particles;          % particles
        
        count = 1;          % data saving
        x_appended = [];    % filtered data
        is_init = "no";
        first_run = 0;
        
        % functions for nonlinear system.
        function_f;
        function_h_ns;
        function_h_nm;
        
        adjacency;          % acjacency matrix
        
        nm = 0;             % neighbor number
        n;                  % this agent number
    end
    methods
        
        function obj = NCLA(ns_, init_state_, A_, B_, Q_, P_, adjacency_, this_agent_num_)
            obj.A = A_;
            obj.B = B_;
            obj.Q = Q_;
            obj.P = P_;
            obj.ns = ns_;
            obj.n = this_agent_num_;
            obj.adjacency = adjacency_;
            obj.x_appended = zeros(size(init_state_, 1), 1);
            obj.x_appended(:,1) = init_state_(:);
            
            obj.particles = zeros(size(init_state_, 1), obj.ns);
            for i = 1:obj.ns
                for j = 1:size(init_state_, 1)
                    obj.particles(j,i) = obj.x_appended(j,1) + normrnd(0, (obj.P(j,j)));
                end
            end
            
            for j = 1:size(obj.adjacency,1)
                if(obj.n ~= j)
                    if(obj.adjacency(obj.n,j) ~= 0)
                        obj.nm = obj.nm + 1;
                    end
                end
            end
            
            % distance calculation matrices
            obj.H1 = zeros(obj.nm * size(init_state_,1), size(init_state_,1));
            obj.H2 = zeros(obj.nm * size(init_state_,1), obj.nm * size(init_state_,1));
            obj.H_after = zeros(obj.nm * 2 ,obj.nm * size(init_state_,1));
            for i = 1:obj.nm
                obj.H1(size(init_state_,1)*(i-1)+1:size(init_state_,1)*(i), 1:size(init_state_,1)) = [1 0 0 0;
                    0 0 0 0;
                    0 0 1 0;
                    0 0 0 0];
            end
            for i = 1:obj.nm
                obj.H2(size(init_state_,1)*(i-1)+1:size(init_state_,1)*(i), size(init_state_,1)*(i-1)+1:size(init_state_,1)*(i)) = [ 1 0 0 0;
                    0 0 0 0;
                    0 0 1 0;
                    0 0 0 0];
            end
            for i = 1:obj.nm
                obj.H_after(2*(i-1)+1:2*(i), 4*(i-1)+1:4*(i)) = [1 0 0 0;
                    0 0 1 0];
            end
            
            obj.count = 2;
            obj.is_init = "ok";
        end
        
        function r = get_distance_from_neighbors_states(obj, neighbors_)
            % arguments
            % neighbors_ (obj.nm * obj.nx) * 1 argumented vector.
            r = zeros(obj.nm, 1);
            h1 = obj.H1 * obj.x_appended(:,end)
            h2 = obj.H2 * neighbors_
            diff = h1 - h2;
            diff = obj.H_after * diff;
            for i = 1:obj.nm
                r(i,1) = norm(diff(2*(i-1)+1:2*(i)))
            end
        end
        
        function NCLA_importance_sampling(obj)
            if obj.is_init == "ok"
                % obj.ns 
                
            end
        end 
        
        function r = KLD_calculation(obj)
            
        end
        
        function r = importance_sampling(obj)
            
        end
        
        function r = residual_systematic_resampling(obj)
            
        end
    end
    
end

function r = KLD_calculation(obj)

end
function r = importance_sampling(obj)

end
function r = resudual_systematic_resampling(obj)

end

