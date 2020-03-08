classdef PF < handle
    properties
        % functions
        function_f;
        function_h;
        
        % particles random point
        P;
        
        % error matrix
        Q;
        R;
        
        % number of particles
        ns;
        % weights
        w;
        % particles
        particles;
        
        % data saving
        count = 1;
        x_appended;
        
        is_init = "no";
        first_run = 0;
        
        fig;
        
        resampling_strategy = 'multinomial_resampling';
        %         resampling_strategy = 'systematic_resampling';
    end
    methods
        %% PF init
        % parameter
        % obj               = PF class
        % ns_               = number of particles
        % init_states_      = init state
        % appended_num_     = size of storing data 
        % P_                = random factor of initial particles
        % Q_ ,R_            = noise covariance
        % function_f_       = system function
        % function_h_       = observation function
        % strategy_         = resampling method
        
        function r = PF_init(obj, ns_, init_state_, appended_num_, P_, Q_, R_, function_f_, function_h_, strategy_)
            
            % resampling strategy
            obj.resampling_strategy = strategy_;
            
            % particles random noise
            obj.P = P_;
            
            % matrix init
            obj.Q = Q_;
            obj.R = R_;
            
            %function init
            obj.function_f = function_f_;   % process
            obj.function_h = function_h_;   % observation
            
            obj.ns = ns_;
            
            disp("init saved datas");
            obj.x_appended = zeros(size(init_state_, 1), appended_num_);
            obj.x_appended(:,1) = init_state_;
            
            disp("init particles");
            obj.particles = zeros(size(init_state_, 1), obj.ns, appended_num_);
            for i = 1:obj.ns
                for j = 1:size(init_state_, 1)
                    obj.particles(j,i,1) = obj.x_appended(j,1) + normrnd(0, (P_(j,j)));
                end
            end
            
            disp("init weights");
            obj.w = zeros(obj.ns, appended_num_);
            obj.w(:,1) = repmat(1/obj.ns, obj.ns, 1);
            
            disp("init ok");
            obj.count = 1;
            obj.is_init = "ok";
            r = obj.is_init;
        end
        
        %% PF run
        % obj = PF class
        % u_  = input
        % z_  = measurement
        function r = PF_run(obj, u_, z_)
            if obj.is_init == "ok"
                obj.count = obj.count + 1;
                %% run
                xk = zeros(size(obj.particles, 1), size(obj.particles, 2));
                wk = zeros(obj.ns, 1);
                
                
                xkm1 = obj.particles(:,:,obj.count-1);
                wkm1 = obj.w(:,obj.count-1);
                
                for i = 1:obj.ns
                    arguments1 = num2cell([xkm1(:,i)' u_']);
                    xk(:, i) = obj.function_f(arguments1{:}) + diag(normrnd(0,obj.Q));
                    arguments2 = num2cell([xk(:,i)' u_']);
                    yk = obj.function_h(arguments2{:});
                    yk = yk + diag(normrnd(0, obj.R));
                    wk(i) = wkm1(i) * (normpdf((z_ - yk), 0, sqrt(1)))'*(normpdf((z_ - yk), 0, sqrt(1)));
                end
                wk = wk./sum(wk);
                Neff = 1/sum(wk.^2);
                percentage = 0.7;
                Nt = percentage * obj.ns;
                
                if(Neff < Nt)
                    fprintf("Resampling Neff = %f, Nt = %f \n", Neff, Nt);
                    [xk, wk] = resample(xk, wk, obj.resampling_strategy);
                end
                %% Compute estimated state
                xhk = zeros(size(xk,1), 1);
                for i = 1 : obj.ns
                    xhk = xhk + wk(i) * xk(:,i);
                end
                obj.w(:,obj.count) = wk;
                obj.particles(:,:,obj.count) = xk;
                
                obj.x_appended(:,obj.count) = xhk;
                r = xhk;
            else
                error("you must init class    : Call (PF_init(obj, ...)");
            end
        end
        function r = draw(obj)
            figure(1);
            for i = 1:obj.ns
                plot(obj.particles(1,i,obj.count), obj.particles(2,i,obj.count) ,"*");
                hold on;
            end
            grid on;
            grid minor;
            pbaspect([1 1 1])
            xlim([-2, 17]);
            ylim([-2, 17]);
            hold off;
        end
    end
end



%% Resampling function
function [xk, wk, idx] = resample(xk, wk, resampling_strategy)

Ns = length(wk);  % Ns = number of particles

% wk = wk./sum(wk); % normalize weight vector (already done)

switch resampling_strategy
    case 'multinomial_resampling'
        with_replacement = true;
        idx = randsample(1:Ns, Ns, with_replacement, wk);
        %{
      THIS IS EQUIVALENT TO:
      edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
      edges(end) = 1;                 % get the upper edge exact
      % this works like the inverse of the empirical distribution and returns
      % the interval where the sample is to be found
      [~, idx] = histc(sort(rand(Ns,1)), edges);
        %}
    case 'systematic_resampling'
        % this is performing latin hypercube sampling on wk
        edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
        edges(end) = 1;                 % get the upper edge exact
        u1 = rand/Ns;
        % this works like the inverse of the empirical distribution and returns
        % the interval where the sample is to be found
        [~, idx] = histc(u1:1/Ns:1, edges);
        % case 'regularized_pf'      TO BE IMPLEMENTED
        % case 'stratified_sampling' TO BE IMPLEMENTED
        % case 'residual_sampling'   TO BE IMPLEMENTED
    otherwise
        error('Resampling strategy not implemented')
end

xk = xk(:,idx);                    % extract new particles
wk = repmat(1/Ns, 1, Ns);          % now all particles have the same weight

end

