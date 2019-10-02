clc
clear all
%% Value Iteration
% gamblers_problem

data.GOAL = 100;
data.STATE = zeros(data.GOAL+1, 1);
data.HEAD_PROB = 0.4;
data.BACK_PROB = 0.6;

data

% state value
state_value = zeros(data.GOAL+1, 1); % adding state 0
state_value(data.GOAL+1) = 1.0; % reward = 1

sweeps_history = zeros(size(state_value));

% value iteration
scnt = 1;
cnt = 1;
while(true)
   old_state_value = state_value;
   sweeps_history(:, cnt) = old_state_value;
   cnt = cnt + 1;
   for state = 1:data.GOAL
%       get_possible actions for current state
        actions = min(state, data.GOAL - state); % possible action number
        action_returns = zeros(actions + 1, 1);
        for action = 1 : actions
            if action == 0
                break;
            end
            p_index = state + action + 1;
            m_index = state - action + 1;
            action_returns(action) = data.HEAD_PROB * state_value(p_index, 1) + (1 - data.HEAD_PROB) * state_value(m_index, 1);
        end
        new_value = max(action_returns);
        state_value(state + 1, 1) = new_value;
        state_value(data.GOAL+1, 1) = 1.0;
   end
   delta = max(abs(state_value - old_state_value));
   if delta < 0.0000000000001
        sweeps_history(:, cnt) = state_value;
        cnt = cnt + 1;
       break;
   end
end

subplot(1,2,1);
cnt
for index = 1 : cnt-1
    plot(sweeps_history(:,index));
    hold on;
end

% compute the optimal policy
policy = zeros(data.GOAL + 1, 1);
for state = 1:data.GOAL
   actions = min(state, data.GOAL - state); % possible action number
   action_returns = zeros(actions + 1, 1);
   for action = 1:actions
       p_index = state + action + 1;
       m_index = state - action + 1;
        action_returns(action) = data.HEAD_PROB * state_value(p_index, 1) + (1 - data.HEAD_PROB) * state_value(m_index, 1);
   end
   action_returns;
   rnd = round(action_returns, 5);
   [value, argmax] = max(rnd);
   policy(state + 1, 1) = argmax;
end


subplot(1,2,2);
plot(policy);
% plot(policy);

