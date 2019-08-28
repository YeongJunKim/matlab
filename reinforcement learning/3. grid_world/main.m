
clc
clear all

WORLD_SIZE = 5;
A_POS = [0, 1];
A_PRIME_POS = [5, 1];
B_POS = [0, 3];
B_PRIME_POS = [3, 3];

DISCOUNT = 0.9;


ACTION_NUM      = 4;
ACTION_LEFT     = [0, -1];
ACTION_UP       = [-1, 0];
ACTION_RIGHT    = [0, 1];
ACTION_DOWN     = [1, 0];

ACTION_P        = 0.25;

value = zeros(WORLD_SIZE, WORLD_SIZE);
while(1)
    new_value = zeros(WORLD_SIZE, WORLD_SIZE);
    for i = 1:WORLD_SIZE
        for j=1:WORLD_SIZE
           for  k=1:ACTION_NUM
                [next_ij, reward] = step([i,j], ACTION_NUM);
% %                 belman ford equation
                new_value(i,j) = new_value(i,j) + ACTION_P * (reward + DISCOUNT * value(next_ij(1),next_ij(2)));
           end
        end
    end
    sub = sum(abs(new_value - value), 'all');
    if sub < 0.1 
        break;
    end
    value = new_value;
end
result_value = new_value;

subplot(1,2,1);
heatmap(result_value);

value = zeros(WORLD_SIZE, WORLD_SIZE);
while(1)
    new_value = zeros(WORLD_SIZE, WORLD_SIZE);
    for i = 1:WORLD_SIZE
        for j = 1:WORLD_SIZE
            values = zeros(1);
            size(values);
            for k = 1:ACTION_NUM    
                [next_ij, reward] = step([i,j], ACTION_NUM);
    % %             value iteration
                values(k) = reward + DISCOUNT * value(next_ij(1), next_ij(2));
            end 
            new_value(i,j) = max(values);
        end
    end
    sub = sum(abs(new_value - value), 'all');
    if sub < 0.001
        break;
    end
    value = new_value
end
result_value = new_value;

subplot(1,2,2);
heatmap(result_value);
% 
% cdata = [45 60 32; 43 54 76; 32 94 68; 23 95 58];
% h = heatmap(cdata);


function [next_state, reward] = step(state, action)
WORLD_SIZE = 5;
A_POS = [1, 2];
A_PRIME_POS = [1, 5];
B_POS = [1, 4];
B_PRIME_POS = [3, 4];
ACTION_LEFT     = [0, -1];
ACTION_UP       = [-1, 0];
ACTION_RIGHT    = [0, 1];
ACTION_DOWN     = [1, 0];

    if state(1) == A_POS(1) && state(2) == A_POS(2)
        next_state = A_PRIME_POS;
        reward = 10;
        return 
    end
    if state(1) == B_POS(1) && state(2) == B_POS(2)
        next_state = B_PRIME_POS;
        reward = 5;
        return
    end
    
    if action == 1
        next_state = state + ACTION_LEFT;
    end
    if action == 2
        next_state = state + ACTION_UP;
    end
    if action == 3
        next_state = state + ACTION_RIGHT;
    end
    if action == 4
        next_state = state + ACTION_DOWN;
    end
    x = next_state(1);
    y = next_state(2);
    if x < 0 || x >= WORLD_SIZE || y < 0 || y >= WORLD_SIZE
        reward = -1.0;
        next_state = state;
    else
        reward = 0;
    end
end
