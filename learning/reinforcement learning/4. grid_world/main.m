clc
clear all;


data.WORLD_SIZE = 4;
data.ACTION_SIZE = 4;
data.ACTIONS = [[0,-1];
                [-1,0];
                [0,1];
                [1,0]];
data.ACTION_PROB = 0.25;

[aa, asycn_iteration] = compute_state_value(1, 1, data);
[values, sync_iteration] = compute_state_value(0, 1, data);

subplot(1,2,1);
h1 = heatmap(aa);
subplot(1,2,2);
h2 = heatmap(values);

asycn_iteration
sync_iteration


function [new_state_values, iteration] = compute_state_value(in_place, discount, data)
    new_state_values = zeros(data.WORLD_SIZE, data.WORLD_SIZE);
    iteration = 0;
    while(true)
       if in_place == 1
           state_values = new_state_values;
       else
           state_values = new_state_values;
       end
       old_state_values = state_values;
       
       for i = 1:data.WORLD_SIZE
          for j = 1:data.WORLD_SIZE
              value = 0;
              for action = 1:data.ACTION_SIZE
                  [next_ij, reward] = step([i,j], action, data);
                  value = value + data.ACTION_PROB * (reward + discount * state_values(next_ij(1), next_ij(2)));
              end
              new_state_values(i,j) = value;
          end
       end
       
       new_state_values
       
       max_delta_value = max((abs(old_state_values - new_state_values)));
       iteration = iteration + 1;  
       if max_delta_value < 0.001
           break;
       end
    end
end

function [next_state, reward] = step(state, action, data)
if is_terminal(state, data) == 1
    next_state = state;
    reward = 0;
    return;
end
    
    next_state = state + data.ACTIONS(action,:);
    
    x = next_state(1);
    y = next_state(2);
    
    if x <= 0 || x > data.WORLD_SIZE || y <= 0 || y > data.WORLD_SIZE
        next_state = state;
    end
    
    reward = -1;
end

function re = is_terminal(state, data)


x = state(1);
y = state(2);

if x == 1 && y == 1
    re = 1;
    return
end
if x == data.WORLD_SIZE && y == data.WORLD_SIZE
    re = 1;
    return
end

if x == data.WORLD_SIZE && y == 1
    re = 1;
    return
end

if x == 1 && y == data.WORLD_SIZE
    re = 1;
    return
end

re = 0;
return
end