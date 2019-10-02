clc
clear all

data.MAX_CARS                   = 20;
data.MAX_MOVE_OF_CARS           = 5;
data.RENTAL_REQUEST_FIRST_LOC   = 3;
data.RENTAL_REQUEST_SECOND_LOC  = 4;
data.RETURN_FIRST_LOC           = 3;
data.RETURN_SECOND_LOC          = 2;
data.DISCOUNT                   = 0.9;
data.RENTAL_CREDIT              = 10;
data.MOVE_CAR_COST              = 2;
data.ACTIONS                    = zeros((data.MAX_MOVE_OF_CARS * 2)+1, 1);
data.POISSON_UPPER_BOUND        = 11;
% data.POISSON_CACHE = dic()


function poisson_cache = poisson_probability(n, lam)

end
function exp_return = expected_return(state, action, state_value, constant_returned_cars)

end