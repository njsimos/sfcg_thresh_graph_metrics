function [posFcg] = pos_fcg_func(fcg)
    % returns input matrix with all negative values set to zero
    posFcg = fcg;
    posFcg(fcg<0) = 0;
end

