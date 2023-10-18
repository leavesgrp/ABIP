function SGM = calculate_SGM(time_list, sh)

    SGM = exp(sum(log(max(1,time_list + sh))) / length(time_list)) - sh;
end