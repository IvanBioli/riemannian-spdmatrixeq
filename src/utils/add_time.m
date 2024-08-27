function info = add_time(info,time_add)
% add_time - Adds time_add to the time field of the struct array info

% Compute the total time
aux = num2cell([info.time] + time_add);
[info.time] = aux{:};
end

