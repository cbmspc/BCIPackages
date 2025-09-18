function timestamps = group_timestamps_into_bounds(timestamps, timebounds)
timestamps = timestamps(:);
for i = 1:size(timebounds,1)
    timestamps(timestamps >= timebounds(i,1) & timestamps <= timebounds(i,2),2) = i;
end
