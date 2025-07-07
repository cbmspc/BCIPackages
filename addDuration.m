% Automatically convert seconds to the appropriate duration
function str = addDuration(seconds)
minutes = 0;
hours = 0;
days = 0;

if seconds >= 60
    minutes = floor(seconds / 60);
    seconds = mod(seconds, 60);
end

if minutes >= 60
    hours = floor(minutes / 60);
    minutes = mod(minutes, 60);
end

if hours >= 24
    days = floor(hours / 24);
    hours = mod(hours, 24);
end

if days > 0
    str = sprintf('%id %.1fh', days, (hours+minutes/60+seconds/3600));
elseif hours > 0
    str = sprintf('%ih %.1fm', hours, (minutes+seconds/60));
elseif minutes > 0
    str = sprintf('%im %.1fs', minutes, seconds);
else
    str = sprintf('%.1fs', seconds);
end

