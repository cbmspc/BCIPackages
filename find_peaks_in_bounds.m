function pklocs_sec = find_peaks_in_bounds(x, Fs, timebounds, MinPeakHeight, MinPeakProminence, MinPeakDistance)
if nargout > 0
    [~,pklocs] = findpeaks(x, 'MinPeakHeight', MinPeakHeight, 'MinPeakProminence', MinPeakProminence, 'MinPeakDistance', MinPeakDistance);
else
    findpeaks(x, 'MinPeakHeight', MinPeakHeight, 'MinPeakProminence', MinPeakProminence, 'MinPeakDistance', MinPeakDistance, 'Annotate', 'extents');
    keyboard
end
pklocs_sec = (pklocs-1)/Fs;
pklocs_sec = group_timestamps_into_bounds(pklocs_sec, timebounds);
pklocs_sec = pklocs_sec(pklocs_sec(:,2)>0,1);
