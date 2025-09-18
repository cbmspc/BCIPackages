% Example usages:
% EventTimeStamps = make_eventtimestamps_from_timestamps_and_tags(t_RightKnee_MinVel, 'Rts', t_LeftKnee_MinVel, 'Lts', t_Right_ToeOff_TerminalSwing_HeelStrike(:,1), 'Rto', t_Left_ToeOff_TerminalSwing_HeelStrike(:,1), 'Lto');
% EventTimeStamps = make_eventtimestamps_from_timestamps_and_tags(t_Right_ToeOff_TerminalSwing_HeelStrike(:,[1 3]), 'Rsw', t_Left_ToeOff_TerminalSwing_HeelStrike(:,[1 3]), 'Lsw');
%
function EventTimeStamps = make_eventtimestamps_from_timestamps_and_tags(varargin)
EventTimeStamps = {};
for i = 1:2:nargin
    for j = 1:size(varargin{i},1)
        if ~any(isnan(varargin{i}(j,:)))
            EventTimeStamps = [EventTimeStamps; {varargin{i}(j,:)} varargin(i+1)];
        end
    end
end
[~,in] = sort(cellfun(@mean,EventTimeStamps(:,1)));
EventTimeStamps = EventTimeStamps(in,:);
