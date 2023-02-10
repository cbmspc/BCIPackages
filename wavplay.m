function wavplay(signal, sampleRate)
% WAVPLAY  Plays an audio signal using the audioplayer object.
%   WAVPLAY(x, Fs) plays the audio signal x at the specified sample rate Fs.
    % Instantiate object to play audio data.
    player = audioplayer(signal, sampleRate);
    
    % Play signal from beginning to end (no overlap with other code possible).
    playblocking(player);
end