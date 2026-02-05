function s_filt = preprocessing_1(s, h, ord, freq_band, lap)
    % remove the last channel
    s = s(:,1:end-1);

    % laplacian filter
    s = s * lap;

    fs = h.SampleRate;
    %   a. Filter the signal in the μ and β bands
        %      Use a Butterworth filter (choose the order and cutoff frequencies)    
    [b_butter, a_butter] = butter(ord, freq_band/(fs/2), 'bandpass');
        %      Apply zero-phase filtering with MATLAB’s filtfilt function
    s_butter = filtfilt(b_butter, a_butter, s);

    %   b. Rectify the signal (square it)
    square_s_filt = s_butter.^2;
    %   c. Apply a moving average using a 1-second window
    window_size = fs*1; 
    b = (1/window_size)*ones(1,window_size);
    a = 1;
    s_filt = filtfilt(b, a, square_s_filt);
end

