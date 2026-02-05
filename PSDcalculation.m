function [selected_PSD,selected_f,h] = PSDcalculation(s_lap,h,wlength,pshift,wshift,mlength,subset_f,winconv)
    % PSD over time
    samplerate = h.SampleRate;
    [PSD, f] = proc_spectrogram(s_lap, wlength, wshift, pshift, samplerate, mlength);

    % Select only a subset of frequency from the PSD
    idx_f_subset = find(f >= subset_f(1) & f <= subset_f(2));
    selected_f = f(idx_f_subset);
    selected_PSD = PSD(:,idx_f_subset,:);

    % Recompute the h.EVENT.POS and .DUR with respect to the PSD windows
    h.EVENT.POS = proc_pos2win(h.EVENT.POS, wshift*samplerate, winconv, wlength*samplerate);
    h.EVENT.DUR = zeros(length(h.EVENT.POS),1);
    for i = 2:length(h.EVENT.POS)
        pre = h.EVENT.POS(i-1);
        current = h.EVENT.POS(i);
        h.EVENT.DUR(i-1) = current-pre;
    end
    h.EVENT.DUR(end) = size(selected_PSD,1) - h.EVENT.POS(end);
    clear pre
    clear current
end

