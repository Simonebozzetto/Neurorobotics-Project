function s_trials = conc_to_trial(s,n_samples,n_channels,n_trials,Tk)
    % Creation of the 3D matrices [n_samples,n_channels,n_trials] starting 
    % from the concatenation signal s and the vector Tk
    s_trials = zeros(n_samples,n_channels,n_trials);
    min_length_trial = size(s_trials,1);
    for i = 1:n_trials
        sig = s(Tk==i,:);
        if size(sig,1) ~= min_length_trial
            sig = sig(1:min_length_trial,:);
        end
        s_trials(:,:,i) = sig;
    end
    clear i
    clear sig
end

