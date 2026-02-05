function [Tk,n_trials,CFbK,Ck] = label_vectors(PSD,EVENT,cf_code,cue_labels)
    n_windows = size(PSD,1);
    % n_frequencies = size(PSD,2);
    % n_channels = size(PSD,3);

    % Tk: trial vector
    Tk = zeros(size(PSD,1),1);
    trial_start = EVENT.POS(EVENT.TYP == 771 | EVENT.TYP == 773);
    trial_end = EVENT.POS(EVENT.TYP == 781) + EVENT.DUR(EVENT.TYP == 781) -1;
    for i = 1: length(trial_start)
        Tk(trial_start(i):trial_end(i)) = i;
    end
    clear i
    clear trial_start
    clear trial_end
    
    trials = unique(Tk);
    trials = trials(trials~=0);
    n_trials = length(trials);
    
    % CFk vector
    CFbK = zeros(size(PSD,1),1);
    cfk_start = EVENT.POS(EVENT.TYP == 781);
    cfk_end = EVENT.POS(EVENT.TYP == 781) + EVENT.DUR(EVENT.TYP == 781) -1;
    for i = 1: length(cfk_start)
        CFbK(cfk_start(i):cfk_end(i)) = cf_code;
    end
    clear i
    clear cfk_start
    clear cfk_end
    
    % Ck vector
    Ck = zeros(n_windows,1);
    trial_start = EVENT.POS(EVENT.TYP == 771 | EVENT.TYP == 773 | EVENT.TYP == 783);
    trial_end = EVENT.POS(EVENT.TYP == 781) + EVENT.DUR(EVENT.TYP == 781) -1;
    idx_cue = find(ismember(EVENT.TYP, cue_labels));
    for i = 1:n_trials
        Ck(trial_start(i):trial_end(i)) = EVENT.TYP(idx_cue(i));
    end
    clear i
end

