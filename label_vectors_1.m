function [Tk,n_trials,min_length_trial,Ck,Fk,fk_start,fk_end,min_length_fix,CFk,min_length_cf] = label_vectors_1(s,EVENT,labels)
    % Extraction of the label vectors
    fixation_cross = labels(1);
    both_hands = labels(2);
    both_feet = labels(3);
    rest = labels(4);
    continuous_feedback = labels(5);
    
    % Tk: trial vector
    Tk = zeros(size(s,1),1);
    trial_start = EVENT.POS(EVENT.TYP == fixation_cross);
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
    
    length_trial = zeros(n_trials,1);
    for i = 1:n_trials
       length_trial(i) = sum(Tk==i);
    end
    clear i
    min_length_trial = min(length_trial);
    clear length_trial
    
    Ck = zeros(n_trials,1);
    cue_labels = [both_hands,both_feet,rest];
    idx_cue = find(ismember(EVENT.TYP, cue_labels));
    for i = 1:n_trials
        cue_pos = EVENT.POS(idx_cue(i));  
        trial_id = Tk(cue_pos);           
    
        if trial_id ~= 0
            Ck(trial_id) = EVENT.TYP(idx_cue(i));  
        end
    end
    clear i
    clear cue_pos
    clear trial_id
    
    % Fk vector
    Fk = zeros(size(s,1),1);
    fk_start = EVENT.POS(EVENT.TYP == fixation_cross);
    fk_end = EVENT.POS(EVENT.TYP == fixation_cross) + EVENT.DUR(EVENT.TYP == fixation_cross) -1;
    length_fix = zeros(n_trials,1);
    for i = 1:n_trials
        Fk(fk_start(i):fk_end(i)) = fixation_cross;
        length_fix(i) = fk_end(i)-fk_start(i);
    end
    min_length_fix = min(length_fix);
    clear i
    clear length_fix
    
    % CFk vector
    cf_code = continuous_feedback;
    CFk = zeros(size(s,1),1);
    cfk_start = EVENT.POS(EVENT.TYP == continuous_feedback);
    cfk_end = EVENT.POS(EVENT.TYP == continuous_feedback) + EVENT.DUR(EVENT.TYP == continuous_feedback) -1;
    length_cf = zeros(n_trials,1);
    for i = 1: length(cfk_start)
        CFk(cfk_start(i):cfk_end(i)) = cf_code;
        length_cf(i) = cfk_end(i)-cfk_start(i);
    end
    min_length_cf = min(length_cf);
    clear i
    clear cfk_start
    clear cfk_end

end

