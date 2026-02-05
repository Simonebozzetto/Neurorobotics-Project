close all
clear all
clc
%% 
n_subjects = 8;
cf_code = 781;
cue_labels = [773, 771];
channels = {'FZ','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPZ','CP2','CP4'};
load chanlocs16.mat;
%% calibration files 
for subj = 1:n_subjects
%subj = 2;
    %% load calibration files
    filename = sprintf('subject_0%d_calibration.mat', subj);
    load(filename);
    
    n_cal_runs = size(cal_s,2);

    Fisher_runs = cell(1,n_cal_runs);
    for i = 1 : n_cal_runs
        [~,~,CFbK,Ck] = label_vectors(cal_PSD{i}, cal_h{i}.EVENT, cf_code, cue_labels);
        [Fisher_runs{i},~] = fisher(cal_PSD{i}, Ck, CFbK, cue_labels);
    end
    
    %% visualization of the calibration runs
    figure;
    for run = 1:n_cal_runs
        subplot(1,n_cal_runs,run)
        imagesc(Fisher_runs{run}')
        xlabel('Frequency [Hz]')
        ylabel('Channel')
        title(sprintf('Calibration run %d', run))
        xticks(1:length(f))
        xticklabels(string(f))
        yticks(1:length(channels))
        yticklabels(channels)
    end
    sgtitle(sprintf('Fisher Score - Subject %d', subj))
    
    %% concatenation of PSDs
    PSD = [];
    for i = 1:n_cal_runs
        PSD = [PSD; cal_PSD{i}];
    end
    
    EVENT.TYP = [];
    EVENT.DUR = [];
    EVENT.POS = [];
    offset = 0;
    for i = 1:n_cal_runs
        EVENT.TYP = [EVENT.TYP; cal_h{i}.EVENT.TYP];
        EVENT.DUR = [EVENT.DUR; cal_h{i}.EVENT.DUR];
        EVENT.POS = [EVENT.POS; cal_h{i}.EVENT.POS + offset];
        offset = offset + size(cal_PSD{i},1);
    end
    
    [Tk,n_trials,CFbK,Ck] = label_vectors(PSD, EVENT, cf_code, cue_labels);
    [Fisher,~] = fisher(PSD, Ck, CFbK, cue_labels);
    
    %% visualization of the concatenated runs
    figure;
    imagesc(Fisher')
    xlabel('Frequency [Hz]')
    ylabel('Channel')
    title(sprintf('Fisher score of concatenated runs - Subject %d', subj))
    xticks(1:length(f))
    xticklabels(string(f))
    yticks(1:length(channels))
    yticklabels(channels)
    
    %% extraction of the features by selecting them manually
    disp('Insert a Nx2 matrix: each row is [frequency channel]');
    pairs = input('Pairs [freq chan]: ');
    
    if size(pairs,2) ~= 2
        error('pairs must be an Nx2 matrix: [freq chan].');
    end
    
    freq_vals = pairs(:,1);
    chan_vals = pairs(:,2);
    
    % Map frequencies to indices in f
    idx_f = zeros(size(freq_vals));
    for k = 1:length(freq_vals)
        [~, ii] = min(abs(f - freq_vals(k)));
        idx_f(k) = ii;
    end
    
    Ns = size(PSD,1);
    Np = size(pairs,1);
    
    PSD_selected = zeros(Ns, Np, 1);
    for k = 1:Np
        PSD_selected(:,k,1) = PSD(:, idx_f(k), chan_vals(k));
    end
    
    %% Train discriminant model
    LabelIdx = CFbK == cf_code;
    model = fitcdiscr(PSD_selected(LabelIdx, :), Ck(LabelIdx), 'DiscrimType','quadratic');
    
    %% Accuracy sul training set
    [Gk, pp] = predict(model, PSD_selected(LabelIdx, :));
    overall_accuracy = sum(Gk == Ck(LabelIdx))/length(Gk);
    accuracy_771 = mean(Gk(Ck(LabelIdx) == 771) == 771);
    accuracy_773 = mean(Gk(Ck(LabelIdx) == 773) == 773);
    
    figure;
    bar([overall_accuracy,accuracy_771,accuracy_773])
    set(gca, 'XTickLabel', {'overall', '771', '773'});
    ylabel('Accuracy');
    title(sprintf('Training accuracy - Subject %d', subj))
    
    %% Save the model
    save(sprintf('sub_0%d_model.mat', subj), 'model', 'pairs');

    %% Single sample accuracy and trial accuracy
    if n_cal_runs >= 2 % due to the fact that the second subject has only one calibration file
        
        % single sample accuracy
        run_id = [];
        for i = 1:n_cal_runs
            run_id = [run_id; i*ones(size(cal_PSD{i},1),1)];
        end
        % selection of only the cue labels
        cue = (Ck == 771) | (Ck == 773);
        idx_use = (CFbK == cf_code) & cue;    
    
        X = PSD_selected(idx_use, :);
        Y = Ck(idx_use);
        T = Tk(idx_use);
        R = run_id(idx_use);
    
        all_pred = zeros(size(Y));
        all_pp   = zeros(size(Y,1), 2);
        
        for i = 1:n_cal_runs
            train_runs = (R ~= i);
            test_run = (R == i);
        
            % Train on (n_cal_runs-1) runs
            model_cv = fitcdiscr(X(train_runs,:), Y(train_runs), 'DiscrimType','quadratic');
        
            % Predict on the left-out run
            [g_cv, pp_cv] = predict(model_cv, X(test_run,:));
        
            all_pred(test_run) = g_cv;
            all_pp(test_run,:) = pp_cv;
        end
    
        offline_sample_acc      = mean(all_pred == Y);
        offline_sample_acc_771  = mean(all_pred(Y==771) == 771);
        offline_sample_acc_773  = mean(all_pred(Y==773) == 773);
    
        figure;
        bar([offline_sample_acc,offline_sample_acc_771,offline_sample_acc_773])
        set(gca, 'XTickLabel', {'overall', '771', '773'});
        ylabel('Accuracy');
        title(sprintf('Offline single sample accuracy - Subject %d', subj))
        
        % trial accuracy
        trial_correct = zeros(n_trials,1);
        trial_true = zeros(n_trials,1);
        trial_hat  = zeros(n_trials,1);
        
        for k = 1 : n_trials
            idx_k = (T == k);

            y_true = mode(Y(idx_k));
        
            y_hat = mode(all_pred(idx_k));
        
            trial_true(k) = y_true;
            trial_hat(k)  = y_hat;
            trial_correct(k) = (y_hat == y_true);
        end
        
        offline_trial_acc = mean(trial_correct);
        offline_trial_acc_771 = mean(trial_correct(trial_true == 771));
        offline_trial_acc_773 = mean(trial_correct(trial_true == 773));

        figure;
        bar([offline_trial_acc,offline_trial_acc_771,offline_trial_acc_773])
        set(gca, 'XTickLabel', {'overall', '771', '773'});
        ylabel('Accuracy');
        title(sprintf('Offline trial accuracy - Subject %d', subj))
    end
    %% save offline matrics 
    if n_cal_runs >= 2
        save(sprintf('sub_0%d_offline_metrics.mat',subj),'offline_sample_acc','offline_sample_acc_771','offline_sample_acc_773','offline_trial_acc','offline_trial_acc_771','offline_trial_acc_773')
    end
end