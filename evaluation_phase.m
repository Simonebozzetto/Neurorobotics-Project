close all
clear all
clc
%% 
n_subjects = 8;
cf_code = 781;
cue_labels = [773,771,783]; % here in the evaluation part we have also the rest coded as 783
% but after an inspection there aren't any
channels = {'FZ','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPZ','CP2','CP4'};
load chanlocs16.mat;
%% evaluation files 
for subj = 1:n_subjects
%subj = 7;
    %% load evaluation files for each subject and concatenate them
    evaluation_file = sprintf('subject_0%d_evaluation.mat', subj);
    load(evaluation_file);
    
    model_file = sprintf('sub_0%d_model.mat', subj);
    load(model_file);

    n_eva_runs = size(eva_s,2);

    PSD = [];
    for i = 1:n_eva_runs
        PSD = [PSD; eva_PSD{i}];
    end
    
    EVENT.TYP = [];
    EVENT.DUR = [];
    EVENT.POS = [];
    offset = 0;
    for i = 1:n_eva_runs
        EVENT.TYP = [EVENT.TYP; eva_h{i}.EVENT.TYP];
        EVENT.DUR = [EVENT.DUR; eva_h{i}.EVENT.DUR];
        EVENT.POS = [EVENT.POS; eva_h{i}.EVENT.POS + offset];
        offset = offset + size(eva_PSD{i},1);
    end
    
    [Tk,n_trials,CFbK,Ck] = label_vectors(PSD, EVENT, cf_code, cue_labels);
    
    %% feature extraction based on the selected ones in the calibration part
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
    
    features = zeros(Ns, Np, 1);
    for k = 1:Np
        features(:,k,1) = PSD(:, idx_f(k), chan_vals(k));
    end

    %features = PSD(:, indices_sel_freq, selected_channels);

    %% evaluate the model created during the calibration phase
    [Gk, pp] = predict(model, features(:, :));

    mi_idx = ismember(Ck, [771 773]);

    overall_accuracy = mean(Gk(mi_idx) == Ck(mi_idx));
    accuracy_773 = sum(Gk(Ck==773) == Ck(Ck==773)) / sum(Ck == 773);
    accuracy_771 = sum(Gk(Ck==771) == Ck(Ck==771)) / sum(Ck == 771);
    
    figure
    bar([overall_accuracy,accuracy_771,accuracy_773])
    set(gca, 'XTickLabel', {'overall', '771', '773'});
    ylabel('Accuracy');
    title('Single sample accuracy on test set');

    %% Exponential accumulation framework on the posterior probabilities 
    n_windows = size(Gk,1);

    classNames = model.ClassNames; 
    i773 = find(classNames == 773, 1);
    i771 = find(classNames == 771, 1);

    pp_both_hands = pp(:,i773);
    pp_both_feet = pp(:,i771);
    D_both_hands = 0.5 * ones(n_windows,1);
    D_both_feet = 0.5 * ones(n_windows,1);
    
    alpha = 0.1;
    fixed_thresholds = [0.2,0.8];
    
    start_trial = EVENT.POS(EVENT.TYP == 781);
    for i = 2:n_windows
        if ismember(i,start_trial)
            D_both_hands(i) = 0.5;
        else
            D_both_hands(i) = alpha * D_both_hands(i-1) + (1-alpha) * pp_both_hands(i);
        end
    end
    
    for i = 2:n_windows
        if ismember(i,start_trial)
            D_both_feet(i) = 0.5;
        else
            D_both_feet(i)  = alpha * D_both_feet(i-1)  + (1-alpha) * pp_both_feet(i);
        end
    end
    %%
    n_plot_trials = 1; % chaning the number we will have more plots  
    for i = 1:n_plot_trials
        id_trial = i;
        sample_trial = Tk==i;
        
        D_trial_h = D_both_hands(sample_trial);
        pp_trial_h = pp_both_hands(sample_trial);

        D_trial_f = D_both_feet(sample_trial);
        pp_trial_f = pp_both_feet(sample_trial);

        start_D_trial_h = find(D_trial_h == 0.5);
        start_D_trial_f = find(D_trial_f == 0.5);
        
        D_trial_h(1:start_D_trial_h-1) = [];
        pp_trial_h(1:start_D_trial_h-1) = [];
        
        D_trial_f(1:start_D_trial_f-1) = [];
        pp_trial_f(1:start_D_trial_f-1) = [];

        samples_h = 1:length(D_trial_h);
        samples_f = 1:length(D_trial_f);

        figure
        subplot(121)
        plot(samples_h,D_trial_h,'k-','LineWidth',1)
        hold on
        plot(samples_h,pp_trial_h,'ok')
        plot(samples_h,0.5*ones(length(samples_h),1),'k--')
        plot(samples_h,fixed_thresholds(1)*ones(length(samples_h),1),'k--')
        plot(samples_h,fixed_thresholds(2)*ones(length(samples_h),1),'k--')
        plot(samples_h,zeros(length(samples_h),1),'r-','LineWidth',3)
        plot(samples_h,ones(length(samples_h),1),'r-','LineWidth',3)
        ylim([0,1])
        xlim([1,length(samples_h)])
        xlabel('sample')
        ylabel('probability/control')
        title(sprintf('Trial %d - Class Both Hands - alpha = 0.1',i))
        legend('integrated prob','raw prob','lower threshold','upper threshold')
        subplot(122)
        plot(samples_f,D_trial_f,'k-','LineWidth',1)
        hold on
        plot(samples_f,pp_trial_f,'ok')
        plot(samples_f,0.5*ones(length(samples_f),1),'k--')
        plot(samples_f,fixed_thresholds(1)*ones(length(samples_f),1),'k--')
        plot(samples_f,fixed_thresholds(2)*ones(length(samples_f),1),'k--')
        plot(samples_f,zeros(length(samples_f),1),'r-','LineWidth',3)
        plot(samples_f,ones(length(samples_f),1),'r-','LineWidth',3)
        ylim([0,1])
        xlim([1,length(samples_h)])
        xlabel('sample')
        ylabel('probability/control')
        title(sprintf('Trial %d - Class Both Feet - alpha = 0.1',i))
        legend('integrated prob','raw prob','lower threshold','upper threshold')
    end
    %% Dynamical control to properly define the system response
    C_both_hands = 0.5 * ones(n_windows,1);
    C_both_feet = 0.5 * ones(n_windows,1);
    
    alpha_control = 0.1; % control integration gain
    
    for i = 2:n_windows
        if ismember(i,start_trial)
            C_both_hands(i) = 0.5;
            C_both_feet(i)  = 0.5;
        else
            % apply dynamic control to the accumulated posterior
            C_both_hands(i) = C_both_hands(i-1) + alpha_control * (D_both_hands(i) - C_both_hands(i-1));
            C_both_feet(i)  = C_both_feet(i-1)  + alpha_control * (D_both_feet(i)  - C_both_feet(i-1));
        end
    end
    %%
    n_plot_trials = 1; % number of trials to plot
    
    for i = 1:n_plot_trials
        id_trial = i;
        sample_trial = Tk==i;
        
        % extract dynamic control and raw posterior for this trial
        C_trial_h = C_both_hands(sample_trial);
        pp_trial_h = pp_both_hands(sample_trial);
    
        C_trial_f = C_both_feet(sample_trial);
        pp_trial_f = pp_both_feet(sample_trial);
    
        % find trial start (where control resets to 0.5)
        start_C_trial_h = find(C_trial_h == 0.5, 1, 'first');
        start_C_trial_f = find(C_trial_f == 0.5, 1, 'first');
        
        % trim data before trial start
        C_trial_h(1:start_C_trial_h-1) = [];
        pp_trial_h(1:start_C_trial_h-1) = [];
        
        C_trial_f(1:start_C_trial_f-1) = [];
        pp_trial_f(1:start_C_trial_f-1) = [];
    
        samples_h = 1:length(C_trial_h);
        samples_f = 1:length(C_trial_f);
    
        figure
        % Plot Both Hands
        subplot(121)
        plot(samples_h,C_trial_h,'k-','LineWidth',1)
        hold on
        plot(samples_h,pp_trial_h,'ok')
        plot(samples_h,0.5*ones(length(samples_h),1),'k--')
        plot(samples_h,fixed_thresholds(1)*ones(length(samples_h),1),'k--')
        plot(samples_h,fixed_thresholds(2)*ones(length(samples_h),1),'k--')
        plot(samples_h,zeros(length(samples_h),1),'r-','LineWidth',3)
        plot(samples_h,ones(length(samples_h),1),'r-','LineWidth',3)
        ylim([0,1])
        xlim([1,length(samples_h)])
        xlabel('sample')
        ylabel('control / probability')
        title(sprintf('Trial %d - Class Both Hands (Dynamic Control)', i))
        legend('control','raw prob','baseline','lower threshold','upper threshold')
    
        % Plot Both Feet
        subplot(122)
        plot(samples_f,C_trial_f,'k-','LineWidth',1)
        hold on
        plot(samples_f,pp_trial_f,'ok')
        plot(samples_f,0.5*ones(length(samples_f),1),'k--')
        plot(samples_f,fixed_thresholds(1)*ones(length(samples_f),1),'k--')
        plot(samples_f,fixed_thresholds(2)*ones(length(samples_f),1),'k--')
        plot(samples_f,zeros(length(samples_f),1),'r-','LineWidth',3)
        plot(samples_f,ones(length(samples_f),1),'r-','LineWidth',3)
        ylim([0,1])
        xlim([1,length(samples_f)])
        xlabel('sample')
        ylabel('control / probability')
        title(sprintf('Trial %d - Class Both Feet (Dynamic Control)', i))
        legend('control','raw prob','baseline','lower threshold','upper threshold')
    end
    % now we have more regular choice of the movement

    %% Online trial accuracy
    upper_th = fixed_thresholds(2);

    cmd_ctrl = zeros(size(Ck));
    cmd_ctrl(C_both_hands >= upper_th) = 773;
    cmd_ctrl(C_both_feet  >= upper_th) = 771;


    wshift = 0.0625; % to perform the conversion from windows to seconds
                     % it's the same used for the preparation part

    trial_true = zeros(n_trials,1);
    trial_decision = zeros(n_trials,1);
    trial_correct = zeros(n_trials,1);
    for i = 1:n_trials
        idx_trial = Tk == i;

        y_true = mode(Ck(idx_trial));
        trial_true(i) = y_true;

        k_dec = find(cmd_ctrl(idx_trial) ~= 0, 1, 'first');

        if isempty(k_dec)
            trial_decision(i) = 0; % no decision
            trial_correct(i) = false;
        else
            cmd_vals = cmd_ctrl(idx_trial);
            y_hat = cmd_vals(k_dec);
            trial_decision(i) = y_hat;
            trial_correct(i) = (y_hat == y_true);
        end
    end
    % drial accuracy only for the trials that decided
    decision = (trial_decision ~= 0);
    online_trial_acc = mean(trial_correct(decision));
    
    % decision rate
    decision_rate = mean(decision);
    
    % trial accuracy for classes
    online_trial_acc_771 = mean(trial_correct(decision & trial_true==771));
    online_trial_acc_773 = mean(trial_correct(decision & trial_true==773));

    figure
    bar([online_trial_acc,online_trial_acc_771,online_trial_acc_773])
    set(gca, 'XTickLabel', {'overall', '771', '773'});
    ylabel('Accuracy');
    title('Online trial accuracy on test set');
    %% Average time to deliver a command
    time_to_cmd = zeros(n_trials,1);
    
    for i = 1:n_trials
        idx_trial = (Tk == i);
    
        k_dec = find(cmd_ctrl(idx_trial) ~= 0, 1, 'first');
    
        if ~isempty(k_dec)
            time_to_cmd(i) = (k_dec - 1) * wshift;
        end
    end
    
    avg_time_to_cmd = mean(time_to_cmd(decision), 'omitnan');
    std_time_to_cmd = std(time_to_cmd(decision), 'omitnan');

    figure
    histogram(time_to_cmd(decision))
    xlabel('Time to command (s)')
    ylabel('Number of trials')
    title('Distribution of Time to Command')
    grid on

    figure
    bar(avg_time_to_cmd)
    hold on
    errorbar(1, avg_time_to_cmd, std_time_to_cmd, 'k', 'LineWidth', 1.5)
    ylabel('Time to command (s)')
    title('Average Time to Deliver Command')
    xticks(1)
    xticklabels({'Subject'})
    grid on
    
    %% save online matrics
    save(sprintf('sub_0%d_online_metrics.mat',subj),'overall_accuracy','accuracy_771','accuracy_773','online_trial_acc','online_trial_acc_771','online_trial_acc_773','decision_rate','avg_time_to_cmd','std_time_to_cmd')
end