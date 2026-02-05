close all
clear all
clc
%%
n_subjects = 8;

% offline sample accuracy
offline_sample_accuracy = nan(n_subjects,1);
offline_sample_accuracy_771 = nan(n_subjects,1);
offline_sample_accuracy_773 = nan(n_subjects,1);
% offline trial accuracy
offline_trial_accuracy  = nan(n_subjects,1);
offline_trial_accuracy_771  = nan(n_subjects,1);
offline_trial_accuracy_773  = nan(n_subjects,1);
% online sample accuracy
online_sample_accuracy = zeros(n_subjects,1);
online_sample_accuracy_771 = zeros(n_subjects,1);
online_sample_accuracy_773 = zeros(n_subjects,1);
% online trial accuracy
online_trial_accuracy  = zeros(n_subjects,1);
online_trial_accuracy_771  = zeros(n_subjects,1);
online_trial_accuracy_773  = zeros(n_subjects,1);
% decision rate
decision_rate  = zeros(n_subjects,1);
% average and std time to command
ttc_avg = zeros(n_subjects,1);
ttc_std  = zeros(n_subjects,1);

for subj = 1:n_subjects
    % offline metrics
    file_name_off = sprintf('sub_%02d_offline_metrics.mat', subj);
    if exist(file_name_off, 'file') % needed to add this for the second subject
        metrics_off = load(file_name_off);
        % offline sample accuracy
        offline_sample_accuracy(subj) = metrics_off.offline_sample_acc;
        offline_sample_accuracy_771(subj) = metrics_off.offline_sample_acc_771;
        offline_sample_accuracy_773(subj) = metrics_off.offline_sample_acc_773;
        % offline trial accuracy
        offline_trial_accuracy(subj) = metrics_off.offline_trial_acc;
        offline_trial_accuracy_771(subj) = metrics_off.offline_trial_acc_771;
        offline_trial_accuracy_773(subj) = metrics_off.offline_trial_acc_773;
    else
        warning('Missing %subj', file_name_off);
    end

    % online metrics
    file_name_on = sprintf('sub_%02d_online_metrics.mat', subj);
    if exist(file_name_on, 'file') 
        metrics_on = load(file_name_on);
        % online sample accuracy
        online_sample_accuracy(subj) = metrics_on.overall_accuracy;
        online_sample_accuracy_771(subj) = metrics_on.accuracy_771;
        online_sample_accuracy_773(subj) = metrics_on.accuracy_773;
        % online trial accuracy
        online_trial_accuracy(subj) = metrics_on.online_trial_acc;
        online_trial_accuracy_771(subj) = metrics_on.online_trial_acc_771;
        online_trial_accuracy_773(subj) = metrics_on.online_trial_acc_773;
        % decision rate
        decision_rate(subj)  = metrics_on.decision_rate;
        % average and std time to command
        ttc_avg(subj) = metrics_on.avg_time_to_cmd;
        ttc_std(subj)  = metrics_on.std_time_to_cmd;
    else
        warning('Missing %subj', file_name_on);
    end
end
clear file_name_on
clear file_name_off
clear subj
clear metrics_on
clear metrics_off

% plots
figure
bar(offline_sample_accuracy)
yline(0.5,'r--') % chance level
xlabel('Subject');
ylabel('Accuracy');
title('Offline single sample accuracy');

figure
bar(offline_trial_accuracy)
yline(0.5,'r--')
xlabel('Subject')
ylabel('Accuracy')
title('Offline Trial Accuracy')

figure
bar(online_sample_accuracy)
yline(0.5,'r--') % chance level
xlabel('Subject');
ylabel('Accuracy');
title('Online single sample accuracy');

figure
bar(online_trial_accuracy)
yline(0.5,'r--')
xlabel('Subject')
ylabel('Accuracy')
title('Online Trial Accuracy')

figure
bar(ttc_avg)
xlabel('Subject')
ylabel('Seconds')
title('Time to deliver a command')


% average of the results
m_ = @(x) mean(x, 'omitnan');
s_ = @(x) std(x,  'omitnan');

fprintf('\nGRAND AVERAGE for the BMI decoding part (mean ± std) \n');
fprintf('Offline sample acc: %.4f ± %.4f\n', m_(offline_sample_accuracy), s_(offline_sample_accuracy));
fprintf('Offline trial  acc: %.4f ± %.4f\n', m_(offline_trial_accuracy),  s_(offline_trial_accuracy));
fprintf('Online  sample acc: %.4f ± %.4f\n', m_(online_sample_accuracy),  s_(online_sample_accuracy));
fprintf('Online  trial  acc: %.4f ± %.4f\n', m_(online_trial_accuracy),   s_(online_trial_accuracy));
fprintf('\nControl metrics:\n');
fprintf('Decision rate      : %.4f ± %.4f\n', m_(decision_rate), s_(decision_rate));
fprintf('Time-to-command (s): %.4f ± %.4f\n', m_(ttc_avg), s_(ttc_avg));

%% 
% boxplot of accuracies
figure
boxplot([offline_sample_accuracy, offline_trial_accuracy, online_sample_accuracy, online_trial_accuracy],'Labels', {'Offline sample','Offline trial','Online sample','Online trial'});
ylabel('Accuracy')
title('Grand Average Accuracies Across Subjects')
grid on

% boxplot of the time to deliver a command
figure
boxplot(ttc_avg)
ylabel('Time to command (s)')
title('Grand Average Time to Command')
grid on

% boxplot of the decision rate
figure
boxplot(decision_rate)
ylabel('Decision rate')
title('Grand Average Decision Rate')
grid on

figure
scatter(ttc_avg, online_trial_accuracy, 'filled')
xlabel('Avg time to command (s)')
ylabel('Online trial accuracy')
title('Accuracy–Latency Trade-off')
grid on