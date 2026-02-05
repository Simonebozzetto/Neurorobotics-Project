close all
clear all
clc
%% parameters common for all the subjects
ord = 4; % butterworth filter order
mu_band = [8,12];
beta_band = [13,30];
load laplacian16.mat; % laplacian mask
%% labels
labels = [786, 773, 771, 783, 781];
%% load the filelist for each subject
n_subjects = 8;
% after done all of this, we decided to remove a calibration file from the
% subject 2, in particular, the first one!
% Simone, also removed the third calibration file of the second subject!
filelist_calibration = {
    {'ai6.20180316.153104.offline.mi.mi_bhbf.gdf', 'ai6.20180316.154006.offline.mi.mi_bhbf.gdf','ai6.20180316.154811.offline.mi.mi_bhbf.gdf'}, ... % sub 01
    {'ai7.20180316.103209.offline.mi.mi_bhbf.gdf'}, ... % sub 02 {'ai7.20180316.102257.offline.mi.mi_bhbf.gdf',...,'ai7.20180316.104023.offline.mi.mi_bhbf.gdf'}
    {'ai8.20180320.112744.offline.mi.mi_bhbf.gdf','ai8.20180320.113734.offline.mi.mi_bhbf.gdf','ai8.20180320.114543.offline.mi.mi_bhbf.gdf'}, ... % sub 03
    {'aj1.20180312.113542.offline.mi.mi_bhbf.gdf','aj1.20180312.114418.offline.mi.mi_bhbf.gdf'}, ... % sub 04
    {'aj3.20180313.113110.offline.mi.mi_bhbf.gdf','aj3.20180313.114118.offline.mi.mi_bhbf.gdf','aj3.20180313.114946.offline.mi.mi_bhbf.gdf'}, ... % sub 05
    {'aj4.20180313.151634.offline.mi.mi_bhbf.gdf','aj4.20180313.152525.offline.mi.mi_bhbf.gdf','aj4.20180313.153339.offline.mi.mi_bhbf.gdf'}, ... % sub 06
    {'aj7.20180323.161608.offline.mi.mi_bhbf.gdf','aj7.20180323.162629.offline.mi.mi_bhbf.gdf'}, ... % sub 07
    {'aj9.20180326.153615.offline.mi.mi_bhbf.gdf','aj9.20180326.154532.offline.mi.mi_bhbf.gdf','aj9.20180326.155323.offline.mi.mi_bhbf.gdf'}  ... % sub 08
};

filelist_evaluation = {
    {'ai6.20180316.160351.online.mi.mi_bhbf.ema.gdf',...
    'ai6.20180316.161026.online.mi.mi_bhbf.dynamic.gdf',...
    'ai6.20180417.164812.online.mi.mi_bhbf.ema.gdf',...
    'ai6.20180417.165259.online.mi.mi_bhbf.dynamic.gdf',...
    'ai6.20180529.151753.online.mi.mi_bhbf.ema.gdf',...
    'ai6.20180529.152240.online.mi.mi_bhbf.dynamic.gdf'}, ... % sub 01
    {'ai7.20180316.110108.online.mi.mi_bhbf.ema.gdf',...
    'ai7.20180316.110748.online.mi.mi_bhbf.dynamic.gdf',...
    'ai7.20180420.101528.online.mi.mi_bhbf.ema.gdf',...
    'ai7.20180420.101934.online.mi.mi_bhbf.dynamic.gdf',...
    'ai7.20180518.101353.online.mi.mi_bhbf.ema.gdf',...
    'ai7.20180518.101711.online.mi.mi_bhbf.dynamic.gdf'}, ... % sub 02
    {'ai8.20180320.120206.online.mi.mi_bhbf.ema.gdf',...
    'ai8.20180320.122511.online.mi.mi_bhbf.dynamic.gdf',...
    'ai8.20180427.152806.online.mi.mi_bhbf.ema.gdf',...
    'ai8.20180427.154229.online.mi.mi_bhbf.dynamic.gdf',...
    'ai8.20180601.153831.online.mi.mi_bhbf.ema.gdf',...
    'ai8.20180601.154427.online.mi.mi_bhbf.dynamic.gdf'}, ... % sub 03
    {'aj1.20180312.121602.online.mi.mi_bhbf.ema.gdf',...
    'aj1.20180312.122731.online.mi.mi_bhbf.dynamic.gdf',...
    'aj1.20180507.112603.online.mi.mi_bhbf.ema.gdf',...
    'aj1.20180507.113054.online.mi.mi_bhbf.dynamic.gdf',...
    'aj1.20180601.115626.online.mi.mi_bhbf.ema.gdf',...
    'aj1.20180601.120114.online.mi.mi_bhbf.dynamic.gdf'}, ... % sub 04
    {'aj3.20180313.120839.online.mi.mi_bhbf.ema.gdf',...
    'aj3.20180313.121507.online.mi.mi_bhbf.dynamic.gdf',...
    'aj3.20180425.111824.online.mi.mi_bhbf.ema.gdf',...
    'aj3.20180425.112235.online.mi.mi_bhbf.dynamic.gdf',...
    'aj3.20180529.101640.online.mi.mi_bhbf.ema.gdf',...
    'aj3.20180529.102142.online.mi.mi_bhbf.dynamic.gdf'}, ... % sub 05
    {'aj4.20180313.154650.online.mi.mi_bhbf.ema.gdf',...
    'aj4.20180313.155139.online.mi.mi_bhbf.dynamic.gdf',...
    'aj4.20180501.110901.online.mi.mi_bhbf.ema.gdf',...
    'aj4.20180501.111938.online.mi.mi_bhbf.ema.gdf',...
    'aj4.20180522.112515.online.mi.mi_bhbf.ema.gdf',...
    'aj4.20180522.113943.online.mi.mi_bhbf.dynamic.gdf'}, ... % sub 06
    {'aj7.20180323.163958.online.mi.mi_bhbf.ema.gdf',...
    'aj7.20180323.165113.online.mi.mi_bhbf.ema.gdf',...
    'aj7.20180420.170712.online.mi.mi_bhbf.ema.gdf',...
    'aj7.20180420.171150.online.mi.mi_bhbf.dynamic.gdf',...
    'aj7.20180518.154751.online.mi.mi_bhbf.ema.gdf',...
    'aj7.20180518.155307.online.mi.mi_bhbf.dynamic.gdf'}, ... % sub 07
    {'aj9.20180326.160922.online.mi.mi_bhbf.ema.gdf',...
    'aj9.20180326.161314.online.mi.mi_bhbf.dynamic.gdf',...
    'aj9.20180430.152635.online.mi.mi_bhbf.ema.gdf',...
    'aj9.20180430.153028.online.mi.mi_bhbf.dynamic.gdf',...
    'aj9.20180528.112241.online.mi.mi_bhbf.ema.gdf',...
    'aj9.20180528.112738.online.mi.mi_bhbf.dynamic.gdf'}  ... % sub 08
};
%% Calibration files
for subj = 1:n_subjects
% subj = 1;
    fprintf('Processing Subject %02d...\n', subj);
    
    % CALIBRATION FILES
    nCal = length(filelist_calibration{subj});
    cal_s_mu = cell(1, nCal);
    cal_s_beta = cell(1, nCal);
    % for temporal visualization
    cal_ERD_mu_hand = cell(1, nCal);
    cal_ERD_beta_hand = cell(1, nCal);
    cal_ERD_mu_feet = cell(1, nCal);
    cal_ERD_beta_feet = cell(1, nCal);
    % for topoplot
    cal_ERD_mu_hand_fix = cell(1, nCal);
    cal_ERD_beta_hand_fix = cell(1, nCal);
    cal_ERD_mu_feet_fix = cell(1, nCal);
    cal_ERD_beta_feet_fix = cell(1, nCal);
    cal_ERD_mu_hand_cf = cell(1, nCal);
    cal_ERD_beta_hand_cf = cell(1, nCal);
    cal_ERD_mu_feet_cf = cell(1, nCal);
    cal_ERD_beta_feet_cf = cell(1, nCal);


    for i = 1:nCal
        [s, h] = sload(filelist_calibration{subj}{i});

        [Tk,n_trials,min_length_trial,Ck,Fk,fk_start,fk_end,min_length_fk,CFk,min_length_cf] = label_vectors_1(s,h.EVENT,labels);
        
        s_mu = preprocessing_1(s, h, ord, mu_band, lap);
        s_beta = preprocessing_1(s, h, ord, beta_band, lap);
        
        n_channels = size(s_mu,2);

        % sure it's possible to adapt the function conc_to_trial to this
        % one
        fixData_mu = zeros(min_length_fk,n_channels,n_trials);
        fixData_beta = zeros(min_length_fk,n_channels,n_trials);
        for j = 1:n_trials
            s_fk_mu = s_mu(fk_start(j):fk_end(j),:);
            if size(s_fk_mu,1) ~= min_length_fk
                s_fk_mu = s_fk_mu(1:min_length_fk,:);
            end
            fixData_mu(:,:,j) = s_fk_mu;
        
            s_fk_beta = s_beta(fk_start(j):fk_end(j),:);
            if size(s_fk_beta,1) ~= min_length_fk
                s_fk_beta = s_fk_beta(1:min_length_fk,:);
            end
            fixData_beta(:,:,j) = s_fk_beta; 
        end
        clear j
        clear s_fk_mu
        clear s_fk_beta
        
        s_mu_trials = conc_to_trial(s_mu,min_length_trial,n_channels,n_trials,Tk);
        s_beta_trials = conc_to_trial(s_beta,min_length_trial,n_channels,n_trials,Tk);

        % ERD
        reference_mu = repmat(mean(fixData_mu), [size(s_mu_trials, 1) 1 1]); 
        ERD_mu = 100 * (s_mu_trials - reference_mu)./ (reference_mu);
        
        reference_beta = repmat(mean(fixData_beta), [size(s_beta_trials, 1) 1 1]); 
        ERD_beta = 100 * (s_beta_trials - reference_beta)./ reference_beta;

        % ERD divided based on the class
        idx_hand = find(Ck == labels(2));
        idx_feet = find(Ck == labels(3));
        idx_rest = find(Ck == labels(4));
        
        n_trials_hand = length(idx_hand);
        n_trials_feet = length(idx_feet);
        n_trials_rest = length(idx_rest);
        
        ERD_mu_hand = ERD_mu(:,:,idx_hand);
        ERD_mu_feet = ERD_mu(:,:,idx_feet);
        if n_trials_rest ~= 0
            ERD_mu_rest = ERD_beta(:,:,idx_rest);
        end
        
        ERD_mu_hand_fix = ERD_mu_hand(1:min_length_fk,:,:);
        ERD_mu_feet_fix = ERD_mu_feet(1:min_length_fk,:,:);
        ERD_mu_hand_cf = ERD_mu_hand(end-min_length_cf:end,:,:);
        ERD_mu_feet_cf = ERD_mu_feet(end-min_length_cf:end,:,:);
        
               
        ERD_beta_hand = ERD_beta(:,:,idx_hand);
        ERD_beta_feet = ERD_beta(:,:,idx_feet);
        if n_trials_rest ~= 0
            ERD_beta_rest = ERD_beta(:,:,idx_rest);
        end

        ERD_beta_hand_fix = ERD_beta_hand(1:min_length_fk,:,:);
        ERD_beta_feet_fix = ERD_beta_feet(1:min_length_fk,:,:);
        ERD_beta_hand_cf = ERD_beta_hand(end-min_length_cf:end,:,:);
        ERD_beta_feet_cf = ERD_beta_feet(end-min_length_cf:end,:,:);

        cal_s_mu{i} = s_mu;
        cal_s_beta{i} = s_beta;
        % for temporal visualization
        cal_ERD_mu_hand{i} = ERD_mu_hand;
        cal_ERD_beta_hand{i} = ERD_beta_hand;
        cal_ERD_mu_feet{i} = ERD_mu_feet;
        cal_ERD_beta_feet{i} = ERD_beta_feet;
        % for topoplot
        cal_ERD_mu_hand_fix{i} = ERD_mu_hand_fix;
        cal_ERD_beta_hand_fix{i} = ERD_beta_hand_fix;
        cal_ERD_mu_feet_fix{i} = ERD_mu_feet_fix;
        cal_ERD_beta_feet_fix{i} = ERD_beta_feet_fix;
        cal_ERD_mu_hand_cf{i} = ERD_mu_hand_cf;
        cal_ERD_beta_hand_cf{i} = ERD_beta_hand_cf;
        cal_ERD_mu_feet_cf{i} = ERD_mu_feet_cf;
        cal_ERD_beta_feet_cf{i} = ERD_beta_feet_cf;
        
        %cal_EVENT{i} = h.EVENT;
    end
     save(sprintf('subject_%02d_cal.mat', subj), 'cal_s_mu', 'cal_s_beta',...
         'cal_ERD_mu_hand', 'cal_ERD_beta_hand','cal_ERD_mu_feet', 'cal_ERD_beta_feet',...
         'cal_ERD_mu_hand_fix', 'cal_ERD_beta_hand_fix', 'cal_ERD_mu_feet_fix','cal_ERD_beta_feet_fix', ...
         "cal_ERD_mu_hand_cf",'cal_ERD_beta_hand_cf','cal_ERD_mu_feet_cf','cal_ERD_beta_feet_cf');
end
%% Evaluation files
for subj = 1:n_subjects
    %subj = 1;
    fprintf('Processing Subject %02d...\n', subj);
    
    nEva = length(filelist_evaluation{subj});
    eva_s_mu = cell(1, nEva);
    eva_s_beta = cell(1, nEva);
    % for temporal visualization
    eva_ERD_mu_hand = cell(1, nEva);
    eva_ERD_beta_hand = cell(1, nEva);
    eva_ERD_mu_feet = cell(1, nEva);
    eva_ERD_beta_feet = cell(1, nEva);
    % for topoplot
    eva_ERD_mu_hand_fix = cell(1, nEva);
    eva_ERD_beta_hand_fix = cell(1, nEva);
    eva_ERD_mu_feet_fix = cell(1, nEva);
    eva_ERD_beta_feet_fix = cell(1, nEva);
    eva_ERD_mu_hand_cf = cell(1, nEva);
    eva_ERD_beta_hand_cf = cell(1, nEva);
    eva_ERD_mu_feet_cf = cell(1, nEva);
    eva_ERD_beta_feet_cf = cell(1, nEva);

    for i = 1:nEva
        [s, h] = sload(filelist_evaluation{subj}{i});

        [Tk,n_trials,min_length_trial,Ck,Fk,fk_start,fk_end,min_length_fk,CFk,min_length_cf] = label_vectors_1(s,h.EVENT,labels);
        
        s_mu = preprocessing_1(s, h, ord, mu_band, lap);
        s_beta = preprocessing_1(s, h, ord, beta_band, lap);
        
        n_channels = size(s_mu,2);

        % sure it's possible to adapt the function conc_to_trial to this
        % one
        fixData_mu = zeros(min_length_fk,n_channels,n_trials);
        fixData_beta = zeros(min_length_fk,n_channels,n_trials);
        for j = 1:n_trials
            s_fk_mu = s_mu(fk_start(j):fk_end(j),:);
            if size(s_fk_mu,1) ~= min_length_fk
                s_fk_mu = s_fk_mu(1:min_length_fk,:);
            end
            fixData_mu(:,:,j) = s_fk_mu;
        
            s_fk_beta = s_beta(fk_start(j):fk_end(j),:);
            if size(s_fk_beta,1) ~= min_length_fk
                s_fk_beta = s_fk_beta(1:min_length_fk,:);
            end
            fixData_beta(:,:,j) = s_fk_beta; 
        end
        clear j
        clear s_fk_mu
        clear s_fk_beta
        
        s_mu_trials = conc_to_trial(s_mu,min_length_trial,n_channels,n_trials,Tk);
        s_beta_trials = conc_to_trial(s_beta,min_length_trial,n_channels,n_trials,Tk);

        % ERD
        reference_mu = repmat(mean(fixData_mu), [size(s_mu_trials, 1) 1 1]); 
        ERD_mu = 100 * (s_mu_trials - reference_mu)./ (reference_mu);
        
        reference_beta = repmat(mean(fixData_beta), [size(s_beta_trials, 1) 1 1]); 
        ERD_beta = 100 * (s_beta_trials - reference_beta)./ reference_beta;
    
        % ERD divided based on the class
        idx_hand = find(Ck == labels(2));
        idx_feet = find(Ck == labels(3));
        idx_rest = find(Ck == labels(4));
        
        n_trials_hand = length(idx_hand);
        n_trials_feet = length(idx_feet);
        n_trials_rest = length(idx_rest);
        
        ERD_mu_hand = ERD_mu(:,:,idx_hand);
        ERD_mu_feet = ERD_mu(:,:,idx_feet);
        if n_trials_rest ~= 0
            ERD_mu_rest = ERD_beta(:,:,idx_rest);
        end
        
        ERD_mu_hand_fix = ERD_mu_hand(1:min_length_fk,:,:);
        ERD_mu_feet_fix = ERD_mu_feet(1:min_length_fk,:,:);
        ERD_mu_hand_cf = ERD_mu_hand(end-min_length_cf:end,:,:);
        ERD_mu_feet_cf = ERD_mu_feet(end-min_length_cf:end,:,:);
        
               
        ERD_beta_hand = ERD_beta(:,:,idx_hand);
        ERD_beta_feet = ERD_beta(:,:,idx_feet);
        if n_trials_rest ~= 0
            ERD_beta_rest = ERD_beta(:,:,idx_rest);
        end

        ERD_beta_hand_fix = ERD_beta_hand(1:min_length_fk,:,:);
        ERD_beta_feet_fix = ERD_beta_feet(1:min_length_fk,:,:);
        ERD_beta_hand_cf = ERD_beta_hand(end-min_length_cf:end,:,:);
        ERD_beta_feet_cf = ERD_beta_feet(end-min_length_cf:end,:,:);

        eva_s_mu{i} = s_mu;
        eva_s_beta{i} = s_beta;
        % for temporal visualization
        eva_ERD_mu_hand{i} = ERD_mu_hand;
        eva_ERD_beta_hand{i} = ERD_beta_hand;
        eva_ERD_mu_feet{i} = ERD_mu_feet;
        eva_ERD_beta_feet{i} = ERD_beta_feet;
        % for topoplot
        eva_ERD_mu_hand_fix{i} = ERD_mu_hand_fix;
        eva_ERD_beta_hand_fix{i} = ERD_beta_hand_fix;
        eva_ERD_mu_feet_fix{i} = ERD_mu_feet_fix;
        eva_ERD_beta_feet_fix{i} = ERD_beta_feet_fix;
        eva_ERD_mu_hand_cf{i} = ERD_mu_hand_cf;
        eva_ERD_beta_hand_cf{i} = ERD_beta_hand_cf;
        eva_ERD_mu_feet_cf{i} = ERD_mu_feet_cf;
        eva_ERD_beta_feet_cf{i} = ERD_beta_feet_cf;    
    end
    save(sprintf('subject_%02d_eva.mat', subj), 'eva_s_mu', 'eva_s_beta',...
         'eva_ERD_mu_hand', 'eva_ERD_beta_hand','eva_ERD_mu_feet', 'eva_ERD_beta_feet',...
         'eva_ERD_mu_hand_fix', 'eva_ERD_beta_hand_fix', 'eva_ERD_mu_feet_fix','eva_ERD_beta_feet_fix', ...
         "eva_ERD_mu_hand_cf",'eva_ERD_beta_hand_cf','eva_ERD_mu_feet_cf','eva_ERD_beta_feet_cf');
end
