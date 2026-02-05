close all
clear all
clc
%% 
n_subjects = 8; 
load 'chanlocs16.mat';
cel_z_cal_ERD_mu_hand = cell(1, n_subjects);
cel_z_cal_ERD_mu_feet = cell(1, n_subjects);
cel_z_cal_ERD_beta_hand = cell(1, n_subjects);
cel_z_cal_ERD_beta_feet = cell(1, n_subjects);

cel_z_cal_ERD_mu_hand_fix = cell(1, n_subjects);
cel_z_cal_ERD_mu_feet_fix = cell(1, n_subjects);
cel_z_cal_ERD_beta_hand_fix = cell(1, n_subjects);
cel_z_cal_ERD_beta_feet_fix = cell(1, n_subjects);

cel_z_cal_ERD_mu_hand_cf = cell(1, n_subjects);
cel_z_cal_ERD_mu_feet_cf = cell(1, n_subjects);
cel_z_cal_ERD_beta_hand_cf = cell(1, n_subjects);
cel_z_cal_ERD_beta_feet_cf = cell(1, n_subjects);

cel_z_eva_ERD_mu_hand = cell(1, n_subjects);
cel_z_eva_ERD_mu_feet = cell(1, n_subjects);
cel_z_eva_ERD_beta_hand = cell(1, n_subjects);
cel_z_eva_ERD_beta_feet = cell(1, n_subjects);

cel_z_eva_ERD_mu_hand_fix = cell(1, n_subjects);
cel_z_eva_ERD_mu_feet_fix = cell(1, n_subjects);
cel_z_eva_ERD_beta_hand_fix = cell(1, n_subjects);
cel_z_eva_ERD_beta_feet_fix = cell(1, n_subjects);

cel_z_eva_ERD_mu_hand_cf = cell(1, n_subjects);
cel_z_eva_ERD_mu_feet_cf = cell(1, n_subjects);
cel_z_eva_ERD_beta_hand_cf = cell(1, n_subjects);
cel_z_eva_ERD_beta_feet_cf = cell(1, n_subjects);
%%
for i = 1:n_subjects
%i = 1; %just to check with one iteration
%% load calibration files
    filename_cal = sprintf('subject_0%d_cal.mat', i);
    load(filename_cal);
    
    % entire trial
    z_cal_ERD_mu_hand     = concat_avg_z_ERD(cal_ERD_mu_hand);
    cel_z_cal_ERD_mu_hand{i} = z_cal_ERD_mu_hand;
    
    z_cal_ERD_beta_hand    = concat_avg_z_ERD(cal_ERD_beta_hand);
    cel_z_cal_ERD_beta_hand{i} = z_cal_ERD_beta_hand;
    
    z_cal_ERD_mu_feet     = concat_avg_z_ERD(cal_ERD_mu_feet);
    cel_z_cal_ERD_mu_feet{i} = z_cal_ERD_mu_feet;
    
    z_cal_ERD_beta_feet    = concat_avg_z_ERD(cal_ERD_beta_feet);
    cel_z_cal_ERD_beta_feet{i} = z_cal_ERD_beta_feet;
    
    % fk only
    z_cal_ERD_mu_hand_fix  = concat_avg_z_ERD(cal_ERD_mu_hand_fix);
    cel_z_cal_ERD_mu_hand_fix{i} = z_cal_ERD_mu_hand_fix;
    
    z_cal_ERD_beta_hand_fix = concat_avg_z_ERD(cal_ERD_beta_hand_fix);
    cel_z_cal_ERD_beta_hand_fix{i} = z_cal_ERD_beta_hand_fix;
    
    z_cal_ERD_mu_feet_fix  = concat_avg_z_ERD(cal_ERD_mu_feet_fix);
    cel_z_cal_ERD_mu_feet_fix{i} = z_cal_ERD_mu_feet_fix;
    
    z_cal_ERD_beta_feet_fix = concat_avg_z_ERD(cal_ERD_beta_feet_fix);
    cel_z_cal_ERD_beta_feet_fix{i} = z_cal_ERD_beta_feet_fix;
    
    % cf only
    z_cal_ERD_mu_hand_cf   = concat_avg_z_ERD(cal_ERD_mu_hand_cf);
    cel_z_cal_ERD_mu_hand_cf{i} = z_cal_ERD_mu_hand_cf;
    
    z_cal_ERD_beta_hand_cf  = concat_avg_z_ERD(cal_ERD_beta_hand_cf);
    cel_z_cal_ERD_beta_hand_cf{i} = z_cal_ERD_beta_hand_cf;
    
    z_cal_ERD_mu_feet_cf   = concat_avg_z_ERD(cal_ERD_mu_feet_cf);
    cel_z_cal_ERD_mu_feet_cf{i} = z_cal_ERD_mu_feet_cf;
    
    z_cal_ERD_beta_feet_cf  = concat_avg_z_ERD(cal_ERD_beta_feet_cf);
    cel_z_cal_ERD_beta_feet_cf{i} = z_cal_ERD_beta_feet_cf;
    
    %% load evaluation files
    filename_eva = sprintf('subject_0%d_eva.mat', i);
    load(filename_eva);
    
    % entire trial
    z_eva_ERD_mu_hand     = concat_avg_z_ERD(eva_ERD_mu_hand);
    cel_z_eva_ERD_mu_hand{i} = z_eva_ERD_mu_hand;
    
    z_eva_ERD_beta_hand    = concat_avg_z_ERD(eva_ERD_beta_hand);
    cel_z_eva_ERD_beta_hand{i} = z_eva_ERD_beta_hand;
    
    z_eva_ERD_mu_feet     = concat_avg_z_ERD(eva_ERD_mu_feet);
    cel_z_eva_ERD_mu_feet{i} = z_eva_ERD_mu_feet;
    
    z_eva_ERD_beta_feet    = concat_avg_z_ERD(eva_ERD_beta_feet);
    cel_z_eva_ERD_beta_feet{i} = z_eva_ERD_beta_feet;
    
    % fk only
    z_eva_ERD_mu_hand_fix  = concat_avg_z_ERD(eva_ERD_mu_hand_fix);
    cel_z_eva_ERD_mu_hand_fix{i} = z_eva_ERD_mu_hand_fix;
    
    z_eva_ERD_beta_hand_fix = concat_avg_z_ERD(eva_ERD_beta_hand_fix);
    cel_z_eva_ERD_beta_hand_fix{i} = z_eva_ERD_beta_hand_fix;
    
    z_eva_ERD_mu_feet_fix  = concat_avg_z_ERD(eva_ERD_mu_feet_fix);
    cel_z_eva_ERD_mu_feet_fix{i} = z_eva_ERD_mu_feet_fix;
    
    z_eva_ERD_beta_feet_fix = concat_avg_z_ERD(eva_ERD_beta_feet_fix);
    cel_z_eva_ERD_beta_feet_fix{i} = z_eva_ERD_beta_feet_fix;
    
    % cf only
    z_eva_ERD_mu_hand_cf   = concat_avg_z_ERD(eva_ERD_mu_hand_cf);
    cel_z_eva_ERD_mu_hand_cf{i} = z_eva_ERD_mu_hand_cf;
    
    z_eva_ERD_beta_hand_cf  = concat_avg_z_ERD(eva_ERD_beta_hand_cf);
    cel_z_eva_ERD_beta_hand_cf{i} = z_eva_ERD_beta_hand_cf;
    
    z_eva_ERD_mu_feet_cf   = concat_avg_z_ERD(eva_ERD_mu_feet_cf);
    cel_z_eva_ERD_mu_feet_cf{i} = z_eva_ERD_mu_feet_cf;
    
    z_eva_ERD_beta_feet_cf  = concat_avg_z_ERD(eva_ERD_beta_feet_cf);
    cel_z_eva_ERD_beta_feet_cf{i} = z_eva_ERD_beta_feet_cf;

    %% Temporal Visualization of the entire trial
    fs = 512; % [Hz]
    time_cal = (0:size(z_cal_ERD_beta_feet,1)-1)/fs; % [s]
    time_eva = (0:size(z_eva_ERD_beta_feet,1)-1)/fs; % [s]
    channels = [7 9 11]; % most appropriate channels
    
    figure(1+(i-1)*6)
    subplot(221)
    plot(time_cal,z_cal_ERD_mu_hand(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in mu band | Calibration phase')
    legend('C3','Cz','C4')
    
    subplot(222)
    plot(time_cal,z_cal_ERD_beta_hand(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in beta band | Calibration phase')
    legend('C3','Cz','C4')
    
    subplot(223)
    plot(time_eva,z_eva_ERD_mu_hand(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in mu band | Evaluation phase')
    legend('C3','Cz','C4')
    
    subplot(224)
    plot(time_eva,z_eva_ERD_beta_hand(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in beta band | Evaluation phase')
    legend('C3','Cz','C4')
    
    sgtitle(sprintf('Subject %d - ERD/ERS - Channels C3,Cz,C4 - Both hands',i))
    %%
    figure(2 +(i-1)*6)
    subplot(221)
    plot(time_cal,z_cal_ERD_mu_feet(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in mu band | Calibration phase')
    legend('C3','Cz','C4')
    
    subplot(222)
    plot(time_cal,z_cal_ERD_beta_feet(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in beta band | Calibration phase')
    legend('C3','Cz','C4')
    
    subplot(223)
    plot(time_eva,z_eva_ERD_mu_feet(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in mu band | Evaluation phase')
    legend('C3','Cz','C4')
    
    subplot(224)
    plot(time_eva,z_eva_ERD_beta_feet(:,channels))
    xlabel('Time [s]')
    ylabel('[ERD/ERS]')
    title('ERD in beta band | Evaluation phase')
    legend('C3','Cz','C4')
    
    sgtitle(sprintf('Subject %d - ERD/ERS - Channels C3,Cz,C4 - Both feet',i))

    %% Spatial Visualization - Topography 
    % ERD beta feet 
    topo_fix_cal = mean(z_cal_ERD_beta_feet_fix, 1);
    topo_cf_cal = mean(z_cal_ERD_beta_feet_cf, 1);
    figure(3 +(i-1)*6)
    subplot(221)
    topoplot(topo_fix_cal, chanlocs16)
    caxis([-25 25])
    colorbar
    title('calibration - fix period')
    
    subplot(223)
    topoplot(topo_cf_cal, chanlocs16)
    caxis([-25 25])
    colorbar
    title('calibration - cf period')
    
    % evaluation ERD beta feet 
    topo_fix_eva = mean(z_eva_ERD_beta_feet_fix, 1);
    topo_cf_eva = mean(z_eva_ERD_beta_feet_cf, 1);
    figure(3 +(i-1)*6)
    subplot(222)
    topoplot(topo_fix_eva, chanlocs16)
    caxis([-25 25])
    colorbar
    title('evaluation - fix period')
    
    subplot(224)
    topoplot(topo_cf_eva, chanlocs16)
    caxis([-25 25])
    colorbar
    title('evaluation - cf period')
    
    sgtitle(sprintf('Subject %d - ERD/ERS - Beta band - Both feet',i))
    %%
    % ERD mu feet 
    topo_fix_cal = mean(z_cal_ERD_mu_feet_fix, 1);
    topo_cf_cal = mean(z_cal_ERD_mu_feet_cf, 1);
    figure(4 +(i-1)*6)
    subplot(221)
    topoplot(topo_fix_cal, chanlocs16)
    caxis([-50 60])
    colorbar
    title('calibration - fix period')
    
    subplot(223)
    topoplot(topo_cf_cal, chanlocs16)
    caxis([-50 60])
    colorbar
    title('calibration - cf period')
    
    % evaluation ERD mu hand 
    topo_fix_eva = mean(z_eva_ERD_mu_feet_fix, 1);
    topo_cf_eva = mean(z_eva_ERD_mu_feet_cf, 1);
    figure(4 +(i-1)*6)
    subplot(222)
    topoplot(topo_fix_eva, chanlocs16)
    caxis([-50 60])
    colorbar
    title('evaluation - fix period')
    
    subplot(224)
    topoplot(topo_cf_eva, chanlocs16)
    caxis([-50 60])
    colorbar
    title('evaluation - cf period')
    
    sgtitle(sprintf('Subject %d - ERD/ERS - Mu band - Both feet',i))
    %%
    % ERD beta hand 
    topo_fix_cal = mean(z_cal_ERD_beta_hand_fix, 1);
    topo_cf_cal = mean(z_cal_ERD_beta_hand_cf, 1);
    figure(5 +(i-1)*6)
    subplot(221)
    topoplot(topo_fix_cal, chanlocs16)
    caxis([-25 25])
    colorbar
    title('calibration - fix period')
    
    subplot(223)
    topoplot(topo_cf_cal, chanlocs16)
    caxis([-25 25])
    colorbar
    title('calibration - cf period')
    
    % evaluation ERD beta hand 
    topo_fix_eva = mean(z_eva_ERD_beta_hand_fix, 1);
    topo_cf_eva = mean(z_eva_ERD_beta_hand_cf, 1);
    figure(5 +(i-1)*6)
    subplot(222)
    topoplot(topo_fix_eva, chanlocs16)
    caxis([-25 25])
    colorbar
    title('evaluation - fix period')
    
    subplot(224)
    topoplot(topo_cf_eva, chanlocs16)
    caxis([-25 25])
    colorbar
    title('evaluation - cf period')
    
    sgtitle(sprintf('Subject %d - ERD/ERS - Beta band - Both hand',i))
    %%
    % ERD mu hand 
    topo_fix_cal = mean(z_cal_ERD_mu_hand_fix, 1);
    topo_cf_cal = mean(z_cal_ERD_mu_hand_cf, 1);
    figure(6+(i-1)*6)
    subplot(221)
    topoplot(topo_fix_cal, chanlocs16)
    caxis([-50 60])
    colorbar
    title('calibration - fix period')
    
    subplot(223)
    topoplot(topo_cf_cal, chanlocs16)
    caxis([-50 60])
    colorbar
    title('calibration - cf period')
    
    % evaluation ERD beta feet 
    topo_fix_eva = mean(z_eva_ERD_mu_hand_fix, 1);
    topo_cf_eva = mean(z_eva_ERD_mu_hand_cf, 1);
    figure(6+(i-1)*6)
    subplot(222)
    topoplot(topo_fix_eva, chanlocs16)
    caxis([-50 60])
    colorbar
    title('evaluation - fix period')
    
    subplot(224)
    topoplot(topo_cf_eva, chanlocs16)
    caxis([-50 60])
    colorbar
    title('evaluation - cf period')
    
    sgtitle(sprintf('Subject %d - ERD/ERS - Mu band - Both hand',i))
end
%% averaging among subjects for each ERD 
% entire trial
subj_avg_cal_ERD_mu_hand = avg_subj_ERD(cel_z_cal_ERD_mu_hand,n_subjects);
subj_avg_cal_ERD_beta_hand = avg_subj_ERD(cel_z_cal_ERD_beta_hand,n_subjects);
subj_avg_eva_ERD_mu_hand = avg_subj_ERD(cel_z_eva_ERD_mu_hand,n_subjects);
subj_avg_eva_ERD_beta_hand = avg_subj_ERD(cel_z_eva_ERD_beta_hand,n_subjects);

subj_avg_cal_ERD_mu_feet = avg_subj_ERD(cel_z_cal_ERD_mu_feet,n_subjects);
subj_avg_cal_ERD_beta_feet = avg_subj_ERD(cel_z_cal_ERD_beta_feet,n_subjects);
subj_avg_eva_ERD_mu_feet = avg_subj_ERD(cel_z_eva_ERD_mu_feet,n_subjects);
subj_avg_eva_ERD_beta_feet = avg_subj_ERD(cel_z_eva_ERD_beta_feet,n_subjects);

% fk trial
subj_avg_cal_ERD_mu_hand_fix = avg_subj_ERD(cel_z_cal_ERD_mu_hand_fix,n_subjects);
subj_avg_cal_ERD_beta_hand_fix = avg_subj_ERD(cel_z_cal_ERD_beta_hand_fix,n_subjects);
subj_avg_eva_ERD_mu_hand_fix = avg_subj_ERD(cel_z_eva_ERD_mu_hand_fix,n_subjects);
subj_avg_eva_ERD_beta_hand_fix = avg_subj_ERD(cel_z_eva_ERD_beta_hand_fix,n_subjects);

subj_avg_cal_ERD_mu_feet_fix = avg_subj_ERD(cel_z_cal_ERD_mu_feet_fix,n_subjects);
subj_avg_cal_ERD_beta_feet_fix = avg_subj_ERD(cel_z_cal_ERD_beta_feet_fix,n_subjects);
subj_avg_eva_ERD_mu_feet_fix = avg_subj_ERD(cel_z_eva_ERD_mu_feet_fix,n_subjects);
subj_avg_eva_ERD_beta_feet_fix = avg_subj_ERD(cel_z_eva_ERD_beta_feet_fix,n_subjects);

% cf trial
subj_avg_cal_ERD_mu_hand_cf = avg_subj_ERD(cel_z_cal_ERD_mu_hand_cf,n_subjects);
subj_avg_cal_ERD_beta_hand_cf = avg_subj_ERD(cel_z_cal_ERD_beta_hand_cf,n_subjects);
subj_avg_eva_ERD_mu_hand_cf = avg_subj_ERD(cel_z_eva_ERD_mu_hand_cf,n_subjects);
subj_avg_eva_ERD_beta_hand_cf = avg_subj_ERD(cel_z_eva_ERD_beta_hand_cf,n_subjects);

subj_avg_cal_ERD_mu_feet_cf = avg_subj_ERD(cel_z_cal_ERD_mu_feet_cf,n_subjects);
subj_avg_cal_ERD_beta_feet_cf = avg_subj_ERD(cel_z_cal_ERD_beta_feet_cf,n_subjects);
subj_avg_eva_ERD_mu_feet_cf = avg_subj_ERD(cel_z_eva_ERD_mu_feet_cf,n_subjects);
subj_avg_eva_ERD_beta_feet_cf = avg_subj_ERD(cel_z_eva_ERD_beta_feet_cf,n_subjects);

%% Temporal Visualization of the entire trial
fs = 512; % [Hz]
time_cal = (0:size(subj_avg_cal_ERD_beta_hand,1)-1)/fs; % [s]
time_eva = (0:size(subj_avg_eva_ERD_beta_hand,1)-1)/fs; % [s]
channels = [7 9 11]; % most appropriate channels

figure(50)
subplot(221)
plot(time_cal,subj_avg_cal_ERD_mu_hand(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('mu band - calibration')
legend('C3','Cz','C4')

subplot(222)
plot(time_cal,subj_avg_cal_ERD_beta_hand(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('beta band | calibration')
legend('C3','Cz','C4')

subplot(223)
plot(time_eva,subj_avg_eva_ERD_mu_hand(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('mu band - evaluation')
legend('C3','Cz','C4')

subplot(224)
plot(time_eva,subj_avg_eva_ERD_beta_hand(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('beta band - evaluation')
legend('C3','Cz','C4')

sgtitle('average ERD/ERS - Channels C3,Cz,C4 - Both hands')

figure(51)
subplot(221)
plot(time_cal,subj_avg_cal_ERD_mu_feet(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('mu band - calibration')
legend('C3','Cz','C4')

subplot(222)
plot(time_cal,subj_avg_cal_ERD_beta_feet(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('beta band - calibration phase')
legend('C3','Cz','C4')

subplot(223)
plot(time_eva,subj_avg_eva_ERD_mu_feet(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('mu band - evaluation')
legend('C3','Cz','C4')

subplot(224)
plot(time_eva,subj_avg_eva_ERD_beta_feet(:,channels))
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('beta band - evaluation')
legend('C3','Cz','C4')

sgtitle('average ERD/ERS - Channels C3,Cz,C4 - Both feet')
%% Spatial Visualization - Topography 
load 'chanlocs16.mat';

% ERD beta feet 
topo_fix_cal = mean(subj_avg_cal_ERD_beta_feet_fix, 1);
topo_cf_cal = mean(subj_avg_cal_ERD_beta_feet_cf, 1);
figure(52)
subplot(221)
topoplot(topo_fix_cal, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – feet – calibration - fix period')

subplot(223)
topoplot(topo_cf_cal, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – feet – calibration - cf period')

% evaluation ERD beta feet 
topo_fix_eva = mean(subj_avg_eva_ERD_beta_feet_fix, 1);
topo_cf_eva = mean(subj_avg_eva_ERD_beta_feet_cf, 1);
figure(52)
subplot(222)
topoplot(topo_fix_eva, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – feet – evaluation - fix period')

subplot(224)
topoplot(topo_cf_eva, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – feet – evaluation - cf period')

sgtitle('ERD/ERS - Beta band - Both feet')
%%
% ERD mu feet 
topo_fix_cal = mean(subj_avg_cal_ERD_mu_feet_fix, 1);
topo_cf_cal = mean(subj_avg_cal_ERD_mu_feet_cf, 1);
figure(53)
subplot(221)
topoplot(topo_fix_cal, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – feet – calibration - fix period')

subplot(223)
topoplot(topo_cf_cal, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – feet – calibration - cf period')

% evaluation ERD beta feet 
topo_fix_eva = mean(subj_avg_eva_ERD_mu_feet_fix, 1);
topo_cf_eva = mean(subj_avg_eva_ERD_mu_feet_cf, 1);
figure(53)
subplot(222)
topoplot(topo_fix_eva, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – feet – evaluation - fix period')

subplot(224)
topoplot(topo_cf_eva, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – feet – evaluation - cf period')

sgtitle('ERD/ERS - Mu band - Both feet')

%%
% ERD beta hand 
topo_fix_cal = mean(subj_avg_cal_ERD_beta_hand_fix, 1);
topo_cf_cal = mean(subj_avg_cal_ERD_beta_hand_cf, 1);
figure(54)
subplot(221)
topoplot(topo_fix_cal, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – hand – calibration - fix period')

subplot(223)
topoplot(topo_cf_cal, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – hand – calibration - cf period')

% evaluation ERD beta hand 
topo_fix_eva = mean(subj_avg_eva_ERD_beta_hand_fix, 1);
topo_cf_eva = mean(subj_avg_eva_ERD_beta_hand_cf, 1);
figure(54)
subplot(222)
topoplot(topo_fix_eva, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – hand – evaluation - fix period')

subplot(224)
topoplot(topo_cf_eva, chanlocs16)
caxis([-25 25])
colorbar
title('ERD beta – hand – evaluation - cf period')

sgtitle('ERD/ERS - Beta band - Both hand')
%%
% ERD mu hand 
topo_fix_cal = mean(subj_avg_cal_ERD_mu_hand_fix, 1);
topo_cf_cal = mean(subj_avg_cal_ERD_mu_hand_cf, 1);
figure(55)
subplot(221)
topoplot(topo_fix_cal, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – hand – calibration - fix period')

subplot(223)
topoplot(topo_cf_cal, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – hand – calibration - cf period')

% evaluation ERD beta feet 
topo_fix_eva = mean(subj_avg_eva_ERD_mu_hand_fix, 1);
topo_cf_eva = mean(subj_avg_eva_ERD_mu_hand_cf, 1);
figure(55)
subplot(222)
topoplot(topo_fix_eva, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – hand – evaluation - fix period')

subplot(224)
topoplot(topo_cf_eva, chanlocs16)
caxis([-50 60])
colorbar
title('ERD mu – hand – evaluation - cf period')

sgtitle('ERD/ERS - Mu band - Both hand')