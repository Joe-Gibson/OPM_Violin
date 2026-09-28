%% Violin Beamformer
clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

% If on new desktop, change onedrive location from C: to D:
if strcmp(getenv('COMPUTERNAME'), 'DUIP78552')
    Onedrive_path = 'D:\OneDrive - The University of Nottingham\';
else
    Onedrive_path = 'C:\Users\ppyjg12\OneDrive - The University of Nottingham\';
end

addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files"))
addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\tools"))
addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\fieldtrip-20231113"))
addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\Beamformer"))
addpath(strcat(Onedrive_path, "Documents\MATLAB_C"))
addpath(strcat(Onedrive_path, "Documents\MATLAB_C\slanCM\slanCM"))

addpath(strcat(Onedrive_path, "Documents\MATLAB_C\Optitrack"))
addpath('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts')
ft_defaults

beta = 1;


hp = 13;
lp = 30;
TFS_nos = 106;
ylim_TFS = [6 35];

run_all_timing = tic;

freq_savename = ['_' num2str(hp) '-' num2str(lp) '_'];

bids_path = 'D:\OneDrive - The University of Nottingham\Documents\Data\Violin\V2\Bids\';


runs = ["run_1"; "run_2"];

load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\MARS\MARS_cortical_Key.mat')


sub_nums = [1:22];
sub_nums(21) = []; % Lots of noisy trials
sub_nums(16) = []; % Exclude sub 16 - wrong triggers
sub_nums(1) = []; % Exclude sub 1 - wrong triggers

n_subjects = length(sub_nums);


TFS_l_playing_avg_sub = zeros(106, 3750, n_subjects, 2);
TFS_r_playing_avg_sub = zeros(106, 3750, n_subjects, 2);
H_VE_l_playing_avg_sub = zeros(3750, n_subjects, 2);
H_VE_r_playing_avg_sub = zeros(3750, n_subjects, 2);
TFS_l_rebound_avg_sub = zeros(106, 3750, n_subjects, 2);
TFS_r_rebound_avg_sub = zeros(106, 3750, n_subjects, 2);
H_VE_l_rebound_avg_sub = zeros(3750, n_subjects, 2);
H_VE_r_rebound_avg_sub = zeros(3750, n_subjects, 2);
TFS_l_rest_avg_sub = zeros(106, 3750, n_subjects, 2);
TFS_r_rest_avg_sub = zeros(106, 3750, n_subjects, 2);
H_VE_l_rest_avg_sub = zeros(3750, n_subjects, 2);
H_VE_r_rest_avg_sub = zeros(3750, n_subjects, 2);


f = 375;
duration = 10;
trig_offset_rest = 10; % Look n seconds after trigger
trial_time_rest = [1:duration*f]./f + trig_offset_rest;
trig_offset_rebound = -5; % Look n seconds after trigger
trial_time_rebound = [1:duration*f]./f + trig_offset_rebound;
trig_offset_playing = 0; % Look n seconds before trigger
trial_time_playing = [1:duration*f]./f + trig_offset_playing;


n_to_load = n_subjects*2;
load_val = 1;
%% Choose files
load([bids_path 'Derivatives\score_orders.mat'])%, "scores_sorted", "trial_score_order");

trial_score_orders = NaN(n_subjects, 60);
for sub = n_subjects:-1:1
    sub_number = sub_nums(sub);
    trial_score_orders(sub, 1:length(trial_score_order.(['sub_' num2str(sub_number)]))) = trial_score_order.(['sub_' num2str(sub_number)]);
    ntrials(sub) = length(trial_score_order.(['sub_' num2str(sub_number)]));
end

ntrials_third = floor(2*ntrials/5);
ylim_vals_env = [-50 75];

file_to_load = ['D:\OneDrive - The University of Nottingham\Documents\Data\Violin\V2\' ...
    'BIDS\Derivatives\BF_VEs\Peak_tstat\Peak_tstat_avg_TFS_H_VE_mat_relative_' num2str(hp) '_' num2str(lp) '_scores_cov.mat'];

if ~isfile(file_to_load)

    % ft_progress('init', 'etf', 'Loading data...')      % ascii progress bar
    wait_fig = waitbar(0,'Please wait...');
    % disp(['Location ' num2str(loc)])

    for sub = 1:n_subjects
        sub_number = sub_nums(sub);
        subject = ['Sub-' num2str(sub_number, '%03.f')];
        % disp(subject)

        TFS_l_playing  = NaN(TFS_nos, duration*f, 60);
        TFS_r_playing = NaN(TFS_nos, duration*f, 60);
        H_VE_l_playing  = NaN(duration*f, 60);
        H_VE_r_playing  = NaN(duration*f, 60);
        TFS_l_rebound = NaN(TFS_nos, duration*f, 60);
        TFS_r_rebound = NaN(TFS_nos, duration*f, 60);
        H_VE_l_rebound = NaN(duration*f, 60);
        H_VE_r_rebound = NaN(duration*f, 60);
        TFS_l_rest = NaN(TFS_nos, duration*f, 60);
        TFS_r_rest = NaN(TFS_nos, duration*f, 60);
        H_VE_l_rest = NaN(duration*f, 60);
        H_VE_r_rest = NaN(duration*f, 60);
        ntrials_run1 = 0;

        for run_no = 1:2
            run = ['run_' num2str(run_no)];
            % disp(run)
            % ft_progress(load_val/n_to_load, 'Loading part %d of %d', load_val, n_to_load)
            waitbar(load_val/n_to_load,wait_fig,sprintf('Loading part %d of %d', load_val, n_to_load));

            load([bids_path,'derivatives\BF_VEs\Peak_tstat\' subject '\' run ...
                '_3trigs_unaveraged_peak_tstat_cov.mat']);

            TFS_l_playing(:, :, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = TFS_l_playing_all;
            TFS_r_playing(:, :, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = TFS_r_playing_all;
            TFS_l_rebound(:, :, ntrials_run1 +(1:size(TFS_r_rebound_all, 3))) = TFS_l_rebound_all;
            TFS_r_rebound(:, :, ntrials_run1 +(1:size(TFS_r_rebound_all, 3))) = TFS_r_rebound_all;
            TFS_l_rest(:, :, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = TFS_l_rest_all;
            TFS_r_rest(:, :, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = TFS_r_rest_all;

            H_VE_l_playing_all =  Hilbert_envelope_l_playing_trials;
            H_VE_r_playing_all =  Hilbert_envelope_r_playing_trials;
            H_VE_l_rebound_all =  Hilbert_envelope_l_rebound_trials;
            H_VE_r_rebound_all =  Hilbert_envelope_r_rebound_trials;
            H_VE_l_rest_all =  Hilbert_envelope_l_rest_trials;
            H_VE_r_rest_all =  Hilbert_envelope_r_rest_trials;


            % Relative
            meanrest_l_H_VE =  repmat(squeeze(mean(H_VE_l_rest_all, 1)), [size(H_VE_l_playing_all, 1), 1]);
            H_VE_l_playing_all = (H_VE_l_playing_all-meanrest_l_H_VE)./meanrest_l_H_VE*100;
            H_VE_l_rebound_all = (H_VE_l_rebound_all - meanrest_l_H_VE)./meanrest_l_H_VE*100;
            H_VE_l_rest_all = (H_VE_l_rest_all - meanrest_l_H_VE)./meanrest_l_H_VE*100;

            meanrest_r_H_VE =  repmat(squeeze(mean(H_VE_r_rest_all, 1)), [size(H_VE_r_playing_all, 1), 1]);
            H_VE_r_playing_all = (H_VE_r_playing_all-meanrest_r_H_VE)./meanrest_r_H_VE*100;
            H_VE_r_rebound_all = (H_VE_r_rebound_all - meanrest_r_H_VE)./meanrest_r_H_VE*100;
            H_VE_r_rest_all = (H_VE_r_rest_all - meanrest_r_H_VE)./meanrest_r_H_VE*100;

            H_VE_l_playing(:, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = H_VE_l_playing_all;
            H_VE_r_playing(:, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = H_VE_r_playing_all;
            H_VE_l_rebound(:, ntrials_run1 +(1:size(TFS_r_rebound_all, 3))) = H_VE_l_rebound_all;
            H_VE_r_rebound(:, ntrials_run1 +(1:size(TFS_r_rebound_all, 3))) = H_VE_r_rebound_all;
            H_VE_l_rest(:, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = H_VE_l_rest_all;
            H_VE_r_rest(:, ntrials_run1 +(1:size(TFS_r_rest_all, 3))) = H_VE_r_rest_all;
            ntrials_run1 = size(TFS_r_rest_all, 3);

            load_val = load_val + 1;
        end


        clear meanrest_l_H_VE meanrest_r_H_VE meanrest_l_TFS meanrest_r_TFS ...
            TFS_r_rebound_all TFS_r_rest_all TFS_r_playing_all TFS_l_rebound_all TFS_l_rest_all TFS_l_playing_all ...
            H_VE_l_rebound_all H_VE_l_rest_all H_VE_l_playing_all H_VE_r_rebound_all H_VE_r_rest_all H_VE_r_playing_all

        TFS_l_playing_avg_sub(:, :, sub, 2) = mean(TFS_l_playing(:,:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 3, 'omitnan');
        TFS_l_playing_avg_sub(:, :, sub, 1) = mean(TFS_l_playing(:,:,trial_score_orders(1:ntrials_third(sub))), 3, 'omitnan');
        TFS_r_playing_avg_sub(:, :, sub, 2) = mean(TFS_r_playing(:,:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 3, 'omitnan');
        TFS_r_playing_avg_sub(:, :, sub, 1) = mean(TFS_r_playing(:,:,trial_score_orders(1:ntrials_third(sub))), 3, 'omitnan');
        TFS_l_rebound_avg_sub(:, :, sub, 2) = mean(TFS_l_rebound(:,:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 3, 'omitnan');
        TFS_l_rebound_avg_sub(:, :, sub, 1) = mean(TFS_l_rebound(:,:,trial_score_orders(1:ntrials_third(sub))), 3, 'omitnan');
        TFS_r_rebound_avg_sub(:, :, sub, 2) = mean(TFS_r_rebound(:,:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 3, 'omitnan');
        TFS_r_rebound_avg_sub(:, :, sub, 1) = mean(TFS_r_rebound(:,:,trial_score_orders(1:ntrials_third(sub))), 3, 'omitnan');
        TFS_l_rest_avg_sub(:, :, sub, 2) = mean(TFS_l_rest(:,:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 3, 'omitnan');
        TFS_l_rest_avg_sub(:, :, sub, 1) = mean(TFS_l_rest(:,:,trial_score_orders(1:ntrials_third(sub))), 3, 'omitnan');
        TFS_r_rest_avg_sub(:, :, sub, 2) = mean(TFS_r_rest(:,:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 3, 'omitnan');
        TFS_r_rest_avg_sub(:, :, sub, 1) = mean(TFS_r_rest(:,:,trial_score_orders(1:ntrials_third(sub))), 3, 'omitnan');

        H_VE_l_playing_avg_sub(:, sub, 2) = mean(H_VE_l_playing(:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 2, 'omitnan');
        H_VE_l_playing_avg_sub(:, sub, 1) = mean(H_VE_l_playing(:,trial_score_orders(1:ntrials_third(sub))), 2, 'omitnan');
        H_VE_r_playing_avg_sub(:, sub, 2) = mean(H_VE_r_playing(:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 2, 'omitnan');
        H_VE_r_playing_avg_sub(:, sub, 1) = mean(H_VE_r_playing(:,trial_score_orders(1:ntrials_third(sub))), 2, 'omitnan');
        H_VE_l_rebound_avg_sub(:, sub, 2) = mean(H_VE_l_rebound(:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 2, 'omitnan');
        H_VE_l_rebound_avg_sub(:, sub, 1) = mean(H_VE_l_rebound(:,trial_score_orders(1:ntrials_third(sub))), 2, 'omitnan');
        H_VE_r_rebound_avg_sub(:, sub, 2) = mean(H_VE_r_rebound(:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 2, 'omitnan');
        H_VE_r_rebound_avg_sub(:, sub, 1) = mean(H_VE_r_rebound(:,trial_score_orders(1:ntrials_third(sub))), 2, 'omitnan');
        H_VE_l_rest_avg_sub(:, sub, 2) = mean(H_VE_l_rest(:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 2, 'omitnan');
        H_VE_l_rest_avg_sub(:, sub, 1) = mean(H_VE_l_rest(:,trial_score_orders(1:ntrials_third(sub))), 2, 'omitnan');
        H_VE_r_rest_avg_sub(:, sub, 2) = mean(H_VE_r_rest(:,trial_score_orders(ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 2, 'omitnan');
        H_VE_r_rest_avg_sub(:, sub, 1) = mean(H_VE_r_rest(:,trial_score_orders(1:ntrials_third(sub))), 2, 'omitnan');

    end
    TFS_l_playing_avg = squeeze(mean(TFS_l_playing_avg_sub, 3, 'omitnan'));
    TFS_r_playing_avg = squeeze(mean(TFS_r_playing_avg_sub, 3, 'omitnan'));
    TFS_l_rebound_avg = squeeze(mean(TFS_l_rebound_avg_sub, 3, 'omitnan'));
    TFS_r_rebound_avg = squeeze(mean(TFS_r_rebound_avg_sub, 3, 'omitnan'));
    TFS_l_rest_avg = squeeze(mean(TFS_l_rest_avg_sub, 3, 'omitnan'));
    TFS_r_rest_avg = squeeze(mean(TFS_r_rest_avg_sub, 3, 'omitnan'));
    H_VE_l_playing_avg = squeeze(mean(H_VE_l_playing_avg_sub, 2, 'omitnan'));
    H_VE_r_playing_avg = squeeze(mean(H_VE_r_playing_avg_sub, 2, 'omitnan'));
    H_VE_l_rebound_avg = squeeze(mean(H_VE_l_rebound_avg_sub, 2, 'omitnan'));
    H_VE_r_rebound_avg = squeeze(mean(H_VE_r_rebound_avg_sub, 2, 'omitnan'));
    H_VE_l_rest_avg = squeeze(mean(H_VE_l_rest_avg_sub, 2, 'omitnan'));
    H_VE_r_rest_avg = squeeze(mean(H_VE_r_rest_avg_sub, 2, 'omitnan'));

    % ft_progress('close')
    close(wait_fig)

    save(file_to_load, "TFS_l_playing_avg_sub", ...
        "TFS_r_playing_avg_sub", "TFS_l_rebound_avg_sub", "TFS_r_rebound_avg_sub", "TFS_l_rest_avg_sub", "TFS_r_rest_avg_sub", "H_VE_l_playing_avg_sub", ...
        "H_VE_r_playing_avg_sub", "H_VE_l_rebound_avg_sub", "H_VE_r_rebound_avg_sub", "H_VE_l_rest_avg_sub", "H_VE_r_rest_avg_sub")

else
    load(file_to_load)
    TFS_l_playing_avg = squeeze(mean(TFS_l_playing_avg_sub, 3, 'omitnan'));
    TFS_r_playing_avg = squeeze(mean(TFS_r_playing_avg_sub, 3, 'omitnan'));
    TFS_l_rebound_avg = squeeze(mean(TFS_l_rebound_avg_sub, 3, 'omitnan'));
    TFS_r_rebound_avg = squeeze(mean(TFS_r_rebound_avg_sub, 3, 'omitnan'));
    TFS_l_rest_avg = squeeze(mean(TFS_l_rest_avg_sub, 3, 'omitnan'));
    TFS_r_rest_avg = squeeze(mean(TFS_r_rest_avg_sub, 3, 'omitnan'));
    H_VE_l_playing_avg = squeeze(mean(H_VE_l_playing_avg_sub, 2, 'omitnan'));
    H_VE_r_playing_avg = squeeze(mean(H_VE_r_playing_avg_sub, 2, 'omitnan'));
    H_VE_l_rebound_avg = squeeze(mean(H_VE_l_rebound_avg_sub, 2, 'omitnan'));
    H_VE_r_rebound_avg = squeeze(mean(H_VE_r_rebound_avg_sub, 2, 'omitnan'));
    H_VE_l_rest_avg = squeeze(mean(H_VE_l_rest_avg_sub, 2, 'omitnan'));
    H_VE_r_rest_avg = squeeze(mean(H_VE_r_rest_avg_sub, 2, 'omitnan'));
end


%%
cbar_lim = 0.3; % colour bar limit (relative change)

if size(TFS_l_playing_avg(:, :, 1), 1) == 110
    highpass = 1:110;
    lowpass = 2.5:111.5;
    fre = highpass + ((lowpass - highpass)./2);
elseif size(TFS_l_playing_avg(:, :, 1), 1) == 106
    highpass = 5:110;
    lowpass = 6.5:111.5;
    fre = highpass + ((lowpass - highpass)./2);
end

figure()
% subplot(323)
tiledlayout(3, 2)
nexttile(3)
TFS_data =[TFS_l_playing_avg(:, :, 1) TFS_l_rebound_avg(:, :, 1) TFS_l_rest_avg(:, :, 1)];
pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp
xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
xtickangle(0)
xlim([0 11250])

clim([-cbar_lim cbar_lim])
title({'';'Left Hemisphere low scores'})
ylabel('Frequency (Hz)')
colormap(slanCM('coolwarm'))
% colorbar()
xline(length(TFS_l_playing_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
xline(length(TFS_l_playing_avg(:, :, 1)) + length(TFS_l_rebound_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
drawnow
ylim(ylim_TFS);


% subplot(324)
nexttile(4)
TFS_data =[TFS_r_playing_avg(:, :, 1) TFS_r_rebound_avg(:, :, 1) TFS_r_rest_avg(:, :, 1)];
pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp

xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
xtickangle(0)
xlim([0 11250])

clim([-cbar_lim cbar_lim])
title({'';'Right Hemisphere Low scores'})
ylabel('Frequency (Hz)')
colormap(slanCM('coolwarm'))
% colorbar()
xline(length(TFS_r_playing_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
xline(length(TFS_r_playing_avg(:, :, 1)) + length(TFS_r_rebound_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
drawnow
ylim(ylim_TFS);

% subplot(325)
nexttile(5)
TFS_data = [TFS_l_playing_avg(:, :, 2) TFS_l_rebound_avg(:, :, 2) TFS_l_rest_avg(:, :, 2)];
pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp

xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
xtickangle(0)
xlim([0 11250])

clim([-cbar_lim cbar_lim])
title(['Left Hemisphere High scores'])
ylabel('Frequency (Hz)')
xlabel('Time (s)');
colormap(slanCM('coolwarm'))
% colorbar()
xline(length(TFS_l_playing_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
xline(length(TFS_l_playing_avg(:, :, 1)) + length(TFS_l_rebound_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
drawnow
ylim(ylim_TFS);

% subplot(326)
nexttile(6)
TFS_data = [TFS_r_playing_avg(:, :, 2) TFS_r_rebound_avg(:, :, 2) TFS_r_rest_avg(:, :, 2)];
pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp

xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
xtickangle(0)
xlim([0 11250])

clim([-cbar_lim cbar_lim])
title(['Right Hemisphere High scores'])
ylabel('Frequency (Hz)')
xlabel('Time (s)');
colormap(slanCM('coolwarm'))
% colorbar()
xline(length(TFS_r_playing_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
xline(length(TFS_r_playing_avg(:, :, 1)) + length(TFS_r_rebound_avg(:, :, 1)), 'LineWidth', 5, 'HandleVisibility', 'off')
drawnow
ylim(ylim_TFS);



%% compare envelopes for before and High scores

% subplot(321)
nexttile(1)
hold on
astd = std(H_VE_l_playing_avg_sub(:,:,1),[],2); % to get std shading
fill([1:length(trial_time_playing) fliplr(1:length(trial_time_playing))],[H_VE_l_playing_avg(:, 1)+astd; flipud(H_VE_l_playing_avg(:, 1) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
astd = std(H_VE_l_rebound_avg_sub(:,:,1),[],2); % to get std shading
fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_l_rebound_avg(:, 1)+astd; flipud(H_VE_l_rebound_avg(:, 1) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

astd = std(H_VE_l_rest_avg_sub(:,:,1),[],2); % to get std shading
fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound)  (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_l_rest_avg(:, 1)+astd; flipud(H_VE_l_rest_avg(:, 1) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

astd = std(H_VE_l_playing_avg_sub(:,:,2),[],2); % to get std shading
fill([1:length(trial_time_playing) length(trial_time_playing):-1:1],[H_VE_l_playing_avg(:, 2)+astd; flipud(H_VE_l_playing_avg(:, 2) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
astd = std(H_VE_l_rebound_avg_sub(:,:,2),[],2); % to get std shading
fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_l_rebound_avg(:, 2)+astd; flipud(H_VE_l_rebound_avg(:, 2) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
astd = std(H_VE_l_rest_avg_sub(:,:,2),[],2); % to get std shading
fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound) (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_l_rest_avg(:, 2)+astd; flipud(H_VE_l_rest_avg(:, 2) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

plot([H_VE_l_playing_avg(:, 1); H_VE_l_rebound_avg(:, 1); H_VE_l_rest_avg(:, 1)],'b', 'DisplayName', 'Low scores')
plot([H_VE_l_playing_avg(:, 2); H_VE_l_rebound_avg(:, 2); H_VE_l_rest_avg(:, 2)],'r', 'DisplayName', 'High scores')

xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
xtickangle(0)
xlim([0 11250])

xline(length(TFS_r_playing_avg(:,:,1)), 'LineWidth', 5, 'HandleVisibility', 'off')
xline(length(TFS_r_playing_avg(:,:,1))+length(TFS_r_rebound_avg(:,:,1)), 'LineWidth', 5, 'HandleVisibility', 'off')
ylabel('Relative change (%)')
drawnow
ylim([ylim_vals_env])
title('Left hemisphere')
legend()


% subplot(322)
nexttile(2)
hold on
astd = std(H_VE_r_playing_avg_sub(:,:,1),[],2); % to get std shading
fill([1:length(trial_time_playing) fliplr(1:length(trial_time_playing))],[H_VE_r_playing_avg(:, 1)+astd; flipud(H_VE_r_playing_avg(:, 1) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
astd = std(H_VE_r_rebound_avg_sub(:,:,1),[],2); % to get std shading
fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_r_rebound_avg(:, 1)+astd; flipud(H_VE_r_rebound_avg(:, 1) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
astd = std(H_VE_r_rest_avg_sub(:,:,1),[],2); % to get std shading
fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound)  (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_r_rest_avg(:, 1)+astd; flipud(H_VE_r_rest_avg(:, 1) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

astd = std(H_VE_r_playing_avg_sub(:,:,2),[],2); % to get std shading
fill([1:length(trial_time_playing) length(trial_time_playing):-1:1],[H_VE_r_playing_avg(:, 2)+astd; flipud(H_VE_r_playing_avg(:, 2) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
astd = std(H_VE_r_rebound_avg_sub(:,:,2),[],2); % to get std shading
fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_r_rebound_avg(:, 2)+astd; flipud(H_VE_r_rebound_avg(:, 2) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
astd = std(H_VE_r_rest_avg_sub(:,:,2),[],2); % to get std shading
fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound) (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_r_rest_avg(:, 2)+astd; flipud(H_VE_r_rest_avg(:, 2) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

plot([H_VE_r_playing_avg(:, 1); H_VE_r_rebound_avg(:, 1); H_VE_r_rest_avg(:, 1)],'b', 'DisplayName', 'Low scores')
plot([H_VE_r_playing_avg(:, 2); H_VE_r_rebound_avg(:, 2); H_VE_r_rest_avg(:, 2)],'r', 'DisplayName', 'High scores')

xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
xtickangle(0)
xlim([0 11250])

xline(length(TFS_r_playing_avg(:,:,1)), 'LineWidth', 5, 'HandleVisibility', 'off')
xline(length(TFS_r_playing_avg(:,:,1))+length(TFS_r_rebound_avg(:,:,1)), 'LineWidth', 5, 'HandleVisibility', 'off')
ylabel('Relative change (%)')
drawnow
ylim([ylim_vals_env])
title('Right hemisphere')
legend()

toc(run_all_timing);

fh = findall(0,'Type','Figure');
set( findall(fh, '-property', 'fontsize'), 'fontsize', 20)