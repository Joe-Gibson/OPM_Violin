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

load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\MARS\MARS_cortical_Key.mat')

motor_inds = [20, 22, 23, 24, 26, 27, 29, 30]; motor_inds = [motor_inds motor_inds+41];

location_names = strrep(strrep(strrep(strrep(strrep(MARS_cortical_Key.names(motor_inds), ' Cortex', ''), 'Left ', ''), 'Right ', ''), '8',  ''), ' ', '_');

runs = ["run_1"; "run_2"];


sub_nums = [1:22];
sub_nums(21) = []; % Lots of noisy trials
sub_nums(16) = []; % Exclude sub 16 - wrong triggers
sub_nums(1) = []; % Exclude sub 1 - wrong triggers

% sub_nums = 2;

n_subjects = length(sub_nums);
modulation_l_calc = zeros(n_subjects, 2);
modulation_r_calc = zeros(n_subjects, 2);

f = 375;
duration = 10;
trig_offset_rest = 10; % Look n seconds after trigger
trial_time_rest = [1:duration*f]./f + trig_offset_rest;
trig_offset_rebound = -5; % Look n seconds after trigger
trial_time_rebound = [1:duration*f]./f + trig_offset_rebound;
trig_offset_playing = 0; % Look n seconds before trigger
trial_time_playing = [1:duration*f]./f + trig_offset_playing;


n_to_load = n_subjects*length(motor_inds)/2*2;
load_val = 1;

ylim_vals_env = [-30 50];


%% Choose files
file_to_load = ['D:\OneDrive - The University of Nottingham\Documents\Data\Violin\V2\' ...
        'BIDS\Derivatives\BF_VEs\MARS\avg_TFS_H_VE_mat_relative_' num2str(hp) '_' num2str(lp) '_lesson.mat'];

if ~isfile(file_to_load)
    % ft_progress('init', 'etf', 'Loading data...')      % ascii progress bar
    wait_fig = waitbar(0,'Please wait...');
    for sub = 1:n_subjects
        sub_number = sub_nums(sub);
        subject = ['Sub-' num2str(sub_number, '%03.f')];
        % disp(subject)
        for run_no = 1:2
            run = ['run_' num2str(run_no)];
            % disp(run)
            for mars_loc = 1:length(motor_inds)/2
                mars_loc_l = motor_inds(mars_loc);
                mars_loc_r = motor_inds(mars_loc) + 41;

                % disp(['Location ' num2str(loc)])

                % ft_progress(load_val/n_to_load, 'Loading part %d of %d', load_val, n_to_load)
                waitbar(load_val/n_to_load,wait_fig,sprintf('Loading part %d of %d', load_val, n_to_load));

                VE = load([bids_path,'derivatives\BF_VEs\MARS\' subject '\' run ...
                    '_3trigs' freq_savename 'unaveraged_location_'  ...
                    num2str(mars_loc_l) '.mat']);

                TFS_l_playing_all.(runs{run_no}).(location_names(mars_loc))(sub, :, :)  = zeros(1, 106, duration*f);
                TFS_r_playing_all.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = zeros(1, 106, duration*f);
                H_VE_l_playing_all.(runs(run_no)).(location_names(mars_loc))(sub, :)  = zeros(1, duration*f);
                H_VE_r_playing_all.(runs(run_no)).(location_names(mars_loc))(sub, :)  = zeros(1, duration*f);
                TFS_l_rebound_all.(runs{run_no}).(location_names(mars_loc))(sub, :, :) = zeros(1, 106, duration*f);
                TFS_r_rebound_all.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = zeros(1, 106, duration*f);
                H_VE_l_rebound_all.(runs(run_no)).(location_names(mars_loc))(sub, :)  = zeros(1, duration*f);
                H_VE_r_rebound_all.(runs(run_no)).(location_names(mars_loc))(sub, :)  = zeros(1, duration*f);
                TFS_l_rest_all.(runs{run_no}).(location_names(mars_loc))(sub, :, :) = zeros(1, 106, duration*f);
                TFS_r_rest_all.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = zeros(1, 106, duration*f);
                H_VE_l_rest_all.(runs(run_no)).(location_names(mars_loc))(sub, :)  = zeros(1, duration*f);
                H_VE_r_rest_all.(runs(run_no)).(location_names(mars_loc))(sub, :)  = zeros(1, duration*f);

                % Average across trials
                TFS_l_playing.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = mean(VE.TFS_l_playing_all,3);
                TFS_r_playing.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = mean(VE.TFS_r_playing_all,3);
                TFS_l_rebound.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = mean(VE.TFS_l_rebound_all,3);
                TFS_r_rebound.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = mean(VE.TFS_r_rebound_all,3);
                TFS_l_rest.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = mean(VE.TFS_l_rest_all,3);
                TFS_r_rest.(runs(run_no)).(location_names(mars_loc))(sub, :, :) = mean(VE.TFS_r_rest_all,3);

                % Envelope from filter
                VE.H_VE_l_playing_all = VE.Hilbert_envelope_l_playing_trials;
                VE.H_VE_l_rebound_all = VE.Hilbert_envelope_l_rebound_trials;
                VE.H_VE_l_rest_all = VE.Hilbert_envelope_l_rest_trials;
                VE.H_VE_r_playing_all = VE.Hilbert_envelope_r_playing_trials;
                VE.H_VE_r_rebound_all = VE.Hilbert_envelope_r_rebound_trials;
                VE.H_VE_r_rest_all = VE.Hilbert_envelope_r_rest_trials;

                % Relative
                meanrest_l_H_VE =  repmat(squeeze(mean(VE.H_VE_l_rest_all, 1)), [size(VE.H_VE_l_playing_all, 1), 1]);
                VE.H_VE_l_playing_all = (VE.H_VE_l_playing_all-meanrest_l_H_VE)./meanrest_l_H_VE*100;
                VE.H_VE_l_rebound_all = (VE.H_VE_l_rebound_all - meanrest_l_H_VE)./meanrest_l_H_VE*100;
                VE.H_VE_l_rest_all = (VE.H_VE_l_rest_all - meanrest_l_H_VE)./meanrest_l_H_VE*100;

                meanrest_r_H_VE =  repmat(squeeze(mean(VE.H_VE_r_rest_all, 1)), [size(VE.H_VE_r_playing_all, 1), 1]);
                VE.H_VE_r_playing_all = (VE.H_VE_r_playing_all-meanrest_r_H_VE)./meanrest_r_H_VE*100;
                VE.H_VE_r_rebound_all = (VE.H_VE_r_rebound_all - meanrest_r_H_VE)./meanrest_r_H_VE*100;
                VE.H_VE_r_rest_all = (VE.H_VE_r_rest_all - meanrest_r_H_VE)./meanrest_r_H_VE*100;

                H_VE_l_playing.(runs(run_no)).(location_names(mars_loc))(sub, :) = mean(VE.H_VE_l_playing_all, 2);
                H_VE_r_playing.(runs(run_no)).(location_names(mars_loc))(sub, :) = mean(VE.H_VE_r_playing_all, 2);
                H_VE_l_rebound.(runs(run_no)).(location_names(mars_loc))(sub, :) = mean(VE.H_VE_l_rebound_all, 2);
                H_VE_r_rebound.(runs(run_no)).(location_names(mars_loc))(sub, :) = mean(VE.H_VE_r_rebound_all, 2);
                H_VE_l_rest.(runs(run_no)).(location_names(mars_loc))(sub, :) = mean(VE.H_VE_l_rest_all, 2);
                H_VE_r_rest.(runs(run_no)).(location_names(mars_loc))(sub, :) = mean(VE.H_VE_r_rest_all, 2);
                load_val = load_val+1;
            end
            clear meanrest_l_H_VE meanrest_r_H_VE meanrest_l_TFS meanrest_r_TFS ...
                VE
        end
    end
    % ft_progress('close')
    close(wait_fig)

    for run_no = 1:2
        for mars_loc = 1:length(motor_inds)/2
            % disp(['Location ' num2str(loc)])
            TFS_l_playing_avg.(runs(run_no)).(location_names(mars_loc)) = squeeze(mean(TFS_l_playing.(runs(run_no)).(location_names(mars_loc)), 1));
            TFS_r_playing_avg.(runs(run_no)).(location_names(mars_loc)) = squeeze(mean(TFS_r_playing.(runs(run_no)).(location_names(mars_loc)), 1));
            H_VE_l_playing_avg.(runs(run_no)).(location_names(mars_loc)) =  mean(H_VE_l_playing.(runs(run_no)).(location_names(mars_loc)), 1);
            H_VE_r_playing_avg.(runs(run_no)).(location_names(mars_loc)) =  mean(H_VE_r_playing.(runs(run_no)).(location_names(mars_loc)), 1);
            TFS_l_rebound_avg.(runs(run_no)).(location_names(mars_loc)) = squeeze(mean(TFS_l_rebound.(runs(run_no)).(location_names(mars_loc)), 1));
            TFS_r_rebound_avg.(runs(run_no)).(location_names(mars_loc)) = squeeze(mean(TFS_r_rebound.(runs(run_no)).(location_names(mars_loc)), 1));
            H_VE_l_rebound_avg.(runs(run_no)).(location_names(mars_loc)) =  mean(H_VE_l_rebound.(runs(run_no)).(location_names(mars_loc)), 1);
            H_VE_r_rebound_avg.(runs(run_no)).(location_names(mars_loc)) =  mean(H_VE_r_rebound.(runs(run_no)).(location_names(mars_loc)), 1);
            TFS_l_rest_avg.(runs(run_no)).(location_names(mars_loc)) = squeeze(mean(TFS_l_rest.(runs(run_no)).(location_names(mars_loc)), 1));
            TFS_r_rest_avg.(runs(run_no)).(location_names(mars_loc)) = squeeze(mean(TFS_r_rest.(runs(run_no)).(location_names(mars_loc)), 1));
            H_VE_l_rest_avg.(runs(run_no)).(location_names(mars_loc)) =  mean(H_VE_l_rest.(runs(run_no)).(location_names(mars_loc)), 1);
            H_VE_r_rest_avg.(runs(run_no)).(location_names(mars_loc)) =  mean(H_VE_r_rest.(runs(run_no)).(location_names(mars_loc)), 1);
        end
    end
    save(file_to_load, "TFS_l_playing_avg", ...
        "TFS_r_playing_avg", "H_VE_l_playing_avg", "H_VE_r_playing_avg", "TFS_l_rebound_avg", "TFS_r_rebound_avg", "H_VE_l_rebound_avg", ...
        "H_VE_r_rebound_avg", "TFS_l_rest_avg", "TFS_r_rest_avg", "H_VE_l_rest_avg", "H_VE_r_rest_avg", "TFS_l_playing", ...
        "TFS_r_playing", "H_VE_l_playing", "H_VE_r_playing", "TFS_l_rebound", "TFS_r_rebound", "H_VE_l_rebound", ...
        "H_VE_r_rebound", "TFS_l_rest", "TFS_r_rest", "H_VE_l_rest", "H_VE_r_rest")
else
    load(file_to_load)
end


%%
for mars_loc = 1:length(motor_inds)/2
    % disp(['Location ' num2str(loc)])

    cbar_lim = 0.3; % colour bar limit (relative change)

    if size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc))(:, :, 1), 1) == 110
        highpass = 1:110;
        lowpass = 2.5:111.5;
        fre = highpass + ((lowpass - highpass)./2);
    elseif size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc))(:, :, 1), 1) == 106
        highpass = 5:110;
        lowpass = 6.5:111.5;
        fre = highpass + ((lowpass - highpass)./2);
    end

    xlim_vals = [-10 25];

    figure(mars_loc)
    tiledlayout(3, 2)
    nexttile(3)
    TFS_data = [TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)) TFS_l_rebound_avg.(runs(1)).(location_names(mars_loc)) TFS_l_rest_avg.(runs(1)).(location_names(mars_loc))];
    pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp
    xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
    xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
    xtickangle(0)
    xlim([0 11250])
    clim([-cbar_lim cbar_lim])
    title({'';'Left Hemisphere before lesson'})
    ylabel('Frequency (Hz)')
    colormap(slanCM('coolwarm'))
    % colorbar()
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2) + size(TFS_l_rebound_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    drawnow
    ylim(ylim_TFS);
   

    figure(mars_loc)
    nexttile(4)
    TFS_data = [TFS_r_playing_avg.(runs(1)).(location_names(mars_loc)) TFS_r_rebound_avg.(runs(1)).(location_names(mars_loc)) TFS_r_rest_avg.(runs(1)).(location_names(mars_loc))];
    pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp
    xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
    xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
    xtickangle(0)
    xlim([0 11250])
    xtickangle(0)
    clim([-cbar_lim cbar_lim])
    title({'';'Right Hemisphere before lesson'})
    ylabel('Frequency (Hz)')
    colormap(slanCM('coolwarm'))
    % colorbar()
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2) + size(TFS_l_rebound_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    drawnow
    ylim(ylim_TFS);

    figure(mars_loc)
    nexttile(5)
    TFS_data = [TFS_l_playing_avg.(runs(2)).(location_names(mars_loc)) TFS_l_rebound_avg.(runs(2)).(location_names(mars_loc)) TFS_l_rest_avg.(runs(2)).(location_names(mars_loc))];
    pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp
    xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
    xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
    xtickangle(0)
    xlim([0 11250])
    xtickangle(0)
    clim([-cbar_lim cbar_lim])
    title(['Left Hemisphere after lesson'])
    ylabel('Frequency (Hz)')
    xlabel('Time (s)');
    colormap(slanCM('coolwarm'))
    % colorbar()
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2) + size(TFS_l_rebound_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    drawnow
    ylim(ylim_TFS);
   

    figure(mars_loc)
    nexttile(6)
    TFS_data = [TFS_r_playing_avg.(runs(2)).(location_names(mars_loc)) TFS_r_rebound_avg.(runs(2)).(location_names(mars_loc)) TFS_r_rest_avg.(runs(2)).(location_names(mars_loc))];
    pcolor(1:size(TFS_data, 2), highpass, TFS_data);shading interp
    xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
    xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
    xtickangle(0)
    xlim([0 11250])
    xtickangle(0)
    clim([-cbar_lim cbar_lim])
    title(['Right Hemisphere after lesson'])
    ylabel('Frequency (Hz)')
    xlabel('Time (s)');
    colormap(slanCM('coolwarm'))
    % colorbar()
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2) + size(TFS_l_rebound_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    drawnow
    ylim(ylim_TFS);
   

    sgtitle(strrep(location_names(mars_loc), '_', ' '))


    %% compare envelopes for before and after lesson
    figure(mars_loc)
    nexttile(1)
    hold on
    astd = std(H_VE_l_playing.(runs(1)).(location_names(mars_loc)),[],1); % to get std shading
    fill([1:length(trial_time_playing) fliplr(1:length(trial_time_playing))],[H_VE_l_playing_avg.(runs(1)).(location_names(mars_loc))+astd fliplr(H_VE_l_playing_avg.(runs(1)).(location_names(mars_loc)) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_l_rebound.(runs(1)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_l_rebound_avg.(runs(1)).(location_names(mars_loc))+astd fliplr(H_VE_l_rebound_avg.(runs(1)).(location_names(mars_loc)) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_l_rest.(runs(1)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound) (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_l_rest_avg.(runs(1)).(location_names(mars_loc))+astd fliplr(H_VE_l_rest_avg.(runs(1)).(location_names(mars_loc)) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

    astd = std(H_VE_l_playing.(runs(2)).(location_names(mars_loc)),[],1); % to get std shading
    fill([1:length(trial_time_playing) length(trial_time_playing):-1:1],[H_VE_l_playing_avg.(runs(2)).(location_names(mars_loc))+astd fliplr(H_VE_l_playing_avg.(runs(2)).(location_names(mars_loc)) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_l_rebound.(runs(2)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_l_rebound_avg.(runs(2)).(location_names(mars_loc))+astd fliplr(H_VE_l_rebound_avg.(runs(2)).(location_names(mars_loc)) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_l_rest.(runs(2)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound) (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_l_rest_avg.(runs(2)).(location_names(mars_loc))+astd fliplr(H_VE_l_rest_avg.(runs(2)).(location_names(mars_loc)) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

    plot([H_VE_l_playing_avg.(runs(1)).(location_names(mars_loc)) H_VE_l_rebound_avg.(runs(1)).(location_names(mars_loc)) H_VE_l_rest_avg.(runs(1)).(location_names(mars_loc))],'b', 'DisplayName', 'Before lesson')
    plot([H_VE_l_playing_avg.(runs(2)).(location_names(mars_loc)) H_VE_l_rebound_avg.(runs(2)).(location_names(mars_loc)) H_VE_l_rest_avg.(runs(2)).(location_names(mars_loc))],'r', 'DisplayName', 'After lesson')

    xlim([0 11250])
    xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
    xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
    xtickangle(0)
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2) + size(TFS_l_rebound_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    ylabel('Relative change (%)')
    drawnow
    ylim([ylim_vals_env])
    title('Left hemisphere')
    legend()



    figure(mars_loc)
    % subplot(322)
    nexttile(2)
    hold on
    astd = std(H_VE_r_playing.(runs(1)).(location_names(mars_loc)),[],1); % to get std shading
    fill([1:length(trial_time_playing) fliplr(1:length(trial_time_playing))],[H_VE_r_playing_avg.(runs(1)).(location_names(mars_loc))+astd fliplr(H_VE_r_playing_avg.(runs(1)).(location_names(mars_loc)) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_r_rebound.(runs(1)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_r_rebound_avg.(runs(1)).(location_names(mars_loc))+astd fliplr(H_VE_r_rebound_avg.(runs(1)).(location_names(mars_loc)) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_r_rest.(runs(1)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound) (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_r_rest_avg.(runs(1)).(location_names(mars_loc))+astd fliplr(H_VE_r_rest_avg.(runs(1)).(location_names(mars_loc)) -astd)],'b', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

    astd = std(H_VE_r_playing.(runs(2)).(location_names(mars_loc)),[],1); % to get std shading
    fill([1:length(trial_time_playing) length(trial_time_playing):-1:1],[H_VE_r_playing_avg.(runs(2)).(location_names(mars_loc))+astd fliplr(H_VE_r_playing_avg.(runs(2)).(location_names(mars_loc)) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_r_rebound.(runs(2)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rebound)) + length(trial_time_playing)  (length(trial_time_rebound):-1:1) + length(trial_time_playing)],[H_VE_r_rebound_avg.(runs(2)).(location_names(mars_loc))+astd fliplr(H_VE_r_rebound_avg.(runs(2)).(location_names(mars_loc)) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');
    astd = std(H_VE_r_rest.(runs(2)).(location_names(mars_loc)),[],1); % to get std shading
    fill([(1:length(trial_time_rest)) + length(trial_time_playing) + length(trial_time_rebound) (length(trial_time_rest):-1:1) + length(trial_time_playing) + length(trial_time_rebound)],[H_VE_r_rest_avg.(runs(2)).(location_names(mars_loc))+astd fliplr(H_VE_r_rest_avg.(runs(2)).(location_names(mars_loc)) -astd)],'r', 'FaceAlpha', 0.5,'linestyle','none', 'HandleVisibility', 'off');

    plot([H_VE_r_playing_avg.(runs(1)).(location_names(mars_loc)) H_VE_r_rebound_avg.(runs(1)).(location_names(mars_loc)) H_VE_r_rest_avg.(runs(1)).(location_names(mars_loc))],'b', 'DisplayName', 'Before lesson')
    plot([H_VE_r_playing_avg.(runs(2)).(location_names(mars_loc)) H_VE_r_rebound_avg.(runs(2)).(location_names(mars_loc)) H_VE_r_rest_avg.(runs(2)).(location_names(mars_loc))],'r', 'DisplayName', 'After lesson')

    xlim([0 11250])
    xticks([250 1875 3500 4000 5625 7250 7750 9375 11000])
    xticklabels(["0" "Playing" "10" "0" "Rebound" "10" "0" "Rest" "10"])
    xtickangle(0)
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    xline(size(TFS_l_playing_avg.(runs(1)).(location_names(mars_loc)), 2) + size(TFS_l_rebound_avg.(runs(1)).(location_names(mars_loc)), 2), 'LineWidth', 5, 'HandleVisibility', 'off')
    ylabel('Relative change (%)')
    drawnow
    ylim([ylim_vals_env])
    title('Right hemisphere')
    legend()

    sgtitle(strrep(location_names(mars_loc), '_', ' '))

end

toc(run_all_timing);

fh = findall(0,'Type','Figure');
set( findall(fh, '-property', 'fontsize'), 'fontsize', 20)