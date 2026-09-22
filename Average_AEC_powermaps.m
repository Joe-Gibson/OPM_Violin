clear all
close all
restoredefaultpath
set(0,'DefaultFigureWindowStyle','docked')

% If on new desktop, change onedrive location from C: to D:
if strcmp(getenv('COMPUTERNAME'), 'DUIP78552')
    Onedrive_path = 'D:\OneDrive - The University of Nottingham\';
else
    Onedrive_path = 'C:\Users\ppyjg12\OneDrive - The University of Nottingham\';
end

addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\Beamformer"))
addpath(genpath("C:\Users\ppyjg12\Documents\Projects\Violin\V2"))
if isempty(which('ft_defaults'))
    addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\fieldtrip-20231113"))
    ft_defaults
end

%%
bids_path = [Onedrive_path 'Documents\Data\Violin\V2\BIDS\'];

Nlocs = 82;
plot_percentile = 0.95;

wins = [13 30];
band_names = "beta";

exp_string = ["before_lesson", "after_lesson"];

load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\BrainPlots\all_voxels_MARS.mat')
load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\MARS\MARS_cortical_Key.mat')
load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\BrainPlots\mars_labels_20484.mat')
vars.MARS_atlas_all_vox = MARS_atlas_all_vox;
vars.MARS_cortical_Key = MARS_cortical_Key;
vars.mars_labels_20484 = mars_labels_20484;

mars_labels = MARS_cortical_Key.label;
location_names = strrep(strrep(strrep(MARS_cortical_Key.names, ' Cortex', ''), '8',  ''), ' ', '_');


% What to average
do_average_AEC = 1;
do_average_powermap = 1;
do_average_AEC_powermap = 1;

% what triggers
use_trial_events = 0;
use_rest_events = 1;
use_rebound_events = 0;

do_mars_positions_plot = 0;
do_averaged_AEC_calc = 1;
significance_thresh = 0.05;
do_save = 0;

if sum([use_trial_events, use_rest_events use_rebound_events]) > 1
    error('Only use one trigger condition')
end

% what to average over
average_all = 0;
average_exp = 1;
average_score = 0;

if sum([average_all, average_exp, average_score]) > 1
    error('You probably only want to use one average at once')
end

sub_nums = [1:22];
sub_nums(21) = []; % Lots of noisy trials
sub_nums(16) = []; % Exclude sub 16 - wrong triggers
sub_nums(1) = []; % Exclude sub 1 - wrong triggers

n_subjects = length(sub_nums);

if use_trial_events
    trigger_string = 'trial';
    trigger_title_string = 'playing';
elseif use_rest_events
    trigger_string = 'rest';
    trigger_title_string = 'rest';
elseif use_rebound_events
    trigger_string = 'rebound';
    trigger_title_string = 'rebound';
end

if do_averaged_AEC_calc
    avg_string = 'Averaged\';
else
    avg_string = '';
end

all_AEC = struct();
all_powermap = struct();
average_node = struct();
for n_band = 1:size(wins, 1)
    if average_all | average_exp
        all_AEC.(band_names{n_band}) = NaN(Nlocs, Nlocs, n_subjects, 2, 30);
        all_powermap.(band_names{n_band}) = NaN(Nlocs,  n_subjects, 2, 30);
    elseif average_score
        all_AEC.(band_names{n_band}) = NaN(Nlocs, Nlocs, n_subjects, 60);
        all_powermap.(band_names{n_band}) = NaN(Nlocs,  n_subjects, 60);
    end
end

%Loop for each frequency band
for n_band = 1:size(wins, 1)
    hp = wins(n_band,1);
    lp = wins(n_band,2);

    for sub = 1:n_subjects
        sub_number = sub_nums(sub);
        subject = ['Sub-' num2str(sub_number, '%03.f')];
        ntrials_run1 = 0;

        for run_no = 1:2
            run = ['run_' num2str(run_no)];

            % clearvars -except All_files n_all_files Onedrive_path vars ...
            % run_all_time bids_path F subject sub run run_no
            % close all

            filename = run;

            path_VEs = [bids_path,'Derivatives\VEs\' avg_string subject '\'];
            path_AEC = [bids_path,'Derivatives\AEC\' avg_string subject '\'];
            files_AEC = [filename,'_AEC'];
            path_figures_AEC = [bids_path,'Derivatives\figures\', ...
                avg_string subject '\AEC\'];
            path_power_map = [bids_path,'Derivatives\power_map\', ...
                avg_string subject,'\'];
            files_power_map = [filename,'_power_map'];
            path_figures_power_map = [bids_path,'Derivatives\figures\', ...
                avg_string subject '\power_map\'];

            % Load power map
            load([path_power_map files_power_map '_' num2str(hp) '_to_'  ...
                num2str(lp) '_' trigger_string '_all_trials.mat'])

            % Load AEC
            load([path_AEC files_AEC '_' num2str(hp) '_to_' num2str(lp) ...
                '_' trigger_string '_all_trials.mat'])


            if average_all | average_exp
                all_AEC.(band_names{n_band})(:, :, sub, run_no, 1:size(AEC_mat, 3)) = (AEC_mat + pagetranspose(AEC_mat))./2;
                all_powermap.(band_names{n_band})(:, sub, run_no, 1:size(power_output_relative_mat, 2)) = power_output_relative_mat;
            elseif average_score
                all_AEC.(band_names{n_band})(:, :, sub, ntrials_run1 +(1:size(AEC_mat, 3))) = (AEC_mat + pagetranspose(AEC_mat))./2;
                all_powermap.(band_names{n_band})(:, sub, ntrials_run1 +(1:size(power_output_relative_mat, 2))) = power_output_relative_mat;
                ntrials_run1 = size(AEC_mat, 3);
            end
        end
    end
end

if average_all
    fnames = fieldnames(all_AEC);
    for i = 1:length(fnames)
        case_names = "All";
        average_AEC_subs.(fnames{i}) = mean(all_AEC.(fnames{i}), [4 5], 'omitnan');
        average_powermap_subs.(fnames{i}) = mean(all_powermap.(fnames{i}), [3 4], 'omitnan');
        average_AEC.(fnames{i}) = mean(average_AEC_subs.(fnames{i}), 3, 'omitnan');
        average_powermap.(fnames{i}) = mean(average_powermap_subs.(fnames{i}), 2, 'omitnan');
    end
elseif average_exp %Averaging to compare before lesson and after lesson
    fnames = fieldnames(all_AEC);
    for i = 1:length(fnames)
        for run_no = 1:2
            average_AEC_subs.(fnames{i})(:,:,:, run_no) = squeeze(mean(all_AEC.(fnames{i})(:,:,:,run_no,:), 5, 'omitnan'));
            average_powermap_subs.(fnames{i})(:,:,run_no) = squeeze(mean(all_powermap.(fnames{i})(:,:,run_no,:), 4, 'omitnan'));
            average_AEC.(fnames{i})(:,:,run_no) = squeeze(mean(average_AEC_subs.(fnames{i})(:,:,:, run_no), 3, 'omitnan'));
            average_powermap.(fnames{i})(:,run_no) = squeeze(mean( average_powermap_subs.(fnames{i})(:,:,run_no), 2, 'omitnan'));
        end
        case_names = char(["Before Lesson"; "After Lesson"]);
        average_save_name = 'Lesson';
    end
elseif average_score   % AVerage to compare high and low scores
    case_names = char(["Low scores"; "High scores"]);
    average_save_name = 'Scores';
    load([bids_path 'Derivatives\score_orders.mat'])%, "scores_sorted", "trial_score_order");
    trial_score_orders = NaN(n_subjects, 60);

    for sub = n_subjects:-1:1
        sub_number = sub_nums(sub);
        trial_score_orders(sub, 1:length(trial_score_order.(['sub_' num2str(sub_number)]))) = trial_score_order.(['sub_' num2str(sub_number)]);
        ntrials(sub) = length(trial_score_order.(['sub_' num2str(sub_number)]));
    end

    ntrials_third = floor(2*ntrials/5);

    fnames = fieldnames(all_AEC);
    for i = 1:length(fnames)

        for sub = 1:n_subjects
            average_AEC_subs.(fnames{i})(:, :, sub, 2) = mean(all_AEC.(fnames{i})(:,:,sub,trial_score_orders(sub, ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 4, 'omitnan');
            average_AEC_subs.(fnames{i})(:, :, sub, 1) = mean(all_AEC.(fnames{i})(:,:,sub,trial_score_orders(sub, 1:ntrials_third(sub))), 4, 'omitnan');
            average_powermap_subs.(fnames{i})(:, sub, 2) = mean(all_powermap.(fnames{i})(:,sub,trial_score_orders(sub, ntrials(sub) - ntrials_third(sub) + 1:ntrials(sub))), 3, 'omitnan');
            average_powermap_subs.(fnames{i})(:, sub, 1) = mean(all_powermap.(fnames{i})(:,sub,trial_score_orders(sub, 1:ntrials_third(sub))), 3, 'omitnan');
        end

        average_AEC.(fnames{i}) = squeeze(mean(average_AEC_subs.(fnames{i}), 3, 'omitnan'));
        average_powermap.(fnames{i}) = squeeze(mean(average_powermap_subs.(fnames{i}), 2, 'omitnan'));
    end

end

Ncases = size(average_AEC.(fnames{1}), 3); % number of 'runs averaged over, if
% average all 1 run, or if everage_Exp  then split into 2, before and after experiment

for run = 1:Ncases
    run_string = case_names(run, :);
    if do_average_powermap
        fnames = fieldnames(average_powermap);
        power_output_max = struct();
        for i = 1:length(fnames)
            % power_output_max.(fnames{i}) = [min(abs(average_powermap.(fnames{i})),[],'all') max(abs(average_powermap.(fnames{i})),[],'all')];

            power_output_max.(fnames{i}) = [0.17 0.24];
        end

        if ~average_all
            Nsubplots = Ncases + 1;
            Nsubplots_AEC = Ncases; % No resiudal plot for AEC
        else
            Nsubplots = Ncases;
            Nsubplots_AEC = Ncases;
        end


        for i = 1:length(fnames)
            powermap_fig = figure(i);
            powermap_fig.Name = ['Average powermap ' fnames{i}];
            subplot(1, Nsubplots, run)
            PaintBrodmannAreas_mars_1view(average_powermap.(fnames{i})(:, run),Nlocs, 256,  ...
                power_output_max.(fnames{i}), '', [], vars)
            cb = colorbar;
            cb.Label.String = 'Relative power (A.U.)';

            title(run_string)
            sgtitle(['Average ' trigger_title_string ' ' strrep(fnames{i}, '_', '  ') ' band power'])
            axis equal
            axis([-68.7943   69.7716 -106.4590   70.1386])

            for region = size(average_powermap_subs.(fnames{i}), 1):-1:1
                powermap_significance.(fnames{i})(region) = signrank(average_powermap_subs.(fnames{i})(region,:, 1), average_powermap_subs.(fnames{i})(region,:, 2));
            end
        end
    end

    if do_average_AEC
        fnames = fieldnames(average_AEC);
        for i = 1:length(fnames)
            AEC_max.(fnames{i}) = max(abs(average_AEC.(fnames{i})),[],'all', 'omitnan');

            thresh = plot_percentile;

            if thresh < 1; thresh = thresh*100; end

            for ii = 1:size(average_AEC.(fnames{i}), 3)
                C = average_AEC.(fnames{i})(:, :, ii);

                C(eye(size(C))==1) = 0;

                limit.(fnames{i})(ii) = prctile(C(triu(ones(size(C)),1)==1),thresh);

                mask = abs(C) >= limit.(fnames{i})(ii);
                C(~mask) = NaN;

                sphereWidths = sum(~isnan(C))*0.5;
                max_sw.(fnames{i})(ii) = max(sphereWidths(:));
            end
            limit.(fnames{i}) = max(limit.(fnames{i}));
            max_sw.(fnames{i}) = max(max_sw.(fnames{i}));
        end

        % AEC_max = 0.055;

        for i = 1:length(fnames)
            AEC_fig = figure(length(fnames)+i);
            AEC_fig.Name = ['Average AEC ' fnames{i}];
            subplot(Nsubplots_AEC, 2, 1 + (run-1)*2)
            imagesc(average_AEC.(fnames{i})(:, :, run));colorbar; clim([-AEC_max.(fnames{i}) AEC_max.(fnames{i})])
            axis square
            title(['AEC matrix for ' strrep(fnames{i}, '_', '  ') ' band ' trigger_title_string ' ' run_string])
            subplot(Nsubplots_AEC, 2, 2 + (run-1)*2)
            go_netviewer_perctl_mars(average_AEC.(fnames{i})(:, :, run), plot_percentile, limit.(fnames{i}), max_sw.(fnames{i}))
            title(run_string)
            axis equal
            axis([-68.7943   69.7716 -106.4590   70.1386])
            drawnow

            for region1 = size(average_AEC_subs.(fnames{i}), 1):-1:1
                for region2 = size(average_AEC_subs.(fnames{i}), 2):-1:1
                    if region1~= region2
                        AEC_significance.(fnames{i})(region1, region2) = signrank(squeeze(average_AEC_subs.(fnames{i})(region1,region2,:, 1)), squeeze(average_AEC_subs.(fnames{i})(region1,region2,:, 2)));
                    else
                        AEC_significance.(fnames{i})(region1, region2) = nan;
                    end
                end
            end

            %% Lukas three view figure
            fh = figure(5.*length(fnames)+i + (run - 1)*length(fnames));
            % set(fh,'WindowStyle','normal')
            fh.Name = ['Average AEC ' fnames{i} ' Lukas'];
            set(fh,'Units','centimeters','Color','w','Renderer','painters');

            fwidth = 30;%for poster
            fheight = 20;%for poster
            % fh.Position([3,4]) = [fwidth,fheight];
            ax_mat = axes;
            imagesc(average_AEC.(fnames{i})(:, :, run)); clim([-AEC_max.(fnames{i}) AEC_max.(fnames{i})])
            cb = colorbar('Location','eastoutside', 'FontSize', 20);
            axis square;
            label_inds = [4, 8, 23, 27, 37, 45, 49, 64, 68, 78];
            yticks(label_inds);
            yticklabels(strrep(strrep(strrep(strrep(MARS_cortical_Key.names(label_inds), ' Cortex', ''), 'Left ', 'L. '), 'Right ', 'R. '), '8',  ''));
            xticks(label_inds);
            xticklabels(strrep(strrep(strrep(strrep(MARS_cortical_Key.names(label_inds), ' Cortex', ''), 'Left ', 'L. '), 'Right ', 'R. '), '8',  ''));
            xtickangle(90)
            ax_mat.Position = [0.48,0.4,0.52,0.52];
            ax_mat.FontSize=6;
            ax_mat.FontSize=20;%for poster
            % = [0.15,0.67];

            drawnow
            ax3 = axes;        ax3.Position = [-0.37,-0.05,1,0.55];
            go_netviewer_perctl_mars(average_AEC.(fnames{i})(:, :, run), plot_percentile, limit.(fnames{i}), max_sw.(fnames{i}))
            view([-180,0]) %FRONT VIEW
            ax2 = axes;        ax2.Position = [0.05,-0.03,0.75,0.53];
            go_netviewer_perctl_mars(average_AEC.(fnames{i})(:, :, run), plot_percentile, limit.(fnames{i}), max_sw.(fnames{i}))
            view([-90,0]) %side view
            ax1 = axes;ax1.Position = [-0.15,0.45,0.65,0.5];
            go_netviewer_perctl_mars(average_AEC.(fnames{i})(:, :, run), plot_percentile, limit.(fnames{i}), max_sw.(fnames{i}))
            st = sgtitle(['AEC matrix for ' strrep(fnames{i}, '_', '  ') ' band ' trigger_title_string ' ' run_string]);
            st.FontSize = 30;
            drawnow
            set(fh,'Renderer','painters')
        end
    end

    if do_average_AEC_powermap
        fnames = fieldnames(average_AEC);
        for i = 1:length(fnames)
            average_AEC_powermap.(fnames{i}) = squeeze(sum(average_AEC.(fnames{i}), 2, 'omitnan'));
            power_output_max_AEC.(fnames{i}) = [min(abs(average_AEC_powermap.(fnames{i})), [], 'all','omitnan') max(abs(average_AEC_powermap.(fnames{i})),[],'all','omitnan')];

            for region = size(average_AEC_subs.(fnames{i}), 1):-1:1
                AEC_powermap_significance.(fnames{i})(region) = signrank(squeeze(sum(average_AEC_subs.(fnames{i})(region,:,:, 1), 2, 'omitnan')), squeeze(sum(average_AEC_subs.(fnames{i})(region,:,:, 2), 2, 'omitnan')));
            end
        end

        for i = 1:length(fnames)
            AEC_powermap_fig = figure(2*length(fnames)+i);
            AEC_powermap_fig.Name = ['Average AEC powermap ' fnames{i}];
            subplot(1, Nsubplots, run)
            PaintBrodmannAreas_mars_1view(average_AEC_powermap.(fnames{i})(:, run),Nlocs, 256,  ...
                power_output_max_AEC.(fnames{i}), '', [], vars)
            cb = colorbar;
            cb.Label.String = 'Connectivity (A.U.)';
            title(run_string)
            sgtitle(['Average power plot from AEC for band ' strrep(fnames{i}, '_', '  ') ' ' trigger_title_string])
            axis equal
            axis([-68.7943   69.7716 -106.4590   70.1386])
        end
    end
end


% Show residuals
if ~average_all && do_average_powermap
    fnames = fieldnames(average_powermap);
    for i = 1:length(fnames)

        power_residual.(fnames{i}) = average_powermap.(fnames{i})(:, 2) - average_powermap.(fnames{i})(:, 1);
        power_residual_clim = [-1*max(abs(power_residual.(fnames{i})), [], 'all') max(abs(power_residual.(fnames{i})), [], 'all')];
        power_residual.(fnames{i})(powermap_significance.(fnames{i})'>significance_thresh) = 0;

        powermap_fig = figure(i);
        subplot(1, Nsubplots, 3)
        PaintBrodmannAreas_mars_1view(power_residual.(fnames{i}),Nlocs, 256,  ...
            power_residual_clim, '', [], vars)
        cb = colorbar;
        cb.Label.String = 'Relative power (A.U.)';
        title('Residual (H - L)')

        axis equal
        axis([-68.7943   69.7716 -106.4590   70.1386])

        fh = findall(0,'Type','Figure');
        set( findall(fh, '-property', 'fontsize'), 'fontsize', 30)

        % Print out signigicant areas
        disp("Power significant areas: ")
        disp([ location_names(powermap_significance.(fnames{i})'<significance_thresh) num2str(find(powermap_significance.(fnames{i})'<significance_thresh))])

        if do_save
            print(gcf,['D:\OneDrive - The University of Nottingham\Documents\Media\Papers\Violin\powermap_' fnames{i} '_' trigger_title_string '_' average_save_name '.jpg'], "-djpeg")
        end
    end
end

% Show residuals
if ~average_all && do_average_AEC_powermap
    for i = 1:length(fnames)
        AEC_power_residual.(fnames{i}) = average_AEC_powermap.(fnames{i})(:, 2) - average_AEC_powermap.(fnames{i})(:, 1);
        AEC_power_residual_clim = [-max(abs(AEC_power_residual.(fnames{i})), [], 'all') max(abs(AEC_power_residual.(fnames{i})), [], 'all')];
        AEC_power_residual.(fnames{i})(AEC_powermap_significance.(fnames{i})'>significance_thresh) = 0;

        AEC_powermap_fig = figure(2*length(fnames)+i);
        subplot(1, Nsubplots, 3)
        PaintBrodmannAreas_mars_1view(AEC_power_residual.(fnames{i}),Nlocs, 256,  ...
            AEC_power_residual_clim, '', [], vars)
        cb = colorbar;
        cb.Label.String = 'Connectivity (A.U.)';
        title('Residual (H - L)')

        axis equal
        axis([-68.7943   69.7716 -106.4590   70.1386])

        fh = findall(0,'Type','Figure');
        set( findall(fh, '-property', 'fontsize'), 'fontsize', 30)

        % Print out signigicant areas
        % disp("Connectivity significant areas: ")
        % disp([location_names(AEC_powermap_significance.(fnames{i})'<significance_thresh) num2str(find(AEC_powermap_significance.(fnames{i})'<significance_thresh))])

        if do_save
            print(gcf,['D:\OneDrive - The University of Nottingham\Documents\Media\Papers\Violin\AEC_powermap_' fnames{i} '_' trigger_title_string '_' average_save_name '.jpg'], "-djpeg")
        end
    end
end
%%
% Statistical test on motor power
if (average_exp)  || (average_score)
    fnames = fieldnames(average_powermap);
    motor_inds = [15, 16, 20, 21, 23, 24, 56, 57, 61, 62, 64, 65];
    visual_inds = [1 2 3 4 5 17 42 43 44 45 46 58];

    if do_mars_positions_plot
        load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\BrainPlots\aalviewer.mat')
        load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\MARS\sourcepos_mars_mni.mat')
        figure;clf
        % Plot the cortical mesh
        set(gcf,'color',[1 1 1]);
        axis off
        p = patch('faces',aalviewer.faces,'vertices',aalviewer.vertices,'edgecolor','none','facecolor','k','facealpha',0.05);
        hold on

        % Plot Spheres
        for ii = 1:length(motor_inds)
            [x y z] = sphere(10);
            sw = 1;
            surf(sw*x+sourcepos_mars(motor_inds(ii),1)*1000,sw*y+sourcepos_mars(motor_inds(ii),2)*1000,sw*z+sourcepos_mars(motor_inds(ii),3)*1000,'edgecolor','none','facecolor','b','facealpha',1);
        end
        text(sourcepos_mars(motor_inds,1)*1000, sourcepos_mars(motor_inds,2)*1000 + 3, sourcepos_mars(motor_inds,3)*1000, string(motor_inds)')
        axis vis3d
        axis equal
        axis off
        hold off
        rotate3d('on')
        set(gcf,'renderer','opengl')
        colorbar off


        figure;clf
        % Plot the cortical mesh
        set(gcf,'color',[1 1 1]);
        axis off
        p = patch('faces',aalviewer.faces,'vertices',aalviewer.vertices,'edgecolor','none','facecolor','k','facealpha',0.05);
        hold on

        % Plot Spheres
        for ii = 1:length(sourcepos_mars)
            [x y z] = sphere(10);
            sw = 1;
            surf(sw*x+sourcepos_mars(ii,1)*1000,sw*y+sourcepos_mars(ii,2)*1000,sw*z+sourcepos_mars(ii,3)*1000,'edgecolor','none','facecolor','b','facealpha',1);
        end
        text(sourcepos_mars(:,1)*1000, sourcepos_mars(:,2)*1000 + 3, sourcepos_mars(:,3)*1000, string(1:82)')
        axis vis3d
        axis equal
        axis off
        hold off
        rotate3d('on')
        set(gcf,'renderer','opengl')
        colorbar off
    end
end

%% Broadman area for stats
if do_mars_positions_plot
    figure()
    stats_powermap = ones(82, 1).*3;
    stats_powermap(motor_inds) = 1;
    PaintBrodmannAreas_mars_1view(stats_powermap,Nlocs, 256,[0 3], '', [], vars)
    axis equal
    colormap('Hot')
    % colorbar


    figure()
    stats_powermap = ones(82, 1).*3;
    stats_powermap(visual_inds) = 1;
    PaintBrodmannAreas_mars_1view(stats_powermap,Nlocs, 256,[0 3], '', [], vars)
    axis equal
end


%% Bar chart of power values with error bars
% if average_exp
%     fnames = fieldnames(average_powermap_subs);
%     for fname_ind = 1:length(fnames)
%         fname = fnames{fname_ind};
%         labels = ["Before Lesson motor" "After Lesson motor" "Before Lesson visual" "After Lesson visual" "Before Lesson global" "After Lesson global"];
%         bar_chart_func(average_powermap_subs.(fname), average_AEC_subs.(fname), labels, motor_inds, visual_inds, [fname ' motor inds'])
%     end
% elseif average_score
fnames = fieldnames(average_powermap_subs);
for fname_ind = 1:length(fnames)
    fname = fnames{fname_ind};

    % regions = [27 68];
    regions = motor_inds;

    regions = [];
    % regions =[30 34 35 39 40]; regions = [regions regions+41];

    % if fname_ind == 1
    %     % Plot number of trials for each score
    %     violin_number_trials_func(scores_sorted, n_subjects, sub_nums)
    % end

    % Bar chart
    labels = ["Lower Score motor" "Upper Score motor" "Lower Score visual" "Upper Score visual" "Lower Score global" "Upper Score global"];
    bar_chart_func(average_powermap_subs.(fname), average_AEC_subs.(fname), labels, motor_inds, visual_inds, [fname ' ' trigger_title_string ' motor_' average_save_name], do_save)

    if ~isempty(regions)
        for region_no = 1:length(regions)+1
            if region_no == length(regions)+1
                region = regions;
                region_name = ['all regions ' num2str(region)];
            else
                region = regions(region_no);
                region_name = strrep([fname ' ' location_names{region}], '_', ' ');
            end

            fig_val = 100*fname_ind + region_no;
            % if fname_ind == 1
            % Region plot
            fig = figure(fig_val);
            fig.Name = [num2str(region) ' region location'];
            subplot(2, 2, [1 3])
            stats_powermap = ones(82, 1).*3;
            stats_powermap(region) = 1;
            PaintBrodmannAreas_mars_1view(stats_powermap,Nlocs, 256,[0 3], '', [], vars)
            axis equal
            colormap('Hot')
            view([180 0])
            sgtitle(region_name)
            % end

            % power/AEC high/low
            violin_plot_func(average_powermap_subs.(fname), average_AEC_subs.(fname), region, region_name, do_average_powermap, do_average_AEC, fig_val)

            % power/AEC vs score
            % violin_plot_scores_func(all_powermap.(fname), all_AEC.(fname), scores_sorted, region,  region_name, n_subjects, sub_nums)
        end
    end
end
% end

%% scatter plot of score against power and connectivity
% if average_score
%     figure
%     hold on
%     for sub = 1:n_subjects
%         sub_number = sub_nums(sub);
%         scatter(scores_sorted.(['sub_' num2str(sub_number)]), squeeze(mean(all_powermap.beta(motor_inds, sub, 1:ntrials(sub)))))
%     end
%     xlabel('Score')
%     ylabel('Beta power')
%     title('Motor region')
%
%     figure
%     hold on
%     for sub = 1:n_subjects
%         sub_number = sub_nums(sub);
%         scatter(scores_sorted.(['sub_' num2str(sub_number)]), squeeze(mean(all_powermap.beta(visual_inds, sub, 1:ntrials(sub)))))
%     end
%     xlabel('Score')
%     ylabel('Beta power')
%     title('Visual region')
%
%
%     figure
%     hold on
%     for sub = 1:n_subjects
%         sub_number = sub_nums(sub);
%         scatter(scores_sorted.(['sub_' num2str(sub_number)]), squeeze(mean(sum(all_AEC.beta(motor_inds, :, sub, 1:ntrials(sub)), 2, 'omitnan'))))
%     end
%     xlabel('Score')
%     ylabel('Beta connectivity')
%     title('Motor region')
%
%     figure
%     hold on
%     for sub = 1:n_subjects
%         sub_number = sub_nums(sub);
%         scatter(scores_sorted.(['sub_' num2str(sub_number)]), squeeze(mean(sum(all_AEC.beta(visual_inds, :, sub, 1:ntrials(sub)), 2, 'omitnan'))))
%     end
%     xlabel('Score')
%     ylabel('Beta connectivity')
%     title('Visual region')
% end

%%
load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\BrainPlots\aalviewer.mat')
load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\MARS\sourcepos_mars_mni.mat')
figure;clf
% Plot the cortical mesh
set(gcf,'color',[1 1 1]);
axis off
p = patch('faces',aalviewer.faces,'vertices',aalviewer.vertices,'edgecolor','none','facecolor','k','facealpha',0.05);
hold on

% Plot Spheres
for ii = 1:length(sourcepos_mars)
    [x y z] = sphere(10);
    sw = 1;
    surf(sw*x+sourcepos_mars(ii,1)*1000,sw*y+sourcepos_mars(ii,2)*1000,sw*z+sourcepos_mars(ii,3)*1000,'edgecolor','none','facecolor','b','facealpha',1);
    text(sourcepos_mars(ii,1)*1000, sourcepos_mars(ii,2)*1000, sourcepos_mars(ii,3)*1000, num2str(ii));
end
% text(sourcepos_mars(:,1)*1000, sourcepos_mars(:,2)*1000, sourcepos_mars(:,3)*1000, strrep(mars_labels, '_', ' '))
axis vis3d
axis equal
axis off
hold off
rotate3d('on')
set(gcf,'renderer','opengl')
colorbar off


%%
load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\BrainPlots\aalviewer.mat')
load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\MARS\sourcepos_mars_mni.mat')
figure;clf
% Plot the cortical mesh
set(gcf,'color',[1 1 1]);
axis off
stats_powermap = 1:82;
PaintBrodmannAreas_mars_1view(stats_powermap,Nlocs, 256,[0 82], '', [], vars, 'prism')

hold on

% Plot Spheres
for ii = 1:length(sourcepos_mars)
    [x y z] = sphere(10);
    sw = 1;
    surf(sw*x+sourcepos_mars(ii,1)*1000,sw*y+sourcepos_mars(ii,2)*1000,sw*z+sourcepos_mars(ii,3)*1000,'edgecolor','none','facecolor','b','facealpha',1);
    text(sourcepos_mars(ii,1)*1000, sourcepos_mars(ii,2)*1000, sourcepos_mars(ii,3)*1000, num2str(ii));
end
% text(sourcepos_mars(:,1)*1000, sourcepos_mars(:,2)*1000, sourcepos_mars(:,3)*1000, strrep(mars_labels, '_', ' '))
axis vis3d
axis equal
axis off
hold off
rotate3d('on')
set(gcf,'renderer','opengl')
colorbar off



fh = findall(0,'Type','Figure');
set( findall(fh, '-property', 'fontsize'), 'fontsize', 30)

%% Funciton
function [p_val] = stats_func(data, region)
if length(size(data)) == 3
    p_val = signrank(reshape(mean(data(region,:,1), 1), 1, []), reshape(mean(data(region,:,2), 1), 1, []));
elseif length(size(data)) == 4
    p_val = signrank(reshape(mean(sum(data(region,:,:,1), 1), 2, 'omitnan'), 1, []),reshape(mean(sum(data(region,:,:,2), 1), 2, 'omitnan'), 1, []));
end
end

function  bar_chart_func(data_power, data_AEC, labels, motor_inds, visual_inds, title_text, do_save)
fig = figure();
fig.Name = 'power bar';
mean_before_motor = mean(mean(data_power(motor_inds,:,1), 'omitnan'), 'omitnan');
mean_after_motor = mean(mean(data_power(motor_inds,:,2), 'omitnan'), 'omitnan');
err_before_motor = std(mean(data_power(motor_inds,:,1), 'omitnan'), 'omitnan')/sqrt(length(mean(data_power(motor_inds,:,1), 'omitnan')));
err_after_motor = std(mean(data_power(motor_inds,:,2), 'omitnan'), 'omitnan')/sqrt(length(mean(data_power(motor_inds,:,2), 'omitnan')));

mean_before_visual = mean(mean(data_power(visual_inds,:,1), 'omitnan'), 'omitnan');
mean_after_visual = mean(mean(data_power(visual_inds,:,2), 'omitnan'), 'omitnan');
err_before_visual = std(mean(data_power(visual_inds,:,1), 'omitnan'), 'omitnan')/sqrt(length(mean(data_power(visual_inds,:,1), 'omitnan')));
err_after_visual = std(mean(data_power(visual_inds,:,2), 'omitnan'), 'omitnan')/sqrt(length(mean(data_power(visual_inds,:,2), 'omitnan')));

mean_before_global = mean(mean(data_power(:,:,1), 'omitnan'), 'omitnan');
mean_after_global = mean(mean(data_power(:,:,2), 'omitnan'), 'omitnan');
err_before_global = std(mean(data_power(:,:,1), 'omitnan'), 'omitnan')/sqrt(length(mean(data_power(:,:,1), 'omitnan')));
err_after_global = std(mean(data_power(:,:,2), 'omitnan'), 'omitnan')/sqrt(length(mean(data_power(:,:,2), 'omitnan')));

y = [mean_before_motor mean_after_motor mean_before_visual mean_after_visual mean_before_global mean_after_global];
y_err = [err_before_motor err_after_motor err_before_visual err_after_visual err_before_global err_after_global];
x = [1 4.5 9 12.5 17 20.5 ];
bar(x,y)
hold on
xticks(x)
xticklabels(labels)
ylim([0.9*min(y-y_err) 1.1*max(y+y_err)])
ylabel('Relative Power (A.U.)')
errorbar(x, y, y_err, y_err, LineStyle="none", LineWidth=3, CapSize=15)
text(0, 1.05*max(y), sprintf('P = %.4f',stats_func(data_power, motor_inds)))
title([title_text ' power'])

if do_save
    fh = findall(0,'Type','Figure');
    set( findall(fh, '-property', 'fontsize'), 'fontsize', 30)
    print(gcf,['D:\OneDrive - The University of Nottingham\Documents\Media\Papers\Violin\' strrep(title_text, ' ', '_') '_power.jpg'], "-djpeg")
end

fig = figure();
fig.Name = 'AEC bar';
mean_before_motor = mean(mean(sum(data_AEC(motor_inds,:,:,1), 2, 'omitnan')));
mean_after_motor = mean(mean(sum(data_AEC(motor_inds,:,:,2), 2, 'omitnan')));
err_before_motor = std(mean(sum(data_AEC(motor_inds,:,:,1), 2, 'omitnan')))/sqrt(length(mean(sum(data_AEC(motor_inds,:,:,1), 2, 'omitnan'))));
err_after_motor = std(mean(sum(data_AEC(motor_inds,:,:,2), 2, 'omitnan')))/sqrt(length(mean(sum(data_AEC(motor_inds,:,:,2), 2, 'omitnan'))));

mean_before_visual = mean(mean(sum(data_AEC(visual_inds,:,:,1), 2, 'omitnan')));
mean_after_visual = mean(mean(sum(data_AEC(visual_inds,:,:,2), 2, 'omitnan')));
err_before_visual = std(mean(sum(data_AEC(visual_inds,:,:,1), 2, 'omitnan')))/sqrt(length(mean(sum(data_AEC(visual_inds,:,:,1), 2, 'omitnan'))));
err_after_visual = std(mean(sum(data_AEC(visual_inds,:,:,2), 2, 'omitnan')))/sqrt(length(mean(sum(data_AEC(visual_inds,:,:,2), 2, 'omitnan'))));

mean_before_global = mean(mean(sum(data_AEC(:,:,:,1), 2, 'omitnan')));
mean_after_global = mean(mean(sum(data_AEC(:,:,:,2), 2, 'omitnan')));
err_before_global = std(mean(sum(data_AEC(:,:,:,1), 2, 'omitnan')))/sqrt(length(mean(sum(data_AEC(:,:,:,1), 2, 'omitnan'))));
err_after_global = std(mean(sum(data_AEC(:,:,:,2), 2, 'omitnan')))/sqrt(length(mean(sum(data_AEC(:,:,:,2), 2, 'omitnan'))));

y = [mean_before_motor mean_after_motor mean_before_visual mean_after_visual mean_before_global mean_after_global];
y_err = [err_before_motor err_after_motor err_before_visual err_after_visual err_before_global err_after_global];
x = [1 4.5 9 12.5 17 20.5 ];
bar(x,y)
hold on
xticks(x)
xticklabels(labels)
ylim([0.9*min(y-y_err) 1.1*max(y+y_err)])
ylabel('Connectivity')
errorbar(x, y, y_err, y_err, LineStyle="none", LineWidth=3, CapSize=15)
text(0, 1.05*max(y), sprintf('P = %.4f',stats_func(data_AEC, motor_inds)))

title([title_text ' AEC'])

if do_save
    fh = findall(0,'Type','Figure');
    set( findall(fh, '-property', 'fontsize'), 'fontsize', 30)
    print(gcf,['D:\OneDrive - The University of Nottingham\Documents\Media\Papers\Violin\' strrep(title_text, ' ', '_') '_AEC.jpg'], "-djpeg")
end

end

function violin_plot_scores_func(data_power, data_AEC, scores_sorted, inds, title_text, n_subjects, sub_nums)
fig = figure();
fig.Name = [num2str(inds) ' power violin scores'];
hold on
ones=NaN(n_subjects, 1);
twos=NaN(n_subjects, 1);
threes=NaN(n_subjects, 1);
fours=NaN(n_subjects, 1);
fives=NaN(n_subjects, 1);
for sub = 1:n_subjects
    sub_number = sub_nums(sub);
    ones(sub) = mean(squeeze(mean(data_power(inds, sub, scores_sorted.(['sub_' num2str(sub_number)])==1))));
    twos(sub) = mean(squeeze(mean(data_power(inds, sub, scores_sorted.(['sub_' num2str(sub_number)])==2))));
    threes(sub) = mean(squeeze(mean(data_power(inds, sub, scores_sorted.(['sub_' num2str(sub_number)])==3))));
    fours(sub) = mean(squeeze(mean(data_power(inds, sub, scores_sorted.(['sub_' num2str(sub_number)])==4))));
    fives(sub) = mean(squeeze(mean(data_power(inds, sub, scores_sorted.(['sub_' num2str(sub_number)])==5))));
end

violinplot([ones twos threes fours fives], [], 'ShowMean',true, 'ShowMedian', true);
xlabel('Scores')
ylabel('Power')
sgtitle([title_text ' power'])


fig = figure();
fig.Name = [num2str(inds) ' AEC violin scores'];
hold on
ones=NaN(n_subjects, 1);
twos=NaN(n_subjects, 1);
threes=NaN(n_subjects, 1);
fours=NaN(n_subjects, 1);
fives=NaN(n_subjects, 1);
for sub = 1:n_subjects
    sub_number = sub_nums(sub);
    ones(sub) = mean(squeeze(mean(sum(data_AEC(inds, :, sub, scores_sorted.(['sub_' num2str(sub_number)])==1), 'omitnan'))));
    twos(sub) = mean(squeeze(mean(sum(data_AEC(inds, :, sub, scores_sorted.(['sub_' num2str(sub_number)])==2), 'omitnan'))));
    threes(sub) = mean(squeeze(mean(sum(data_AEC(inds, :, sub, scores_sorted.(['sub_' num2str(sub_number)])==3), 'omitnan'))));
    fours(sub) = mean(squeeze(mean(sum(data_AEC(inds, :, sub, scores_sorted.(['sub_' num2str(sub_number)])==4), 'omitnan'))));
    fives(sub) = mean(squeeze(mean(sum(data_AEC(inds, :, sub, scores_sorted.(['sub_' num2str(sub_number)])==5), 'omitnan'))));
end

violinplot([ones twos threes fours fives], [], 'ShowMean',true, 'ShowMedian', true);
xlabel('Scores')
ylabel('Connectivity')
sgtitle([title_text ' AEC'])

end

function violin_number_trials_func(scores_sorted, n_subjects, sub_nums)
fig = figure();
fig.Name = 'N trials';
hold on
ones=NaN(n_subjects, 1);
twos=NaN(n_subjects, 1);
threes=NaN(n_subjects, 1);
fours=NaN(n_subjects, 1);
fives=NaN(n_subjects, 1);
for sub = 1:n_subjects
    sub_number = sub_nums(sub);
    ones(sub) = sum(scores_sorted.(['sub_' num2str(sub_number)])==1);
    twos(sub) = sum(scores_sorted.(['sub_' num2str(sub_number)])==2);
    threes(sub) = sum(scores_sorted.(['sub_' num2str(sub_number)])==3);
    fours(sub) = sum(scores_sorted.(['sub_' num2str(sub_number)])==4);
    fives(sub) = sum(scores_sorted.(['sub_' num2str(sub_number)])==5);
end

violinplot([ones twos threes fours fives], [], 'ShowMean',true, 'ShowMedian', true);
xlabel('Scores')
ylabel('Number trials')
sgtitle('N trials for each score')
end

function violin_plot_func(data_power, data_AEC, inds, title_text, do_average_powermap, do_average_AEC, fig_val)
% Violin plot
if do_average_powermap
    figure(fig_val);
    % fig.Name = [num2str(inds) ' power violin low high'];
    subplot(2, 2, 2)
    violinplot([reshape(data_power(inds,:,1), 1, []); reshape(data_power(inds,:,2), 1, [])]', ["Low"; "High"], 'ShowMean',true, 'ShowMedian', true);
    pval = stats_func(data_power, inds);
    % title([title_text ' power p=' num2str(pval)])
    title(['Power p=' num2str(pval)])
end

if do_average_AEC
    figure(fig_val);
    % fig.Name = [num2str(inds) ' AEC violin low high'];
    subplot(2, 2, 4)
    violinplot([reshape(sum(data_AEC(inds,:,:,1), 2, 'omitnan'), 1, []); reshape(sum(data_AEC(inds,:,:,2), 2, 'omitnan'), 1, [])]', ["Low"; "High"], 'ShowMean',true, 'ShowMedian', true);
    pval = stats_func(data_AEC, inds);
    % title([title_text ' AEC p=' num2str(pval)])
    title(['AEC p=' num2str(pval)])
end
end

function [cmap] = RdBu()


ncols = 256;
pn = 'div';
deep = 0;

% ncols = varargin{1};
% pn = 'div';
% deep = 0;
%
% ncols = varargin{1};
% if sum(strcmp('type',varargin));
%     pn = varargin{find(strcmp('type',varargin))+1};
% else
%     pn = 'div';
% end
% if sum(strcmp('deep',varargin));
%     deep = 1;
% end


if rem(ncols,2)~=0;
    error('Can only accept even numbers');
end

ncols = ncols./2;


% colours
lo        = [5 48 97] / 255;
bottom    = [5 113 176] / 255;
botmiddle = [146 197 222] / 255;
middle    = [247 247 247] / 255;
topmiddle = [244 165 130] / 255;
top       = [202   0  32] / 255;
hi        = [103 0 31] / 255;

% Find ratio of negative to positive
if strncmp(pn,'div',3) || strncmp(pn,'neg',3)


    % Just negative
    if deep
        neg = [lo; bottom; botmiddle; middle];
    else
        neg = [bottom; botmiddle; middle];
    end
    len = length(neg);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, ncols);
    neg128 = zeros(ncols, 3);

    for i=1:3
        % Interpolate over RGB spaces of colormap
        neg128(:,i) = min(max(interp1(oldsteps, neg(:,i), newsteps)', 0), 1);
    end

    cmap = neg128;

end

if strncmp(pn,'div',3) || strncmp(pn,'pos',3)
    % Just positive
    if deep
        pos = [middle; topmiddle; top; hi];
    else
        pos = [middle; topmiddle; top];
    end
    len = length(pos);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, ncols);
    pos128 = zeros(ncols, 3);

    for i=1:3
        % Interpolate over RGB spaces of colormap
        pos128(:,i) = min(max(interp1(oldsteps, pos(:,i), newsteps)', 0), 1);
    end
    cmap = pos128;
end

if strmatch(pn,'div')
    % And put 'em together
    cmap = [neg128; pos128];
end
end