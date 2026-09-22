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
addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files"))
addpath(genpath('scripts'))

if isempty(which('ft_defaults'))
    addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\fieldtrip-20231113"))
    ft_defaults
end

vars.do_save = 1;
vars.plot_psd = 0;
vars.do_notch = 1;
vars.do_HFC = 0;
vars.do_tSSS = 1;
vars.do_AEC_analysis = 1;

vars.MARS = 1;

if vars.MARS
    pos_string = '';
end

bids_path = 'D:\OneDrive - The University of Nottingham\Documents\Data\Violin\V2\Bids\';

run_all_time = tic;

%% Choose files
% for trig_to_use = 1:3
% 
%     % trig_to_use = 3;
% 
%     if trig_to_use == 1
%         vars.use_playing_events = 1;
%         vars.use_rest_events = 0;
%         vars.use_rebound_events = 0;
%     elseif trig_to_use == 2
%         vars.use_playing_events = 0;
%         vars.use_rest_events = 1;
%         vars.use_rebound_events = 0;
%     elseif trig_to_use == 3
%         vars.use_playing_events = 0;
%         vars.use_rest_events = 0;
%         vars.use_rebound_events = 1;
%     end

        vars.use_playing_events = 0;
        vars.use_rest_events = 0;
        vars.use_rebound_events = 0;


    for sub = 2:22
        subject = ['Sub-' num2str(sub, '%03.f')];
        disp(subject)

        for run_no = 1:2
            run = ['run_' num2str(run_no)];
            disp(run)

            clearvars -except All_files n_all_files Onedrive_path vars ...
                run_all_time bids_path F subject sub run run_no trig_to_use pos_string
            close all

            F.OPM_data.path = [bids_path 'MEG\' subject '\'];
            F.OPM_data.name = run;
            F.OPM_data.ext = '.lvm';


            load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\BrainPlots\all_voxels_MARS.mat')
            load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\MARS\MARS_cortical_Key.mat')
            load('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts\BrainPlots\mars_labels_20484.mat')

            brodmann_vars.MARS_atlas_all_vox = MARS_atlas_all_vox;
            brodmann_vars.MARS_cortical_Key = MARS_cortical_Key;
            brodmann_vars.mars_labels_20484 = mars_labels_20484;

            % Broadband and bands of interest
            wins = [1 150; 13 30];
            % wins = [1 150; 55 100];

            if sum([vars.use_playing_events, vars.use_rest_events vars.use_rebound_events]) ~= 1
                error('Only use one trigger condition')
            end


            plot_percentile = 0.90;

            if vars.use_playing_events
                trig_string = 'trial';
            elseif vars.use_rest_events
                trig_string = 'rest';
            elseif vars.use_rebound_events
                trig_string = 'rebound';
            end

            filename = run;

            path_VEs = [bids_path,'Derivatives\VEs\Averaged\' pos_string '\' subject '\'];
            if ~exist(path_VEs,'dir'); mkdir(path_VEs);end
            files_VEs = [filename,'_VE'];

            if vars.do_AEC_analysis
                path_AEC = [bids_path,'Derivatives\AEC\Averaged\' pos_string '\' subject '\'];
                if ~exist(path_AEC,'dir'); mkdir(path_AEC);end
                files_AEC = [filename,'_AEC'];

                path_global_connect_vals = [bids_path,'Derivatives\', ...
                    'global_connect_vals\Averaged\' pos_string '\' subject '\'];
                if ~exist(path_global_connect_vals,'dir'); mkdir(path_global_connect_vals);end
                files_global_connect_vals = [filename,'_global_connect'];

                path_figures_AEC = [bids_path,'Derivatives\figures\Averaged\' ...
                    pos_string '\' subject '\AEC\'];
                if ~exist(path_figures_AEC,'dir'); mkdir(path_figures_AEC);end
            end

            path_power_map = [bids_path,'Derivatives\power_map\Averaged\' ...
                pos_string '\' subject,'\'];
            files_power_map = [filename,'_power_map'];
            if ~exist(path_power_map,'dir'); mkdir(path_power_map);end

            path_figures_power_map = [bids_path,'Derivatives\figures\Averaged\' ...
                pos_string '\' subject '\power_map\'];
            if ~exist(path_figures_power_map,'dir'); mkdir(path_figures_power_map);end
            %% Load MEG data
            if ~exist([F.OPM_data.path F.OPM_data.name '.mat'], 'file')
                disp('Reading Data')
                QZFM_data = read_N1lvm_inscript([F.OPM_data.path F.OPM_data.name F.OPM_data.ext]);
                disp('Saving Data')
                save([F.OPM_data.path F.OPM_data.name '.mat'],"QZFM_data")
                disp('Done')
            else
                disp('Loading Data')
                load([F.OPM_data.path F.OPM_data.name '.mat'])
                disp('Done')
            end

            %% Re-order data
            Data.data = QZFM_data(:,2:end)';
            time = QZFM_data(:,1)';
            Data.samp_frequency = round(1./mean(diff(time)));

            Data.Channel_Info.name = cellstr([[repmat([[repmat('A',8,1);repmat('B',8,1);repmat('C',8,1);...
                repmat('D',8,1);repmat('E',8,1);repmat('F',8,1);...
                repmat('G',8,1);repmat('H',8,1)] repmat(num2str([1:8]'),8,1)],3,1)],...
                [repmat(' [X]',64,1);repmat(' [Y]',64,1);repmat(' [Z]',64,1)]]);


            Layout_Info = tdfread([bids_path 'Helmet_config\' subject '\' run '.tsv']);

            last_line_idx = strmatch('Helmet:',Layout_Info.Name);
            field_names = fieldnames(Layout_Info);
            for n = 1:size(field_names,1)
                Layout_Info.(field_names{n})(last_line_idx,:) = [];
                if max(strcmpi(field_names{n},{'Px';'Py';'Pz';'Ox';'Oy';'Oz';'Layx';'Layz'}))
                    try
                        Layout_Info.(field_names{n}) = str2num(Layout_Info.(field_names{n}));
                    end
                end
            end
            Data.Layout_Info.SlotSensorPairs = [cellstr(Layout_Info.Name) cellstr(Layout_Info.Sensor)];
            CI_name = Data.Channel_Info.name;
            for n = 1:size(CI_name,1)
                CI_name{n} = CI_name{n}(isstrprop(CI_name{n},'alphanum'));
            end
            for n = 1:size(Data.Layout_Info.SlotSensorPairs,1)
                sens_name = Data.Layout_Info.SlotSensorPairs{n,2};
                sens_name = sens_name(isstrprop(sens_name,'alphanum'));

                data_idx = find(strcmpi(sens_name,CI_name));
                if isempty(data_idx)
                    Data.Layout_Info.SlotSensorPairs{n,3} = 0;
                else
                    Data.Layout_Info.SlotSensorPairs{n,3} = data_idx;
                end
            end

            Data.Layout_Info.Position = [Layout_Info.Px,Layout_Info.Py,Layout_Info.Pz];
            Data.Layout_Info.Orientation = [Layout_Info.Ox,Layout_Info.Oy,Layout_Info.Oz];

            % Read sens transform opts
            opts = delimitedTextImportOptions("NumVariables", 4);
            % Specify range and delimiter
            opts.DataLines = [1, Inf];
            opts.Delimiter = "\t";
            % Specify column names and types
            opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4"];
            opts.VariableTypes = ["double", "double", "double", "double"];
            % Specify file level properties
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";

            if isfile([bids_path 'Anat\' subject '\' subject '.nii'])
                try
                    SensorTransform = readtable([bids_path 'Einscan\' subject '\' run '.txt'], opts);
                    SensorTransform = table2array(SensorTransform);
                    disp('Loaded transform')
                catch
                    % call coreg_func with real MRI
                    SensorTransform = coreg_func(bids_path, subject, run, F, ...
                        [bids_path 'Anat\' subject '\' subject '.ply' ], opts, 0);
                end
            else
                try
                    SensorTransform = readtable([bids_path 'Einscan\' subject '\' run '_MNI.txt'], opts);
                    SensorTransform = table2array(SensorTransform);
                    disp('Loaded MNI brain transform')
                catch
                    % call coreg_func with MNI

                    SensorTransform = coreg_func(bids_path, subject, run, F, "D:\\OneDrive - " + ...
                        "The University of Nottingham\\Documents\\Data\\R_backup\\MRI" + ...
                        "\\BabyMRIs\\templates\\Adult\\adultMNI152_T1_1mm_3000.ply", opts, 1);
                end
            end


            Data.Layout_Info.Position_crg = [SensorTransform(1:3,1:3)*Data.Layout_Info.Position' + SensorTransform(1:3,4)/1000]';
            Data.Layout_Info.Orientation_crg = [SensorTransform(1:3,1:3)*Data.Layout_Info.Orientation']';
            Data.Layout_Info.Lay = [Layout_Info.Layx,Layout_Info.Layy];

            data_idx = cell2mat(Data.Layout_Info.SlotSensorPairs(:,3));
            empty_slot_idx = find(data_idx == 0);
            data_idx(empty_slot_idx) = [];
            sensornamesinuse = Data.Channel_Info.name(data_idx);
            f = Data.samp_frequency; % sampling rate

            Sens_pos = Data.Layout_Info.Position_crg;
            Sens_pos(empty_slot_idx,:) = [];
            Sens_ors = Data.Layout_Info.Orientation_crg;
            Sens_ors(empty_slot_idx,:) = [];
            Sens_names = Data.Layout_Info.SlotSensorPairs(:,2);
            Sens_names(empty_slot_idx,:) = [];
            Lay.pos = Data.Layout_Info.Lay;
            Lay.pos(empty_slot_idx,:) = [];

            %% Notch filter
            OPM_tmp = Data.data(data_idx,:).*1e3; % In fT
            OPM_tmp = OPM_tmp - mean(OPM_tmp,2); % Mean Correct

            if vars.do_notch
                disp('Notch Filter')
                if sub == 1
                    % Notch out optitrack peaks at 37 and 83Hz in sub 1
                    notch_freqs = [37 50 83 100];
                    ab = [-2 -22 -2 -6];
                else
                    notch_freqs = [50 100];
                    ab = [-30 -6];
                end

                for nn = 1:length(notch_freqs)
                    Wo =  notch_freqs(nn)/(f/2);  BW = Wo/380;
                    [b,a] = iirnotch(Wo,BW,ab(nn));
                    disp('Applying Notch filter')
                    OPM_tmp = filter(b,a,OPM_tmp,[],2);
                end

            end

            %% Plot PSD
            if vars.plot_psd
                [po1,fs1] = get_PSD_violin(OPM_tmp,f);
                f_psd2 = figure;
                f_psd2.Name = 'PSD';
                p1 = semilogy(fs1,po1);
                hold on
                p2 = plot(0:600,ones(1,601)*15,'--k');
                p2.LineWidth=2;
                xlabel('Frequency (Hz)')
                ylabel('Magnitude (fT/rtHz)')
                % xlim([1,20]);
                ylim([1,1e5])
                xlim([0 100])
                grid on
                legend([cell2mat(Sens_names) repmat(' ',size(Sens_names)) num2str([1:size(Sens_names,1)]')])
                drawnow
            end
            %% Bad Channels
            good_ch = 1:size(OPM_tmp,1);
            % If any of the channels used are undetermined, autofind the bad channels.
            % Otherwise, use the bad channels already defined in the channels tsv
            if ~exist([F.OPM_data.path F.OPM_data.name '_bad_chans.mat'], 'file')
                disp('Finding bad channels')
                ch_mean_psd = mean(po1,2);
                ch_std_psd = std(po1,[],2);

                % Set upper and lower limits for: 1) the range to assess the noise floor
                % (e.g., lo_range = 60 and up_range = 80 will find the mean noise for each
                % sensor between 60 and 80 Hz), and 2) the cut-off floors (e.g., up_lim = 25
                % and low_lim = 7 will set any channels with noise above or below these
                % values as bad)
                lo_range = 60; % Hz
                up_range = 80; % Hz
                range_inds = find(fs1>lo_range & fs1<up_range);

                up_lim = 30; % fT/sqrt(Hz)
                low_lim = 7; % fT/sqrt(Hz)

                bad_ch4 = find(mean(po1(range_inds,:),1)>up_lim | ...
                    mean(po1(range_inds,:),1)<low_lim);
                bad_chans = [Sens_names,cellstr(repmat('Good',size(Sens_names)))];
                bad_chans(bad_ch4,2) = {'Bad'};
                save([F.OPM_data.path F.OPM_data.name '_bad_chans.mat'],'bad_chans')
            else
                load([F.OPM_data.path F.OPM_data.name '_bad_chans.mat'])
            end

            bad_ch = find(contains(bad_chans(:,2),'Bad'));

            if vars.plot_psd
                for n = 1:length(bad_ch)
                    p1(bad_ch(n)).Color = 'r';
                    p1(bad_ch(n)).LineWidth = 5;
                end
                drawnow
            end
            good_ch(bad_ch) = [];

            disp([num2str(length(bad_ch)) ' bad channels found'])

            OPM_data = OPM_tmp(good_ch,:);
            Sens_names = Sens_names(good_ch,:);
            sensornamesinuse = sensornamesinuse(good_ch,:);
            Sens_pos = Sens_pos(good_ch,:);
            Sens_ors = Sens_ors(good_ch,:);
            Lay.pos = Lay.pos(good_ch,:);

            if vars.do_HFC == 1
                if vars.plot_psd
                    disp('Running PSD')
                    [po,fs] = get_PSD_violin(OPM_data,f);
                    po_b4 = po;
                    mean_psd_b4 = mean(po_b4(31:1001,:),2);
                    median_psd_b4 = median(mean_psd_b4);
                    f_psd1 = figure;
                    f_psd1.Name = 'PSD Before HFC';
                    p1 = semilogy(fs,po);
                    hold on
                    p2 = plot(0:600,ones(1,601)*15,'--k');
                    legend(Sens_names)
                    p2.LineWidth=2;
                    title(['Median between 3 and 100Hz = ',num2str(median_psd_b4),'fT/rt(Hz)'])
                    xlabel('Frequency (Hz)')
                    ylabel('Magnitude (fT/rtHz)')
                    % xlim([1,20]);
                    ylim([1,1e5])
                    xlim([0 100])
                end
                disp('HFC')
                N = Sens_ors;
                M = eye(length(N)) - N*pinv(N);
                Sdata = M*OPM_data;

                OPM_data = Sdata;

                if vars.plot_psd
                    disp('Running PSD')
                    [po,fs] = get_PSD_violin(Sdata,f);
                    mean_psd = mean(po(31:1001,:),2);
                    median_psd = median(mean_psd);
                    f_psd2 = figure;
                    f_psd2.Name = 'PSD After MFC';
                    p1 = semilogy(fs,po);
                    hold on
                    p2 = plot(0:600,ones(1,601)*15,'--k');
                    p2.LineWidth=2;
                    title(['Median between 3 and 100Hz = ',num2str(median_psd),'fT/rt(Hz)'])
                    xlabel('Frequency (Hz)')
                    ylabel('Magnitude (fT/rtHz)')
                    % xlim([1,20]);
                    ylim([1,1e5])
                    xlim([0 100])
                    grid on
                    legend([cell2mat(Sens_names) repmat(' ',size(Sens_names)) num2str([1:size(Sens_names,1)]')])
                    disp('Done!')

                    f_psd3 = figure;
                    f_psd3.Name = 'Average PSD Before and After HFC';
                    semilogy(fs,mean(po_b4,2),'r')
                    hold on
                    semilogy(fs,mean(po,2),'b')
                    p2 = plot(0:600,ones(1,601)*15,'--k');
                    p2.LineWidth=2;
                    xlabel('Frequency (Hz)')
                    ylabel('Magnitude (fT/rtHz)')
                    % xlim([1,20]);
                    ylim([1,1e5])
                    xlim([0 100])
                    grid on

                    drawnow
                end

            end

            if vars.do_tSSS
                if vars.plot_psd
                    disp('Running PSD')
                    [po,fs] = get_PSD_violin(OPM_data,f);
                    po_b4 = po;
                    mean_psd_b4 = mean(po_b4(31:1001,:),2);
                    median_psd_b4 = median(mean_psd_b4);
                    f_psd1 = figure;
                    f_psd1.Name = 'PSD Before tSSS';
                    p1 = semilogy(fs,po);
                    hold on
                    p2 = plot(0:600,ones(1,601)*15,'--k');
                    legend(Sens_names)
                    p2.LineWidth=2;
                    title(['Median between 3 and 100Hz = ',num2str(median_psd_b4),'fT/rt(Hz)'])
                    xlabel('Frequency (Hz)')
                    ylabel('Magnitude (fT/rtHz)')
                    % xlim([1,20]);
                    ylim([1,1e5])
                    xlim([0 100])
                    grid on
                end

                disp('tSSS')
                Lin = 10;
                Lout = 3;

                if exist([F.OPM_data.path F.OPM_data.name '_tSSS.mat'], 'file')
                    load([F.OPM_data.path F.OPM_data.name '_tSSS.mat'])
                    disp('Loaded tSSS data')
                else
                    segment_length = 5*f;
                    data_tSSS=[];
                    ft_progress('init', 'text', 'Performing tSSS...')      % ascii progress bar
                    for n = 1:floor(size(OPM_data,2)/segment_length)
                        [sss_outs]=tSSS_opm(OPM_data(:,((n-1)*segment_length) + ...
                            1:n*segment_length),Sens_pos,Sens_ors,Lin,Lout,1,0.95);
                        data_tSSS=[data_tSSS,sss_outs.tSSS];
                        % disp(['Processing tSSS for segment ',num2str(n),'/', ...
                        %     num2str(floor(size(OPM_data,2)/segment_length))])
                        ft_progress(n/floor(size(OPM_data,2)/segment_length), 'Processing segment %d of %d', n, floor(size(OPM_data,2)/segment_length))
                    end
                    ft_progress('close')

                    disp('Saved tSSS data')
                    save([F.OPM_data.path F.OPM_data.name '_tSSS.mat'], 'data_tSSS' ...
                        , 'sss_outs')
                end
                if vars.plot_psd
                    disp('Running PSD')
                    [po,fs] = get_PSD_violin(data_tSSS,f);
                    mean_psd = mean(po(31:1001,:),2);
                    median_psd = median(mean_psd);
                    f_psd2 = figure;
                    f_psd2.Name = 'PSD After tSSS';
                    p1 = semilogy(fs,po);
                    hold on
                    p2 = plot(0:600,ones(1,601)*15,'--k');
                    p2.LineWidth=2;
                    legend(Sens_names)
                    title(['Median between 3 and 100Hz = ',num2str(median_psd),'fT/rt(Hz)'])
                    xlabel('Frequency (Hz)')
                    ylabel('Magnitude (fT/rtHz)')
                    % xlim([1,20]);
                    ylim([1,1e5])
                    xlim([0 100])
                    grid on
                    disp('Done!')

                    f_psd3 = figure;
                    f_psd3.Name = 'Average PSD Before and After tSSS';
                    semilogy(fs,mean(po_b4,2),'r')
                    hold on
                    semilogy(fs,mean(po,2),'b')
                    p2 = plot(0:600,ones(1,601)*15,'--k');
                    p2.LineWidth=2;
                    xlabel('Frequency (Hz)')
                    ylabel('Magnitude (fT/rtHz)')
                    % xlim([1,20]);
                    ylim([1,1e5])
                    xlim([0 100])
                    grid on

                    drawnow
                end

                OPM_data = data_tSSS;
            end

            % bb filter
            [b,a] = butter(4,2*[1 150]/f);
            OPM_dataf_bb = [filtfilt(b,a,OPM_data')]';
            drawnow
            Nchans = size(OPM_data, 1);


            %% Events information
            % Three triggers, one for playing, rebound, rest

            % rest from optitrack - taken at last movement of any rigid body
            trig_chan = 9;
            trig_offset_rest = 10; % Look n seconds after trigger
            duration = 10;
            trial_time_rest = [1:duration*f]./f + trig_offset_rest;

            load([F.OPM_data.path F.OPM_data.name '_opti_events.mat']);

            rest_trigger_matching_inds = find(ismember(events_list_opti.Trig_ON_Chan,trig_chan));

            rest_on_inds = events_list_opti.Trig_ON_Sample(rest_trigger_matching_inds,:)...
                + trig_offset_rest.*f;
            rest_trial_status = cellstr(events_list_opti.Status(rest_trigger_matching_inds,:));
            Ntrials_rest = length(rest_on_inds);


            % rebound from optitrack - taken at last movement of any rigid body
            trig_chan = 9;
            trig_offset_rebound = -5; % Look n seconds after trigger
            duration = 10;
            trial_time_rebound = [1:duration*f]./f + trig_offset_rebound;
            events_list_rebound = load([F.OPM_data.path F.OPM_data.name '_rebound_events.mat']);

            % Rebound events list taken from opti events - to avoid same name
            % as rest events file load in and overwrite variable
            try
                events_list_rebound = events_list_rebound.events_list_rebound;
            catch
                events_list_rebound = events_list_rebound.events_list_opti;
            end

            rebound_trigger_matching_inds = find(ismember(events_list_rebound.Trig_ON_Chan,trig_chan));
            rebound_on_inds = events_list_rebound.Trig_ON_Sample(rebound_trigger_matching_inds,:)...
                + trig_offset_rebound.*f;
            rebound_trial_status = cellstr(events_list_rebound.Status(rebound_trigger_matching_inds,:));
            Ntrials_rebound = length(rebound_on_inds);

            % playing from audio - taken from playing window
            trig_chan = [9];
            trig_offset_playing = 0; % Look n seconds before trigger
            trial_time_playing = [1:duration*f]./f + trig_offset_playing;
            load([F.OPM_data.path F.OPM_data.name '_audio_events.mat']);

            playing_trigger_matching_inds = find(ismember(events_list_audio.Trig_ON_Chan,trig_chan));

            playing_on_inds = events_list_audio.Trig_ON_Sample(playing_trigger_matching_inds,:)...
                + trig_offset_playing.*f;

            % Some events list exported without statuses included
            try
                playing_trial_status = cellstr(events_list_audio.Status(playing_trigger_matching_inds,:));
            catch
                playing_trial_status = repmat("undetermined", length(events_list_audio.Trig_ON_Chan), 1);
            end
            Ntrials_playing = length(playing_on_inds);

            duration_all = 50;
            trial_time_all = [1:duration_all*f]./f + trig_offset_playing;
            %% Do bad trials and combine from both triggers
            if Ntrials_rest ~= Ntrials_playing
                warning('Different number of events in each trigger channel')

                % sort trigger indices and find where triggers aren't
                % sequential
                all_trigs = [playing_on_inds; rebound_on_inds; rest_on_inds];
                trig_chans = [ones(size(playing_on_inds)); ones(size(rebound_on_inds)).*2; ones(size(rest_on_inds)).*3];
                trig_matching_inds = [playing_trigger_matching_inds; rebound_trigger_matching_inds; rest_trigger_matching_inds];
                trial_statuses = [playing_trial_status; rebound_trial_status; rest_trial_status];
                [all_trigs_sorted, order] = sort(all_trigs);
                trig_order = trig_chans(order);
                trig_matching_inds_sorted = trig_matching_inds(order);
                trial_statuses_sorted = trial_statuses(order);

                % First trig should be a 1 (playing)
                if trig_order(1) ~= 1
                    all_trigs_sorted(1) = [];
                    trig_matching_inds_sorted(1) = [];
                    trial_statuses_sorted(1) = [];
                    trig_order(1) = [];
                end

                trig_order_correct = [1; 2; 3];
                start = 0;
                while start<length(all_trigs_sorted)
                    if any(trig_order(1+start :3+start) ~= trig_order_correct)
                        all_trigs_sorted(start + 1) = [];
                        trig_matching_inds_sorted(start + 1) = [];
                        trial_statuses_sorted(start + 1) = [];
                        trig_order(start + 1) = [];
                    else
                        start = start + 3;
                    end
                end

                playing_on_inds = all_trigs_sorted(1:3:end);
                rebound_on_inds= all_trigs_sorted(2:3:end);
                rest_on_inds = all_trigs_sorted(3:3:end);
                playing_trigger_matching_inds = trig_matching_inds_sorted(1:3:end);
                rebound_trigger_matching_inds = trig_matching_inds_sorted(2:3:end);
                rest_trigger_matching_inds = trig_matching_inds_sorted(3:3:end);
                playing_trial_status = trial_statuses_sorted(1:3:end);
                rebound_trial_status = trial_statuses_sorted(2:3:end);
                rest_trial_status = trial_statuses_sorted(3:3:end);
                Ntrials_playing = length(playing_on_inds);
                Ntrials_rebound = length(rebound_on_inds);
                Ntrials_rest = length(rest_on_inds);
            end

            % Remove any triggers too close to the end of the experiment
            too_long = find(rest_on_inds+duration*f - 1 > length(OPM_data));
            playing_on_inds(too_long) = [];
            rebound_on_inds(too_long) = [];
            rest_on_inds(too_long) = [];
            playing_trigger_matching_inds(too_long) = [];
            rebound_trigger_matching_inds(too_long) = [];
            rest_trigger_matching_inds(too_long) = [];
            playing_trial_status(too_long) = [];
            rebound_trial_status(too_long) = [];
            rest_trial_status(too_long) = [];
            Ntrials_playing = length(playing_on_inds);
            Ntrials_rebound = length(rebound_on_inds);
            Ntrials_rest = length(rest_on_inds);

            % Do bad if any undetermined or any don't match between palying and
            % rest
            do_bad_trials = any(ismember(rest_trial_status, ...
                'undetermined')) || any(ismember(playing_trial_status, ...
                'undetermined')) || any(ismember(rebound_trial_status, ...
                'undetermined')) || any(~cellfun(@strcmp, playing_trial_status, ...
                rest_trial_status)) || any(~cellfun(@strcmp, playing_trial_status, ...
                rebound_trial_status)) || any(~cellfun(@strcmp, rebound_trial_status, ...
                rest_trial_status));
            % do_bad_trials = 1;

            if do_bad_trials
                % Manual bad trials on remaining unclassified trials
                disp('Bad Trials')
                gb_trial = playing_trial_status;
                f_bt = figure;
                f_bt.Name = 'Bad Trials';
                res_val = 5e3;
                to_determine = unique([find(~cellfun(@strcmp, playing_trial_status, ...
                    rest_trial_status)); find(~cellfun(@strcmp, playing_trial_status, ...
                    rebound_trial_status)); find(~cellfun(@strcmp, rebound_trial_status, ...
                    rest_trial_status))]);

                % to_determine = 1:Ntrials_playing;
                trial = 1;
                while trial < length(to_determine) + 1
                    n = to_determine(trial);
                    ax1 = subplot(131);
                    playing = OPM_dataf_bb(:,playing_on_inds(n):playing_on_inds(n)+duration*f - 1);
                    plot(trial_time_playing,[linspace(1,size(playing,1).*res_val, ...
                        size(playing,1))]'+playing)
                    f_bt.CurrentAxes.XLim = [trig_offset_playing duration+trig_offset_playing];
                    f_bt.CurrentAxes.YLim = [-res_val (size(playing,1)+1).*res_val];
                    xlabel('Trial Time (s)')
                    title('Playing')

                    ax2 = subplot(132);
                    rebound = OPM_dataf_bb(:,rebound_on_inds(n):rebound_on_inds(n)+duration*f - 1);
                    plot(trial_time_rebound,[linspace(1,size(rebound,1).*res_val, ...
                        size(rebound,1))]'+rebound)
                    f_bt.CurrentAxes.XLim = [trig_offset_rebound duration+trig_offset_rebound];
                    f_bt.CurrentAxes.YLim = [-res_val (size(rebound,1)+1).*res_val];
                    xlabel('Trial Time (s)')
                    title('Rebound')

                    ax3 = subplot(133);
                    rest = OPM_dataf_bb(:,rest_on_inds(n):rest_on_inds(n)+duration*f - 1);
                    plot(trial_time_rest,[linspace(1,size(rest,1).*res_val, ...
                        size(rest,1))]'+rest)
                    f_bt.CurrentAxes.XLim = [trig_offset_rest duration+trig_offset_rest];
                    f_bt.CurrentAxes.YLim = [-res_val (size(rest,1)+1).*res_val];
                    xlabel('Trial Time (s)')
                    title('Rest')

                    sgtitle(['Trial ' num2str(n) '/' num2str(Ntrials_rest)])
                    gb_trial{n} = questdlg('Is this a good or bad trial?', ...
                        'Good or bad trials','Good','Bad','Go back','Good');
                    if strcmpi(gb_trial{n},'Go back')
                        trial = trial-1;
                    elseif isempty(gb_trial{n})
                        error('Cancelled')
                    else
                        trial = trial + 1;
                    end
                    clear playing rest
                end

                % Update TSVs
                trial_status = gb_trial';
                events_list_opti.Status(rest_trigger_matching_inds) = gb_trial';
                save([F.OPM_data.path F.OPM_data.name '_opti_events.mat'],...
                    'events_list_opti');
                events_list_rebound.Status(rebound_trigger_matching_inds) = gb_trial';
                save([F.OPM_data.path F.OPM_data.name '_rebound_events.mat'],...
                    'events_list_rebound');
                events_list_audio.Status(playing_trigger_matching_inds) = gb_trial';
                save([F.OPM_data.path F.OPM_data.name '_audio_events.mat'],...
                    'events_list_audio');
                close()
            else
                trial_status = playing_trial_status;
            end

            if vars.use_playing_events
                on_inds = playing_on_inds;
                trial_status = playing_trial_status;
                trial_time = trial_time_playing;
                Ntrials = Ntrials_playing;
            elseif vars.use_rest_events
                on_inds = rest_on_inds;
                trial_status = rest_trial_status;
                trial_time = trial_time_rest;
                Ntrials = Ntrials_rest;
            elseif vars.use_rebound_events
                on_inds = rebound_on_inds;
                trial_status = rebound_trial_status;
                trial_time = trial_time_rebound;
                Ntrials = Ntrials_rebound;
            end

            clear   too_long playing_on_inds rest_on_inds playing_trigger_matching_inds ...
                rest_trigger_matching_inds playing_trial_status rest_trial_status ...
                Ntrials_playing Ntrials_rest


            %% Remove bad trials from triggers
            bad_trials = find(strcmpi(trial_status,'bad'));
            on_inds(bad_trials) = [];
            Ntrials = length(on_inds);
            clear OPM_trial_f

            %% MRI
            % Load MRI
            if isfile([bids_path 'Anat\' subject '\' subject '.nii'])
                F.mri.path = [bids_path 'Anat\' subject '\'];
                F.mri.name = subject;
                F.mri.ext = '.nii';
            else
                F.mri.path = [bids_path 'Anat\MNI\'];
                F.mri.name = 'adultMNI152_T1_1mm';
                F.mri.ext = '.nii';
            end

            mri = ft_read_mri([F.mri.path F.mri.name F.mri.ext]);
            if ~exist([F.mri.path 'meshes.mat'], 'file')
                if ~exist([F.mri.path 'segmentedmri.mat'], 'file')
                    disp('Segmenting MRI...')
                    cfg = [];
                    cfg.output    = {'brain' 'scalp' 'skull'};
                    cfg.scalpthreshold = 0.02; % 0.1 by default
                    segmentedmri  = ft_volumesegment(cfg, mri);
                    save([F.mri.path 'segmentedmri.mat'],'segmentedmri')
                    disp('Done!')
                else
                    disp('Loading Segmented MRI...')
                    load([F.mri.path 'segmentedmri.mat'])
                    disp('Done!')
                end

                cfg = [];
                cfg.tissue = {'brain' 'scalp' 'skull'};
                cfg.numvertices = [5000 5000 5000];

                mesh2 = ft_prepare_mesh(cfg,segmentedmri);
                mesh1 = ft_convert_units(mesh2,'m');
                for n = size(mesh1,2):-1:1
                    meshes(n).pnt = mesh1(n).pos;
                    meshes(n).tri = mesh1(n).tri;
                    meshes(n).unit = mesh1(n).unit;
                    meshes(n).name = cfg.tissue{n};
                end
                save([F.mri.path 'meshes.mat'],'meshes','segmentedmri')
                disp('Mesh saved')
            else
                load([F.mri.path 'meshes.mat'])
            end


            if vars.MARS
                if ~exist([F.mri.path 'centroids_MARS.mat'], 'file')
                    disp('Getting MARS regions')
                    if ~exist([F.mri.path 'MARS_regions_' F.mri.name '-FLIRT_fnirt' F.mri.ext])
                        disp('Unzipping MARS regions MRI')
                        gunzip([F.mri.path 'MARS_regions_' F.mri.name '-FLIRT_fnirt' F.mri.ext '.gz'])
                    end
                    MARS_locations = ft_read_mri([F.mri.path 'MARS_regions_' F.mri.name '-FLIRT_fnirt' F.mri.ext]);
                    % change units to metres
                    MARS_locations = ft_convert_units(MARS_locations,'m');

                    %%
                    S.mri_file = [F.mri.path F.mri.name F.mri.ext];
                    S.meshes_file = [F.mri.path 'meshes.mat'];
                    S.sensor_info.pos = Sens_pos;
                    S.sensor_info.ors = Sens_ors;
                    S.meshes_file = [F.mri.path 'meshes.mat'];
                    % find the location indices
                    load('MARS_cortical_Key.mat')
                    [sourcepos_vox] = get_MARS_coords(MARS_locations,S,MARS_cortical_Key);

                    % apply transform to find locations in metres
                    sourcepos = ft_warp_apply(MARS_locations.transform,sourcepos_vox);

                    save([F.mri.path,'centroids_MARS.mat'],'sourcepos')
                    disp('MARS regions extracted')
                    clear S
                else
                    load([F.mri.path,'centroids_MARS.mat'])
                end
            end


            %Loop for each frequency band
            for n_band = 1:length(wins)
                hp = wins(n_band,1);
                lp = wins(n_band,2);
                % Filter data to band of interest
                [b,a] = butter(4,2*[hp lp]/f);
                OPM_data_f = [filtfilt(b,a,OPM_data')]';
                for ii = length(on_inds):-1:1
                    OPM_data_f_T_mat(:, : ,ii) = OPM_data_f(:, on_inds(ii):on_inds(ii) + duration*f -1);
                end
                OPM_data_f_T = reshape(OPM_data_f_T_mat, size(OPM_data_f_T_mat, 1), size(OPM_data_f_T_mat, 2)*Ntrials);

                if n_band == 1
                    % Covariance matrix
                    mu = 0.001;
                    % calculate cov for each trial then average otherwise get steps
                    for ii = length(on_inds):-1:1
                        C_1(:, :, ii) = cov(OPM_data_f_T_mat(:, :, ii)');
                    end
                    C = mean(C_1,3);

                    maxEV = max(svd(C));
                    Noise_Cr = min(svd(C)).*eye(size(C));
                    C = C + mu.*maxEV.*eye(size(C));
                    condnr = cond(C);
                    Cinv = inv(C);
                    % f_cov = figure;
                    % f_cov.Name = 'Covariance';
                    % imagesc(C)
                    % axis square
                    % title(sprintf('\n mu = %.3f, cond(C) = %.2f',...
                    %     mu,condnr))

                    % Beamforming
                    S.mri_file = [F.mri.path F.mri.name F.mri.ext];
                    S.meshes_file = [F.mri.path 'meshes.mat'];
                    S.sensor_info.pos = Sens_pos;
                    S.sensor_info.ors = Sens_ors;
                    S.C = C;
                    S.Cinv = Cinv;
                    if vars.do_HFC
                        S.M = M;
                    end

                    if vars.do_tSSS
                        S.Lin = Lin;
                        S.Lout = Lout;

                        S.sss_ins.S=sss_outs.S;
                        S.sss_ins.Sin=sss_outs.Sin;
                        S.sss_ins.Sout=sss_outs.Sout;
                        S.sss_ins.Sn=sss_outs.Sn;
                    end

                    [bf_outs] = run_beamformer_v3_SSS_quick('shell',sourcepos,S,1,'N/A',vars.do_HFC, vars.do_tSSS);

                    close
                    % close
                end
                %% Powermap
                Nlocs = 82;
                % - bb filtered VE
                VE = zeros(length(OPM_data_f_T),size(sourcepos, 1));
                VE_mat = zeros(size(OPM_data_f_T_mat, 2), size(OPM_data_f_T_mat, 3), size(sourcepos, 1));
                for n_source = 1:size(sourcepos, 1)
                    w_aal_region = bf_outs.Weights(:,n_source);
                    VE(:,n_source) = transpose((w_aal_region'*OPM_data_f_T)./ ...
                        sqrt(w_aal_region'*w_aal_region));
                    % Unnormalised
                    % VE(:,n_source) = (w_aal_region'*OPM_data_f_T)';
                    for n_trial = 1:Ntrials
                        VE_mat(:,n_trial,n_source) = (w_aal_region'*squeeze(OPM_data_f_T_mat(:, :, n_trial)))';
                    end
                end

                if vars.do_save
                    save([path_VEs files_VEs '_' num2str(hp) '_to_' ...
                        num2str(lp) '_' trig_string '.mat'],'VE')
                end
                duration_in_s = size(OPM_data_f_T,2)./f;

                % power averaged over trials
                for n_trial = Ntrials:-1:1
                    [power_values,fs_band] = pwelch(VE_mat(:,n_trial, :),f*5,[],[],f);
                    fs_band_interval = fs_band(2);
                    hp_index = round(hp./fs_band_interval + 1);
                    lp_index = round(lp./fs_band_interval + 1);
                    integral_power_mat(:, n_trial) = trapz(fs_band(hp_index:lp_index), ...
                        power_values(hp_index:lp_index,:));
                end
                integral_power = mean(integral_power_mat, 2);

                if n_band == 1
                    power_output_bb = integral_power;
                    power_output_bb_mat = integral_power_mat;
                else
                    power_output_relative =  integral_power./power_output_bb;
                    power_output_relative_mat =  integral_power_mat./power_output_bb_mat;

                    if vars.do_save
                        save([path_power_map files_power_map '_' num2str(hp) '_to_'  ...
                            num2str(lp) '_' trig_string '.mat'],'power_output_relative')
                        save([path_power_map files_power_map '_' num2str(hp) '_to_'  ...
                            num2str(lp) '_' trig_string '_all_trials.mat'],'power_output_relative_mat')
                    end
                    power_output_max = max(power_output_relative);

                    if vars.MARS
                        figure()
                        PaintBrodmannAreas_mars_1view(power_output_relative,Nlocs, 256,  ...
                            [0 power_output_max], '', [], brodmann_vars)
                        colorbar
                        title(['Power plot for band ' num2str(hp) '-' num2str(lp)  ...
                            'Hz sub ' subject ' ' run ' ' trig_string ' PSD'])
                        axis square
                        if vars.do_save
                            savefig([path_figures_power_map files_power_map '_' num2str(hp)  ...
                                '_to_' num2str(lp) '_' trig_string '_fig'])
                        end
                    end
                end

                %%
                if vars.do_AEC_analysis
                    if n_band>= 2 % Anything not broad band
                        down_f = 1;
                        %% amplitude envelope correlation - Connectivity
                        AEC = zeros(Nlocs,Nlocs);
                        tic
                        ft_progress('init', 'text', 'Calculating AEC...')      % ascii progress bar
                        for seed = 1:Nlocs
                            for n_trial = Ntrials:-1:1
                                % caclculate connectivity over each trial
                                % and then average together
                                Xsig = squeeze(VE_mat(:,n_trial, seed));
                                for test = 1:Nlocs
                                    if seed == test
                                        AEC_mat(seed,test, n_trial) = NaN;
                                    else
                                        Ysig = squeeze(VE_mat(:, n_trial, test));
                                        X_win = (Xsig - mean(Xsig));
                                        Y_win = (Ysig - mean(Ysig));
                                        %%regress leakage
                                        beta_leak = (pinv(X_win)*Y_win);
                                        Y_win_cor = Y_win - X_win*beta_leak;
                                        %%calculate envelopes
                                        H_X = abs(hilbert(X_win));
                                        H_X_d = H_X;
                                        % H_X_d = mean(reshape(H_X,f/down_f,
                                        % round(duration_in_s*down_f),1));

                                        %%calculate envelopes
                                        H_Y = abs(hilbert(Y_win_cor));
                                        H_Y_d = H_Y;
                                        % H_Y_d = mean(reshape(H_Y,f/down_f,
                                        % round(duration_in_s*down_f),1));
                                        AEC_mat(seed,test, n_trial) = corr(H_X_d,H_Y_d);
                                    end
                                    AEC(seed, test) = mean(AEC_mat(seed, test, :), 3, 'omitnan');
                                end
                            end
                            ft_progress(seed/Nlocs, 'Processing region %d of %d', seed, Nlocs)

                            %if seed >1;fprintf(repmat('\b',1,fl));end
                            %fl=fprintf('Getting AEC Freq %2d/%2d (%3d-%3d Hz) | AEC ...
                            % conn. region %2d\n',f_i,length(hpfs),hp,lp,seed);
                        end
                        ft_progress('close')
                        toc
                        AEC = 0.5*(AEC + AEC');
                        if vars.do_save
                            save([path_AEC files_AEC '_' num2str(hp) '_to_' ...
                                num2str(lp) '_' trig_string ...
                                '.mat'],'AEC')
                            save([path_AEC files_AEC '_' num2str(hp) '_to_' ...
                                num2str(lp) '_' trig_string ...
                                '_all_trials.mat'],'AEC_mat')
                        end
                        %%
                        AEC_max = max(abs(AEC),[],'all');
                        figure()
                        subplot(121)
                        imagesc(AEC);colorbar; clim([-AEC_max AEC_max])
                        title(['AEC matrix for band ' num2str(hp) '-' num2str(lp)  ...
                            'Hz sub ' subject ' ' run ' ' trig_string])
                        subplot(122)
                        go_netviewer_perctl_mars(AEC,plot_percentile)
                        drawnow
                        if vars.do_save
                            savefig([path_figures_AEC files_AEC '_' num2str(hp) '_to_'  ...
                                num2str(lp) '_' trig_string '_fig'])
                        end
                        %% - takes mean of only above diagonal elements
                        x_mat_row = 1:Nlocs;
                        y_mat_col = [1:Nlocs]';
                        x_mat = zeros(Nlocs,Nlocs);
                        y_mat = zeros(Nlocs,Nlocs);
                        for n = 1:Nlocs
                            x_mat(n,:) = x_mat_row;
                            y_mat(:,n) = y_mat_col;
                        end
                        x_y_mat = x_mat./y_mat;
                        x_y_mat(x_y_mat<=1) = 0;
                        x_y_mat(x_y_mat>1) = 1;

                        global_connect_mat = AEC.*x_y_mat;
                        global_connect = mean(global_connect_mat,'all','omitnan');
                        if vars.do_save
                            save([path_global_connect_vals files_global_connect_vals '_'  ...
                                num2str(hp) '_to_' num2str(lp) '_' trig_string ...
                                '.mat'],'global_connect')
                        end
                    end
                end
            end
        end
    end
% end

disp(['Run whole script in ' num2str(toc(run_all_time)) 's']);

%% Functions
function [po, fs] = get_PSD_violin(data,f)
S = [];
eD = [];
S.trialength = 10*f;
for n = 1:floor(length(data)/S.trialength)
    eD(:,:,n) = data(:,((n-1)*S.trialength) + 1:n*S.trialength);
end
% set flattop window
nepochs = size(eD,3);
F = [0:.1:600]';
pow = zeros(size(F,1),size(eD,1),nepochs);
wind = window(@flattopwin,size(eD,2));
% loop over windows dc correcting along the way
for j = 1:nepochs
    Btemp = eD(:,:,j)';
    mu = mean(Btemp);
    zf = bsxfun(@minus,Btemp,mu);
    nchan = size(zf,2);
    fzf = zf;
    [pxx,fs] = pwelch(fzf,wind, 0, F,f);
    pow(:,:,j) = sqrt(pxx);
end
Nchans = size(data,1);
chans = [1:Nchans];
po = median(pow(:,chans,:),3);
end
