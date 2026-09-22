%% Violin Beamformer
clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

% If on new desktop, change onedrive location from C: to D:
if strcmp(getenv('COMPUTERNAME'), 'DUIP78552')
    Onedrive_path = 'D:\OneDrive - The University of Nottingham\';
else
    Onedrive_path = 'C:\Users\ppyjg12\OneDrive - The University of Nottingham \';
end

% addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files"))
% addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\tools"))
% addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\Beamformer"))
% addpath(strcat(Onedrive_path, "Documents\MATLAB_C"))
% addpath(strcat(Onedrive_path, "Documents\MATLAB_C\Optitrack"))
% addpath(genpath('C:\Users\ppyjg12\Documents\Projects\Violin\V2\scripts'))

addpath('C:\Users\ppyjg12\Documents\Projects\Violin\publication\dependencies')
if isempty(which('ft_defaults'))
    addpath(strcat(Onedrive_path, "Documents\MATLAB\Matlab_files\fieldtrip-20231113"))
    ft_defaults
end
hp = 13;
lp = 30;

run_all_timing = tic;

% only do one of HFC or tSSS
do_HFC = 0;
do_tSSS = 1;

do_notch = 1;
do_regression = 0;

% Only do one of MARs / peak location / whole head at once
% leave both to 0 to look at whole head
MARS = 0;
peak_location = 1;  % Peak location taken from averaged tstat maps

save_vars = 0;
save_unaveraged_vars = 0;
save_overlay = 0;
plot_psd = 0;


bids_path = 'D:\OneDrive - The University of Nottingham\Documents\Data\Violin\V2\Bids\';
%% Choose files
% subject 1 had issues with triggers so does not have the correct events
% files
for sub = 2:22

    subject = ['Sub-' num2str(sub, '%03.f')];
    disp(subject)

    for run_no = 1:2
        run = ['run_' num2str(run_no)];
        disp(run)

        clearvars -except hp lp run_all_timing do_HFC do_tSSS bids_path sub ...
            subject run_no Onedrive_path plot_psd MARS save_overlay save_vars ...
            do_notch save_unaveraged_vars do_regression run ...
            peak_location
        close all

        F.OPM_dat.path = [bids_path 'MEG\' subject '\'];
        F.OPM_dat.name = run;
        F.OPM_dat.ext = '.lvm';
        if ~exist([F.OPM_dat.path F.OPM_dat.name '.mat'], 'file')
            disp('Reading Data')
            QZFM_data = read_N1lvm_inscript([F.OPM_dat.path F.OPM_dat.name F.OPM_dat.ext]);
            disp('Saving Data')
            save([F.OPM_dat.path F.OPM_dat.name '.mat'],"QZFM_data")
            disp('Done')
        else
            disp('Loading Data')
            load([F.OPM_dat.path F.OPM_dat.name '.mat'])
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

        % Data.Layout_Info.Position_crg = [Layout_Info.Px,Layout_Info.Py,Layout_Info.Pz];
        % Data.Layout_Info.Orientation_crg = [Layout_Info.Ox,Layout_Info.Oy,Layout_Info.Oz];
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

        if do_notch
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
                OPM_tmp = filter(b,a,OPM_tmp,[],2);
            end
        end

        %% Plot PSD
        if plot_psd
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
        if ~exist([F.OPM_dat.path F.OPM_dat.name '_bad_chans.mat'], 'file')
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
            save([F.OPM_dat.path F.OPM_dat.name '_bad_chans.mat'],'bad_chans')
        else
            load([F.OPM_dat.path F.OPM_dat.name '_bad_chans.mat'])
        end

        bad_ch = find(contains(bad_chans(:,2),'Bad'));

        if plot_psd
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

        if do_HFC == 1
            if plot_psd
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

            if plot_psd
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

        if do_tSSS
            if plot_psd
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

            if exist([F.OPM_dat.path F.OPM_dat.name '_tSSS.mat'], 'file')
                load([F.OPM_dat.path F.OPM_dat.name '_tSSS.mat'])
                disp('Loaded tSSS data')
            else
                segment_length = 5*f;
                data_tSSS=[];
                ft_progress('init', 'text', 'Performing tSSS...')      % ascii progress bar
                for n = 1:floor(size(OPM_data,2)/segment_length)
                    [sss_outs]=tSSS_opm(OPM_data(:,((n-1)*segment_length) + ...
                        1:n*segment_length),Sens_pos,Sens_ors,Lin,Lout,1,0.95);
                    data_tSSS=[data_tSSS,sss_outs.tSSS];
                    ft_progress(n/floor(size(OPM_data,2)/segment_length), 'Processing segment %d of %d', n, floor(size(OPM_data,2)/segment_length))
                end
                ft_progress('close')

                disp('Saved tSSS data')
                save([F.OPM_dat.path F.OPM_dat.name '_tSSS.mat'], 'data_tSSS' ...
                    , 'sss_outs')
            end
            if plot_psd
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

        %% Clear OPM_tmp and any other large vars
        clear OPM_tmp OPM_tmp_f Sdata Sdata_f

        %%  Find trigger data
        find_trigs = 0;
        if find_trigs

            trig_std.num = 194:226;
            trig_std.name = string(trig_std.num);
            trig_nums = [15 6 13 4 18 11 2 9]; %Add one because time is first column
            trig_names = string({'Trig_1', 'Trig_2', 'Trig_3', 'Trig_4', ...
                'Trig_5', 'Trig_6', 'Trig_7', 'Trig_8'});
            for i =1:8
                trig_std.name(trig_nums(i)) = trig_names(i);
            end
            trig_std.std = std(QZFM_data(:, 194:226), 0, 1);
            clear trig_nums trig_names i

            med = median(trig_std.std, 'omitnan');
            trig_std.name(trig_std.std<med*2) = [];
            trig_std.num(trig_std.std<med*2) = [];
            trig_std.std(trig_std.std<med*2) = [];
            trig_std.name(trig_std.std>med*100) = [];
            trig_std.num(trig_std.std>med*100) = [];
            trig_std.std(trig_std.std>med*100) = [];
            trig_std.name(isnan(trig_std.std)) = [];
            trig_std.num(isnan(trig_std.std)) = [];
            trig_std.std(isnan(trig_std.std)) = [];

            figure('Name','All trigger channels')
            for i = 1:length(trig_std.name)
                plot(time, QZFM_data(:, trig_std.num(i)), 'DisplayName', trig_std.name(i));
                hold on
            end

            legend()
        end

        % clear trig_std
        %% Events information
        % Three triggers, one for playing, rebound, rest

        % rest from optitrack - taken at last movement of any rigid body
        trig_chan = 9;
        trig_offset_rest = 10; % Look n seconds after trigger
        duration = 10;
        trial_time_rest = [1:duration*f]./f + trig_offset_rest;

        load([F.OPM_dat.path F.OPM_dat.name '_opti_events.mat']);

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
        events_list_rebound = load([F.OPM_dat.path F.OPM_dat.name '_rebound_events.mat']);

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
        load([F.OPM_dat.path F.OPM_dat.name '_audio_events.mat']);

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
            save([F.OPM_dat.path F.OPM_dat.name '_opti_events.mat'],...
                'events_list_opti');
            events_list_rebound.Status(rebound_trigger_matching_inds) = gb_trial';
            save([F.OPM_dat.path F.OPM_dat.name '_rebound_events.mat'],...
                'events_list_rebound');
            events_list_audio.Status(playing_trigger_matching_inds) = gb_trial';
            save([F.OPM_dat.path F.OPM_dat.name '_audio_events.mat'],...
                'events_list_audio');
            close()
        else
            trial_status = playing_trial_status;
        end

        %% Remove bad trials from triggers
        bad_trials = find(strcmpi(trial_status,'bad'));
        playing_on_inds(bad_trials) = [];
        rebound_on_inds(bad_trials) = [];
        rest_on_inds(bad_trials) = [];
        Ntrials = length(playing_on_inds);

        disp([num2str(length(bad_trials)) ' bad trials removed'])

        %% Filter data to band of interest
        [b,a] = butter(4,2*[hp lp]/f);
        OPM_data_f = [filtfilt(b,a,OPM_data')]';

        %% Create matrix of trials
        for n = Ntrials:-1:1
            data_window_playing = playing_on_inds(n):playing_on_inds(n) + duration*f - 1;
            OPM_trialsf_playing(:,:,n) = OPM_data_f(:,data_window_playing).*1e-15;
            OPM_trialsbb_playing(:,:,n) = OPM_dataf_bb(:,data_window_playing).*1e-15;

            data_window_rebound = rebound_on_inds(n):rebound_on_inds(n) + duration*f - 1;
            OPM_trialsf_rebound(:,:,n) = OPM_data_f(:,data_window_rebound).*1e-15;
            OPM_trialsbb_rebound(:,:,n) = OPM_dataf_bb(:,data_window_rebound).*1e-15;

            data_window_rest = rest_on_inds(n):rest_on_inds(n) + duration*f - 1;
            OPM_trialsf_rest(:,:,n) = OPM_data_f(:,data_window_rest).*1e-15;
            OPM_trialsbb_rest(:,:,n) = OPM_dataf_bb(:,data_window_rest).*1e-15;
        end
        OPM_data_f_T_playing = reshape(OPM_trialsf_playing,Nchans,round(duration*f)*Ntrials);
        OPM_data_f_T_rebound = reshape(OPM_trialsf_rebound,Nchans,round(duration*f)*Ntrials);
        OPM_data_f_T_rest = reshape(OPM_trialsf_rest,Nchans,round(duration*f)*Ntrials);

        OPM_data_bb_T_playing = reshape(OPM_trialsbb_playing,Nchans,round(duration*f)*Ntrials);
        OPM_data_bb_T_rebound = reshape(OPM_trialsbb_rebound,Nchans,round(duration*f)*Ntrials);
        OPM_data_bb_T_rest = reshape(OPM_trialsbb_rest,Nchans,round(duration*f)*Ntrials);

        %% TFS
        cbar_lim = 0.5; % colour bar limit (relative change)
        highpass = 5:110;
        lowpass = 6.5:111.5;
        fre = highpass + ((lowpass - highpass)./2);

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
                cfg.scalpthreshold = 0.05; % 0.1 by default - change per particiapnt
                cfg.write = 'yes';
                cfg.name = 'Segmentedmri_output.nii';
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


        if MARS
            % if ~exist([F.mri.path 'centroids_MARS.mat'], 'file')
                disp('Getting MARS regions')
                if ~exist([F.mri.path 'MARS_regions_' F.mri.name F.mri.ext])
                    disp('Unzipping MARS regions MRI')
                    gunzip([F.mri.path 'MARS_regions_' F.mri.name F.mri.ext '.gz'])
                end
                MARS_locations = ft_read_mri([F.mri.path 'MARS_regions_' F.mri.name F.mri.ext]);
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

                % save([F.mri.path,'centroids_MARS.mat'],'sourcepos')
                disp('MARS regions extracted')
                clear S
            % else
            %     load([F.mri.path,'centroids_MARS.mat'])
            % end
        
        elseif peak_location
            if ~exist([F.mri.path 'peak_location_avg_tstat_cov.mat'], 'file')
                disp('Getting peak region')
                if ~exist([strrep(F.mri.path, 'Anat', 'Derivatives\Tstat') 'mean_t_stat_13-30_2trigs_cov_2MRI' F.mri.ext], 'file')
                    disp('Unzipping average tstat')
                    gunzip([strrep(F.mri.path, 'Anat', 'Derivatives\Tstat') 'mean_t_stat_13-30_2trigs_cov_2MRI' F.mri.ext '.gz'])
                end
                avg_tstat = ft_read_mri([strrep(F.mri.path, 'Anat', 'Derivatives\Tstat') 'mean_t_stat_13-30_2trigs_cov_2MRI' F.mri.ext]);

                %% do without separate script
                voxID = find(avg_tstat.anatomy);
                [x, y, z] = ind2sub(avg_tstat.dim(1:3),voxID);

                sourcepos_vox = [x y z];
                sourcepos1 = ft_warp_apply(avg_tstat.transform,sourcepos_vox);
                sourcepos = sourcepos1/1000; % convert to metres

                tstat = avg_tstat.anatomy(voxID);

                tstat_LH = tstat(sourcepos(:, 1)<0, :);
                sourcepos_LH = sourcepos(sourcepos(:, 1)<0, :);
                [val, LH_peak_ind] = min(tstat_LH);
                sourcepos_LH = sourcepos_LH(LH_peak_ind, :);

                tstat_RH = tstat(sourcepos(:, 1)>0, :);
                sourcepos_RH = sourcepos(sourcepos(:, 1)>0, :);
                [val, RH_peak_ind] = min(tstat_RH);
                sourcepos_RH = sourcepos_RH(RH_peak_ind, :);
              
                save([F.mri.path,'peak_location_avg_tstat_cov.mat'],'sourcepos_LH', 'sourcepos_RH')
                disp('Peak region extracted')
            else
                load([F.mri.path,'peak_location_avg_tstat_cov.mat'])
            end

            %Plot to check peak LH centroid
            f_mri_LH_motor = figure;
            f_mri_LH_motor.Name = 'Scalp surface and peak LH voxel';
            ft_plot_mesh(meshes(1),'facecolor',[0.5 0.5 0.5],'facealpha',0.3,'edgecolor','none')
            hold on
            scatter3(sourcepos_LH(:, 1),sourcepos_LH(:, 2),sourcepos_LH(:, 3), 'b')
            scatter3(sourcepos_RH(:, 1),sourcepos_RH(:, 2),sourcepos_RH(:, 3), 'r')
            view([180,0])
            f_mri_LH_motor.Color = [1,1,1];
            axis equal
            hold on
            camlight
            drawnow
            view([0 0])

            sourcepos = [sourcepos_LH; sourcepos_RH];
        else
            cfg = [];
            cfg.downsample = 4;
            brainDS = ft_volumedownsample(cfg,segmentedmri);
            voxID = find(brainDS.brain);
            [x, y, z] = ind2sub(brainDS.dim(1:3),voxID);
            sourcepos1 = ft_warp_apply(brainDS.transform,[x y z]);
            sourcepos = sourcepos1/1000; % convert to metres
        end
        %% Covariance full array

        % Covariance matrix
        mu = 0.001;
        for n = size(OPM_trialsf_rest,3):-1:1
            C_1(:,:,n) = cov(OPM_trialsf_playing(:,:,n)') + cov(OPM_trialsf_rebound(:,:,n)') + cov(OPM_trialsf_rest(:,:,n)');
        end
        C_1 = C_1./3; %take mean of covariances over each of the three regions 
        C = mean(C_1,3);
        maxEV = max(svd(C));
        Noise_Cr = min(svd(C)).*eye(size(C));
        C = C + mu.*maxEV.*eye(size(C));
        condnr = cond(C);
        Cinv = inv(C);
        f_cov = figure;
        f_cov.Name = 'Covariance';
        imagesc(C)
        axis square
        title(sprintf('\n mu = %.3f, cond(C) = %.2f',...
            mu,condnr))

        T_type = 'Min'; 

        % Take cov of each trial and average
        for n = size(OPM_trialsf_playing, 3):-1:1
            OPM_data_T_A_COV(:,:,n) = cov(OPM_trialsf_playing(:,:,n)');
            % OPM_data_T_C_COV(:,:,n) =  cov(OPM_trialsf_rebound(:,:,n)');
            OPM_data_T_C_COV(:,:,n) =  cov(OPM_trialsf_rest(:,:,n)');
        end
        Ca = mean(OPM_data_T_A_COV, 3);
        Cc = mean(OPM_data_T_C_COV, 3);

        % Beamforming
        S.mri_file = [F.mri.path F.mri.name F.mri.ext];

        S.meshes_file = [F.mri.path 'meshes.mat'];
        S.sensor_info.pos = Sens_pos;
        S.sensor_info.ors = Sens_ors;
        S.C = C;
        S.Cinv = Cinv;
        S.Noise_Cr = Noise_Cr;
        S.Ca = Ca;
        S.Cc = Cc;
        if do_HFC
            S.M = M;
        end

        if do_tSSS
            S.Lin = Lin;
            S.Lout = Lout;

            S.sss_ins.S=sss_outs.S;
            S.sss_ins.Sin=sss_outs.Sin;
            S.sss_ins.Sout=sss_outs.Sout;
            S.sss_ins.Sn=sss_outs.Sn;
        end

        [bf_outs_f] = run_beamformer_v3_SSS_quick('shell',sourcepos,S,1,T_type,do_HFC, do_tSSS);

        %%
        % Covariance matrix
        mu = 0.001;
        for n = size(OPM_trialsf_rest,3):-1:1
            C_1(:,:,n) = cov(OPM_trialsbb_playing(:,:,n)') + cov(OPM_trialsbb_rebound(:,:,n)') + cov(OPM_trialsbb_rest(:,:,n)');
        end
        C_1 = C_1./3; %take mean of covariances over each of the three regions 
        C = mean(C_1,3);
        maxEV = max(svd(C));
        Noise_Cr = min(svd(C)).*eye(size(C));
        C = C + mu.*maxEV.*eye(size(C));
        condnr = cond(C);
        Cinv = inv(C);
        f_cov = figure;
        f_cov.Name = 'Covariance';
        imagesc(C)
        axis square
        title(sprintf('\n mu = %.3f, cond(C) = %.2f',...
            mu,condnr))

        T_type = 'Min'; 

        % Take cov of each trial and average
        for n = size(OPM_trialsf_playing, 3):-1:1
            OPM_data_T_A_COV(:,:,n) = cov(OPM_trialsf_playing(:,:,n)');
            % OPM_data_T_C_COV(:,:,n) =  cov(OPM_trialsf_rebound(:,:,n)');
            OPM_data_T_C_COV(:,:,n) =  cov(OPM_trialsf_rest(:,:,n)');
        end
        Ca = mean(OPM_data_T_A_COV, 3);
        Cc = mean(OPM_data_T_C_COV, 3);

        % Beamforming
        S.mri_file = [F.mri.path F.mri.name F.mri.ext];

        S.meshes_file = [F.mri.path 'meshes.mat'];
        S.sensor_info.pos = Sens_pos;
        S.sensor_info.ors = Sens_ors;
        S.C = C;
        S.Cinv = Cinv;
        S.Noise_Cr = Noise_Cr;
        S.Ca = Ca;
        S.Cc = Cc;
        if do_HFC
            S.M = M;
        end

        if do_tSSS
            S.Lin = Lin;
            S.Lout = Lout;

            S.sss_ins.S=sss_outs.S;
            S.sss_ins.Sin=sss_outs.Sin;
            S.sss_ins.Sout=sss_outs.Sout;
            S.sss_ins.Sn=sss_outs.Sn;
        end

        [bf_outs] = run_beamformer_v3_SSS_quick('shell',sourcepos,S,1,T_type,do_HFC, do_tSSS);


        clear Cc Ca Cinv C Noise_Cr  f_cov maxEV

        %% Virtual Electrode
        if MARS
            % Do for each centroid
            motor_inds = [20, 22, 23, 24, 26, 27, 61, 63, 64, 65, 67, 68];
            
            % motor_inds = [29, 30]; motor_inds = [motor_inds motor_inds+41];
            % motor_inds = [30 34 35 39 40]; motor_inds = [motor_inds motor_inds+41];

            for mars_loc = motor_inds(1:end/2)
                mars_loc_l = mars_loc;
                mars_loc_r = mars_loc + 41;

                [TFS_l_playing, TFS_r_playing, TFS_l_rebound, TFS_r_rebound, ...
                    TFS_l_rest, TFS_r_rest, H_VE_l_playing, H_VE_r_playing, ...
                    H_VE_l_rebound, H_VE_r_rebound, H_VE_l_rest, H_VE_r_rest, ...
                    TFS_r_playing_all, TFS_l_playing_all, TFS_r_rebound_all, ...
                    TFS_l_rebound_all, TFS_r_rest_all, TFS_l_rest_all, ...
                    Hilbert_envelope_r_playing_trials, Hilbert_envelope_l_playing_trials, ...
                    Hilbert_envelope_r_rebound_trials, Hilbert_envelope_l_rebound_trials, ...
                    Hilbert_envelope_r_rest_trials, Hilbert_envelope_l_rest_trials] = ...
                    show_VE_3_TFS(OPM_data, OPM_data_f_T_playing, OPM_data_f_T_rebound, ...
                    OPM_data_f_T_rest, highpass, lowpass, fre, playing_on_inds, ...
                    rebound_on_inds, rest_on_inds, duration, f, Ntrials, ...
                    trial_time_playing, trig_offset_playing, trig_offset_rebound, ...
                    trig_offset_rest, cbar_lim, meshes, sourcepos, mars_loc_l, ...
                    mars_loc_r, bf_outs, do_regression);


                %%
                if ~isfolder([bids_path,'derivatives\BF_VEs\MARS\' subject '\'])
                    mkdir([bids_path,'derivatives\BF_VEs\MARS\' subject '\'])
                end
                if do_regression
                    regression_save_name = '_regression';
                else
                    regression_save_name ='';
                end

                freq_savename = ['_' num2str(hp) '-' num2str(lp) '_'];

                if save_unaveraged_vars
                    save([bids_path,'derivatives\BF_VEs\MARS\' subject '\' ...
                        run regression_save_name '_3trigs' freq_savename 'unaveraged_location_' ...
                        num2str(mars_loc) '.mat'], 'TFS_r_playing_all', ...
                        'TFS_l_playing_all', 'Hilbert_envelope_r_playing_trials', ...
                        'Hilbert_envelope_l_playing_trials', 'TFS_r_rebound_all', ...
                        'TFS_l_rebound_all', 'Hilbert_envelope_r_rebound_trials', ...
                        'Hilbert_envelope_l_rebound_trials', 'TFS_r_rest_all', ...
                        'TFS_l_rest_all', 'Hilbert_envelope_r_rest_trials', ...
                        'Hilbert_envelope_l_rest_trials');
                end
                if save_vars
                    save([bids_path,'derivatives\BF_VEs\MARS\' subject '\' run ...
                        regression_save_name '_3trigs' freq_savename 'location_'  ...
                        num2str(mars_loc) '.mat'], 'TFS_r_playing', 'TFS_l_playing', ...
                        'H_VE_l_playing', 'H_VE_r_playing', 'TFS_r_rebound', ...
                        'TFS_l_rebound', 'H_VE_l_rebound', 'H_VE_r_rebound', 'TFS_r_rest', ...
                        'TFS_l_rest', 'H_VE_l_rest', 'H_VE_r_rest');
                end
            end

        elseif peak_location
            [TFS_l_playing, TFS_r_playing, TFS_l_rebound, TFS_r_rebound, ...
                TFS_l_rest, TFS_r_rest, H_VE_l_playing, H_VE_r_playing, ...
                H_VE_l_rebound, H_VE_r_rebound, H_VE_l_rest, H_VE_r_rest, ...
                TFS_r_playing_all, TFS_l_playing_all, TFS_r_rebound_all, ...
                TFS_l_rebound_all, TFS_r_rest_all, TFS_l_rest_all, ...
                Hilbert_envelope_r_playing_trials, Hilbert_envelope_l_playing_trials, ...
                Hilbert_envelope_r_rebound_trials, Hilbert_envelope_l_rebound_trials, ...
                Hilbert_envelope_r_rest_trials, Hilbert_envelope_l_rest_trials] = ...
                show_VE_3_TFS(OPM_data, OPM_data_f_T_playing, OPM_data_f_T_rebound, ...
                OPM_data_f_T_rest, highpass, lowpass, fre, playing_on_inds, ...
                rebound_on_inds, rest_on_inds, duration, f, Ntrials, ...
                trial_time_playing, trig_offset_playing, trig_offset_rebound, ...
                trig_offset_rest, cbar_lim, meshes, sourcepos, 1, ...
                2, bf_outs, do_regression);


            %%
            if ~isfolder([bids_path,'derivatives\BF_VEs\Peak_tstat\' subject '\'])
                mkdir([bids_path,'derivatives\BF_VEs\Peak_tstat\' subject '\'])
            end

            if do_regression
                regression_save_name = '_regression';
            else
                regression_save_name ='';
            end

            if save_unaveraged_vars
                save([bids_path,'derivatives\BF_VEs\Peak_tstat\' subject '\' run ...
                    regression_save_name '_3trigs_unaveraged_peak_tstat_cov.mat'], ...
                    'TFS_r_playing_all', 'TFS_l_playing_all', 'Hilbert_envelope_r_playing_trials', ...
                    'Hilbert_envelope_l_playing_trials', 'TFS_r_rebound_all', ...
                    'TFS_l_rebound_all', 'Hilbert_envelope_r_rebound_trials', ...
                    'Hilbert_envelope_l_rebound_trials', 'TFS_r_rest_all', ...
                    'TFS_l_rest_all', 'Hilbert_envelope_r_rest_trials', ...
                    'Hilbert_envelope_l_rest_trials');
            end
            if save_vars
                save([bids_path,'derivatives\BF_VEs\Peak_tstat\' subject '\' run ...
                    regression_save_name '_3trigs_peak_tstat.mat'], 'TFS_r_playing', 'TFS_l_playing', ...
                    'H_VE_l_playing', 'H_VE_r_playing', 'TFS_r_rebound', ...
                    'TFS_l_rebound', 'H_VE_l_rebound', 'H_VE_r_rebound', 'TFS_r_rest', ...
                    'TFS_l_rest', 'H_VE_l_rest', 'H_VE_r_rest');
            end

        else
            left_idx = sourcepos(:,1) < 0;
            right_idx = sourcepos(:,1) > 0;
            if strcmpi(T_type,'Min')
                [v_l, l_l] = min(cell2mat(bf_outs.T_stat).*left_idx');
                [v_r, l_r] = min(cell2mat(bf_outs.T_stat).*right_idx');
            elseif strcmpi(T_type,'Max')
                [v_l, l_l] = max(cell2mat(bf_outs.T_stat).*left_idx');
                [v_r, l_r] = max(cell2mat(bf_outs.T_stat).*right_idx');
            end

          
            [TFS_l_playing, TFS_r_playing, TFS_l_rebound, TFS_r_rebound, ...
                TFS_l_rest, TFS_r_rest, H_VE_l_playing, H_VE_r_playing, ...
                H_VE_l_rebound, H_VE_r_rebound, H_VE_l_rest, H_VE_r_rest, ...
                TFS_r_playing_all, TFS_l_playing_all, TFS_r_rebound_all, ...
                TFS_l_rebound_all, TFS_r_rest_all, TFS_l_rest_all, ...
                Hilbert_envelope_r_playing_trials, Hilbert_envelope_l_playing_trials, ...
                Hilbert_envelope_r_rebound_trials, Hilbert_envelope_l_rebound_trials, ...
                Hilbert_envelope_r_rest_trials, Hilbert_envelope_l_rest_trials] = ...
                show_VE_3_TFS(OPM_data, OPM_data_f_T_playing, OPM_data_f_T_rebound, ...
                OPM_data_f_T_rest, highpass, lowpass, fre, playing_on_inds, ...
                rebound_on_inds, rest_on_inds, duration, f, Ntrials, ...
                trial_time_playing, trig_offset_playing, trig_offset_rebound, ...
                trig_offset_rest, cbar_lim, meshes, sourcepos, l_l, ...
                l_r, bf_outs, do_regression);


            %%
            if ~isfolder([bids_path,'derivatives\BF_VEs\' subject '\'])
                mkdir([bids_path,'derivatives\BF_VEs\' subject '\'])
            end

            if do_regression
                regression_save_name = '_regression';
            else
                regression_save_name ='';
            end

            if save_unaveraged_vars
                save([bids_path,'derivatives\BF_VEs\' subject '\' run ...
                    regression_save_name '_3trigs_unaveraged_MARS_motor_LH.mat'], ...
                    'TFS_r_playing_all', 'TFS_l_playing_all', 'Hilbert_envelope_r_playing_trials', ...
                    'Hilbert_envelope_l_playing_trials', 'TFS_r_rebound_all', ...
                    'TFS_l_rebound_all', 'Hilbert_envelope_r_rebound_trials', ...
                    'Hilbert_envelope_l_rebound_trials', 'TFS_r_rest_all', ...
                    'TFS_l_rest_all', 'Hilbert_envelope_r_rest_trials', ...
                    'Hilbert_envelope_l_rest_trials');
            end
            if save_vars
                save([bids_path,'derivatives\BF_VEs\' subject '\' run ...
                    regression_save_name '_3trigs_central.mat'], 'TFS_r_playing', 'TFS_l_playing', ...
                    'H_VE_l_playing', 'H_VE_r_playing', 'TFS_r_rebound', ...
                    'TFS_l_rebound', 'H_VE_l_rebound', 'H_VE_r_rebound', 'TFS_r_rest', ...
                    'TFS_l_rest', 'H_VE_l_rest', 'H_VE_r_rest');
            end

            %% Save Tstat as a nifti overlay
            if save_overlay
                % set directory to save overlay to
                overlay_filePath = [bids_path 'Derivatives\Tstat\' subject '\'];

                if exist(overlay_filePath, 'dir') == 0
                    mkdir(overlay_filePath);
                end

                image = brainDS;
                image.anatomy = zeros(brainDS.dim(1:3));
                image.anatomy(voxID) = cell2mat(bf_outs_f.T_stat);

                ft_write_mri([overlay_filePath run '_' num2str(hp) '-' num2str(lp) '_2trigs_cov.nii'],image.anatomy,...
                    'dataformat','nifti','transform',image.transform)

                disp('Saved overlay');
            end
        end
    end
end
toc(run_all_timing);
%% Functions
function SensorTransform = coreg_func(bids_path, subject, run, F, mri_meshpath, opts, MNI)
do_coreg = input('Do you want to perform a coreg? 1 for yes, 0 for no: ');
if do_coreg
    % Path to the Python executable (adjust as necessary)
    pythonPath = 'D:\Users\ppyjg12\AppData\Local\anaconda3\envs\meshlab\python.exe';  % E.g., 'C:\Python39\python.exe'

    % Path to the Python script
    scriptPath = "D:\OneDrive - The University of Nottingham\Documents\Python Scripts\run_meshlab.py";  % Adjust this path

    einscan_meshpath = [bids_path 'Einscan\' subject '\' run '.ply'];

    % Run the Python script
    command = sprintf('"%s" "%s"  "%s" "%s', pythonPath, scriptPath, einscan_meshpath, mri_meshpath);
    status = system(command);

    if status == 0
        disp('MeshLab opened successfully.');
    else
        disp('Error opening MeshLab.');
    end

    SensorTransform = zeros(4, 4);

    for i = 1:16
        [row, col] = ind2sub([4, 4], i);
        SensorTransform(row, col) = input(['Enter element ' num2str(row) ', ' num2str(col) ': ']);
    end

    SensorTransform
    while ~input("If transform is correct enter 1, else press 0 to input again: ")
        for i = 1:16
            [row, col] = ind2sub([4, 4], i);
            SensorTransform(row, col) = input(['Enter element ' num2str(row) ', ' num2str(col) ': ']);
        end
        SensorTransform
    end

    if MNI
        einscan_meshpath = strrep(einscan_meshpath, '.ply', '_MNI.ply');
    end

    writematrix(SensorTransform, strrep(einscan_meshpath, 'ply', 'txt'), 'Delimiter', '\t');
    disp('Sensor Transform saved')

    % Command to close MeshLab by its process name
    status = system('taskkill /IM meshlab.exe /F');

    if status == 0
        disp('MeshLab closed successfully.');
    else
        disp('Failed to close MeshLab. Make sure it is running.');
    end
    SensorTransform = readtable(strrep(einscan_meshpath, 'ply', 'txt'), opts);
    SensorTransform = table2array(SensorTransform);
else
    SensorTransform = eye(4, 4);
end
end

function data = read_N1lvm_inscript(filename, dataLines)
%IMPORTFILE Import data from a text file
%  EPILEPSY1RIGID = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  EPILEPSY1RIGID = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  Epilepsy1rigid = importfile("C:\Users\ppzrh\The University of Nottingham\Matthew Brookes (staff) - Brussels\Epi01\CercaHelmet\Epilepsy_1_rigid.lvm", [24, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 11-Jul-2023 15:37:30

% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [24, Inf];
end

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 227);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["LabVIEWMeasurement", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56", "VarName57", "VarName58", "VarName59", "VarName60", "VarName61", "VarName62", "VarName63", "VarName64", "VarName65", "VarName66", "VarName67", "VarName68", "VarName69", "VarName70", "VarName71", "VarName72", "VarName73", "VarName74", "VarName75", "VarName76", "VarName77", "VarName78", "VarName79", "VarName80", "VarName81", "VarName82", "VarName83", "VarName84", "VarName85", "VarName86", "VarName87", "VarName88", "VarName89", "VarName90", "VarName91", "VarName92", "VarName93", "VarName94", "VarName95", "VarName96", "VarName97", "VarName98", "VarName99", "VarName100", "VarName101", "VarName102", "VarName103", "VarName104", "VarName105", "VarName106", "VarName107", "VarName108", "VarName109", "VarName110", "VarName111", "VarName112", "VarName113", "VarName114", "VarName115", "VarName116", "VarName117", "VarName118", "VarName119", "VarName120", "VarName121", "VarName122", "VarName123", "VarName124", "VarName125", "VarName126", "VarName127", "VarName128", "VarName129", "VarName130", "VarName131", "VarName132", "VarName133", "VarName134", "VarName135", "VarName136", "VarName137", "VarName138", "VarName139", "VarName140", "VarName141", "VarName142", "VarName143", "VarName144", "VarName145", "VarName146", "VarName147", "VarName148", "VarName149", "VarName150", "VarName151", "VarName152", "VarName153", "VarName154", "VarName155", "VarName156", "VarName157", "VarName158", "VarName159", "VarName160", "VarName161", "VarName162", "VarName163", "VarName164", "VarName165", "VarName166", "VarName167", "VarName168", "VarName169", "VarName170", "VarName171", "VarName172", "VarName173", "VarName174", "VarName175", "VarName176", "VarName177", "VarName178", "VarName179", "VarName180", "VarName181", "VarName182", "VarName183", "VarName184", "VarName185", "VarName186", "VarName187", "VarName188", "VarName189", "VarName190", "VarName191", "VarName192", "VarName193", "VarName194", "VarName195", "VarName196", "VarName197", "VarName198", "VarName199", "VarName200", "VarName201", "VarName202", "VarName203", "VarName204", "VarName205", "VarName206", "VarName207", "VarName208", "VarName209", "VarName210", "VarName211", "VarName212", "VarName213", "VarName214", "VarName215", "VarName216", "VarName217", "VarName218", "VarName219", "VarName220", "VarName221", "VarName222", "VarName223", "VarName224", "VarName225", "VarName226", "VarName227"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName227", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName227", "EmptyFieldRule", "auto");

% Import the data
tic
disp('Reading table')
data0 = readtable(filename, opts);
toc

% Turning into array
tic
data = [];
for n = size(data0,2)-1:-1:2
    %     disp(n)
    data(:,n) = data0.(['VarName' num2str(n)]);
end
data(:,1) = data0.LabVIEWMeasurement;
toc
end

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
    fzf = zf;
    [pxx,fs] = periodogram(fzf,wind,F,f);
    pow(:,:,j) = sqrt(pxx);
end
Nchans = size(data,1);
chans = 1:Nchans;
po = median(pow(:,chans,:),3);
end

function [TFS_r_playing, TFS_r_rest, H_VE_r_playing, H_VE_r_rest, ...
    TFS_r_playing_all, TFS_r_rest_all, ...
    Hilbert_envelope_r_playing_trials, Hilbert_envelope_r_rest_trials] = ...
    show_VE_2_TFS_RH(OPM_data, OPM_data_f_T_playing, OPM_data_f_T_rest, highpass, lowpass, fre, playing_on_inds, rest_on_inds, duration, ...
    f, Ntrials, trial_time, trig_offset_playing, trig_offset_rest, cbar_lim, meshes, sourcepos, l_r, bf_outs, do_regression)

f_VE_All = figure;
f_VE_All.Name = 'VE';

subplot(221)
trisurf(meshes(3).tri,meshes(3).pnt(:,1),meshes(3).pnt(:,2),meshes(3).pnt(:,3),...
    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7)
hold on
trisurf(meshes(1).tri,meshes(1).pnt(:,1),meshes(1).pnt(:,2),meshes(1).pnt(:,3),...
    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5)
axis equal
plot3(sourcepos(l_r,1),sourcepos(l_r,2),sourcepos(l_r,3),'ro', 'MarkerSize',5,...
    'MarkerFaceColor','r')

W1_r = bf_outs.Weights(:,l_r);
VE_f_r_playing = W1_r'*OPM_data_f_T_playing*1e-9;
VE_f_r_rest = W1_r'*OPM_data_f_T_rest*1e-9;

OPM_ch_r = W1_r'*OPM_data;
OPM_ch_fb_r = zeros(length(fre),length(OPM_ch_r));

if do_regression
    % regression of muscle artefacts
    for fb = length(highpass):-1:56
        noise_tmp_r(fb, :) = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4, ...
            highpass(fb),lowpass(fb),f,1)'))';
    end
    noise_r = mean(noise_tmp_r(56:end, :))';

    ft_progress('init', 'text', 'Calculating TFS bands...')      % ascii progress bar
    for fb = 1:length(highpass)
        ft_progress(fb/length(highpass), 'Processing band %d of %d', fb, length(highpass))
        sig = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4,highpass(fb),lowpass(fb),f,1)));
        beta = pinv(noise_r-mean(noise_r)) * (sig - mean(sig));
        remove = noise_r*beta;
        remove = remove - mean(remove);
        OPM_ch_fb_r(fb,:) = (sig - remove)';
    end
    ft_progress('close')
    regression_string = '- regression';
else
    ft_progress('init', 'text', 'Calculating TFS bands...')      % ascii progress bar
    for fb = 1:length(highpass)
        ft_progress(fb/length(highpass), 'Processing band %d of %d', fb, length(highpass))
        filt_VE_r = nut_filter3(OPM_ch_r','butter','bp',4,highpass(fb),lowpass(fb),f,1)';
        OPM_ch_fb_r(fb,:) = abs(hilbert(filt_VE_r));
    end
    ft_progress('close')
    regression_string = '';
end

OPM_ch_fb_trials_r_playing = zeros(length(fre),duration.*f,Ntrials);
OPM_ch_fb_trials_r_rest = zeros(length(fre),duration.*f,Ntrials);

for fb = 1:length(highpass)
    % Chop data
    for i = 1:Ntrials
        OPM_ch_fb_trials_r_playing(fb,:,i) = OPM_ch_fb_r(fb,playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_r_rest(fb,:,i) = OPM_ch_fb_r(fb,rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
    end

    meanrest_r_TFS = permute(repmat(squeeze(mean(OPM_ch_fb_trials_r_rest, 2)), [1, 1, size(OPM_ch_fb_trials_r_rest, 2)]), [1 3 2]);
    TFS_r_rest_all = (OPM_ch_fb_trials_r_rest-meanrest_r_TFS)./meanrest_r_TFS;
    TFS_r_playing_all = (OPM_ch_fb_trials_r_playing-meanrest_r_TFS)./meanrest_r_TFS;

end

TFS_r_playing = mean(TFS_r_playing_all,3);
TFS_r_rest = mean(TFS_r_rest_all,3);

subplot(224)
pcolor([TFS_r_playing TFS_r_rest]);shading interp
xticks(linspace(0, length([TFS_r_playing TFS_r_rest]), 2*duration+1))
xticklabels([(0:9)+trig_offset_playing (0:10)+trig_offset_rest])
clim([-cbar_lim cbar_lim])
title(['Right Hemisphere ' regression_string])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% axis square
c = colorbar;
c.Label.String = 'Relative Change';
drawnow
xline(length(TFS_r_rest), 'LineWidth', 5)
if do_regression
    % regression for envelope
    sig = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4,13,30,f,1)));

    beta = pinv(noise_r-mean(noise_r)) * (sig - mean(sig));
    remove = noise_r*beta;
    remove = remove - mean(remove);
    beta_regressed_r = (sig - remove)';


    for i = Ntrials:-1:1
        beta_regressed_trials_r_playing(:,i) = beta_regressed_r(playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        beta_regressed_trials_r_rest(:,i) = beta_regressed_r(rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
    end

    % Average across trials
    beta_regressed_mean_r_playing = mean(beta_regressed_trials_r_playing,2);
    beta_regressed_mean_r_rest = mean(beta_regressed_trials_r_rest,2);

    H_VE_r_playing = (beta_regressed_mean_r_playing-mean(beta_regressed_mean_r_playing))./mean(beta_regressed_mean_r_playing)*100;
    H_VE_r_rest = (beta_regressed_mean_r_rest-mean(beta_regressed_mean_r_rest))./mean(beta_regressed_mean_r_rest)*100;

    Hilbert_envelope_r_playing_trials = beta_regressed_trials_r_playing;
    Hilbert_envelope_r_rest_trials = beta_regressed_trials_r_rest;
else


    H_VE_r_rest_orig = mean(reshape(abs(hilbert(VE_f_r_rest)),duration.*f,Ntrials),2);
    H_VE_r_rest = (H_VE_r_rest_orig - mean(H_VE_r_rest_orig))./mean(H_VE_r_rest_orig)*100;

    Hilbert_envelope_r_rest_trials = reshape(abs(hilbert(VE_f_r_rest)),duration.*f,Ntrials);

    H_VE_r_playing = mean(reshape(abs(hilbert(VE_f_r_playing)),duration.*f,Ntrials),2);
    H_VE_r_playing = (H_VE_r_playing - mean(H_VE_r_rest_orig))./mean(H_VE_r_rest_orig)*100;

    Hilbert_envelope_r_playing_trials = reshape(abs(hilbert(VE_f_r_playing)),duration.*f,Ntrials);

end

subplot(222)
hold on
plot([trial_time trial_time+max(trial_time)],[H_VE_r_playing' H_VE_r_rest'],'r', 'DisplayName', 'Right Hemisphere')
xlabel('Time (s)');ylabel('Relative change (%)')
xticks(linspace(0, max([trial_time trial_time+max(trial_time)]), 2*duration+1))
xticklabels([(0:9)+trig_offset_playing (0:10)+trig_offset_rest])
axis square
legend()
drawnow
xline(max(trial_time), 'LineWidth', 5)
end


function [TFS_playing, TFS_rebound, TFS_rest, H_VE_playing, H_VE_rebound, H_VE_rest, ...
    TFS_playing_all, TFS_rebound_all, TFS_rest_all, ...
    Hilbert_envelope_playing_trials, Hilbert_envelope_rebound_trials, Hilbert_envelope_rest_trials] = ...
    show_VE_3_TFS_single(OPM_data, OPM_data_f_T_playing, OPM_data_f_T_rebound, OPM_data_f_T_rest, highpass, lowpass, fre, playing_on_inds, rebound_on_inds, rest_on_inds, duration, ...
    f, Ntrials, trial_time, trig_offset_playing, trig_offset_rebound, trig_offset_rest, cbar_lim, meshes, sourcepos, l_r, bf_outs, do_regression)

f_VE_All = figure;
f_VE_All.Name = 'VE';

subplot(221)
trisurf(meshes(3).tri,meshes(3).pnt(:,1),meshes(3).pnt(:,2),meshes(3).pnt(:,3),...
    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7)
hold on
trisurf(meshes(1).tri,meshes(1).pnt(:,1),meshes(1).pnt(:,2),meshes(1).pnt(:,3),...
    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5)
axis equal
plot3(sourcepos(l_r,1),sourcepos(l_r,2),sourcepos(l_r,3),'ro', 'MarkerSize',5,...
    'MarkerFaceColor','r')

W1_r = bf_outs.Weights(:,l_r);
VE_f_playing = W1_r'*OPM_data_f_T_playing*1e-9;
VE_f_rebound = W1_r'*OPM_data_f_T_rebound*1e-9;
VE_f_rest = W1_r'*OPM_data_f_T_rest*1e-9;

OPM_ch_r = W1_r'*OPM_data;
OPM_ch_fb_r = zeros(length(fre),length(OPM_ch_r));

if do_regression
    % regression of muscle artefacts
    for fb = length(highpass):-1:56
        noise_tmp_r(fb, :) = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4, ...
            highpass(fb),lowpass(fb),f,1)'))';
    end
    noise_r = mean(noise_tmp_r(56:end, :))';

    ft_progress('init', 'text', 'Calculating TFS bands...')      % ascii progress bar
    for fb = 1:length(highpass)
        ft_progress(fb/length(highpass), 'Processing band %d of %d', fb, length(highpass))
        sig = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4,highpass(fb),lowpass(fb),f,1)));
        beta = pinv(noise_r-mean(noise_r)) * (sig - mean(sig));
        remove = noise_r*beta;
        remove = remove - mean(remove);
        OPM_ch_fb_r(fb,:) = (sig - remove)';
    end
    ft_progress('close')
    regression_string = '- regression';
else
    ft_progress('init', 'text', 'Calculating TFS bands...')      % ascii progress bar
    for fb = 1:length(highpass)
        ft_progress(fb/length(highpass), 'Processing band %d of %d', fb, length(highpass))
        filt_VE_r = nut_filter3(OPM_ch_r','butter','bp',4,highpass(fb),lowpass(fb),f,1)';
        OPM_ch_fb_r(fb,:) = abs(hilbert(filt_VE_r));
    end
    ft_progress('close')
    regression_string = '';
end

OPM_ch_fb_trials_playing = zeros(length(fre),duration.*f,Ntrials);
OPM_ch_fb_trials_rebound = zeros(length(fre),duration.*f,Ntrials);
OPM_ch_fb_trials_rest = zeros(length(fre),duration.*f,Ntrials);

for fb = 1:length(highpass)
    % Chop data
    for i = 1:Ntrials
        OPM_ch_fb_trials_playing(fb,:,i) = OPM_ch_fb_r(fb,playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_rebound(fb,:,i) = OPM_ch_fb_r(fb,rebound_on_inds(i):rebound_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_rest(fb,:,i) = OPM_ch_fb_r(fb,rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
    end


    meanrest_TFS = permute(repmat(squeeze(mean(OPM_ch_fb_trials_rest, 2)), [1, 1, size(OPM_ch_fb_trials_rest, 2)]), [1 3 2]);
    TFS_rebound_all = (OPM_ch_fb_trials_rebound-meanrest_TFS)./meanrest_TFS;
    TFS_rest_all = (OPM_ch_fb_trials_rest-meanrest_TFS)./meanrest_TFS;
    TFS_playing_all = (OPM_ch_fb_trials_playing-meanrest_TFS)./meanrest_TFS;

end

TFS_playing = mean(TFS_playing_all,3);
TFS_rebound = mean(TFS_rebound_all,3);
TFS_rest = mean(TFS_rest_all,3);


subplot(224)
pcolor([TFS_playing TFS_rebound TFS_rest]);shading interp
xticks(linspace(0, length([TFS_playing TFS_rebound TFS_rest]), 3*duration+1))
xticklabels([(0:9)+trig_offset_playing (0:9)+trig_offset_rebound (0:10)+trig_offset_rest])
clim([-cbar_lim cbar_lim])
title(['Right Hemisphere ' regression_string])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% axis square
c = colorbar;
c.Label.String = 'Relative Change';
drawnow
xline(length(TFS_playing), 'LineWidth', 5)
xline(length(TFS_playing) + length(TFS_rebound), 'LineWidth', 5)
if do_regression
    % regression for envelope
    sig = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4,13,30,f,1)));

    beta = pinv(noise_r-mean(noise_r)) * (sig - mean(sig));
    remove = noise_r*beta;
    remove = remove - mean(remove);
    beta_regressed_r = (sig - remove)';


    for i = Ntrials:-1:1
        beta_regressed_trials_playing(:,i) = beta_regressed_r(playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        beta_regressed_trials_rebound(:,i) = beta_regressed_r(rebound_on_inds(i):rebound_on_inds(i) + duration*f - 1);
        beta_regressed_trials_rest(:,i) = beta_regressed_r(rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
    end

    % Average across trials
    beta_regressed_mean_playing = mean(beta_regressed_trials_playing,2);
    beta_regressed_mean_rebound = mean(beta_regressed_trials_rebound,2);
    beta_regressed_mean_rest = mean(beta_regressed_trials_rest,2);

    H_VE_playing = (beta_regressed_mean_playing-mean(beta_regressed_mean_playing))./mean(beta_regressed_mean_playing)*100;
    H_VE_rebound = (beta_regressed_mean_rebound-mean(beta_regressed_mean_rest))./mean(beta_regressed_mean_rest)*100;
    H_VE_rest = (beta_regressed_mean_rest-mean(beta_regressed_mean_rest))./mean(beta_regressed_mean_rest)*100;

    Hilbert_envelope_playing_trials = beta_regressed_trials_playing;
    Hilbert_envelope_rebound_trials = beta_regressed_trials_rebound;
    Hilbert_envelope_rest_trials = beta_regressed_trials_rest;
else


    H_VE_rest_orig = mean(reshape(abs(hilbert(VE_f_rest)),duration.*f,Ntrials),2);
    H_VE_rest = (H_VE_rest_orig - mean(H_VE_rest_orig))./mean(H_VE_rest_orig)*100;

    Hilbert_envelope_rest_trials = reshape(abs(hilbert(VE_f_rest)),duration.*f,Ntrials);

    H_VE_rebound = mean(reshape(abs(hilbert(VE_f_rebound)),duration.*f,Ntrials),2);
    H_VE_rebound = (H_VE_rebound - mean(H_VE_rest_orig))./mean(H_VE_rest_orig)*100;

    Hilbert_envelope_rebound_trials = reshape(abs(hilbert(VE_f_rebound)),duration.*f,Ntrials);

    H_VE_playing = mean(reshape(abs(hilbert(VE_f_playing)),duration.*f,Ntrials),2);
    H_VE_playing = (H_VE_playing - mean(H_VE_rest_orig))./mean(H_VE_rest_orig)*100;

    Hilbert_envelope_playing_trials = reshape(abs(hilbert(VE_f_playing)),duration.*f,Ntrials);

end

subplot(222)
hold on
plot([trial_time trial_time+max(trial_time) trial_time+2*max(trial_time)],[H_VE_playing' H_VE_rebound' H_VE_rest'],'r', 'DisplayName', 'Right Hemisphere')
xlabel('Time (s)');ylabel('Relative change (%)')
xticks(linspace(0, max([trial_time+2*max(trial_time)]), 3*duration+1))
xticklabels([(0:9)+trig_offset_playing (0:9)+trig_offset_rebound (0:10)+trig_offset_rest])
axis square
legend()
drawnow
xline(max(trial_time), 'LineWidth', 5)
xline(2*max(trial_time), 'LineWidth', 5)
end


function[TFS_l_playing, TFS_r_playing, TFS_l_rebound, TFS_r_rebound, ...
    TFS_l_rest, TFS_r_rest, H_VE_l_playing, H_VE_r_playing, ...
    H_VE_l_rebound, H_VE_r_rebound, H_VE_l_rest, H_VE_r_rest, ...
    TFS_r_playing_all, TFS_l_playing_all, TFS_r_rebound_all, ...
    TFS_l_rebound_all, TFS_r_rest_all, TFS_l_rest_all, ...
    Hilbert_envelope_r_playing_trials, Hilbert_envelope_l_playing_trials, ...
    Hilbert_envelope_r_rebound_trials, Hilbert_envelope_l_rebound_trials, ...
    Hilbert_envelope_r_rest_trials, Hilbert_envelope_l_rest_trials] =...
    show_VE_3_TFS(OPM_data, OPM_data_f_T_playing, OPM_data_f_T_rebound, ...
    OPM_data_f_T_rest, highpass, lowpass, fre, playing_on_inds, ...
    rebound_on_inds, rest_on_inds, duration, f, Ntrials, ...
    trial_time, trig_offset_playing, trig_offset_rebound, ...
    trig_offset_rest, cbar_lim, meshes, sourcepos, l_l, ...
    l_r, bf_outs, do_regression)
f_VE_All = figure;
f_VE_All.Name = 'VE';

subplot(221)
trisurf(meshes(3).tri,meshes(3).pnt(:,1),meshes(3).pnt(:,2),meshes(3).pnt(:,3),...
    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.7)
hold on
trisurf(meshes(1).tri,meshes(1).pnt(:,1),meshes(1).pnt(:,2),meshes(1).pnt(:,3),...
    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5)
axis equal
plot3(sourcepos(l_l,1),sourcepos(l_l,2),sourcepos(l_l,3),'bo', 'MarkerSize',5,...
    'MarkerFaceColor','b')
plot3(sourcepos(l_r,1),sourcepos(l_r,2),sourcepos(l_r,3),'ro', 'MarkerSize',5,...
    'MarkerFaceColor','r')

W1_l = bf_outs.Weights(:,l_l);
W1_r = bf_outs.Weights(:,l_r);
VE_f_l_playing = W1_l'*OPM_data_f_T_playing*1e-9;
VE_f_r_playing = W1_r'*OPM_data_f_T_playing*1e-9;
VE_f_l_rebound = W1_l'*OPM_data_f_T_rebound*1e-9;
VE_f_r_rebound = W1_r'*OPM_data_f_T_rebound*1e-9;
VE_f_l_rest = W1_l'*OPM_data_f_T_rest*1e-9;
VE_f_r_rest = W1_r'*OPM_data_f_T_rest*1e-9;

OPM_ch_l = W1_l'*OPM_data;
OPM_ch_fb_l = zeros(length(fre),length(OPM_ch_l));
OPM_ch_r = W1_r'*OPM_data;
OPM_ch_fb_r = zeros(length(fre),length(OPM_ch_r));

if do_regression
    error('Regression not working yet!')
    % regression of muscle artefacts
    for fb = length(highpass):-1:56
        noise_tmp_r(fb, :) = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4, ...
            highpass(fb),lowpass(fb),f,1)'))';
        noise_tmp_l(fb, :) = abs(hilbert(nut_filter3(OPM_ch_l','butter','bp',4, ...
            highpass(fb),lowpass(fb),f,1)'))';
    end
    noise_r = mean(noise_tmp_r(56:end, :))';
    noise_l = mean(noise_tmp_l(56:end, :))';

    ft_progress('init', 'text', 'Calculating TFS bands...')      % ascii progress bar
    for fb = 1:length(highpass)
        ft_progress(fb/length(highpass), 'Processing band %d of %d', fb, length(highpass))
        sig = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4,highpass(fb),lowpass(fb),f,1)));
        beta = pinv(noise_r-mean(noise_r)) * (sig - mean(sig));
        remove = noise_r*beta;
        remove = remove - mean(remove);
        OPM_ch_fb_r(fb,:) = (sig - remove)';

        sig = abs(hilbert(nut_filter3(OPM_ch_l','butter','bp',4,highpass(fb),lowpass(fb),f,1)));
        beta = pinv(noise_l-mean(noise_l)) * (sig - mean(sig));
        remove = noise_l*beta;
        remove = remove - mean(remove);
        OPM_ch_fb_l(fb,:) = (sig - remove)';
    end
    ft_progress('close')
    regression_string = '- regression';
else
    ft_progress('init', 'text', 'Calculating TFS bands...')      % ascii progress bar
    for fb = 1:length(highpass)
        ft_progress(fb/length(highpass), 'Processing band %d of %d', fb, length(highpass))
        filt_VE_r = nut_filter3(OPM_ch_r','butter','bp',4,highpass(fb),lowpass(fb),f,1)';
        OPM_ch_fb_r(fb,:) = abs(hilbert(filt_VE_r));
        filt_VE_l = nut_filter3(OPM_ch_l','butter','bp',4,highpass(fb),lowpass(fb),f,1)';
        OPM_ch_fb_l(fb,:) = abs(hilbert(filt_VE_l));
    end
    ft_progress('close')
    regression_string = '';
end

OPM_ch_fb_trials_r_playing = zeros(length(fre), duration.*f,Ntrials);
OPM_ch_fb_trials_l_playing = zeros(length(fre), duration.*f,Ntrials);
OPM_ch_fb_trials_r_rebound = zeros(length(fre), duration.*f,Ntrials);
OPM_ch_fb_trials_l_rebound = zeros(length(fre), duration.*f,Ntrials);
OPM_ch_fb_trials_r_rest = zeros(length(fre), duration.*f,Ntrials);
OPM_ch_fb_trials_l_rest = zeros(length(fre), duration.*f,Ntrials);

for fb = 1:length(highpass)
    % Chop data
    for i = 1:Ntrials
        OPM_ch_fb_trials_r_playing(fb,:,i) = OPM_ch_fb_r(fb,playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_l_playing(fb,:,i) = OPM_ch_fb_l(fb,playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_r_rebound(fb,:,i) = OPM_ch_fb_r(fb,rebound_on_inds(i):rebound_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_l_rebound(fb,:,i) = OPM_ch_fb_l(fb,rebound_on_inds(i):rebound_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_r_rest(fb,:,i) = OPM_ch_fb_r(fb,rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
        OPM_ch_fb_trials_l_rest(fb,:,i) = OPM_ch_fb_l(fb,rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
    end

    meanrest_r_TFS = permute(repmat(squeeze(mean(OPM_ch_fb_trials_r_rest, 2)), [1, 1, size(OPM_ch_fb_trials_r_rest, 2)]), [1 3 2]);
    TFS_r_rest_all = (OPM_ch_fb_trials_r_rest-meanrest_r_TFS)./meanrest_r_TFS;
    TFS_r_rebound_all = (OPM_ch_fb_trials_r_rebound-meanrest_r_TFS)./meanrest_r_TFS;
    TFS_r_playing_all = (OPM_ch_fb_trials_r_playing-meanrest_r_TFS)./meanrest_r_TFS;

    meanrest_l_TFS = permute(repmat(squeeze(mean(OPM_ch_fb_trials_l_rest, 2)), [1, 1, size(OPM_ch_fb_trials_l_playing, 2)]), [1 3 2]);
    TFS_l_rest_all = (OPM_ch_fb_trials_l_rest-meanrest_l_TFS)./meanrest_l_TFS;
    TFS_l_rebound_all = (OPM_ch_fb_trials_l_rebound-meanrest_l_TFS)./meanrest_l_TFS;
    TFS_l_playing_all = (OPM_ch_fb_trials_l_playing-meanrest_l_TFS)./meanrest_l_TFS;
end

TFS_l_playing = mean(TFS_l_playing_all,3);
TFS_l_rebound = mean(TFS_l_rebound_all,3);
TFS_l_rest = mean(TFS_l_rest_all,3);

subplot(223)
pcolor([TFS_l_playing TFS_l_rebound TFS_l_rest]);shading interp
xticks(linspace(0, length([TFS_l_playing TFS_l_rebound TFS_l_rest]), 3*duration+1))
xticklabels([(0:9)+trig_offset_playing (0:9)+trig_offset_rebound (0:10)+trig_offset_rest])
clim([-cbar_lim cbar_lim])
title(['Left Hemisphere ' regression_string])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% axis square
c = colorbar;
c.Label.String = 'Relative Change';
drawnow
xline(length(TFS_l_playing), 'LineWidth', 5)
xline(length(TFS_l_playing) + length(TFS_l_rebound), 'LineWidth', 5)


TFS_r_playing = mean(TFS_r_playing_all,3);
TFS_r_rebound = mean(TFS_r_rebound_all,3);
TFS_r_rest = mean(TFS_r_rest_all,3);

subplot(224)
pcolor([TFS_r_playing TFS_r_rebound TFS_r_rest]);shading interp
xticks(linspace(0, length([TFS_r_playing TFS_r_rebound TFS_r_rest]), 3*duration+1))
xticklabels([(0:9)+trig_offset_playing (0:9)+trig_offset_rebound (0:10)+trig_offset_rest])
clim([-cbar_lim cbar_lim])
title(['Right Hemisphere ' regression_string])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% axis square
c = colorbar;
c.Label.String = 'Relative Change';
drawnow
xline(length(TFS_r_playing), 'LineWidth', 5)
xline(length(TFS_r_playing) + length(TFS_r_rebound), 'LineWidth', 5)


if do_regression
    % regression for envelope
    sig = abs(hilbert(nut_filter3(OPM_ch_r','butter','bp',4,13,30,f,1)));

    beta = pinv(noise_r-mean(noise_r)) * (sig - mean(sig));
    remove = noise_r*beta;
    remove = remove - mean(remove);
    beta_regressed_r = (sig - remove)';

    sig = abs(hilbert(nut_filter3(OPM_ch_l','butter','bp',4,13,30,f,1)));

    beta = pinv(noise_l-mean(noise_l)) * (sig - mean(sig));
    remove = noise_l*beta;
    remove = remove - mean(remove);
    beta_regressed_l = (sig - remove)';

    for i = Ntrials:-1:1
        beta_regressed_trials_r_playing(:,i) = beta_regressed_r(playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        beta_regressed_trials_r_rest(:,i) = beta_regressed_r(rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
        beta_regressed_trials_l_playing(:,i) = beta_regressed_l(playing_on_inds(i):playing_on_inds(i) + duration*f - 1);
        beta_regressed_trials_l_rest(:,i) = beta_regressed_l(rest_on_inds(i):rest_on_inds(i) + duration*f - 1);
    end

    % Average across trials
    beta_regressed_mean_r_playing = mean(beta_regressed_trials_r_playing,2);
    beta_regressed_mean_r_rest = mean(beta_regressed_trials_r_rest,2);
    beta_regressed_mean_l_playing = mean(beta_regressed_trials_l_playing,2);
    beta_regressed_mean_l_rest = mean(beta_regressed_trials_l_rest,2);

    H_VE_r_rest = (beta_regressed_mean_r_rest-mean(beta_regressed_mean_r_rest))./mean(beta_regressed_mean_r_rest)*100;
    H_VE_r_playing = (beta_regressed_mean_r_playing-mean(beta_regressed_mean_r_rest))./mean(beta_regressed_mean_r_rest)*100;
    H_VE_l_rest = (beta_regressed_mean_l_rest-mean(beta_regressed_mean_l_rest))./mean(beta_regressed_mean_l_rest)*100;
    H_VE_l_playing = (beta_regressed_mean_l_playing-mean(beta_regressed_mean_l_rest))./mean(beta_regressed_mean_l_rest)*100;

    Hilbert_envelope_r_playing_trials = beta_regressed_trials_r_playing;
    Hilbert_envelope_r_rest_trials = beta_regressed_trials_r_rest;
    Hilbert_envelope_l_playing_trials = beta_regressed_trials_l_playing;
    Hilbert_envelope_l_rest_trials = beta_regressed_trials_l_rest;
else
    H_VE_l_rest_orig = mean(reshape(abs(hilbert(VE_f_l_rest)),duration.*f,Ntrials),2);
    H_VE_l_rest = (H_VE_l_rest_orig - mean(H_VE_l_rest_orig))./mean(H_VE_l_rest_orig)*100;
    H_VE_r_rest_orig = mean(reshape(abs(hilbert(VE_f_r_rest)),duration.*f,Ntrials),2);
    H_VE_r_rest = (H_VE_r_rest_orig - mean(H_VE_r_rest_orig))./mean(H_VE_r_rest_orig)*100;

    Hilbert_envelope_r_rest_trials = reshape(abs(hilbert(VE_f_r_rest)),duration.*f,Ntrials);
    Hilbert_envelope_l_rest_trials = reshape(abs(hilbert(VE_f_l_rest)),duration.*f,Ntrials);

    H_VE_l_playing = mean(reshape(abs(hilbert(VE_f_l_playing)),duration.*f,Ntrials),2);
    H_VE_l_playing = (H_VE_l_playing - mean(H_VE_l_rest_orig))./mean(H_VE_l_rest_orig)*100;
    H_VE_r_playing = mean(reshape(abs(hilbert(VE_f_r_playing)),duration.*f,Ntrials),2);
    H_VE_r_playing = (H_VE_r_playing - mean(H_VE_r_rest_orig))./mean(H_VE_r_rest_orig)*100;

    Hilbert_envelope_l_playing_trials = reshape(abs(hilbert(VE_f_l_playing)),duration.*f,Ntrials);
    Hilbert_envelope_r_playing_trials = reshape(abs(hilbert(VE_f_r_playing)),duration.*f,Ntrials);

    H_VE_l_rebound = mean(reshape(abs(hilbert(VE_f_l_rebound)),duration.*f,Ntrials),2);
    H_VE_l_rebound = (H_VE_l_rebound - mean(H_VE_l_rest_orig))./mean(H_VE_l_rest_orig)*100;
    H_VE_r_rebound = mean(reshape(abs(hilbert(VE_f_r_rebound)),duration.*f,Ntrials),2);
    H_VE_r_rebound = (H_VE_r_rebound - mean(H_VE_r_rest_orig))./mean(H_VE_r_rest_orig)*100;

    Hilbert_envelope_l_rebound_trials = reshape(abs(hilbert(VE_f_l_rebound)),duration.*f,Ntrials);
    Hilbert_envelope_r_rebound_trials = reshape(abs(hilbert(VE_f_r_rebound)),duration.*f,Ntrials);


end

subplot(222)
hold on
plot([trial_time trial_time+max(trial_time) trial_time+2*max(trial_time)],[H_VE_l_playing' H_VE_l_rebound' H_VE_l_rest'],'b', 'DisplayName', 'Left Hemisphere')
plot([trial_time trial_time+max(trial_time) trial_time+2*max(trial_time)],[H_VE_r_playing' H_VE_r_rebound' H_VE_r_rest'],'r', 'DisplayName', 'Right Hemisphere')
xlabel('Time (s)');ylabel('Relative change (%)')

xticks(linspace(0, max([trial_time+2*max(trial_time)]), 3*duration+1))
xticklabels([(0:9)+trig_offset_playing (0:9)+trig_offset_rebound (0:10)+trig_offset_rest])

legend()
drawnow
xline(max(trial_time), 'LineWidth', 5)
xline(2*max(trial_time), 'LineWidth', 5)
end

