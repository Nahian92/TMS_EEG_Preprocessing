%% Set Directories
wpms.plugins = 'Insert EEGLAB and FieldTrip Directory'
wpms.chanlocs = 'Insert channel locations directory'
wpms.data = 'Insert data directory with each subject ID identifier starting with 0'

%% Add all plugins
addpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1'])
addpath(genpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\functions']));
addpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\TESA_v1.1.0'])
addpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\FastICA_25'])
addpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\xdfimport1.16']);
addpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\bva-io1.5.13']);
addpath(genpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\ARfitStudio0.41']));
addpath(genpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\mcolonFolder']));
addpath(genpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\SIFT1.52']));
addpath(genpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\PrepPipeline0.55.3']));
addpath(genpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\firfilt']));
addpath(genpath([wpms.plugins filesep 'EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\plugins\grandaverage0.9']));
addpath(([wpms.plugins filesep 'EEGLAB_FieldTrip\fieldtrip-20200215']));
addpath([wpms.plugins filesep 'EEGLAB_FieldTrip\fieldtrip-20200215\external\eeglab'])
addpath([wpms.plugins '\EEGLAB_FieldTrip\eeglab_current\eeglab2019_1\functions\sigprocfunc'])
%% Load File

[EEG, ~] = pop_loadbv([wpms.data filsep 'Insert Subject Data Folder', 'Insert Subject Data File.vhdr'); 
load([wpms.chanlocs filsep 'chanlocs.mat'])
EEG.chanlocs = chanlocs;

%% Inspect EEG Data for Bad channels
pop_eegplot(EEG, 1, 1, 1)

%% Remove Bad Channels
EEG = pop_select( EEG,'nochannel',[Insert Bad Channel Numbers]); 

%% Use arfit Studio to remove and interpolate TMS pulse artefact
arfitStudio;

%% Downsample to 1000Hz
EEG = pop_resample( EEG, 1000);

%% Epoch the data, baseline correct, then remove Bad trials

tms_trigger_values = {'S 14', 'S 13','S 11'}; %name of TMS trigger for each condition
conditions = {'Prepain', 'Pain','Postpain'};
subject_folder_identifier = '0';
subject_folders = dir(fullfile(wpms.data,[subject_folder_identifier '*']));

subject_index = %Insert Subject ID
current_subject = subject_folders(subject_index).name;
current_subject_folder = [subject_folders(subject_index).folder '\' current_subject];
for cond = 1:length(conditions);
    EEG = pop_loadset('filename', [current_subject '.set'], 'filepath', current_subject_folder);
    if subject_index == 1
        for i = 1:length(EEG.event);
            this_event = string(EEG.event(i).type);
            if  EEG.event(i).type == "Pre Pain";
                EEG.event(i).type = 'S 14';
            elseif EEG.event(i).type == "Pain";
                EEG.event(i).type = 'S 13';
            elseif EEG.event(i).type == "Post Pain";
                EEG.event(i).type = 'S 11';
            end
        end
    end
    EEG = pop_epoch(EEG, tms_trigger_values(cond), [-1, 1], 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-1000 -5]);
    EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,5,5,0,0);
    pop_rejmenu(EEG,1);
    pause_script = input('Highlight bad trials, update marks and then press enter');
    EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
    EEG = pop_rejepoch( EEG, EEG.BadTr ,0);
    Data{1, cond} =  EEG;
end
save(fullfile(current_subject_folder ,'\badchan_arfit_downsample_epoch_badtrials'), 'Data')

%% Use ICA, SOUND and Filtering to clean the data

close all
subject_index =  % change this to the participant number you want to analyze
current_subject = subject_folders(subject_index).name;
current_subject_folder = [subject_folders(subject_index).folder '\' current_subject];
load([current_subject_folder filesep 'badchan_arfit_downsample_epoch_badtrials.mat'])
for cond = 1:length(conditions);
    EEG = Data{cond}
    EEG = pop_epoch(EEG, tms_trigger_values(cond), [-1, 1], 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-1000 -5]);
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
    EEG = pop_tesa_compselect( EEG,'compCheck','off','comps',[],'figSize','large',...
        'plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh',...
        8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',...
        2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','on','moveThresh',...
        2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','elecNoise','off',...
        'elecNoiseThresh',4,'elecNoiseFeedback','off' );

    EEG = pop_tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 5, 'leadfieldInFile', [], 'leadfieldChansFile', [], 'replaceChans', [], 'multipleConds', [])

    EEG = pop_tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' );
    EEG = pop_tesa_filtbutter( EEG, 48, 52, 4, 'bandstop' );
    EEG = pop_interp( EEG, chanlocs );

    save(fullfile(current_subject_folder ,'\final_TEPs'), 'EEG')

end

