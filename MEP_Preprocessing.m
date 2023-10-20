clear all;
%% Set-up Workspace %%
cwd = [pwd filesep]; % Store current working directory
wpms = []; % pre-allocate a workspace variables in struct

% Locate folders of interest - Change if necessary
wpms.DATAIN     = fullfile(cwd, [filesep 'Data'], filesep); % /../ skips back a folder
wpms.DATAOUT    = [cwd 'MEP_Data_Output' filesep];

% Structure and store subject codes
subjlist = dir([wpms.DATAIN '0*']);

% Prepare OUTPUT folders - Do this only once to set-up. Comment out afterwards
mkdir(wpms.DATAOUT); % Create Output directory
for subout = 1:length(subjlist)
    mkdir([wpms.DATAOUT subjlist(subout).name]) % Create folder for each px
end

% Experiment variables of interest - structure and store
exp           = [];
exp.blocks  = {'prepain','pain','postpain'};

%% Run script %%

% Run script looping through each subject code. 'px' is used becuase it
% stands out more within the loop amoung the variable names

for px = 1:length(subjlist) %Change the 1 here if you want start at a specific participant e.g. change 1 to 16
    
    % print out subject codes to review
    fprintf(['\n Analysing participant: ' subjlist(px).name '\n\n']);
    
    % Look for subject data files
    spikefiles  = dir([wpms.DATAIN subjlist(px).name filesep '*pain.mat']);
    
    % Loop through sessions, making sure the sessions match and there is
    % not more than one file for each session. 'block' is used becuase it
    % stands out more within the loop amoung variable names
    
    for block = 1:length(spikefiles) %Change the 1 here if you want start at a specific session e.g. change 1 to 3 for the 3rd session
        
        % print out session codes to review
        fprintf(['\n Analysing session: ' exp.blocks{block} '\n\n']);
        
        % Check this data has not already been analysed. If so skip this
        % participant
        if exist([wpms.DATAOUT subjlist(px).name filesep subjlist(px).name '_' exp.blocks{block} '_TrialData.mat']) % TODO: Change to what the final data will be called.
            fprintf(['\n This participants ' exp.blocks{block} ' data has already been processed. Moving on... \n\n']);
            
        else
            % Print out exact files to be loaded for QC
            fprintf(['Working on the following spike file: ' spikefiles(block).name '\n'])           
            Spike_data = load([spikefiles(block).folder filesep spikefiles(block).name]);

            %% Analysis Settings
            
            % The following variables are hard coded and can be adjusted at a later
            % date if so desired
            plot_out_lims = [0 0 0.75 0.75];
            plot_range = [0.12, 0.12];
            
            %% Organize the Data
            Spike_data   = struct2cell(Spike_data);
            n_trials     = length(Spike_data{2}.times);
            tms_times    = round(Spike_data{2}.times / 0.0005 ) * 0.0005; %
            Data         = [];
            Data(:,1)    = Spike_data{1}.values;
            Data(:,2)    = Spike_data{1}.times;
                        
            %% Determine offsets and onsets and background_emg periods trial by trial
            
            fprintf('Starting timing definitions... \n')
            if exist([wpms.DATAOUT subjlist(px).name filesep subjlist(px).name '_' exp.blocks{block} '_TempCursor.mat']) % TODO: Change to what the final data will be called.
                load([wpms.DATAOUT subjlist(px).name filesep subjlist(px).name '_' exp.blocks{block} '_TempCursor.mat']);
                first_trial = length(info_struct) +1;
                fprintf(['\n This participant has existing trial cursors selected. Starting from next trial \n\n']);
            else
                first_trial = 1;
            end
            
            for this_trial = first_trial:n_trials
                fig = figure('units','normalized','outerposition',[0 0 0.75 0.75]);
                % save the start / stop data within the figure,
                % indice 1: 1 = bkg, 2 = mep
                
                setappdata(fig,'startStop',nan(5,2));
                setappdata(fig,'rectified',false);
                
                % marking completed flag to be turn true when enter is pressed?
                setappdata(fig,'markingCompleted',false);
                
                % set the current period (if only one point is selected
                setappdata(fig,'currentPlot',[]);
                setappdata(fig,'currentIndex',[]);
                
                mep = Data((Data(:,2)>(tms_times(this_trial)-plot_range(1))) & (Data(:,2)<(tms_times(this_trial)+plot_range(2))),:);
                mep(:,2) = mep(:,2)-tms_times(this_trial);
                
                dataPlot = line(mep(:,2), mep(:,1), 'LineWidth',1,'Color','black','Tag','data');
                
                % plot the 3 periods as NaN, so it doesn't show up on the plot
                region = 5; % 
                regionName = {'BKGPlot','MEPPlot'};
                lineDescription = {'Data','Background','Magnetically Evoked Potential'};
                regionColour = {[0 1 0],[0 1 1],[0 0.5 0.5],[1 0 1],[0.5 0 0.5]}; %green, light blue, blue, light red, red
                
                for jj = 1:region
                    regionPlot(jj) = line(mep(:,2), nan(length(mep),1), 'LineWidth',1,'Color',regionColour{jj},'Tag',regionName{jj},'marker','o');
                end
                
             
                title(['Trial ' num2str(this_trial)])
                xline(0, '--r');
                upper_limit = max(mep((mep(:,2)>(0.007)),:))*4;
                lower_limit = min(mep((mep(:,2)>(0.007)),:))*4;
                xlabel('Time relative to TMS (s)');
                ylabel('Amplitude (V)');
                ylim([lower_limit(1) upper_limit(1)]);
                xlim([-0.08 0.08])
                ax = gca;
                ax.Toolbar.Visible = 'off';
             
                %     ax.Interactions = [zoomInteraction regionZoomInteraction];
                vertical_cursors(fig,regionName);
                legend([dataPlot,regionPlot], lineDescription)
                waitfor(fig,'ToolBar','none')

                if ishghandle('fig')
                    % todo: if figure get prematurely closed stop gracefully
                    startStop = getappdata(fig,'startStop');
                else
                    startStop = getappdata(fig,'startStop');
                end   
                
                
                for kk = 1:size(startStop,1)
                    if all(~isnan(startStop(kk,:)))
                        info_struct(this_trial,kk).Position = mep(startStop(kk,:),2)';
                        info_struct(this_trial,kk).Index = startStop(kk,:);
                    else
                        info_struct(this_trial,kk).Position = [];
                        info_struct(this_trial,kk).Index = [];
                    end
                end
                
                
                
                if all(~isnan(startStop(1,:)))
                    info_struct(this_trial,1).Position = mep(startStop(1,:),2)';
                else
                    info_struct(this_trial,1).Position = [];
                end
                
                if all(~isnan(startStop(2,:)))
                    info_struct(this_trial,2).Position = mep(startStop(2,:),2)';
                else
                    info_struct(this_trial,2).Position = [];
                end
                
                if all(~isnan(startStop(3,:)))
                    info_struct(this_trial,3).Position = mep(startStop(3,:),2)';
                else
                    info_struct(this_trial,3).Position = [];
                end
                
                
                if all(~isnan(startStop(4,:)))
                    info_struct(this_trial,4).Position = mep(startStop(4,:),2)';
                else
                    info_struct(this_trial,4).Position = [];
                end
                
                
                if all(~isnan(startStop(5,:)))
                    info_struct(this_trial,5).Position = mep(startStop(5,:),2)';
                else
                    info_struct(this_trial,5).Position = [];
                end
                
                save([wpms.DATAOUT subjlist(px).name filesep subjlist(px).name '_' exp.blocks{block} '_TempCursor.mat'], 'info_struct');% 
                switch fig.Tag
                    case 'Next'
                        this_trial = this_trial + 1;
                    case 'Prev'
                        this_trial = this_trial - 1;
                    case 'Select'
                        promptMessage = sprintf(strcat('Select trial',' (',string(1),'-',string(n_trials),')'));
                        this_trial = str2double(inputdlg(promptMessage, 'Trial Selection'));
                end
                close(fig);                               
            end
                    
            fprintf('Finished timing definitions. \n')
            
                     
            %% For trials with no MEPs, replace with mean on set and offsets of all other trials
            
            fprintf('replacing trials with no MEPs with defined windows... \n')
            
            %First determine the mean onsets and offsets that we chose
            for this_set = 1:n_trials
                 if length(info_struct(this_set,2).Position) ~= 0;
                    chosen_onsets_offsets(this_set,:) = [info_struct(this_set,2).Position(1), info_struct(this_set,2).Position(2)];
                 end   
            end
            chosen_onsets_offsets = sort(chosen_onsets_offsets,2);
            chosen_onsets = chosen_onsets_offsets(:,1);
            chosen_offsets = chosen_onsets_offsets(:,2);
            mean_chosen_onset = mean(chosen_onsets(chosen_onsets>0));
            mean_chosen_offset = mean(chosen_offsets(chosen_offsets>0));
            
            %Then for the trials that didnt have MEPs, use the mean of all
            %other onsets and offsets
            for this_set = 1:n_trials
                if length(info_struct(this_set,2).Position) == 0;
                    info_struct(this_set,2).Position(1) = mean_chosen_onset;
                    info_struct(this_set,2).Position(2) = mean_chosen_offset;
                end
            end
                    
            %% Replace background EMG activity cursors with default background
            fprintf('replacing trials where background EMG cursors were not chosen with default EMG window... \n')
            
            for this_set = 1:n_trials
                if length(info_struct(this_set,1).Position) == 0 | length(info_struct(this_set,1).Position) == 1;
                    info_struct(this_set,1).Position(1) = -0.005;
                    info_struct(this_set,1).Position(2) = -0.055;
                end
            end
            
            fprintf('Finished data processing. \n')
            
            %% Store all onsets and offsets
            fprintf('Storing data... \n')
            
            onsets_offsets = [];
            for this_set = 1:n_trials
                onsets_offsets(this_set,:) = [info_struct(this_set,2).Position(2), info_struct(this_set,2).Position(1)];
            end
            onsets_offsets = sort(onsets_offsets,2); %make sure onsets come before offsets
            
           
            %% Store all background_emg periods
            
            background_emg = [];
            for this_set = 1:n_trials
                background_emg(this_set,:) = [info_struct(this_set,1).Position(2), info_struct(this_set,1).Position(1)];
            end
            background_emg = sort(background_emg,2); %make sure onsets come before offsets         
            fprintf('Finshed storing data. \n')
                    
            %% Calculate rms of background_emg
            fprintf('Starting rms calculations... \n')
            rms_background_emg = [];       
            for this_trial = 1:n_trials
                background_emg_window = Data((Data(:,2)>(tms_times(this_trial)+background_emg(this_trial,1))) & (Data(:,2)<(tms_times(this_trial)+background_emg(this_trial,2))),:);
                background_emg_window = abs( background_emg_window);
                rms_background_emg(this_trial,1) = sqrt(mean((background_emg_window(:,1)).^2));
            end
        
            %% Calculate RMS of each MEP window and subtract from background activity
            rms_mep = [];
            for this_mep = 1:n_trials
                mep = Data((Data(:,2)>(tms_times(this_mep)-0.05)) & (Data(:,2)<(tms_times(this_mep)+0.05)),:);
                mep(:,2) = mep(:,2)-tms_times(this_mep);
                mep(:,1) = abs(mep(:,1)); %rectify the signal
                trial_window = mep((mep(:,2)>=onsets_offsets(this_mep,1)) & (mep(:,2)<=onsets_offsets(this_mep,2)),:);
                trial_window_emg_column = trial_window(:,1);
                rms_mep(this_mep,1) = sqrt(mean((trial_window_emg_column).^2));
            end
            
            rms = rms_mep - rms_background_emg;
            rms(rms<0) = 0; %code all negative values as zero
            
            %% Store all data
            Final_data = [];
            Final_data = [rms(:,1), rms_mep(:,1), rms_background_emg(:,1), onsets_offsets(:,1), onsets_offsets(:,2)];
          
            %% Save all data
            fprintf('Saving trial data... \n')
            Stored_Final_data = array2table(Final_data, 'VariableNames', {'Corrected_Amplitude','MEP_Amplitude','Background_Amplitude','MEP_Onset','MEP_Offset'});

            save([wpms.DATAOUT subjlist(px).name filesep subjlist(px).name '_' exp.blocks{block} '_TrialData.mat'], 'Stored_Final_data');% 

            fprintf('Trial Data saved. \n')
        end
            
    end % Sesssion number for loop
    
end % Particpant for loop

%% Below is the vertical cursors functions which allows for selection of MEP, Cortical Silent Period and Burst Windows  

function vertical_cursors(varargin)

if isempty(varargin)
    object = gcf;
else
    object = varargin{1};
    regionNameOrder = varargin{2};
end

mDefaultXLim = xlim;
mDefaultYLim = ylim;

hAx = gca;
pt = get(hAx, 'CurrentPoint');

% Set up cursor text
allLines = findobj(object, 'type', 'line', 'tag','data');
hTextX = nan(1, length(allLines));
hTextY = nan(1, length(allLines));
for id = 1:length(allLines)
    hTextX(id) = text(NaN, NaN, '', ...
        'Parent', get(allLines(id), 'Parent'), ...
        'BackgroundColor', 'yellow', ...
        'Color', get(allLines(id), 'Color'));
    hTextY(id) = text(NaN, NaN, '', ...
        'Parent', get(allLines(id), 'Parent'), ...
        'BackgroundColor', 'yellow', ...
        'Color', get(allLines(id), 'Color'));
end
% Set up cursor lines
allAxes = findobj(object, 'Type', 'axes');
hCur = nan(1, length(allAxes));
for id = 1:length(allAxes)
    hCur(id) = line([NaN NaN], ylim(allAxes(id)), ...
        'Color', 'black', 'Parent', allAxes(id));
end

set(object, ...
    'WindowButtonMotionFcn', @dragFcn,...
    'WindowButtonDownFcn', @clickFcn,...
    'WindowScrollWheel', @WindowScrollWheel, ...
    'WindowKeyPressFcn', @WindowKeyPressCallback);

    function WindowScrollWheel(src,evnt)
        xrange = range(xlim);
%         yrange = range(ylim);
        xrangeNew = xrange * (100 + evnt.VerticalScrollCount * evnt.VerticalScrollAmount)/100;
        set(hAx,'xlim',[pt(1)-xrangeNew/2,pt(1)+xrangeNew/2],'XLimMode','Manual');
    end

    function WindowKeyPressCallback(src, evnt)
        % To do:spacebar to add / modify segments, 
        % Done : enter to go to next trial
        % when enter is pressed check for 2 points per period, otherwise
        % prompt on screen
        % spacebar to switch CSP and MEP
        
        % get existiing startStop points
        startStopLoc = getappdata(object,'startStop');
        dataLine = findobj(object,'type','line','tag','data');
        xdata = get(dataLine,'XData');
        
        pressedPosition=pt;
        
        [~,closestLoc] = min(abs(xdata - pt(1)));
        
        switch evnt.Key
            case 'backspace'
                object.Tag='Prev';
            case 'tab'
                object.Tag='Select';
            case 'return'
                % no NaN = all periods filled out, even NaN means a whole
                % period is left out
                %if sum(isnan(startStopLoc),'all') == 0 || mod(sum(isnan(startStopLoc),'all'),2) == 0
                    object.Tag='Next';
                %else
                %    periodsTxt = {'bkg','mep'};
                %    notSelectedPeriods = find(sum(isnan(startStopLoc),2) ~= 0);
                %    msgbox(['Please select the start/end of  ' strjoin(periodsTxt(notSelectedPeriods),', ') '.']);
                %end
            case {'2','3','4','5'}
                withinPeriod = find((startStopLoc(:,1)-closestLoc <=0) & (startStopLoc(:,2)-closestLoc >=0));
                if ~isempty(withinPeriod)
                    if withinPeriod > 1
                        startStopLocTemp = startStopLoc(str2double(evnt.Key),:);
                        startStopLoc(str2double(evnt.Key),:) = startStopLoc(withinPeriod,:);
                        startStopLoc(withinPeriod,:) = startStopLocTemp;
                        setappdata(object,'startStop',startStopLoc);
                        
                        withinLine = findobj(object,'type','line','tag',regionNameOrder{withinPeriod});
                        targetLine = findobj(object,'type','line','tag',regionNameOrder{str2double(evnt.Key)});
                        
                        ydataWithin = get(withinLine,'YData');
                        ydataTarget = get(targetLine,'YData');
                        
                        set(withinLine,'ydata',ydataTarget);
                        set(targetLine,'ydata',ydataWithin);
                    end
                end                
            case 'd'
                withinPeriod = find((startStopLoc(:,1)-closestLoc <=0) & (startStopLoc(:,2)-closestLoc >=0));
                if ~isempty(withinPeriod)
                    startStopLoc(withinPeriod,:)=nan;
                    setappdata(object,'startStop',startStopLoc);
                    
                    deleteLine = findobj(object,'type','line','tag',regionNameOrder{withinPeriod});
                    ydataDelete = get(deleteLine,'YData');
                    ydataDelete(:) = nan;
                    set(deleteLine,'ydata',ydataDelete);
                end
            case 'o'
                set(hAx,'xlim',mDefaultXLim,'XLimMode','Manual');
            case 'r'
                % rectify signals
                dataLine = findobj(object,'type','line','-regexp','tag','data.*');
                if ~getappdata(object,'rectified')
                    % add an app property where the ydata from the data stream is negative
                    setappdata(object,dataLine.Tag,dataLine.YData);
                    tempLines = findobj(object,'type','line','-regexp','tag','.+');
                    for jj = 1:length(tempLines)
                        set(tempLines(jj),'YData',abs(tempLines(jj).YData));
                    end
                    set(hAx,'ylim',[0,mDefaultYLim(2)],'YLimMode','Manual');
                    setappdata(object,'rectified',true);
                else
                    tempLines = findobj(object,'type','line','-regexp','tag','.+');
                    tempYData = getappdata(object,dataLine.Tag);
                    for jj = 1:length(tempLines)
                        tempNAN = ~isnan(tempLines(jj).YData);
                        tempNAN2=abs(tempNAN);
                        tempNAN2(~tempNAN) = nan;
                        set(tempLines(jj),'YData',...
                            tempNAN2 .* tempYData);
                    end
                    set(hAx,'ylim',mDefaultYLim,'YLimMode','Manual');
                    setappdata(object,dataLine.Tag,[]);
                    setappdata(object,'rsed',true);
                    setappdata(object,'rectified',false);
                end
            case 's'
                % rectify signals
                dataLine = findobj(object,'type','line','-regexp','tag','data.*');
                if ~getappdata(object,'rsed')
                    if ~getappdata(object,'rectified')
                        % add an app property saving the ydata for the data
                        setappdata(object,dataLine.Tag,...
                            dataLine.YData);
                    end
                    tempLines = findobj(object,'type','line','-regexp','tag','.+');
                    for jj = 1:length(tempLines)
                        tempNAN = ~isnan(tempLines(jj).YData);
                        tempNAN2=abs(tempNAN);
                        tempNAN2(~tempNAN) = nan;
                        % 5 local points = 5ms window size for smoothing
                        localPoints = 10;
                        set(tempLines(jj),'YData',...
                            tempNAN2 .* movmean(abs(tempLines(jj).YData),localPoints,'omitnan'));
                    end
                    set(hAx,'ylim',[0,mDefaultYLim(2)],'YLimMode','Manual');
                    setappdata(object,'rectified',true);
                    setappdata(object,'rsed',true);
                else
                    tempLines = findobj(object,'type','line','-regexp','tag','.+');
                    tempYData = getappdata(object,dataLine.Tag);
                    
                    for jj = 1:length(tempLines)
                        % restore original ydata
                        tempNAN = ~isnan(tempLines(jj).YData);
                        tempNAN2=abs(tempNAN);
                        tempNAN2(~tempNAN) = nan;
                        set(tempLines(jj),'YData',...
                            tempNAN2 .* tempYData);
                    end
                    set(hAx,'ylim',mDefaultYLim,'YLimMode','Manual');
                    setappdata(object,dataLine.Tag,[]);
                    setappdata(object,'rectified',false);
                    setappdata(object,'rsed',false);
                end
            case 'space'
                %bound range within plot
                xrange = xlim;
                yrange = ylim;
                
                if xrange(1) < pressedPosition(1,1) && xrange(2) > pressedPosition(1,1) &&...
                        yrange(1) < pressedPosition(1,2) && yrange(2) > pressedPosition(1,2)
                    
                    dataLine = findobj(object,'type','line','tag','data');
                     
                    [~,tempLoc]=min(abs(startStopLoc-closestLoc),[],'all','linear');
                    
                    colTemp = floor((tempLoc-1) / size(startStopLoc,1))+1;
                    rowTemp =  mod(tempLoc-1,size(startStopLoc,1))+1;
                    %startStopLoc = startStopLoc([1,3,2],:);
                    %setappdata(object,'startStop',startStopLoc);
                    tempLine = findobj(object,'type','line','tag',regionNameOrder{rowTemp}); 
                    
                    ydataTemp = get(tempLine,'YData');
                    
                    startStopLoc(rowTemp,colTemp) = closestLoc;
                    ydataTemp(1:end) = nan;
                    if ~isnan(startStopLoc(rowTemp,2))
                    ydataTemp(startStopLoc(rowTemp,1):startStopLoc(rowTemp,2)) = dataLine.YData(startStopLoc(rowTemp,1):startStopLoc(rowTemp,2));
                    else
                        ydataTemp(startStopLoc(rowTemp,1)) = dataLine.YData(startStopLoc(rowTemp,1));
                    end
                    
                    set(tempLine,'ydata',ydataTemp)
                    setappdata(object,'startStop',startStopLoc);
                end
        end
    end
    function clickFcn(src, evnt)
        % AC look at points in pairs? 
        % change properties
        
        xrange = xlim;
        yrange = ylim;
        pressPosition = pt;
        % check if current point is within plotted range
        if pressPosition(1,1) >=  xrange(1) && pressPosition(1,1) <= xrange(2) &&...
                pressPosition(1,2) >=  yrange(1) && pressPosition(1,2) <= yrange(2)
            dataLine = findobj(object,'type','line','tag','data');
            
%             region=findobj(object,'type','line','-regexp','tag','(Plot)+');
%             regionNameTemp = {region.Tag};
            startStopLoc = getappdata(object,'startStop');
            currentPlot = getappdata(object,'currentPlot');
            currentIndex = getappdata(object,'currentIndex');
            
            if isempty(currentPlot)
                if pressPosition(1) <= 0
                    % background if pt is earlier than 0
                    saveIndex = 1;
                else
                    if ~isempty(currentIndex)
                        % find next unfilled region
                        whereUnfilled = find(isnan(startStopLoc(:,1)));
                        if  length(whereUnfilled) > 1
                            if whereUnfilled(1) <= 1
                                whereUnfilled(1) = [];
                            end
                            saveIndex = whereUnfilled(1);
                        else
                            if whereUnfilled(1) <= 1
                                saveIndex = [];
                            else
                                saveIndex = whereUnfilled(1);
                            end
                        end
                    else
                        saveIndex = 2;
                    end
                end
                
                if ~isempty(saveIndex)
                    selectPlot = regionNameOrder{saveIndex};
                else
                    selectPlot = [];
                end
                
                currentPlot = selectPlot;
                currentIndex = saveIndex;
            else
                selectPlot = currentPlot;
                saveIndex = currentIndex;
            end
            
            if ~isempty(currentIndex)
                selectLine = findobj(object,'type','line','tag',selectPlot);
                xdata = get(selectLine,'XData');
                ydata = get(selectLine,'YData');
                [~,closestLoc] = min(abs(xdata - pressPosition(1)));
                
                if sum(~isnan(selectLine.YData)) == 0
                    % flip a single point to the respective colour if there are no
                    % existing marks
                    ydata(closestLoc) = dataLine.YData(closestLoc);
                    startStopLoc(saveIndex,1) = closestLoc;
                    currentPlot = selectPlot;
                elseif sum(~isnan(selectLine.YData)) == 1
                    % flip a section to the respective colour if there is one
                    % existing marked points
                    
                    yRangeTemp = min(find(~isnan(selectLine.YData)),closestLoc):...
                        max(find(~isnan(selectLine.YData)),closestLoc);
                    if length(yRangeTemp) > 1
                        ydata(yRangeTemp) = dataLine.YData(yRangeTemp);
                        startStopLoc(saveIndex,2) = closestLoc;
                        if startStopLoc(saveIndex,1) > startStopLoc(saveIndex,2)
                            tempLoc = startStopLoc(saveIndex,2);
                            startStopLoc(saveIndex,2) = startStopLoc(saveIndex,1);
                            startStopLoc(saveIndex,1) = tempLoc;
                        end
                        if startStopLoc(saveIndex,1) == startStopLoc(saveIndex,2)
                            startStopLoc(saveIndex,2) = [];
                        end
                        if strcmp(currentPlot,'MEPPlot')
                            selectCSP = findobj(object,'type','line','tag','CSPPlot');
                            ydataCSP = get(selectCSP,'YData');
                            % set start of CSP to be the same as end of MEP
                            ydataCSP(closestLoc) = dataLine.YData(closestLoc);
                            startStopLoc(3,1) = startStopLoc(saveIndex,2);
                            currentPlot = 'CSPPlot';
                            currentIndex = 3;
                            set(selectCSP,'YData',ydataCSP);
                        else
                            currentPlot = [];nge
                        end
                    else
                        % overlapping dots selection delete dot
                        ydata(yRangeTemp) = nan;
                        startStopLoc(saveIndex,1) = nan;
                        currentPlot = [];
                    end
                else
                    currentPlot = [];
                end
                
                setappdata(object,'startStop',startStopLoc);
                setappdata(object,'currentPlot',currentPlot);
                setappdata(object,'currentIndex',currentIndex);
                set(selectLine,'YData',ydata);
            end
        end
    end
    function dragFcn(src,evnt)
        % Get mouse location
        pt = get(gca, 'CurrentPoint');
        % Update cursor line position
        % Update cursor text
        for idx = 1:length(allLines)
            xdata = get(allLines(idx), 'XData');
            ydata = get(allLines(idx), 'YData');
            if pt(1) >= xdata(1) && pt(1) <= xdata(end)
                set(hCur, 'XData', [pt(1), pt(1)]);
                y = interp1(xdata, ydata, pt(1));
                yPos = ylim;
                xPos = xlim;
                set(hTextX(idx), 'Position', [pt(1), yPos(1)], ...
                    'String', sprintf('(%0.3f)', pt(1)));
                set(hTextY(idx), 'Position', [xPos(1), y], ...
                    'String', sprintf('(%0.3f)', y));
            else
                set(hTextX(idx), 'Position', [NaN NaN]);
                set(hTextY(idx), 'Position', [NaN NaN]);
            end
        end
    end
end


