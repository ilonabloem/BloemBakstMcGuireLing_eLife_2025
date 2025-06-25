
function design = attWindow_taskData(p)

saveName    = sprintf('AllBehavioralData_%s.mat', p.Task);
saveDir     = fullfile(p.dataDir, 'Data', 'behavioralData');

if ~exist(saveDir, 'dir'), mkdir(saveDir); end

if exist(fullfile(saveDir, saveName), 'file') > 0

    fprintf('Loading design information... \n')
    load(fullfile(saveDir, saveName), 'design');

    
else
    
    fprintf('Generating design information... \n')
    
    design  = struct('trueMean', cell(1,p.nSubs), ...
                     'trueWidth', cell(1,p.nSubs), ...  
                     'trialAccuracy', cell(1,p.nSubs));

    for subj = 1:p.nSubs
                 
                 
        widthConditions = []; attCueCenter      = [];
        Accuracy        = []; ReactionTime      = [];
        count           = 0;
        % Load trial information
        for run = 1:p.numRuns(subj)
            count   = count + 1;
            switch p.Task
                case 'CM'
                    list    = dir(fullfile(p.maindir, p.SubjNames{subj}, 'BehavioralData', sprintf('run_%02d.csv', count)));
                case 'PCM'
                    list    = dir(fullfile(p.maindir, p.SubjNames{subj}, 'BehavioralData', sprintf('loc_run_%02d.csv', count)));
            end
            if isempty(list) % one participant has a skipped run
                count   = count + 1;
            end

            switch p.Task
                case 'CM'
                    trialEvents = readtable(fullfile(p.maindir, p.SubjNames{subj}, 'BehavioralData', sprintf('run_%02d.csv', count)));
                    widthConditions = cat(1, widthConditions, [NaN(5,1); trialEvents.width; NaN(5,1)]);
                    attCueCenter    = cat(1, attCueCenter, [NaN(5,1); p.attLocations(trialEvents.winCtrLetterPos+1)'; NaN(5,1)]);
                    Accuracy        = cat(1, Accuracy, [NaN(5,1); trialEvents.Correct; NaN(5,1)]);
                    ReactionTime    = cat(1, ReactionTime, trialEvents.keyTimes-trialEvents.changeTimes);

                case 'PCM'
                    trialEvents = readtable(fullfile(p.maindir, p.SubjNames{subj}, 'BehavioralData', sprintf('loc_run_%02d.csv', count)));
                    widthConditions = cat(1, widthConditions, [NaN(1,1); table2array(trialEvents(1,:))'; NaN(1,1)]);
                    attCueCenter    = cat(1, attCueCenter, [NaN(1,1); p.attLocations(table2array(trialEvents(2,:))+1)'; NaN(1,1)]);
                    Accuracy        = cat(1, Accuracy, [NaN(1,1); NaN(length(table2array(trialEvents(2,:))),1); NaN(1,1)]);
            end

            clear trialEvents;
        end

        tmp_widthConditions = repmat(widthConditions, [1 p.numTRreps])';
        tmp_attCueCenter    = repmat(attCueCenter, [1 p.numTRreps])';
        tmp_accuracy        = repmat(Accuracy, [1 p.numTRreps])';
        widthConditions     = tmp_widthConditions(:);
        attCueCenter        = tmp_attCueCenter(:);
        trialAccuracy       = tmp_accuracy(:);
        trueMean            = [NaN(p.timeLag,1); attCueCenter(1:end-p.timeLag)];
        trueWidth           = [NaN(p.timeLag,1); widthConditions(1:end-p.timeLag)];
        trialAccuracy       = [NaN(p.timeLag,1); trialAccuracy(1:end-p.timeLag)];


        design(subj).trueMean       = trueMean;
        design(subj).trueWidth      = trueWidth;
        design(subj).trialAccuracy  = trialAccuracy;
        
    end
    
    % save file
    save(fullfile(saveDir, saveName), 'design');
end