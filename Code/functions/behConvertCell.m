function [allData] = behConvertCell(allData,nFiles,nSubs)
%% Convert cell data types into usable numbers and strings

for sub=1:nSubs
    for run=1:nFiles(sub)
        for trial=1:allData{sub}.eye{run}.nTrials
            buff    = cell2mat(allData{sub}.data{run}.windowCenter(trial));
            buff    = strsplit(buff,',');
            buff1   = cell2mat(buff(1));
            buff2   = cell2mat(buff(2));
            allData{sub}.data{run}.winCtrX(trial) = str2double(buff1(2:end));
            allData{sub}.data{run}.winCtrY(trial) = str2double(buff2(1:end-1));
            tempKey = cell2mat(allData{sub}.data{run}.keyEvenOdd(trial));
            if isempty(tempKey)
                allData{sub}.data{run}.keyEO(trial) = NaN;
            else
                if tempKey(3)=='r'
                    allData{sub}.data{run}.keyEO(trial) = 0;
                elseif tempKey(3)=='l'
                    allData{sub}.data{run}.keyEO(trial) = 1;
                else
                    allData{sub}.data{run}.keyEO(trial) = str2double(tempKey(3));
                end
            end
        end
    end
end

