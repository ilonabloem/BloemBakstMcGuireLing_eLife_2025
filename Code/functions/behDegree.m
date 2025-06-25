function [allData] = behDegree(allData,nFiles,pxdeg)
%% Transform data from pixels into degrees

nSubs = length(allData);
for sub=1:nSubs
    for run=1:nFiles(sub)
        allData{sub}.data{run}.winCtrX  = allData{sub}.data{run}.winCtrX/pxdeg;
        allData{sub}.data{run}.winCtrY  = allData{sub}.data{run}.winCtrY/pxdeg;
        allData{sub}.eye{run}.xpos      = allData{sub}.eye{run}.xpos/pxdeg;
        allData{sub}.eye{run}.xposT     = allData{sub}.eye{run}.xposT/pxdeg;
        allData{sub}.eye{run}.ypos      = allData{sub}.eye{run}.ypos/pxdeg;
        allData{sub}.eye{run}.yposT     = allData{sub}.eye{run}.yposT/pxdeg;
        allData{sub}.fix{run}.xavg      = allData{sub}.fix{run}.xavg/pxdeg;
        allData{sub}.fix{run}.yavg      = allData{sub}.fix{run}.yavg/pxdeg;
        allData{sub}.sacc{run}.xstart   = allData{sub}.sacc{run}.xstart/pxdeg;
        allData{sub}.sacc{run}.xend     = allData{sub}.sacc{run}.xend/pxdeg;
        allData{sub}.sacc{run}.ystart   = allData{sub}.sacc{run}.ystart/pxdeg;
        allData{sub}.sacc{run}.yend     = allData{sub}.sacc{run}.yend/pxdeg;
    end
end

