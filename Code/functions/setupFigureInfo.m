% All parameters needed for visualization

function p = setupFigureInfo(p)

    p.histEdgesAng      = linspace(-pi,pi,21); % 9 degrees bins 
    p.histCentersAng    = (p.histEdgesAng(1:end-1)+p.histEdgesAng(2:end))/2;
    
    p.histEdgesWidth    = linspace(0, 2*pi,21); % Spaces in 4.5 degrees bins - matches distance between stimuli 
    p.histCentersWidth  = (p.histEdgesWidth(1:end-1)+p.histEdgesWidth(2:end))/2;
    
    p.histEdgesResp     = linspace(-1,4,21); % 
    p.histCentersResp   = (p.histEdgesResp(1:end-1)+p.histEdgesResp(2:end))/2;
    
    p.histEdgesAmpl     = linspace(0,10,21); % 
    p.histCentersAmpl   = (p.histEdgesAmpl(1:end-1)+p.histEdgesAmpl(2:end))/2;
    
    p.histEdgesBase     = linspace(-5,5,21);
    p.histCentersBase   = (p.histEdgesBase(1:end-1)+p.histEdgesBase(2:end))/2;
    
    p.histEdgesBeta     = linspace(1.8,50,21);
    p.histCentersBeta   = (p.histEdgesBeta(1:end-1)+p.histEdgesBeta(2:end))/2;
    
    p.histEdgesStd      = linspace(deg2rad(6),pi,21);
    p.histCentersStd    = (p.histEdgesStd(1:end-1)+p.histEdgesStd(2:end))/2;
    
    p.R2_cutoff         = 0.8; % proportion of R2 to use
    
    
    customColor         = [ones(1,256/2) linspace(1,0,256/2); ...
                           linspace(1,0,256); ...
                           0.5 * ones(1,256)]';
    switch p.Task
        case 'PCM'
            p.CondColors    = customColor([65 90 120 150 210],:);
            exmpleSubj      = '021';
        case 'CM'
            Colors          = customColor([65 90 120 150 210],:);
            p.CondColors    = Colors([1:3 5],:);
            exmpleSubj      = '019';
    end
    
    p.exampleSubj       = find(strcmp(p.SubjNames, exmpleSubj)); %{'022','012','021','003','025'}
end