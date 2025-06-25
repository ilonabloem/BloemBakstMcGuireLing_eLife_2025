
%% simulation

model           = fitGenGauss([], {'initialize'});
startPrm        = model.init;
startPrm(1)     = deg2rad(60);
startPrm(5)     = -1;

wrapxdata       = linspace(-pi, 3*pi, 120)'; % 2x the space to account for wrap around
xdata           = wrapxdata(wrapxdata >= 0 & wrapxdata < 2*pi);

trueMu          = startPrm(1);
trueAmpl        = 5;
trueBase        = -4;
trueSD          = startPrm(2);
trueBeta        = 10.^startPrm(3);

% simulate attentional field
attField        = exp(-(abs(wrapxdata-trueMu)/trueSD).^trueBeta);
    
% Wrap
attSpace        = sum(cat(2, attField(wrapxdata >= 0 & wrapxdata < 2*pi), ...
                    [attField(wrapxdata >= 2*pi); attField(wrapxdata < 0)]),2);          
attSpace        = normalize(attSpace, 'range', [0 1]);                     
            
ydata           = (attSpace * trueAmpl) + trueBase;

% add some noise
ydata           = ydata + ((trueAmpl/10) * randn(size(ydata)));

% find best params
fixedParams     = cell(3,1);
fixedParams{1}  = 'optimize';
fixedParams{2}  = wrapxdata;
fixedParams{3}  = ydata(:);
[est_params, SSE, exitflag] = ...
    fmincon(@(x) fitGenGauss(x, fixedParams), startPrm, [], [], [], [], model.lb, model.ub, [],  model.opt);
 

fixedParams{1}  = 'prediction';
genGausResults = fitGenGauss(est_params, fixedParams);


% visualize
figure, 
plot(xdata, ydata)
hold on, 
plot(xdata, genGausResults.y_est)

text(2*pi-0.5, 1, {sprintf('Amplitude: %d ', trueAmpl), ...
    sprintf('Baseline: %d', trueBase), ...
    sprintf('Max response: %d',round(max(ydata)))})
sgtitle('Simulation amplitude and baseline parameters attentional field')
box off




