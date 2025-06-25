function out = fitGenGauss(est_params, fixed_params)

mode    = fixed_params{1};

switch mode
    
    case 'initialize'       
        
        %% setup generalized gaussian model 
        % init = [mu, SD, Beta, ampl, offset]

        init    = [pi       pi/8        log10(4)        1       0];
        lb      = [-pi      deg2rad(6)  log10(1.8)      0       -Inf];
        ub      = [3*pi     pi          log10(50)       20      +Inf];
        opts    = optimset( 'display', 'off');

        out.init    = init;
        out.lb      = lb;
        out.ub      = ub;
        out.opt     = opts;
        out.labels  = {'mu', 'SD', 'Beta', 'ampl', 'offset'};
        
    case {'optimize', 'prediction'}
        
        % pull out free parameters
        mu_est      = wrapTo2Pi(est_params(1));
        sd_est      = est_params(2);
        beta_est    = 10.^est_params(3);
        ampl_est    = est_params(4);
        base_est    = est_params(5);

        % pull out fixed parameters
        sp_data     = fixed_params{2}(:);
        data        = fixed_params{3}(:);
        
        %% Simulate attentional window over space  

        R_est       = exp(-(abs(sp_data-mu_est)/sd_est).^beta_est);
    
        % Wrap
        R_space     = sum(cat(2, R_est(sp_data >= 0 & sp_data < 2*pi), ...
                        [R_est(sp_data >= 2*pi); R_est(sp_data < 0)]),2);

        R_space     = normalize(R_space, 'range', [0 1]);            

        widthModel  = (R_space * ampl_est) + base_est;                 

        if strcmp(mode, 'optimize')
            
            % Get goodness of fit for this set of parameters
            % value to minimize.  
            SSE         = nansum((widthModel - data).^2); 
            out         = SSE;

        elseif strcmp(mode, 'prediction')
            
            if ~isempty(data)
                out.SSE     = nansum((widthModel - data).^2); 
                out.R2      = 1 - (nansum((widthModel - data).^2) / nansum((data - nanmean(data)).^2)); 
            end
            out.y_est   = widthModel;

        end
end