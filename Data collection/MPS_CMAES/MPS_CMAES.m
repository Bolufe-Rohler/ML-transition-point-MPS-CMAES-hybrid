function [fval, xmin] = MPS_CMAES( FUN, bound, DIM, fEvalsPerD)
    Max_FEs = DIM*fEvalsPerD;
    k = 70; % percent for eval in MPS
    Max_FEs_MPS = (Max_FEs * k)/100;
    max_fun_evals = Max_FEs - Max_FEs_MPS;
    
    %[x0 fval] = MPS(FUN, DIM, Max_FEs_MPS, ubound)
    
    %parameters for CMAES
    cmaesopts = {};
    cmaesopts.LBounds = -bound;
    cmaesopts.UBounds = bound;
    cmaesopts.PopSize = '(4 + floor(3 * log(N)))';
    cmaesopts.ParentNumber = 'floor(popsize / 2)';
    cmaesopts.DispFinal = 'off';
    cmaesopts.DispModulo = '0';
    cmaesopts.SaveVariables = 'off';
    cmaesopts.LogModulo = '0';
    cmaesopts.LogTime = '0';
    cmaesopts.LogPlot = 'off';
    cmaesopts.stopOnStagnation = 'off';
    cmaesopts.MaxFunEvals = max_fun_evals;
    
%     cmaes('fun', x0, ((ubound - lbound) / 3), cmaesopts);

    [inits_sols, fvalues] = MPS_Mod(FUN, DIM, Max_FEs_MPS, bound);
    
    fvals = inf(1,DIM);
    xmins = zeros(DIM, DIM);
    
    i=1;
    while i<=DIM && cmaesopts.MaxFunEvals > 0   %Using the standard MPS version with poopsize=DIM
        [xmin, fmin, counteval] = cmaes(FUN, inits_sols(i,:), ((2*bound) / 100), cmaesopts);        
        cmaesopts.MaxFunEvals = cmaesopts.MaxFunEvals - counteval;
        fvals(i) = fmin;
        xmins(i,:) = xmin';
        i=i+1;
    end
    
    [fval, idx] = min(fvals);
    xmin = xmins(idx,:);

end
