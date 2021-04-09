function [fmin, xmin] = CMA_ES(FUN, bound, DIM, fEvalsPerD, fNum, run)
%function [fmin, xmin, transition] = CMA_ES(FUN, bound, DIM, fEvalsPerD, fNum, run)
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
    cmaesopts.stopOnStagnation = 'on';
    cmaesopts.MaxFunEvals = DIM*fEvalsPerD;

    
    
    [xmin, fmin] = cmaes(FUN, unifrnd(-bound, bound,DIM,1), (bound / 3), cmaesopts);   
    %transition = -1;
       
end