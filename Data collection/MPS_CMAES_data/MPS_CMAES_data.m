function [fval, xmin] = MPS_CMAES_data(FUN, bound, DIM, fEvalsPerD, Fnum, run)
   
    maxFEs = DIM*fEvalsPerD;


    hyb_FEs = maxFEs/20;
    hyb_counter = 1;
        
    
    % MPS Parameters
    alpha = 0.3;
    gamma = 3;
    popsize = DIM;   
    d = sqrt(DIM)*2*bound;  % Search Space Diagonal
    
    xmins = zeros(popsize, DIM);
    dataRun = zeros(20,125); %Data collected, 19 hybridization points + the final result; hybrid result + 124 attributes
    best30Fval = zeros(30,1); 
    median30Fval = zeros(30,1); 
    worse30Fval = zeros(30,1); 
    newSolsCount30= zeros(30,1); 
    
    Population = zeros(2*popsize, DIM);
    f_Pop = 1e+50*ones(2*popsize,1);
    
    % Initial Population
    Population(1:popsize, :) =  (bound/2)*(((unifrnd(0,1,popsize,DIM)<0.5)*2)-1); 
    f_Pop(1:popsize,:) = feval(FUN, Population(1:popsize,:)');
    FEs = popsize;
        
    while  (FEs < maxFEs) 
        
        [~, indexes] = sort(f_Pop);
        
        % Updating threshold   
        min_step =  max(alpha*d* ((maxFEs-FEs)/maxFEs)^gamma, 1e-05);
        max_step = 2*min_step;         
        
        % Population Centroid
        centroid = repmat(sum(Population(indexes(1:popsize),:))/popsize, popsize, 1);
        
        % Difference Vectors
        dif = normr(centroid - Population(indexes(1:popsize),:));
        
        % Difference Vector Scaling Factor
        F = unifrnd(-max_step, max_step, popsize,1);
        
        % Orthogonal Vectors
        orth = normr(normrnd(0,1,popsize,DIM));
        orth = normr(orth - repmat(dot(orth',dif')',1,DIM).*dif);
        
        % Orthogonal Step Scaling Factor
        min_perp = sqrt(max(min_step^2-abs(F).^2,0));
        max_perp = sqrt(max(max_step^2-abs(F).^2,0));         
        FO = unifrnd(min_perp, max_perp)'; 
         
        % New Solutions & Clamping & Evaluation
        Population(indexes(popsize+1:2*popsize),:) =  max( min(Population(indexes(1:popsize),:) + ...           % Current Population
                                                                repmat(F,1,DIM).*dif + ...                  % Difference Vector Step
                                                                repmat(FO',1,DIM).*orth, bound), -bound);   % Orthogonal Step
        f_Pop(indexes(popsize+1:2*popsize)) = feval(FUN, Population(indexes(popsize+1:2*popsize),:)');   
        FEs = FEs + popsize;
        
        %%% Moving the values to keep vector updated  -- THE LAST ONES ARE
        %%% THE MOST RECENT ONES
        for k = 1: 29
            best30Fval(k) = best30Fval(k+1);
            median30Fval(k) = median30Fval(k+1); 
            worse30Fval(k) = worse30Fval(k+1); 
            newSolsCount30(k) = newSolsCount30(k+1); 
        end
        
        best30Fval(30) = min(f_Pop);
        median30Fval(30) = median(f_Pop); 
        worse30Fval(30) = max(f_Pop); 
        newSolsCount30(30) = sum(median(f_Pop(indexes(1:popsize)))>f_Pop(indexes(popsize+1:2*popsize)));
        
        %Hybridization point!
        if FEs > hyb_counter*hyb_FEs    %Moments for hybrid
       
            [sorted, indexes] = sort(f_Pop);
            for i = 1:popsize
                xmins(i,:) = Population(indexes(i),:);
            end
            
            fval = CMAE_ES_Hyb(FUN, xmins, bound, maxFEs - FEs, DIM);
            dataRun(hyb_counter, 1) = fval; %first value is the hybrid result
            dataRun(hyb_counter, 2) = FEs; 
            dataRun(hyb_counter, 3) = maxFEs - FEs; 
            dataRun(hyb_counter, 4) = min_step; 
            dataRun(hyb_counter, 5) = max_step; 
            dataRun(hyb_counter, 6:35) = best30Fval; 
            dataRun(hyb_counter, 36:65) = median30Fval; 
            dataRun(hyb_counter, 66:95) = worse30Fval; 
            dataRun(hyb_counter, 96:125) = newSolsCount30; 
            
            hyb_counter =  hyb_counter +1;
        end
        
    end
    
    % Final Result
    [sorted, indexes] = sort(f_Pop);
    xmin = Population(indexes(1),:);
    fval = sorted(1);
    dataRun(hyb_counter, 1) = fval; %The final result of MPS
    save(sprintf('data_run_F%i_%i',Fnum, run), 'dataRun');
end

function fval = CMAE_ES_Hyb(FUN, inits_sols, bound, max_fun_evals, DIM)
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


    fvals = inf(1,DIM);

    i=1;
    while i<=DIM && cmaesopts.MaxFunEvals > 0   %Using the standard MPS version with poopsize=DIM
        [xmin, fmin, counteval] = cmaes(FUN, inits_sols(i,:), ((2*bound) / 50), cmaesopts);        
        cmaesopts.MaxFunEvals = cmaesopts.MaxFunEvals - counteval;
        fvals(i) = fmin;
        i=i+1;
    end

    [fval, ~] = min(fvals);
end

function n = normr(m)
    [~,mc]=size(m);
    if (mc == 1)
      n = m ./ abs(m);
    else
        n=sqrt(ones./(sum((m.*m)')))'*ones(1,mc).*m;
    end
end