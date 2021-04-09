function [fmin, xmin] = MPS(FUN, bound, DIM, fEvalsPerD, fNum, run)
%function [fmin, xmin, transition] = MPS(FUN, bound, DIM, fEvalsPerD, fNum, run)
   
    maxFEs = DIM*fEvalsPerD;
     
    
    % MPS Parameters
    alpha = 0.3;
    gamma = 3;
    popsize = DIM;   
    d = sqrt(DIM)*2*bound;  % Search Space Diagonal
    
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
            
    end
    
    % Final Result if no hybridization 
    [sorted, indexes] = sort(f_Pop);
    xmin = Population(indexes(1),:);
    fmin = sorted(1);
    %transition = -1;
end

function n = normr(m)
    [~,mc]=size(m);
    if (mc == 1)
      n = m ./ abs(m);
    else
        n=sqrt(ones./(sum((m.*m)')))'*ones(1,mc).*m;
    end
end