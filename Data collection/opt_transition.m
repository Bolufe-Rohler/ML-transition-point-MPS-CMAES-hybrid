clear

fDeltas = [-1400, -1300, -1200, -1100, -1000, -900, -800, -700, ...
           -600, -500, -400, -300, -200, -100, 100, 200, 300, ...
           400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400];
result = zeros(1, 28);
for fun=1:28
    for run=1:100
        load(sprintf('data_run_F%i_%i.mat',fun, run)); 
        result(fun) = result(fun) + min(dataRun(:,1)); 
    end
    result(fun) = result(fun)/100 - fDeltas(fun);
    if result(fun) < 1e-8
        result(fun) = 0;
    end
end

save opt_transition result