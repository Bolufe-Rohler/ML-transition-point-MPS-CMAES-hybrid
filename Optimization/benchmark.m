function benchmark(fNums, D, fEvalsPerD, runs, algid, algfunc, algargs)

    t0 = clock;

    fprintf('%s started at: %s\n\n', algid, datestr(t0, 'yyyy-mm-dd HH:MM:SS'));

    fDeltas = [-1400, -1300, -1200, -1100, -1000, -900, -800, -700, ...
               -600, -500, -400, -300, -200, -100, 100, 200, 300, ...
               400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400];

    for fNum = fNums
        csvfile = sprintf('%s_f%d_D%d.csv', algid, fNum, D);

        if ~(exist(sprintf('%s/%s', pwd(), csvfile), 'file') == 2)

            fErrors = zeros(12, runs);
            %fErrors = zeros(11, runs);

            f = @(x) cec13_func(x, fNum);

            for run = 1:runs
                rand('state', sum(100 * clock));
                %gbestCost = feval(algfunc, f, 100, D, fEvalsPerD, algargs{:});
                [gbestCost, xmin, transition] = feval(algfunc, f, 100, D, fEvalsPerD, fNum, run);
                %gbestCost = feval(algfunc, f, 100, D, fEvalsPerD, fNum, run);
                
                % Currently only reporting the results after D*fEvalsPerD FEs.
                %fValues = [NaN(1,10), gbestCost];
                fValues = [NaN(1,10), gbestCost, transition];

                fErrors(1:11,run) = fValues(1,1:11) - fDeltas(fNum);
                fErrors(12,run) = fValues(1,12);
                for fErrorIndex = 1:11
                    if fErrors(fErrorIndex,run) <= 1e-8
                        fErrors(fErrorIndex,run) = 0;
                    end
                end

                fprintf('f%0.2d in %d-D, error: %.4e, cumul time[h]: %.2f\n', ...
                        fNum, D, fErrors(11,run), etime(clock, t0) / 60 / 60);
            end

            fprintf('\n');

            % The format of CSV corresponds to the description of the 
            % plain-text output files in the CEC-2013 technical report.
            csvwrite(csvfile, fErrors);
        end
    end

    fprintf('%s ended at: %s, total time[h]: %.2f\n\n', ...
            algid, datestr(clock, 'yyyy-mm-dd HH:MM:SS'), ...
            etime(clock, t0) / 60 / 60);
end
