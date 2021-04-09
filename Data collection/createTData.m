data = -1*ones(100*19,125);
row_counter = 1;

for run=1:100
    for Fnum =1:28
        load(sprintf('data_run_F%i_%i',Fnum, run));
        for rows =1:19
            data(row_counter, 1:124) = dataRun(rows,2:125); %here to choose the attributes 
            if dataRun(rows, 1) < dataRun(20, 1)
                data(row_counter, 125) = 1;
            else
                data(row_counter, 125) = 0;
            end
            
            row_counter = row_counter + 1;
        end
    end
end

save('ClassifierData', 'data');