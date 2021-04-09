
D = 30;
fNums = 27:27;
fEvalsPerD = 10000;
runs = 1;

benchmark(fNums, D, fEvalsPerD, runs, 'MPS_CMAES_data', 'MPS_CMAES_data');
%benchmark(fNums, D, fEvalsPerD, runs, 'MPS_CMAES', 'MPS_CMAES', {});

