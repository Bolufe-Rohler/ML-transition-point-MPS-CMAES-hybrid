D = 30;
fNums = 1:28;
fEvalsPerD = 10000;
runs = 51;

benchmark(fNums, D, fEvalsPerD, runs, 'MLHybrid', 'MLHybrid');