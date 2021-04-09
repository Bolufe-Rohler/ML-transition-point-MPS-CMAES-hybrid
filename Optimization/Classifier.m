function [result] = Classifier(data, trainedModel)
    result = trainedModel.predictFcn(data);    
    
    % For fixed transitions use this code:
    %result = data(1) >+ 285000;   %95%
    %result = data(1) >+ 255000;   %85%
    %result = data(1) >+ 225000;   %75%
    %result = data(1) >+ 150000;  %50%
    %result = data(1) >+ 75000;   %25%
    %result = data(1) >+ 45000;   %15%
    %result = data(1) >+ 15000;   %5%
end

