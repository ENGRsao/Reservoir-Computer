%% clearing data
clc;
%clear;
%close all;

%% code implementation
inputfilename = 'NCEQ_TASKIN_RS2.csv';
outputfilename = 'NCEQ_TASKOUT_RS2.csv';
dataIN = readmatrix(inputfilename);
dataOUT = readmatrix(outputfilename);
X_train = dataIN(1:7000,:);
Y_train = dataOUT(1:7000);
X_train = X_train./max(X_train,[],'all');
%% Ridge regression Algorithm
ridge_coefficient = 0.000005; %0: 1e-5: 5e-3;
scaling_factor = 0;
model = ridge(Y_train,X_train,ridge_coefficient,scaling_factor);


%% test data 

% inputfilename = 'NARMA_TASKIN_RS.csv';
% outputfilename = 'NARMA_TASKOUT_RS.csv';
% dataIN = readmatrix(inputfilename);
% dataOUT = readmatrix(outputfilename);
X_test = dataIN(7001:10000,:);
Y_test = round(dataOUT(7001:10000),2);
X_test = X_test./max(X_test,[],'all');
Y_Symbol = zeros(length(Y_predict),1);
error = zeros(length(Y_predict),1);
if(scaling_factor == 1)
    Y_predict = (X_test * model(1:end));
else
    Y_predict = (model(1) + X_test * model(2:end));    
end
Y_predict = round(Y_predict,2);
for i = 1:length(Y_predict)
    if ((Y_predict(i) >= 0) && (Y_predict(i) < 0.67))
        Y_Symbol(i) = 0.33;
    elseif ((Y_predict(i) < 0) && (Y_predict(i) > -0.67))
        Y_Symbol(i) = -0.33;
    elseif ((Y_predict(i) >= 0.67))
        Y_Symbol(i) = 1;
    elseif ((Y_predict(i) <= -0.67))
        Y_Symbol(i) = -1;
    end

end
for k = 1:length(Y_predict)
    error(k) = isequal(Y_test(k),Y_Symbol(k));
end

SER = (length(Y_predict) - sum(error))/length(Y_predict);
    %Compute cost function
% MSE = sum((Y_test - Y_predict).^2)/length(Y_predict);
% variance_y = var(Y_predict);
% NMSE_ridge = MSE / variance_y;
% figure; plot((1:length(Y_predict)),Y_predict,'*-'); hold on; plot((1:length(Y_predict)),Y_test,'o-');
