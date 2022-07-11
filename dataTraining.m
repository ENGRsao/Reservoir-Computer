%% clearing data
clc;
%clear;
%close all;

%% code implementation
inputfilename = 'NARMA_TASKIN_RS.csv';
outputfilename = 'NARMA_TASKOUT_RS.csv';
dataIN = readmatrix(inputfilename);
dataOUT = readmatrix(outputfilename);
X_train = dataIN(1:2500,:);
Y_train = dataOUT(1:2500);
X_train = X_train./max(X_train,[],'all');
%% Ridge regression Algorithm
ridge_coefficient = 0.000005; %0: 1e-5: 5e-3;
scaling_factor = 0;
model = ridge(Y_train,X_train,ridge_coefficient,scaling_factor);


%% test data 

inputfilename = 'NARMA_TASKIN_RS2.csv';
outputfilename = 'NARMA_TASKOUT_RS2.csv';
dataIN = readmatrix(inputfilename);
dataOUT = readmatrix(outputfilename);
X_test = dataIN(1:1000,:);
Y_test = dataOUT(1:1000);
X_test = X_test./max(X_test,[],'all');

if(scaling_factor == 1)
    Y_predict = (X_test * model(1:end));
else
    Y_predict = (model(1) + X_test * model(2:end));    
end
    %Compute cost function
MSE = sum((Y_test - Y_predict).^2)/length(Y_predict);
variance_y = var(Y_predict);
NMSE_ridge = MSE / variance_y;
%figure; plot((1:length(Y_predict)),Y_predict,'*'); hold on; plot((1:length(Y_predict)),Y_test,'o');

%% Lasso Regression Algorithm
% alpha = 0.1;
% fold_validation = 20;
% [model,FitInfo] = lasso(X_train,Y_train,'Alpha',alpha,'CV',fold_validation);
% idxLambda1SE = FitInfo.Index1SE;
% Weight = model(:,idxLambda1SE);
% constant_coef = FitInfo.Intercept(idxLambda1SE);
% Y_predict = X_test * Weight + constant_coef;
% Compute cost function
% MSE = sum((Y_test - Y_predict).^2)/length(Y_predict);
% variance_y = sum((Y_predict - mean(Y_predict)).^2)/length(Y_predict);
% NMSE_lasso = MSE / variance_y;
% figure;plot((1:length(Y_predict)),Y_predict,'*'); hold on; plot((1:length(Y_predict)),Y_test,'o');

