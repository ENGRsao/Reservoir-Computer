%% clearing data
clc;
clear;
close all;

%% code implementation
filename = 'NARMA10_task_GVD_NL_midpoint2.csv';
data = readmatrix(filename);
X_train = data(1:350,1:4096);
Y_train = data(1:350,4098);
X_test = data(351:500,1:4096);
Y_test = data(351:500,4098);

%% Ridge regression Algorithm
ridge_coefficient = 0.005; %0: 1e-5: 5e-3;
scaling_factor = 0;
model = ridge(Y_train,X_train,ridge_coefficient,scaling_factor);

if(scaling_factor == 1)
    Y_predict = X_test * model(1:end);
else
    Y_predict = model(1) + X_test * model(2:end);    
end
    %Compute cost function
MSE = sum((Y_test - Y_predict).^2)/length(Y_predict);
variance_y = sum((Y_predict - mean(Y_predict)).^2)/length(Y_predict);
NMSE_ridge = MSE / variance_y;
figure; plot((1:length(Y_predict)),Y_predict,'*'); hold on; plot((1:length(Y_predict)),Y_test,'o');

%% Lasso Regression Algorithm
alpha = 0.75;
fold_validation = 20;
[model,FitInfo] = lasso(X_train,Y_train,'Alpha',alpha,'CV',fold_validation);
idxLambda1SE = FitInfo.Index1SE;
Weight = model(:,idxLambda1SE);
constant_coef = FitInfo.Intercept(idxLambda1SE);
Y_predict = X_test * Weight + constant_coef;
%Compute cost function
MSE = sum((Y_test - Y_predict).^2)/length(Y_predict);
variance_y = sum((Y_predict - mean(Y_predict)).^2)/length(Y_predict);
NMSE_lasso = MSE / variance_y;
figure;plot((1:length(Y_predict)),Y_predict,'*'); hold on; plot((1:length(Y_predict)),Y_test,'o');

