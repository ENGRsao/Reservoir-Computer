function [input,output] = loadEQ_Dataset()
    %% Load Dataset base on the required Task
    % DATASET for Nonlinear Channel Equalizer Dataset
    input_filename = 'EqualizerInput.csv';    
    input = readmatrix(input_filename);
    output_filename = 'EqualizerExpected.csv';
    output = readmatrix(output_filename);
end