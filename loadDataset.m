function [input,output] = loadDataset()
    %% Load Dataset base on the required Task
    % val - Input valu determining the type of dataset
    %   val = 1  - NARMA task
    input_filename = 'NARMA10_input.csv';    
    input = readmatrix(input_filename);
    output_filename = 'NARMA10_output.csv';
    output = readmatrix(output_filename);
end