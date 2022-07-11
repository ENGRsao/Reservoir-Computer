function [temp,temp_f] = splitStepMethod(signal, soliton_order, fiber_length, dispersion_constant, division_factor, frequencyGrid)
    % Split step fourier transform algorithm 
    % Input Data
    % fiber length : Normalize to the Dispersion length

    no_of_steps = round(fiber_length * division_factor * soliton_order^2);
    deltaZ = fiber_length/no_of_steps;

    dispersion_deltaz = exp(1i * 0.5 * dispersion_constant * frequencyGrid.^2 * deltaZ);

    [temp,~] = NonLinearEffect(signal,soliton_order,deltaZ/2);
    for n=1:no_of_steps
        f_temp = ifft(temp) .* dispersion_deltaz;
        uu = fft(f_temp);
        [temp,~] = NonLinearEffect(uu,soliton_order,deltaZ);
    end
    [temp,~] = NonLinearEffect(temp,soliton_order,deltaZ/2);%Final Field
    temp_f = (ifft(uu)); %.*(fft_point*delta_tau)/sqrt(2*pi);  
end