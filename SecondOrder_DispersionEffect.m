function [disperse_signal_t, disperse_signal] = SecondOrder_DispersionEffect(signal,dispersion_constant, fiberlength, frequencyGrid)
    signal_spectrum = fft(abs(signal));
    Second_order_effect = 1i * 0.5 * dispersion_constant * frequencyGrid.^2;
    disperse_signal = signal_spectrum.* exp ((Second_order_effect) * fiberlength); %dispersion effect 2nd Order
    disperse_signal_t = ifft((disperse_signal));% convert back to time domain
end