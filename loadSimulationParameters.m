function [Pulse,Fiber] = loadSimulationParameters()
    Pulse = struct();
    Fiber = struct();
    
    %% Pulse Parameters
    Pulse.fft_point = 4096;   %FFT points   N = 2^10
    Pulse.time_span = 32;     %Time period of the desire signal
    Pulse.delta_tau = (2*Pulse.time_span)/Pulse.fft_point;  % smallest division of time
    Pulse.tau = (-Pulse.fft_point/2:Pulse.fft_point/2-1)*Pulse.delta_tau; %time duration spread out of the signal
    Pulse.omega = [(0:Pulse.fft_point/2-1) (-Pulse.fft_point/2:-1)] .*(pi/Pulse.time_span);
    Pulse.normalize_omega = fftshift(Pulse.omega)/(2*pi); %frequency
    Pulse.T_FWHM = 10; %ps
    Pulse.T_nat = Pulse.T_FWHM/1.665; %Input signal Pulse width in ps; for sech => T_nat = T_FWHM/1.763;
    Pulse.Normalize_tau = Pulse.tau/Pulse.T_nat;
    Pulse.guassian_order = 3;
    Pulse.chirp = 0;
    Pulse.speed_of_light = 3 * 10^5; %in nm/ps
    Pulse.wavelength = 1500; % in nm
    Pulse.Peak_power = 0.90; %in Watt
    
    %% Fiber Parameters
        %Dispersion parameters
    Fiber.D = 17; %GVD parameter in ps/nm km
    Fiber.Beta2 = 1 * (Pulse.wavelength^2 * Fiber.D)/(2*pi*Pulse.speed_of_light); % ps^2/Km ; Dispersion Constant
    Fiber.dispersion_length = Pulse.T_nat^2 / abs(Fiber.Beta2); %Km
    Fiber.length = 1.2 * Fiber.dispersion_length;  %Km
    
        %NonLinear Parameters 
    Fiber.gamma = 1.3; %W^-1.km^-1 for standard single mode fiber; Nolinear Coefficeint, 1.76
    Fiber.alpha = 0.2;  % fiber loss db/km
    Fiber.L_effective = (1 - exp(-1 * Fiber.alpha * Fiber.length))/Fiber.alpha; %Km effective length
    Fiber.NonLinear_length = 1 / (Fiber.gamma * Pulse.Peak_power);   %Km
    
        %Numerical Simulation parameter
    Pulse.soliton_order  = round(sqrt(Fiber.dispersion_length/Fiber.NonLinear_length));  %constant
    Fiber.Division_factor = 20;
        
        %Input signal type
    Pulse.sech = sech(Pulse.Normalize_tau) .* exp(-1i*0.5*Pulse.chirp.*Pulse.Normalize_tau.^2); %sech signal
    Pulse.guassian = exp(-0.5 * (1+1i*Pulse.chirp) .* Pulse.Normalize_tau.^2);  %define guassian pulse
    Pulse.super_guassian = exp(-0.5 * (1+1i*Pulse.chirp) .* Pulse.Normalize_tau.^(2*Pulse.guassian_order));  %super Guassian
end