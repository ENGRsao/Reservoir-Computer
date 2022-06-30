clc;
close all;
clear;

%% Simulation parameters
fft_point = 4096;   %FFT points
time_span = 64;     %Time period of the desire signal
delta_tau = (2*time_span)/fft_point;  % smallest division of time
tau = (-fft_point/2:fft_point/2-1)*delta_tau; %time duration spread out of the signal
omega = [(0:fft_point/2-1) (-fft_point/2:-1)] .*(pi/time_span);
normalize_omega = fftshift(omega)/(2*pi); %frequency
T_FWHM = 10; %ps
T_nat = T_FWHM/1.655; %Input signal Pulse width in ps; for sech => T_nat = T_FWHM/1.763;
Normalize_tau = tau/T_nat;
guassian_order = 3;
chirp = 0;
speed_of_light = 3 * 10^5; %in nm/ps
pulse_wavelength = 1550; % in nm

% Dispersion parameters
D = 17;%GVD parameter in ps/nm km
dispersion_constant = -1 * (pulse_wavelength^2 * D)/(2*pi*speed_of_light); % ps^2/Km
dispersion_length = T_nat^2 / abs(dispersion_constant); %Km
fiber_length = 1.5 * dispersion_length;  %Km

% NonLinear Parameters 
gamma = 1.76; %W^-1.km^-1 for standard single mode fiber
alpha = 0.2; %fiber loss db/km
L_effective = (1 - exp(-1 * alpha * fiber_length))/alpha; %Km
Peak_power = 0.2; %in Watt
NonLinear_length = 1 / (gamma * Peak_power);   %Km
soliton_order  = round(sqrt(dispersion_length/NonLinear_length));  %constant
fiber_division_factor = 4000;
%% Input signal
%u_t = sech(Normalize_tau) .* exp(-1i*0.5*chirp.*Normalize_tau.^2); %sech signal
u_t_guassian = exp(-0.5 * (1+1i*chirp) .* Normalize_tau.^2);  %define guassian pulse

%% Reservoir Concept Verification
u_n = ones(1,50);
u_t = u_n(1) .* u_t_guassian;
readout_t = zeros(length(u_n),fft_point);
readout_f = zeros(length(u_n),fft_point);
No_of_recurrent = length(u_n);

figure;     plot(Normalize_tau,u_t); hold on;
for i = 2:No_of_recurrent     
     %xlabel ('Normalize time'); ylabel('Intensity'); title('Reservoir Response');
     %[NL_signal,NL_signal_f] = NonLinearEffectStudy(u_t,NonLinear_length,fiber_length);
     [NL_signal,NL_signal_f] = splitStepMethod(u_t,soliton_order,fiber_length, dispersion_constant, fiber_division_factor, omega );
     %NL_signal_f =  NL_signal_f .*(sqrt(2*pi)/(fft_point*delta_tau));
     readout_t(i-1,:) = abs(NL_signal).^2; 
     %readout_f(i-1,:) = abs(NL_signal_f).^2; 
     figure;
     plot(Normalize_tau,abs(NL_signal).^2);  hold on; 
     u_t = (u_n(i) .* u_t_guassian) + NL_signal;
     plot(Normalize_tau,abs(u_t).^2);     
end
[NL_signal,NL_signal_f] = splitStepMethod(u_t,soliton_order,fiber_length, dispersion_constant, fiber_division_factor, omega );
%NL_signal_f =  NL_signal_f .*(sqrt(2*pi)/(fft_point*delta_tau));
readout_t(No_of_recurrent,:) = abs(NL_signal).^2; 
%readout_f(No_of_recurrent,:) = abs(NL_signal_f).^2; 
%     figure;
%     for i = 1:No_of_recurrent
%         subplot(No_of_recurrent,1,i);
%         plot(Normalize_tau,readout_t(i,:));
%         xlabel ('Normalize time'); ylabel('Intensity'); title('Reservoir Response');
%      end    

%     %% Masking Input data
% input_mask = [0.689435100924171	0.0504550031686729	0.184434115699354	0.0456583234773829	0.885041502041251	0.839794385798114	0.118155245153682	0.410414853697837	0.120228589990650];
% input_mask = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
% 
% mask_len = length(input_mask);
% signal_len = length(u_n);
% 
% reservoirmemory = zeros(mask_len,fft_point);
% readout_t = zeros(signal_len,fft_point);
% readout_f = zeros(signal_len,fft_point);
% %% run test
% for i = 1:signal_len
%     time_multiplexed_signal = input_mask .* u_n(i);    
%     for k = 1:length(time_multiplexed_signal)
%         modulated_signal = time_multiplexed_signal(k) .* u_t_guassian;
%         modulated_signal = modulated_signal + reservoirmemory(k,:);
%         [NL_signal,NL_signal_f] = splitStepMethod(modulated_signal,soliton_order,fiber_length, dispersion_constant, fiber_division_factor, omega ); %Non Linear Node effect on input data generating reservoir state i
%         NL_signal_f =  NL_signal_f .* (sqrt(2*pi)/(fft_point*delta_tau));        
%         reservoirmemory(k,:) = NL_signal;    
%     end
%     readout_t(i,:) = abs(NL_signal).^2;
%     figure; plot(Normalize_tau,readout_t(i,:));
% end  