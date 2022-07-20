clc;
close all;
clear;

[Input,Output] = loadEQ_Dataset();
[Pulse,Fiber] = loadSimulationParameters();

%% test split step function
[SSFM_t,SSFM_f] = splitStepMethod(Pulse.guassian,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega );
plot(Pulse.Normalize_tau,abs(Pulse.guassian).^2);hold on;
plot(Pulse.Normalize_tau,abs(SSFM_t).^2);

%% Reservoir Network application to NARMA10 task 
%input_mask = [0.7 0.2 0.4 0.5 0.9 0.8 0.1 0.4 0.3 0.6 0.95 0.3 0.99	0.35 0.2 0.7 0.8 0.6 0.1 0.4 0.2 0.2 0.1 0.9 0.8 0.4 0.8 0.75 0.6	0.4	0.35 0.15 0.3 0.55 0.45	0.65 0.5	0.8	0.1	0.4	0.6	0.2	0.6	0.1	0.8	0.7	0.3	0.4	1.0	0.3];
input_mask = [0.2 0.5 0.9 0.4 0.8 0.4 0.0 0.2 0.9 0.7 0.5 0.2 0.6 1.0 0.7 0.4 0.4 0.8 0.5 0.2 0.8 0.5 0.3 0.6 0.0 0.7 0.0 0.5 0.3 0.4 0.7 0.8 0.2 0.7 0.4 0.3 0.2 0.6 0.7 0.1 0.4 0.6 0.9 0.8 0.7 0.3 0.9 0.1 0.2 0.4];
u_n = rescale(Input(1:10000)); %Normalize input between 0 and 1
Y_n = normalize(Output(1:10000), 'range' , [-1 1]); %Normalize between -1 and 1

mask_len = length(input_mask);
signal_len = length(u_n);
    
reservoirmemory = zeros(mask_len,Pulse.fft_point);
readout_t = zeros(signal_len,Pulse.fft_point);
readout_f = zeros(signal_len,Pulse.fft_point);
readout_RS = zeros(signal_len, mask_len);
%% Algorithm test, Single delay loop with T = (N+1)Theta
% for i = 1:signal_len
%     time_multiplexed_signal = input_mask .* u_n(i);   % time multiplexing of input signal    
%     for k = 1:length(time_multiplexed_signal)-1   %T = (N+1)Theta
%         modulated_signal = time_multiplexed_signal(k) .* Pulse.guassian; %convert data to guassian input
%         MZM_output = sqrt(Pulse.Peak_power) * sin((pi/(2*Fiber.MZM_Vpi)) .* modulated_signal); %modulation of input data
%         memory_modulated_signal = (MZM_output .* sqrt(1 - Fiber.Inputcoupling)) + (1i * sqrt(Fiber.Inputcoupling) .* reservoirmemory(k,:)); %fiber coupling
%         [SSFM_t,SSFM_f]  = NonLinearEffect(memory_modulated_signal ,Pulse.soliton_order,Fiber.length);
%         %[SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
%         %SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
%         reservoirmemory(k+1,:) = SSFM_t .* sqrt(Fiber.Outputcoupling);  
%         readout_RS(i,k) = mean(abs(SSFM_t(1024:3068).*sqrt(1 - Fiber.Outputcoupling)).^2);
%         %figure; plot(Pulse.Normalize_tau,abs(reservoirmemory(k,:)).^2);
%     end
%     %%last data
%     modulated_signal = time_multiplexed_signal(length(time_multiplexed_signal)); % .* Pulse.guassian;
%     MZM_output = sqrt(Pulse.Peak_power) * sin((pi/(2*Fiber.MZM_Vpi)) .* modulated_signal);
%     memory_modulated_signal = (MZM_output .* sqrt(1 - Fiber.Inputcoupling)) + (1i * sqrt(Fiber.Inputcoupling) .* reservoirmemory(k,:)); %fiber coupling
%     %[SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
%     [SSFM_t,SSFM_f]  = NonLinearEffect(memory_modulated_signal ,Pulse.soliton_order,Fiber.length);
%     %SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
%     reservoirmemory(1,:) = SSFM_t .* sqrt(Fiber.Outputcoupling);  
%     readout_t(i,:) = abs(reservoirmemory(round(k/2),:)).^2;
%     readout_RS(i,mask_len) = mean(abs(SSFM_t(1024:3068).*sqrt(Fiber.Outputcoupling)).^2);
%     %figure; plot(Pulse.Normalize_tau,readout_t(i,:));
% end  

%% Multiple Delay algorithm according to paper
input_mask2 = [0.91	0.13	0.91	0.63	0.10	0.28	0.55	0.96	0.96	0.16	0.97	0.96	0.49	0.80	0.14	0.42	0.92	0.79	0.96	0.66	0.04	0.85	0.93	0.68	0.76	0.74	0.39	0.66	0.17	0.71	0.03	0.28	0.05	0.10	0.82	0.69	0.32	0.95	0.03	0.44	0.38	0.77	0.8	0.2	0.49	0.45	0.65	0.7	0.8	0.3];
delayVal = 9;
delayed_u_n = [zeros(delayVal,1);u_n(1:(signal_len-delayVal))];
%Algorithm test
for i = delayVal+1:signal_len
    time_multiplexed_signal = (input_mask .* u_n(i)) + (input_mask2 .* u_n(i-delayVal));    % time multiplexing of input signal    
    for k = 1:length(time_multiplexed_signal)-1   %T = (N+1)Theta
        modulated_signal = time_multiplexed_signal(k) .* Pulse.guassian; %convert data to guassian input
        MZM_output = sqrt(Pulse.Peak_power) * sin((pi/(2*Fiber.MZM_Vpi)) .* modulated_signal); %modulation of input data
        memory_modulated_signal = (MZM_output .* sqrt(1 - Fiber.Inputcoupling)) + (1i * sqrt(Fiber.Inputcoupling) .* reservoirmemory(k,:)); %fiber coupling
        [SSFM_t,SSFM_f]  = NonLinearEffect(memory_modulated_signal ,Pulse.soliton_order,Fiber.length);
        %[SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
        %SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
        reservoirmemory(k+1,:) = SSFM_t .* sqrt(Fiber.Outputcoupling);  
        readout_RS(i,k) = mean(abs(SSFM_t(1024:3068).*sqrt(1 - Fiber.Outputcoupling)).^2);
        %figure; plot(Pulse.Normalize_tau,abs(reservoirmemory(k,:)).^2);
    end
    %%last data
    modulated_signal = time_multiplexed_signal(length(time_multiplexed_signal)); % .* Pulse.guassian;
    MZM_output = sqrt(Pulse.Peak_power) * sin((pi/(2*Fiber.MZM_Vpi)) .* modulated_signal);
    memory_modulated_signal = (MZM_output .* sqrt(1 - Fiber.Inputcoupling)) + (1i * sqrt(Fiber.Inputcoupling) .* reservoirmemory(k,:)); %fiber coupling
    %[SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
    [SSFM_t,SSFM_f]  = NonLinearEffect(memory_modulated_signal ,Pulse.soliton_order,Fiber.length);
    %SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
    reservoirmemory(1,:) = SSFM_t .* sqrt(Fiber.Outputcoupling);  
    readout_t(i,:) = abs(reservoirmemory(round(k/2),:)).^2;
    readout_RS(i,mask_len) = mean(abs(SSFM_t(1024:3068).*sqrt(Fiber.Outputcoupling)).^2);
    %figure; plot(Pulse.Normalize_tau,readout_t(i,:));
end  

%% output data to file
%writematrix(readout_t,'NARMA_TASKIN_INIT2.csv');
writematrix(readout_RS,'NCEQ_TASKIN_RS2.csv');
writematrix(Y_n,'NCEQ_TASKOUT_RS2.csv');