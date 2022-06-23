function [NonLinearResponse_t,NonLinearResponse_f] = NonLinearEffect(signal,solitonOrder,FiberLength)
    Phase_factor = (1i * solitonOrder^2 * FiberLength)  ;  %phase factor
    NonLinearResponse_t =  signal .* exp( Phase_factor.* abs(signal).^2); % Non linear effect over the span of the fiber on the original signal
    NonLinearResponse_f = fftshift(fft(NonLinearResponse_t));
end
