% calculate minimum phase component of impulse response
% based on the invFIR function - look there for explanations
% https://uk.mathworks.com/matlabcentral/fileexchange/19294-inverse-fir-filter
% alternatively use [y,ym] = rceps(x)
function [h_min] = minph(h)
    [~,h_min] = rceps(h);
%     n = length(h);
%     h_cep = real(ifft(log(abs(fft(h(:,1))))));
%     odd = fix(rem(n,2));
%     wn = [1; 2*ones((n+odd)/2-1,1) ; ones(1-rem(n,2),1); zeros((n+odd)/2-1,1)];
%     h_min = zeros(size(h(:,1)));
%     h_min(:) = real(ifft(exp(fft(wn.*h_cep(:)))));
end