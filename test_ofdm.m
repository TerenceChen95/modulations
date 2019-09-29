% Script for computing the Bit Error probability using OFDM modulation

nFFT = 16; % fft size
nDSC = 52; % number of data subcarriers
nBitPerSym = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nSym = 10^4; % number of symbols

EbN0dB = 6.7895; % SNR of QPSK with 10^-3 BER
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % converting to symbol to noise ratio

% Transmitter
ipBit = (2*(rand(1,nBitPerSym*nSym)>0.5)-1) + 1i*(2*(rand(1,nBitPerSym*nSym)>0.5)-1); 
ipBit = reshape(ipBit,nBitPerSym,nSym).'; % grouping into multiple symbolsa

% Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
xF = [zeros(nSym,6) ipBit(:,[1:nBitPerSym/2]) zeros(nSym,1) ipBit(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,5)] ;

% Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

% Appending cylic prefix
xt = [xt(:,[49:64]) xt];

% Concatenating multiple symbols to form a long vector
xt = reshape(xt.',1,nSym*80);

%Channel 
nTap=3;
E=0.5^2+0.25^2+0.125^2;
ht=[0.5/(sqrt(E)),0.25/sqrt(E),0.125/sqrt(E)];
L  = length(ht);
chanOut = conv(xt,ht);  
nt = 1/sqrt(2)*[randn(1,nSym*80+L-1) + 1j*randn(1,nSym*80+L-1)];

% Noise addition,the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
yt = sqrt(80/64)*chanOut + 10^(-EbN0dB/20)*nt; % additive white gaussian noise
yt = yt(:,[1:800000]);

% Receiver
yt = reshape(yt.',80,nSym).'; % formatting the received vector into symbols
yt = yt(:,[17:80]); % removing cyclic prefix

% converting to frequency domain
yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).'; 
yMod = yF(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]);

% QPSK demodulation
y_re = real(yMod); % real
y_im = imag(yMod); % imaginary
yMod(find(y_re < 0 & y_im < 0)) = -1 + -1*1i;
yMod(find(y_re >= 0 & y_im > 0)) = 1 + 1*1i;
yMod(find(y_re < 0 & y_im >= 0)) = -1 + 1*1i;
yMod(find(y_re >= 0 & y_im < 0)) = 1 - 1*1i;


% counting the errors
nErr = size(find(yMod - ipBit),2);
simBer = nErr/(nSym*nBitPerSym);
%{
%calculate PAR
par=zeros(nSym);
for i=1:nSym
    x=transpose(yMod(i,:));
    x_abs=abs(x);
    par(i)=10*log(max(x_abs.^2)/mean(x_abs.^2));
end
figure
hist(par);
%}
figure
semilogy(EbN0dB,simBer,'mx-','LineWidth',2);
grid on
legend('simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for QPSK using OFDM')
