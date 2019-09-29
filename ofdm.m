% Script for computing the Bit Error probability using OFDM modulation

nFFT = 16; % fft size
nDSC = 52; % number of data subcarriers
nBitPerSym = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nSym = 10^4; % number of symbols

EbN0dB = 6.7895; % SNR of QPSK with 10^-3 BER
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % converting to symbol to noise ratio


% Transmitter
ipBit = (2*(rand(1,nBitPerSym*nSym)>0.5)-1) + 1i*(2*(rand(1,nBitPerSym*nSym)>0.5)-1); 
%ipMod = real(ipBit); 
ipMod = reshape(ipBit,nBitPerSym,nSym).'; % grouping into multiple symbolsa

% Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
xF = [zeros(nSym,6) ipMod(:,[1:nBitPerSym/2]) zeros(nSym,1) ipMod(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,5)] ;

% Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

% Appending cylic prefix
xt = [xt(:,[49:64]) xt];

% Concatenating multiple symbols to form a long vector
xt = reshape(xt.',1,nSym*80);


% Gaussian noise of unit variance, 0 mean
nt = 1/sqrt(2)*[randn(1,nSym*80) + j*randn(1,nSym*80)];

% Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
yt = sqrt(80/64)*xt + 10^(-EsN0dB/20)*nt;

% Receiver
yt = reshape(yt.',80,nSym).'; % formatting the received vector into symbols
yt = yt(:,[17:80]); % removing cyclic prefix

% converting to frequency domain
yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).'; 
yMod = yF(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]); 


% QPSK demodulation
ipModHat = 2*floor(real(yMod/2)) + 1;
y_re = real(ipModHat); % real
y_im = imag(ipModHat); % imaginary
ipModHat(find(y_re < 0 & y_im < 0)) = -1 + -1*1i;
ipModHat(find(y_re >= 0 & y_im > 0)) = 1 + 1*1i;
ipModHat(find(y_re < 0 & y_im >= 0)) = -1 + 1*1i;
ipModHat(find(y_re >= 0 & y_im < 0)) = 1 - 1*1i;

% converting modulated values into bits
ipBitHat = (ipModHat+1)/2;
ipBitHat = reshape(ipBitHat.',nBitPerSym*nSym,1).';

   % counting the errors
   nErr = size(find(ipBitHat - ipBit),2);

simBer = nErr/(nSym*nBitPerSym);
theoryBer = (1/2)*erfc(sqrt(10.^(EbN0dB/10)));

close all; figure
semilogy(EbN0dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(EbN0dB,simBer,'mx-','LineWidth',2);
axis([0 10 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for BPSK using OFDM')

