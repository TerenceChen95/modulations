out=0.1;
QPAM_x=zeros(3,1);
QPAM_y=zeros(3,1);
for i=1:3
    QPAM_F=qfuncinv(out/0.75)^2/0.8;
    QPAM_x(i)=10*log10(QPAM_F);
    QPAM_y(i)=out;
    out=out/10;
end
Eb_N0_dB=QPAM_x;
N = 10^5; % number of symbols
alpha4pam = [-3 -1 1 3]; % 4-PAM alphabets
ipHat = zeros(1,N);

for i = 1:length(Eb_N0_dB)
    ip = randsrc(1,N,alpha4pam);
    s = (1/sqrt(5))*ip; % normalization of energy to 1
    n = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)]; % white guassian noise, 0dB variance
    y = s + 10^(-Eb_N0_dB(i)/20)*n; % additive white gaussian noise
    % demodulation
    r = real(y); % taking only the real part
    ipHat(find(r< -2/sqrt(5))) = -3;
    ipHat(find(r>= 2/sqrt(5))) = 3;
    ipHat(find(r>=-2/sqrt(5) & r<0)) = -1;
    ipHat(find(r>=0 & r<2/sqrt(5))) = 1;
    nErr(i) = size(find([ip- ipHat]),2); % couting the number of errors
end
simBer = nErr/N;
figure
semilogy(QPAM_x,QPAM_y,'b.-');
hold on
semilogy(QPAM_x,simBer,'mx-');
grid on
legend('theory', 'simulation');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 4-PAM modulation')