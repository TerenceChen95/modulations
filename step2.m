
SNRdB=[-4.5080,3.2406,6.1720];           %SNR in dB
SNR=10.^(SNRdB./10);                     %SNR in linear scale   
info_word_length=100000;                 %No. of information words   
n=7;k=4;                                 %Parameters of hamming code   
ber=zeros(length(SNR),1);                %Simulated BER   
info_word=floor(2*rand(k,info_word_length));    %Generation of 0 and 1 for infromation bits

%use 0,1,2,3 as 4-PAM values
%apply gray code
A=[0 0 0;0 0 1;0 1 1;0 1 0]; 
G=[eye(k) A];    %Generator matrix
H=[A' eye(n-k)]; %Parity-check matrix

code_bit5=xor(info_word(1,:),xor(info_word(2,:),info_word(3,:)));   %First Parity Bit
code_bit6=xor(info_word(1,:),xor(info_word(3,:),info_word(4,:)));   %Second Parity Bit
code_bit7=xor(info_word(1,:),xor(info_word(2,:),info_word(4,:)));   %Third Parity Bit

code_word=[info_word;code_bit5;code_bit6;code_bit7];       %Coded information Word with parity bits

code_word(code_word==0)=-1;              %Converting 0 bits to 1  

decoded_block=zeros(n,info_word_length);           %SOFT Decoding Output

C=de2bi((0:2^(k)-1));                       %All bits of length k(Stored in valid code words matrix 'C')
C(1:16,5)=xor(C(:,1),xor(C(:,2),C(:,3)));   %First Parity Bit
C(1:16,6)=xor(C(:,1),xor(C(:,3),C(:,4)));   %Second Parity Bit
C(1:16,7)=xor(C(:,1),xor(C(:,2),C(:,4)));   %Third Parity Bit
distance=zeros(1,2^k);
for i=1:length(SNR)
    y=(sqrt(SNR(i))*code_word)+randn(n,info_word_length);     %Received Codes
     
    %Decoding Received Codes into valid codewords
    for l=1:info_word_length
        %SOFT Decoding
        for m=1:(k^2)           %Tacking distance of each column of the received word to a valid codeword
            distance(m)=norm(y(:,l)-C(m,:)');
        end
        [minval,minind]=min(distance);       %Finding index of the minimum distance valid codeword
        
        decoded_block(:,l)=C(minind,:);      %Decoding as the min distance codeword 
    end
    ber(i,1)=length(find(decoded_block(1:4,:)~=info_word));   %BER in BLOCK Detection
    
end
ber=ber/(k*info_word_length);
%{
out=0.1;
BPSK_x=zeros(3,1);
BPSK_y=zeros(3,1);
for i=1:3
    BPSK_F=(qfuncinv(out/0.5)^2)/2;
    BPSK_x(i)=10*log10(BPSK_F);
    BPSK_y(i)=out;
    out=out/10;
end
Eb_N0_dB=BPSK_x;
%BPSK
ipHat = zeros(1,info_word_length);
for i=1:length(Eb_N0_dB)
    input=info_word;
    s=2*input-1;    %BPSK modulation
    noise = 1/sqrt(2)*[randn(1,info_word_length) + 1i*randn(1,info_word_length)]; % white guassian noise, 0dB variance
    y_bpsk = s + 10^(-Eb_N0_dB(i)/20)*noise; % additive white gaussian noise

    ipHat = real(y_bpsk)>0;
    nErr(i) = size(find([input- ipHat]),2); % couting the number of errors   
end
simBer = nErr/info_word_length;
%}
open('bpsk.fig');
hold on
semilogy(SNRdB,ber(:,1),'m-<','linewidth',2.0)    %Simulated BER in Hamming code

title('BLOCK Detection for (7,4) Hamming Code');xlabel('SNR(dB)');ylabel('BER');
legend('BPSK(uncoded)','Simulated BER(Hamming code)');
axis tight
grid on



