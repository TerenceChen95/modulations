%AWGN Channels
%step1
%BPSK: Eb=E/log2(2)=E
BPSK_eqn1= '0.5*qfunc(sqrt(2*SNR))==0.1';
BPSK_eqn2= '0.5*qfunc(sqrt(2*SNR))==0.01';
BPSK_eqn3= '0.5*qfunc(sqrt(2*SNR))==0.001';
out=0.1;
BPSK_x=zeros(3,1);
BPSK_y=zeros(3,1);
for i=1:3
    BPSK_F=(qfuncinv(out/0.5)^2)/2;
    fprintf('BPSK:The SNR to achieve %f BER is %d\n',out,10*log10(BPSK_F));
    BPSK_x(i)=10*log10(BPSK_F);
    BPSK_y(i)=out;
    out=out/10;
end
%QPSK: Eb=E/log2(4)=E/2
QPSK_eqn1='qfunc(sqrt(2*SNR)==0.1)';
QPSK_eqn2='qfunc(sqrt(2*SNR)==0.01)';
QPSK_eqn3='qfunc(sqrt(2*SNR)==0.001';
out=0.1;
QPSK_x=zeros(3,0);
QPSK_y=zeros(3,0);
for i=1:3
    QPSK_F=qfuncinv(out)^2/2;
    fprintf('QPSK:The SNR to achieve %3f BER is %d\n',out,10*log10(QPSK_F));
    QPSK_x(i)=10*log10(QPSK_F);
    QPSK_y(i)=out;
    out=out/10;
end
%4-PAM: Pb=0.75*qfunc(sqrt(0.8*SNR))
out=0.1;
QPAM_x=zeros(3,1);
QPAM_y=zeros(3,1);
for i=1:3
    QPAM_F=qfuncinv(out/0.75)^2/0.8;
    fprintf('4-PAM:The SNR to achieve %f BER is %d\n',out,10*log10(QPAM_F));
    QPAM_x(i)=10*log10(QPAM_F);
    QPAM_y(i)=out;
    out=out/10;
end
%step2
%BPSK

openfig('ber1'); %simulated plot with BER analysis
hold on
plot(BPSK_x,BPSK_y,'ro')
hold off
xlabel('E_b/N_0(dB)')
ylabel('Bit Error Rate(10^y)')
title('BPSK')

%{
%QPSK
openfig('ber2'); %simulated result
hold on
plot(QPSK_x,QPSK_y,'ro');
hold off
xlabel('E_b/N_0(dB)')
ylabel('Bit Error Rate(10^y)')
title('QPSK')
%}
%{
%QPAM
openfig('ber3'); %simulated result
hold on
plot(QPAM_x,QPAM_y,'ro');
hold off
xlabel('E_b/N_0(dB)')
ylabel('Bit Error Rate(10^y)')
title('QPAM')
%}



