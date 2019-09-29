%use 0,1,2,3 as 4-PAM values
%apply gray code
n=7;
k=4;
A=[0 0 0;0 0 1;0 1 1;0 1 0]; 
G=[eye(k) A];    %Generator matrix
H=[A' eye(n-k)]; %Parity-check matrix
%Encoder
msg=[1 1 1 0];
code = mod(msg*G,2);  
% CHANNEL ERROR(add one error to code)%
%code(1)= ~code(1);
%code(2)= ~code(2);
%code(3)= ~code(3);
code(4)= ~code(4);%Pick one,comment out others
%code(5)= ~code(5);
%code(6)= ~code(6);
%code(7)= ~code(7);
recd = code; %Received codeword with error
% DECODER%
syndrome = mod(recd * H',2);
%Find position of the error in codeword (index)
find = 0;
for ii = 1:n
    if ~find
        errvect = zeros(1,n);
        errvect(ii) = 1;
        search = mod(errvect * H',2);
        if search == syndrome
            find = 1;
            index = ii;
        end
    end
end
disp(['Position of error in codeword=',num2str(index)]);
correctedcode = recd;
correctedcode(index) = mod(recd(index)+1,2); %Corrected codeword
%Strip off parity bits
msg_decoded=correctedcode;
msg_decoded=msg_decoded(1:4);

