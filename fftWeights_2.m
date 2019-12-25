clc;clear;
nfft = 1024;
w = {};
pos = 1;
M = 1;
F1 = fimath ('RoundingMethod','Nearest', 'OverflowAction', 'Saturate', ...
                 'ProductMode', 'SpecifyPrecision' , ...
                 'ProductFractionLength', 10,...
                 'ProductWordLength', 16, ...
                 'SumMode', 'SpecifyPrecision',...
                 'SumFractionLength', 10,...
                 'SumWordLength', 16);
             
%              
% while M < nfft
%     w{pos} = exp(-1i*2*pi*(0:(M-1))/(2*M));
%     M = M*2;
%     pos = pos+1;
%     
%     
% end

layer = 1;
while M < nfft
    wM = exp(-1i*2*pi*(0:(M-1))/(2*M));
    %wM = fi(wM,1,8,6,F1);

    wM = fi(wM,1,12,10,F1);

    Istep = M*2;
    
    for mcnt = 1:M
        for icnt = mcnt:Istep:nfft 
            j = icnt + M;
            Temp = wM(mcnt);
            strp = sprintf('w%d', (mcnt-1)*nfft/M/2 );
            %strn = sprintf('-w%d',2*(mcnt-1));
            wReal(j,layer) =  real(-1*Temp);
            wImag(j,layer) = imag(-1*Temp);
            wReal(icnt,layer) = real(Temp);
            wImag(icnt,layer) = imag(Temp);
            
            coeffStr{j,layer} = strp;
            %coeffStr{icnt,layer} = strp;
        end
    end
    
    M = M*2;
    layer = layer+1;
end

weights = exp(1i*2*pi*(0:(nfft-1))/(nfft));

% w is a 11 element cell array. Each cell is a variable length vector.