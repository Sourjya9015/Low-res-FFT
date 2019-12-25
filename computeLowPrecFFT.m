function [ Af ] = computeLowPrecFFT( At, nfft, nprec)
%COMPUTELOWPRECFFT Summary of this function goes here
%   Detailed explanation goes here

if (size(At,1) ~= nfft)
    error('Input vector length should be consistant with number of FFT points');
end

 Af = zeros(size(At)) + 1i*zeros(size(At));
 
 F = fimath ('RoundingMethod','Nearest', 'OverflowAction', 'Saturate', ...
            'ProductMode', 'SpecifyPrecision' , ...
            'ProductWordLength', 16, ...
            'ProductFractionLength', 16-4,  ...
            'SumMode', 'SpecifyPrecision', ...
            'SumWordLength', 16, 'SumFractionLength', 16-4);
        
% F = fimath ('RoundingMethod','Nearest', 'OverflowAction', 'Saturate', ...
%                 'ProductMode', 'FullPrecision' , 'SumMode', 'FullPrecision');

 indx = 0:(nfft-1);
 binIndx = dec2bin(indx);
 binIndx = fliplr(binIndx);
 indx = bin2dec(binIndx) + 1;
 At = At(indx,:);
 
 %At = fi(At, 1, nprec, (nprec/2), F);
 %Af = fi (Af, 1, nprec, (nprec/2), F);
 weight = cell(log2(nfft));
 for iWt=1:log2(nfft)
     w0 = exp(-1i*2*pi*(0:(2^iWt/2-1))/(2^iWt));
     if (iWt < 3)
        weight{iWt} = fi(w0, 1, 2, 0, F);
%      elseif (iWt >= 3 && iWt <=8)
%         weight{iWt} = fi(w0, 1, 6, 4, F);
     else
        weight{iWt} = fi(w0, 1, nprec, nprec-2, F);
     end
 end
 
 % loop over the matrix
 for i=1:size(At,2)
     vec = fi(At(:,i), 1, 16, 14);
     M = 1;
     % Butterfly loop
     while nfft > M
         Istep = bitshift(M, 1);
         %w0 = exp(-1i*2*pi*(0:(Istep/2-1))/Istep);
         %prec = (Istep < nprec)*Istep + (Istep >= nprec)*nprec;
         %prec = (prec < 4) * 4 + (prec >=4)*prec;
         %w = fi(w0, 1, prec, (prec-2), F);
         
         w = weight{log2(Istep)};
         for mcnt = 1:M
             for icnt = mcnt:Istep:nfft 
                 j = icnt + M;
                 Temp = w(mcnt)* vec(j);
                 vec(j) = vec(icnt) - Temp;
                 vec(icnt) = vec(icnt) + Temp;
             end
         end
         M = Istep;
     end
     
     Af(:,i) = vec;
 end

end

