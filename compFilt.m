% FUNCTION FOR IMPLEMENTING COMPETITIVE FILTER

function [win,rt] = compFilt(const,parms,fminparms,vin,W)

% vin  = max(zeros(1,const.nItems),vin);
% vout = W * vin';
% for t = 1:parms.maxc
% 
%     if (max(vout) > parms.tau)
%         [~,win] = max(vout);
%         rt = t;
%     break
% 
%     else
%         vin  = max(zeros(1,const.nItems),vout');
%         vout = W*vin' + fminparms.noiseSL .* randn(1,const.nItems)';
%     continue
%     
%     end
%     
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% UNCOMMENT THIS CODE (AND COMMENT THE ABOVE) TO RUN THE COMPETITIVE FILTER 
% FOR A FIXED NUMBER OF ITERATIONS AT EACH POSITION BY SETTING PARMS.MAXC TO
% HOWEVER MANY ITERATIONS YOU WANT THE FILTER TO OPERATE FOR

vin  = max(zeros(1,const.nItems),vin);
vout = W * vin';
for t = 1:parms.maxc
    
    vin  = max(zeros(1,const.nItems),vout'); 
    vout = W*vin' + fminparms.noiseSL .* randn(1,const.nItems)';
%     vin   = max(zeros(1,const.nItems),vout');
%     noise = fminparms.noiseSL .* randn(1,const.nItems); noise(5) = 0; 
%     vout = W*vin' + noise';
       
end
[~,win] = max(vout);
rt = t;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~