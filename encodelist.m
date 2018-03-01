% ENCODE INPUT LIST

function C = encodelist(const,parms,context,list,items)

% Create item-context weight matrix
C = zeros(const.nItems,parms.n); % This should never be reset surely??? 

for p = 1:const.ll
    
    % Calculate learning rate for current position
    eta = parms.theta * parms.gamma.^(p-1);
        
    % Form association between current context marker and current item
    C = C + (items(list(p),:)' * context(p,:)) .* eta;

end

C = C./norm(C);