% CREATE SET OF DISTRIBUTED CONTEXT VECTORS WITH DEFINED INTER-SIMILARITY

function [context,cosines] = createContext(const,parms,fminparms)

detas       = ones(1,parms.n);
walshmatrix = walsh(parms.ndim)'; 
context     = zeros(parms.n,parms.ndim); 
% Create context vectors
for i = 1:parms.n
    for j = 1:parms.ndim
        if i == 1
            context(i,j) = walshmatrix(i,j);
        else
            for k = 1:i
                context(i,j) = context(i,j) + detas(k) * walshmatrix(k,j);
            end
        end
    end
    context(i,:) = normvec(context(i,:)')' .* sqrt(parms.ndim);
    for j = 1:i
        detas(j) = detas(j) * fminparms.consim;
    end
    detas(i+1) = sqrt(1 - (fminparms.consim * fminparms.consim));
end

% Normalize context vectors
%for i = 1:const.ll
%    context(i,:) = context(i,:) ./ norm(context(i,:));
%end

% Compute cosine similarity between c(i) and c(j)
cosines = zeros(parms.n);
for i = 1:parms.n
    for j = 1:parms.n
        cosines(i,j) = cosim(context(:,i), context(:,j));
    end
end