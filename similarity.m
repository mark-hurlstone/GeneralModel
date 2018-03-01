% FUNCTION FOR CALCULATING INTER-SIMILARITY BETWEEN ITEMS

function [cb_Mat sim_Mat] = similarity(const,parms,fminparms)

% Spatial coordinates
xy = [1 1; 2 1; 3 1; 1 2; 2 2; 3 2; 1 3; 2 3; 3 3];

% Find the Manhattan distances using pdist
cb = pdist(xy,'cityblock');

% Convert to matrix form
cb_Mat = squareform(cb);

% Calculate the matrix of similarity values
sim_Mat = zeros(const.nItems,const.nItems);
for i = 1:const.nItems
    for j = 1:const.nItems
        sim_Mat(i,j) = exp((-fminparms.c*abs(cb_Mat(i,i)-cb_Mat(i,j))));
    end
end