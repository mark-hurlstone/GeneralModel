% FUNCTION FOR CREATING LATERAL INHIBITION WEIGHT MATRIX

function W = weightMatrix(const,parms)

W = zeros(const.nItems,const.nItems);
for i=1:const.nItems
	for j=1:const.nItems
		if (i==j)
			W(i,j) = parms.alpha;
		else
			W(i,j) = parms.beta;
		end
	end
end