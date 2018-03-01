% GENERATE INPUT LISTS

function lists = genList(const)

lists = zeros(const.nTrials,const.nItems);
n = 0;
while n < const.nTrials    
    list = randperm(const.nItems);
    
    if (sum(list(1:const.ll) ~= 5) == const.ll)
    
    if (n == 0)
        n = n+1;
        lists(n,:) = list;
    else
        if (sum(list == lists(n,:)) == 0) 
            n = n+1;
            lists(n,:) = list;
        end
    end
    
    end
end
lists = lists(:,1:const.ll);

% while n < const.nTrials    
%     list = randperm(const.nItems);
%     list = list(1:const.ll);
%     
%     if (sum(list(1:const.ll) ~= 5) == const.ll)
%     
%     if (n == 0)
%         n = n+1;
%         lists(n,:) = list;
%     else
%         if (sum(list == lists(n,:)) == 0) 
%         %if (sum(ismember(list,lists(n,:))) == 0) 
%             n = n+1;
%             lists(n,:) = list;
%         end
%     end
%     
%     end
% end