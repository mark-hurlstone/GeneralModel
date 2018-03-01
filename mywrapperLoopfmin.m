% FUNCTION FOR IMPLEMENTING FMINSEARCH

function [bestx,bestFval,bestdummy,bestoutput] = mywrapperLoopfmin(parmarray,obs,const,parms) 

bestFval = realmax;
bestx = parmarray.*0;
bestdummy = 0;
bestoutput = 0;

for p1 = .02                     % fminparms.noiseAL                          
    for p2 = .01                 % fminparms.noiseSL          
        for p3 = .75             % fminparms.consim
            for p4 = .85         % fminparms.theta
                for p5 = .04     % fminparms.c
                    for p6 = .25 % fminparms.r                         
                        tic
                        [x,fval,dummy,output] = fminsearchbnd(@bof,...
                            [p1 p2 p3 p4 p5 p6],parms.LB,parms.UB);
                        toc
                        if fval < bestFval
                            bestFval = fval;
                            bestx = x;
                            bestdummy = dummy;
                            bestoutput = output;
                        end
                    end
                end
            end
        end
    end
end

% Nested function inherits data from wrapper function
    function chi2 = bof(fminparms)
        preds = cq(fminparms,const,parms);
        pred  = [reshape(preds.accspc',1,8) reshape(preds.repspc',1,8)...
            preds.transdist reshape(preds.spatldist',1,8) ... 
            preds.fillin preds.protrusions];    
        pred(isnan(pred))=0; obs(isnan(obs))=0;
        pred(pred<eps) = eps; 
        chi2 = sum(100*sum(((obs-pred).^2)./pred));               
    end
end
