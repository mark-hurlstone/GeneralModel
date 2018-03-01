% FUNCTION FOR SCORING MODEL RESPONSES

function scoring(const,parms,fminparms,lists,allRes,allRts) 

global score

% Recode responses by input position
% -------------------------------------------------------------------------
allRes2 = zeros(const.nTrials,const.ll);
for i=1:const.ll     
    for j=1:const.ll
        match = allRes(:,i)==lists(:,j); 
        allRes2(match,i) = j; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% SERIAL POSITION CURVES %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accuracy SPC
% -------------------------------------------------------------------------
acctmp = allRes2 == repmat(1:const.ll,const.nTrials,1);
score.accspc = sum(acctmp)./const.nTrials; 

% Latency SPC
% -------------------------------------------------------------------------
incorrect = allRts == 0;
allRts(incorrect) = 0; 
score.rtspc = sum(allRts)./sum(acctmp);

% Intrusion & protrusion SPCs
% -------------------------------------------------------------------------
protrusions = zeros(const.nTrials,const.ll);
intrusions = zeros(const.nTrials,const.ll);
for i=2:const.nTrials
    for j=1:const.ll
        if (allRes(i,j)==lists(i-1,j) ...
                && sum(ismember(allRes(i,j),lists(i,:)))==0) 
            protrusions(i,j) = 1;
        elseif (ismember(allRes(i,j),lists(i-1,:)) ... 
                && sum(ismember(allRes(i,j),lists(i,:)))==0)
            intrusions(i,j) = 1;
        end
    end
end
score.prospc = sum(protrusions)./const.nTrials;
score.intspc = sum(intrusions)./const.nTrials;
protrusions = protrusions>0; intrusions = intrusions>0;   
allRes2(protrusions) = NaN; allRes2(intrusions) = NaN;   

% Repetition SPC
% -------------------------------------------------------------------------
repetitions1 = zeros(const.nTrials,const.ll);
repetitions2 = zeros(const.nTrials,const.ll);
for i=1:const.ll-1
    for j=i+1:const.ll
        reps1 = allRes2(:,i) == allRes2(:,j);
        reps2 = allRes2(:,j) == allRes2(:,i);
        repetitions1(:,j) = repetitions1(:,j)+reps2;
        repetitions2(:,i) = repetitions2(:,i)+reps1;
        repetitions2(:,i) = repetitions2(:,i)+reps2;
    end
end
%---- Exclude correct recalls ----
repetitions1 = repetitions1>0; 
correctRepeats = allRes2==repmat(1:const.ll,const.nTrials,1);
repetitions1(correctRepeats) = 0;
score.repspc = sum(repetitions1)./const.nTrials;
allRes2(repetitions1) = NaN;

% Transposition SPC
% -------------------------------------------------------------------------
transpositions = allRes2 ~= repmat(1:const.ll,const.nTrials,1);
removeItemErrors = isnan(allRes2); transpositions(removeItemErrors) = 0;
score.transpc = sum(transpositions)./const.nTrials;


%%%%%%%%%%%%%%%%%%%
% ERROR GRADIENTS %
%%%%%%%%%%%%%%%%%%%

% 1. Spatial error gradients
% -------------------------------------------------------------------------
[cb_Mat,~] = similarity(const,parms,fminparms);
distances  = zeros(const.nTrials,max(max(cb_Mat))); 
for i=1:const.nTrials
    for j=1:const.ll
        correct = lists(i,j); 
        response = allRes(i,j);        
        if (isnan(correct)==0 && isnan(response)==0)
            distances(i,j) = cb_Mat(correct,response);
        else
            distances(i,j) = NaN;
        end
    end
end
sdistance = 1:max(max(cb_Mat));
spatldist = zeros(1,max(max(cb_Mat)));
spatldistop = zeros(1,max(max(cb_Mat)),const.ll); 
for i=1:const.ll
    for j=1:max(max(cb_Mat))
        conftmp = distances(:,i) == sdistance(j); 
        spatldistop(:,j,i) = sum(conftmp');
    end
    spatldist = spatldist+spatldistop(:,:,i);
end
score.spatldistop = spatldistop ./repmat(sum(spatldistop),1,max(max(cb_Mat)));
score.spatldist = spatldist ./sum(spatldist);

% Calculate chance levels
spatldistchance = zeros(1,const.ll-1);
for i=1:max(max(cb_Mat))
    spatldistchance(i) = sum(sum(cb_Mat==i));
end
score.spatldistchance = spatldistchance./sum(spatldistchance);


% 2. Transposition error gradients
% -------------------------------------------------------------------------

% Remove residual zeros in allRes2 (extra-list intrusions need removing)
censor = allRes2 == 0; allRes2(censor) = NaN;

tdistance = 1:const.ll-1; 
transdist = zeros(const.nTrials,const.ll-1);
for i=1:const.ll-1
    transtmp = abs((repmat(1:const.ll,const.nTrials,1)-allRes2)==tdistance(i)); 
	transdist(:,i) = sum(transtmp');
end
score.transdist = sum(transdist)./sum(sum(transdist)); 

% Calculate chance levels
transdistchance = zeros(1,const.ll-1); 
for i=1:const.ll-1
    if (i==1) 
        transdistchance(i) = const.ll*2-2;
    else
        transdistchance(i) = transdistchance(i-1)-2;
    end
end
score.transdistchance = transdistchance./sum(transdistchance);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL-IN & INFILL ERRORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Immediate anticipations + postponements
% -------------------------------------------------------------------------
minus1 = (repmat(1:const.ll,const.nTrials,1)-allRes2(:,1:const.ll))==-1;
plus1 = (repmat(1:const.ll,const.nTrials,1)-allRes2(:,1:const.ll))==1;

libfillint = minus1(:,1:(const.ll-2)) & plus1(:,2:(const.ll-1));
score.libfillin = sum(libfillint);

libinfillt = minus1(:,1:(const.ll-2)) & minus1(:,2:(const.ll-1));
score.libinfill = sum(libinfillt);
