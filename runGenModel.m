 
clear all
close all

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%                                             %
%                                             %
% Competitive queuing model for simulating    %
% the monkey data of Botvinick et al. (2009)  %
%                                             %
% VERSION OF 24-07-16 			              %
%                                             %
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%                                     %
% Mark Hurlstone (2016)               %
% School of Psychology                %
% University of Western Australia     %
% E-mail: mark.hurlstone@uwa.edu.au   %
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

global score

%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%

const.ll          = 4;           % List length (*MAX = 9*)
const.nItems      = 9;           % N items in spatial array (*FIXED*)
const.nTrials     = 2500;        % Number of simulation trials
fminparms.noiseAL = .02;         % Noise for activation layer
fminparms.noiseSL = .01;         % Noise for confusion/selection layer

%%%%%%%%%%%%%%%%%%%%
% Model parameters %
%%%%%%%%%%%%%%%%%%%%

parms.n          = const.ll*4-2; % n nodes in context vector
parms.na         = const.ll*2;   % n 'active' nodes in context vector
parms.s          = 2;            % n nodes moved with each update of context C = 
parms.n          = const.ll;     % Number of item markers
parms.ndim       = 16;           % Dimensionality of item markers
fminparms.consim = .65;          % Similarity of item markers
fminparms.theta  = .85;          % Slope parameter for primacy gradient 
fminparms.c      = .08;          % Sensitivity parameter for spatial confusions
parms.inh        = -1;           % Item inhibition level
fminparms.r      = .25;          % Item recovery rate 
parms.phi        = .04;          % Output interference
parms.tau        = .8;           % Response threshold for competitive filter
parms.maxc       = 20;           % Maximum cycles for competitive filter
parms.alpha      = 1.1;          % Self-excitation for competitive filter  
parms.beta       = -.1;          % Mutual inhibition for competitive filter 

%%%%%%%%%%%%%%%%%%%%%%%
% Randomisation stuff %
%%%%%%%%%%%%%%%%%%%%%%%

parms.seed = 111211;           

%%%%%%%%%%%%%%%%%
% Fit The Model %
%%%%%%%%%%%%%%%%%

% Data for fitting
data.accspc      = [.97 .88 .62 NaN .97 .86 .54 .27];
data.rtspc       = [449 491 544 615];
data.repspc      = [0 0.005 .03 NaN 0 0.005 .02 .1];
data.transdist   = [.88 .11 .01];
data.spatldist   = [.73 .19 .08 .003 .42 .33 .22 .04];
data.spatldistop = [.66 .25 .09 .002 .385 .35 .23 .05 .39 .35 .23 .04];
data.fillin      = 0.82; % This is 46/(46+10)
data.protrusions = .33;  % This is 51/157

% Combine into a single data vector for fitting
obs = [data.accspc data.repspc data.transdist data.spatldist...
     data.fillin data.protrusions]; %data.spatldistop

% Converts parameters for fminsearch
temp = struct2cell(fminparms);
for i = 1:length(temp)
	parmarray(i) = temp{i};
end

% Parameter boundaries for fminsearch
parms.UB = ones(1,length(parmarray));
parms.LB = zeros(1,length(parmarray));

% Set up fminsearch
x = parmarray; 
nfuncevals = 750;
tolerance = .1;
defopts = optimset ('fminsearch');
options = optimset (defopts,'Display','iter','TolFun',tolerance,...
    'MaxFunEvals',nfuncevals);

% Commence fitting
[x,fval,dummy,output] = mywrapperLoopfmin(parmarray,obs,const,parms);

% Generate final predictions
preds = cq(x,const,parms);

% Calculate chi2, lnL, AIC, and BIC scores 
pred  = [reshape(preds.accspc',1,8) reshape(preds.repspc',1,8)... 
        preds.transdist reshape(preds.spatldist',1,8) ... %reshape(preds.spatldistop(:,:,2:4),1,12) ...
        preds.fillin preds.protrusions];
pred(isnan(pred))=0; obs(isnan(obs))=0;
pred(pred<eps) = eps; 
chi2 = sum(100*sum(((obs-pred).^2)./pred)); 
n    = length(obs) - sum(isnan(obs-pred));
ssd  = nansum((pred-obs).^2);        
lnL  = n * log(ssd/n);

% Here we have added an extra parameter for response suppression
AIC  = lnL + 2*(length(parmarray)+1);
BIC  = lnL + (length(parmarray)+1)*log(n);

% Print predictions
preds.accspc                                   
preds.repspc
preds.transdist
preds.spatldist
preds.spatldistop
preds.libfillin        
preds.libinfill        
preds.fillin              
preds.protrusions         

% Plot key results
% 1. Serial position curve
subplot(2,2,1)
plot(1:const.ll,preds.accspc(1,:),'Marker','s',...
    'MarkerSize',5,'Color','b','DisplayName',...
    'Accuracy (Length 3)')
hold on
plot(1:const.ll,preds.accspc(2,:),'Marker','o',...
    'MarkerSize',5,'Color','b','DisplayName',...
    'Accuracy (Length 4)')
plot(1:const.ll,preds.repspc(1,:),'Marker','s',...
    'MarkerSize',5,'Color','g','LineStyle','--',...
    'DisplayName','Repetitions (Length 3)')  
plot(1:const.ll,preds.repspc(2,:),'Marker','o',...
    'MarkerSize',5,'Color','g','LineStyle','--',...
    'DisplayName','Repetitions (Length 4)')     
xlim([.5 const.ll+.5])
ylim([0 1]);
title('Accuracy + Repetiton SPCs')
xlabel('Serial Position')
ylabel('Proportion')
legend('show');

% 2. Aggregate transposition error gradient
subplot(2,2,2)
plot(1:const.ll-1,preds.transdist,'Marker','o',...
    'MarkerSize',5,'Color','b','DisplayName',...
    'Observed')  
hold on
plot(1:const.ll-1,score.transdistchance,'Marker','o',... 
    'MarkerSize',5,'Color','r','LineStyle','--',...
    'DisplayName','Chance')
xlim([.5 const.ll-.5])
ylim([0 1]);
title('Transposition Errors')
xlabel('Transposition Distance')
ylabel('Proportion')
legend('show');

% 3. Aggregate spatial error gradient
subplot(2,2,3)
plot(1:4,preds.spatldist(1,:),'Marker','o',...
    'MarkerSize',5,'Color','b','DisplayName',...
    'Length 3')
hold on
plot(1:4,preds.spatldist(2,:),'Marker','s',...
    'MarkerSize',5,'Color','b','DisplayName',...
    'Length 4')
hold on
plot(1:4,score.spatldistchance,'Marker','o',...
    'MarkerSize',5,'Color','r','LineStyle','--',...
    'DisplayName','Chance')
xlim([.5 4.5])
ylim([0 1]);
title('Aggregate Spatial Errors')
xlabel('Spatial Distance')
ylabel('Proportion')
legend('show');

% 4. Spatial error gradients by output position
subplot(2,2,4)
col = {'b','g','c'}; mrk = {'o','s','^'};
lbl = {'Position 2','Position 3','Position 4'};
for i=2:const.ll
    plot(1:4,preds.spatldistop(:,:,i),'Marker',... 
        mrk{i-1},'Markersize',5,'Color',col{i-1},...
        'DisplayName',lbl{i-1})
    hold on
end
hold on
plot(1:4,score.spatldistchance,'Marker','o',...
    'MarkerSize',5,'Color','r','LineStyle','--',...
    'DisplayName','Chance')
xlim([.5 4.5])
ylim([0 1]);
title('Spatial Errors By Output Position')
xlabel('Spatial Distance')
ylabel('Proportion')
legend('show');
