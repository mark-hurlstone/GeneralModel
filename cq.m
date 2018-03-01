% MAIN FUNCTION FOR CQ MODEL

function preds = cq(parmarray,const,parms)

global score

% Retrieve parameter values from fminsearch
fminparms.noiseAL = parmarray(1);
fminparms.noiseSL = parmarray(2);
fminparms.consim  = parmarray(3);
fminparms.theta   = parmarray(4);
fminparms.c       = parmarray(5);
fminparms.r       = parmarray(6);

% Randomisation stuff
randn('state',parms.seed); 
rand( 'state',parms.seed);

% Initialisation stuff
items       = createItems(const);                % Retrieve orthogonal item environment
C           = zeros(const.nItems,parms.ndim);    % Create item-context weight matrix
lists       = genList(const);                    % Generate study lists
[~,sim_Mat] = similarity(const,parms,fminparms); % Create similarity weight matrix
W           = weightMatrix(const,parms);         % Competitive filter weight matrix

ll = [3 4];
%:::::::::: Run the model
for l = 1:2
    
    const.ll    = ll(l);
    [context,~] = createContext(const,parms,fminparms); % Retrieve distributed context vectors
    allRes      = zeros(const.nTrials,const.ll); 
    allRts      = zeros(const.nTrials,const.ll); 
    
    for t = 1:const.nTrials    
        
        %:::::::::::::::::::::::::: ENCODING ::::::::::::::::::::::::::::::
        
        % Retrieve study list 
        list = lists(t,:); 

        for e = 1:const.ll

            % Calculate learning rate for current position
            if (e <= 2)
                eta = fminparms.theta.^(e-1);
            elseif(e == 3)
                eta = fminparms.theta.^(e+2-1);
            else
                eta = fminparms.theta.^(e+4-1);
            end
            
            % Form association between current context marker and current item
            C = C + eta * (items(list(e),:)' * context(e,:));
            
        end

        %:::::::::::::::::::::::::: RETRIEVAL :::::::::::::::::::::::::::::
        
        % Netinput to activation layer at previous position
        netintmin1 = zeros(1,const.nItems);
        
        for r = 1:const.ll

            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Calculate activations in activation layer
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % Cue item layer with context marker for current position       
            f = C * context(r,:)';

            % Rescale the activations to sum to 1 
            f = f ./ sum(f);

            % Add some selection noise
            f = f' + fminparms.noiseAL .* randn(1,const.nItems); 

            % Now pass the activations in f through the activation function 
            negvals = netintmin1 < 0; f(negvals) = f(negvals) + netintmin1(negvals);
            
            % Select the winning item
            [a_f,win_f] = max(f);
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Compute activations in confusion/selection layer
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % Forward winning item to confusion/selection layer
            g = a_f * sim_Mat(win_f,:);
            
            % Add selection noise
            g = g + fminparms.noiseSL .* randn(1,const.nItems);

            % Enter netinputs into competitive filter
            [win_g,rt_g] = compFilt(const,parms,fminparms,g,W);
            
            % Store the winning item and its RT in allRes and allRts
            allRes(t,r) = win_g;
            allRts(t,r) = rt_g; 
            
            % Add output interference to C           
            C = C + randn(size(C)) .* parms.phi;

            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Inhibit recalled items
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            % Copy current activation profile over f into netintmin1
            netintmin1 = f;

            % Suppress the item recalled in the activation layer
            netintmin1(win_f) = parms.inh;
            
            % Suppress the item recalled in the selection layer
            netintmin1(win_g) = parms.inh;
            
            % Implement recovery from inhibition before the next output position
            netintmin1(netintmin1 < 0) = netintmin1(netintmin1 < 0) * exp(-fminparms.r);
                
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Normalise context-item weights
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        for i = 1:const.nItems, C(i,:) = C(i,:)/norm(C(i,:)); end
                
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Implement scoring routine
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scoring(const,parms,fminparms,lists,allRes,allRts);
    if l==1,score.accspc(end+1)=NaN; score.repspc(end+1)=NaN; score.rtspc(end+1)=NaN; end  
    if l==1,score.libfillin(end+1)=NaN; score.libinfill(end+1)=NaN; end
    preds.accspc(l,:)         = score.accspc;                         
    preds.repspc(l,:)         = score.repspc; 
    preds.spatldist(l,:)      = score.spatldist;                                          
    preds.libfillin(l,:)      = score.libfillin; 
    preds.libinfill(l,:)      = score.libinfill; 
    preds.fillin              = sum(score.libfillin)...
        /(sum(score.libfillin)+sum(score.libinfill));
    preds.protrusions         = sum(score.prospc.*const.nTrials)...
        / (sum(score.prospc.*const.nTrials)+sum(score.intspc.*const.nTrials));
    if l==2,preds.rtspc       = score.rtspc; end
    if l==2,preds.transdist   = score.transdist; end
    if l==2,preds.spatldistop = score.spatldistop; end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
end
