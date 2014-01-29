function [Q,Qp, I2, I2up, I2down] = DataSetInconsistency(s,sd)
% this function provide Q and I^2 of two studies, whose parameters are in s, with standard deviation in sd
% The Q provides a measurement of heterogenuity between the parameters
% the I^2 provdes a measurement of the inconsistency between the paremeters
% s - p*1 vector of the parameter values from differenct data sets
% sd - p*1 vetor of the standard deviation of the parameters
% Q - the heterogenuity measurements
% Qp - the p-value of the Q using the chi square one-tail probabilty
% I2 - the inconsistency measurement. Negative I2 is put at 0 so I2 lies
% between 0 and 1. A value of 0 indicate no observed heterogeneity, and
% larger values show increasing heterogeneity.

% prepare
df = length(s)-1; % degree of freedom

% Q computation
wt = 1./(sd.^2); % the weight is 1/sigma^2
nbar = sum(wt.*s)/sum(wt);
% nsigma = 1/sqrt(w);
Q = sum(((s-nbar).^2).*wt);
varQ = var(((s-nbar).^2).*wt);
Qp = 1-chi2cdf(Q,df);

% I^2
% Old I2
% if Q > (df+1),
%     B = 0.5*((log(Q)-log(df))/(sqrt(2*Q)-sqrt((2*df)-1)));
% else
%     B = sqrt(1/(2*(df-1)*(1-(1/(3*((df-1)^2))))));
% end
% 
% L = exp(0.5*log(Q/df)-1.96*B);
% U = exp(0.5*log(Q/df)+1.96*B);
% 
% I2_LL = (L^2-1)/(L^2);
% I2_UL = (U^2-1)/(U^2);
% 
% I2up=I2_UL;
% I2down=I2_LL;



% Q & I2 calculation from [Higgins, 2002], eqns (2), (8.5) (9) (10)
% with CI calc from [Biggerstaff, 1997] (eq. 7)
% S1 = sum(wt);
% S2 = sum(wt.^2);
% S3 = sum(wt.^3);

if (Q==0),
    I2=0;I2up=0;I2down=0;
else
    cur_I2 = (Q-df)/Q;
    
    %tau2 = max(0,(Q-df)/(S1-(S2/S1)));
    
    % [Biggerstaff and Tweedie Eq (7)]
    %tmp
%     varQ = 2*df+...
%         4*(S1-(S2/S1))*tau2+...
%         2*(S2-2*(S3/S1)+((S2^2)/(S1^2)))*(tau2^2);
    
    % H2 = Q/df; H2 from [Higgins A1 (1554)]
    Qup = (Q+1.96*sqrt(varQ));
    Qdown = (Q-1.96*sqrt(varQ));
        
    relQ= abs(Qup-Qdown)/Q;
    relI2 = relQ;
    
    I2up = cur_I2*(1+(relI2/2));
    I2up = min(1,I2up);
    I2down = cur_I2*(1-(relI2/2));
    I2down = max(0,I2down);
    
    I2 = max(0,cur_I2);
    I2 = min(1,cur_I2);
%     
%    

end
