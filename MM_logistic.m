% This Matlab code provides a function that uses the Newton-Raphson algorithm
% to calculate ML estimates of a simple logistic regression. Most of the
% code comes from Anders Swensen, "Non-linear regression." There are two
% elements in the beta vector, which we wish to estimate.
%
function [beta,llhd] = MM_logistic(data,cg,beta_start,do_plot)
x=data(:,1); % x is first column of data
y=data(:,2); % y is second column of data

beta = beta_start; % initial values
% survival/complication time
f2 = ~cellfun('isempty',{cg.mGrp.mDateComp}); % patients with no complication date
f3 = ~cellfun('isempty',{cg.mGrp.mDateLastFollowup}); % patients with no last follow up date
compdate = inf(cg.mNumInGrp,1);
lastfollowup = inf(cg.mNumInGrp,1);
compdate(f2) = ([cg.mGrp(f2).mDateComp] - [cg.mGrp(f2).mDateBaseline])' / 30;
lastfollowup(f3) = ([cg.mGrp(f3).mDateLastFollowup] - [cg.mGrp(f3).mDateBaseline])' / 30;
compdate = min( lastfollowup, compdate );
flgcensor = [cg.mGrp.mFlgCensor]';

% fit with just complication datesj only
compdates = compdate(~flgcensor);
mu_sig = lognfit(sort(compdates));

if do_plot
    fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';
    screen_size=get(0,'ScreenSize');
    cur_fig=figure(1000);clf reset;hold on;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    x_axis = linspace(0,max(compdate),500);
   
    [ax,h1,h2]=plotyy(x_axis,lognpdf(x_axis,mu_sig(1),mu_sig(2)),...
        x_axis,logncdf(x_axis,mu_sig(1),mu_sig(2)));hold on;
    set(h1,'LineWidth',2);
    set(h2,'LineWidth',2);
    %set(ax(1),'YLim',[0.001 1]);
    set(get(ax(1),'Ylabel'),'String','P','FontSize',15);    

    %set(ax(2),'YLim',[1 10]);
    set(get(ax(2),'Ylabel'),'String','P','FontSize',15);
    %semilogy(x_axis,repmat(0.05,length(x_axis)),'r--','LineWidth',1);
    set(ax(1),'XLim',[min(x_axis) max(x_axis)]);
    set(ax(2),'XLim',[min(x_axis) max(x_axis)]);
    set(ax,'FontSize',14);
    set(get(ax(2),'XLabel'),'String','Complication Date','FontSize',15);
    title('Log-Normal fits to complication data','FontSize',15);
    lgnd=legend([h1 h2],'PDF','CDF','Location','Best');
    set(lgnd,'FontSize',15);
    
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'log_normal_fits'],'-pdf');
    
end
%[sorted_dates, dates_idx] = sort(compdate);

f_t = lognpdf(compdate,mu_sig(1),mu_sig(2));
F_t = logncdf(compdate,mu_sig(1),mu_sig(2));


%%tmp
f_t(f_t==0)=0.0001;
F_t(F_t==0)=0.0001;

        
%tmp_p = exp(beta_start(1)+beta_start(2)*x)./(1+exp(beta_start(1)+beta_start(2)*x));
%ll = sum(y.*log(tmp_p)+(1-y).*log(1-tmp_p));
%ll = sum(y.*log(tmp_p.*f_t)+(1-y).*log((1-(tmp_p.*F_t))));

b1 = linspace(beta_start(1)-2.5,beta_start(1)+2.5,100);
b2 = linspace(beta_start(2)-3,beta_start(2)+3,100);
simpost = zeros(100,100);
for i=1:length(b1)
    for j=1:length(b2)
        tmp_p = exp(b1(i)+b2(j)*x)./(1+exp(b1(i)+b2(j)*x));
        %MLE
        %ll = sum(y.*log(tmp_p)+(1-y).*log(1-tmp_p));
        %MM
        ll = sum(y.*log(tmp_p.*f_t)+(1-y).*log((1-(tmp_p.*F_t))));
        %ll = sum(y.*log(tmp_p)+(1-y).*log((1-(tmp_p.*F_t))));
        simpost(i,j)=ll;
    end;
end;

[~,position] = max(simpost(:));
[b1_idx,b2_idx] = ind2sub(size(simpost),position);
beta(1) = b1(b1_idx);
beta(2) = b2(b2_idx);
llhd = simpost(b1_idx,b2_idx);
end


%         
%         
% while diff>0.0001 % convergence criterion
% beta_old = beta;
% p = exp(beta(1)+beta(2)*x)./(1+exp(beta(1)+beta(2)*x));
% p_org = p;
% % (1 - (p.*F_t)) = (1-p) + p*(1-F_t);
% p = (p.*f_t).*y + (1-(p.*F_t)).*(1-y);
% 
% %l_org = sum(y.*log(p)+(1-y).*log((1-p)))
% %MM
% %l = sum(y.*log(p)+(1-y).*log((1-p) + p.*prob_no_obs))
% 
% s = [sum(y-p); % scoring function
% sum((y-p).*x)];
% J_bar = [sum(p.*(1-p)) sum(p.*(1-p).*x); % information matrix
% sum(p.*(1-p).*x) sum(p.*(1-p).*x.*x)]
% beta = beta_old + J_bar\s % new value of beta
% diff = sum(abs(beta-beta_old)) % sum of absolute differences
% end


% get probability of no observation
% from KM
% flgcensor = sa.mFlgCensor;
% comp_times = sa.mSurvivalTime{1};
% prob_no_obs_sorted = sa.mSurvivalCurve{1};
% prob_obs_sorted = 1-prob_no_obs_sorted;
% norm_prob_obs_sorted = prob_obs_sorted./max(prob_obs_sorted);
% prob_no_obs_sorted = 1-norm_prob_obs_sorted;
% 
% comp_times_sorted = sa.mSurvivalTimeSorted{1};
% 
% prob_no_obs = -inf(size(comp_times));
% for i=1:length(comp_times)
%     cur_time = comp_times(i);
%     cur_times = find(cur_time<comp_times_sorted);
%     if ~isempty(cur_times)
%         cur_time_idx = cur_times(1);
%         prob_no_obs(i) = prob_no_obs_sorted(cur_time_idx);
%     else
%         prob_no_obs(i) = 0;
%     end
% end


%       %% mixed-model
%         mm_pr = exp(B0+B1*euds(:,k));
%         mm_pr = mm_pr./(1+mm_pr); % logistic probability
%         %mm_pr(flgcensor) = ((1-mm_pr(flgcensor)).*(mm_pr(flgcensor).*prob_no_obs(flgcensor)));
%         mm_pr(flgcensor) = ((1-mm_pr(flgcensor)) + (mm_pr(flgcensor).*prob_no_obs(flgcensor)));
%         mm_pr = log(mm_pr);
%         mm_loglikelihood(k) = sum(mm_pr);