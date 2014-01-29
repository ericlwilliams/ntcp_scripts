function [fig, pvx, cox_hr]=fPlotKaplanMeierCurve_VDx(cur_CGobj,dose_val,split_val,num_fx)
        screen_size=get(0,'ScreenSize');
        if ~exist('num_fx','var')
            num_fx=-1;
        end
        
        disp([num2str(dose_val),', ',num2str(split_val)]);
        
        sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
        
        CG = cur_CGobj{1};
        % survival/complication time
        f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
        f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
        compdate = inf(CG.mNumInGrp,1);
        lastfollowup = inf(CG.mNumInGrp,1);
        compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
        lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
        compdate = min( lastfollowup, compdate );
        flgcensor = [CG.mGrp.mFlgCensor]';
        
        flgfx = ones(size(flgcensor));
        if num_fx>=0
            tmp_fx = [CG.mGrp.mFxNum];
           flgfx = tmp_fx==num_fx; 
        end
        
        % V33
        %median split
        [~,fdose_val] = min(abs(CG.mBinsDose - dose_val));
        f = cellfun(@(x) strcmpi('VDx',x),CG.mLogRank(:,1));
        
        
        vd=zeros(CG.mNumInGrp,1); % volume v at dose d
        numstart=CG.mLogRankMinSize;
        pvx=1;
        for d=fdose_val:fdose_val
            % volume under d
            vd(:)=0;
            for k=1:CG.mNumInGrp
                vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(d) );
            end
            g=find(vd>-1); % non-zeros volume cases
            g = intersect(g,find(flgfx'));
            % median volume
            if split_val<0
                vol = median(vd);
            else
                vol = split_val;
            end
            [~,fvol] = min(abs(CG.mBinsVol - vol));
            disp(['Median value of V_{',num2str(dose_val),'} = ',num2str(vol)]);
            
            
            % (di,vj)
            numend=length(g)-numstart;
            for v=fvol:fvol
                % check smaple size
                flg_volbelow1=vd(g)<=CG.mBinsVol(v); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
                if f<numstart || f>numend % one group has too less patients, or the volume at the dose is too small, skip it
                    continue;
                end
                cox_beta=coxphfit(~flg_volbelow1,compdate(g),'baseline',0,'censoring',flgcensor(g));
                cox_hr = exp(cox_beta);
                disp(['Hazard Ratio: ',num2str(cox_hr)]);
                % assign properties of object sa
                survivedate={compdate(g(flg_volbelow1)); compdate(g(~flg_volbelow1))}; % survive time of each group
                fcensor={flgcensor(g(flg_volbelow1)); flgcensor(g(~flg_volbelow1))}; % censor flag for each group
                sa.mSurvivalTime=survivedate;
                sa.mFlgCensor=fcensor;
                sa.mHR = cox_hr;
                disp(['V_{',num2str(dose_val),'} split at ',num2str(split_val),...
                    ' has HR: ',num2str(cox_hr)]);
                % compute survival curves and compare them
                sa=sa.fCalculateSurvivalCurve();
                sa=sa.fCombineSurvivalTime();
                sa=sa.fCompareSurvivalByLogrank();
                pvx=sa.mpValue;
            end
        end
        disp(['Log-Rank p-value of V_{',num2str(dose_val),...
                '} split at ',num2str(vol),...
                ' = ',num2str(pvx)]);
            
        [~,twoyr_val] = min(abs(sa.mSurvivalTimeSorted{1} - 24));
        disp(['Actuarial survival rate below median at 2 yrs: ',...
            num2str(1-sa.mSurvivalCurve{1}(twoyr_val))]);
        
        
        fig=figure;  clf reset; hold on; % grid on;
        set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
        sa_times_below=sa.mSurvivalTimeSorted{1};
        sa_curve_below=1-sa.mSurvivalCurve{1};
        sa_times_below(sa_times_below<0)=0;%set negative onset to t=0
        
        sa_times_above=sa.mSurvivalTimeSorted{2};
        sa_curve_above=1-sa.mSurvivalCurve{2};
        sa_times_above(sa_times_above<0)=0;%set negative onset to t=0
        
        h_km(1)=stairs(sa_times_below./12,sa_curve_below,'LineWidth',2);
        plot(sa_times_below(sa.mCensorStatistics{1}(:,1))./12,...
     sa_curve_below(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',20);
             
    
        h_km(2)=stairs(sa_times_above./12,sa_curve_above,'r','LineWidth',2);
        plot(sa_times_above(sa.mCensorStatistics{2}(:,1))./12,...
            sa_curve_above(sa.mCensorStatistics{2}(:,1)),'r+','MarkerSize',20);
                    
        str_pval2 = ['Log-Rank p-value = ',num2str(pvx,'%3.1e\n'),10,...
            'HR = ',num2str(cox_hr,'%3.1f')];
        %text(38,0.25,str_pval2,'FontSize',18,'Location',);
        
        lgnd=legend(h_km,...
            strcat('V$_{',num2str(dose_val,'%3.1f'),'\rm{Gy}} \leq',num2str(vol,'%3.1f'),'$~cc'),...
            strcat('V$_{',num2str(dose_val,'%3.1f'),'\rm{Gy}} >',num2str(vol,'%3.1f'),'$~cc'),...
            'Location','Best');
        set(lgnd,'FontSize',18);
        h_lgnd=legend;
        set(h_lgnd,'interpreter','latex');
        
        textbp(str_pval2,'FontSize',18);
        
        %tmp
        if num_fx>0
            ylim([0 0.51]);
            xlim([0 70]);
        end
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'FontSize',18);
        xlabel(['Years'],'fontsize',20);
        ylabel(['Probability of Complication'],'fontsize',20);

        %if num_fx<0
        %    title(strcat('V_{',num2str(dose_val),'}, ',num2str(vol,4),' cc threshold'),'fontsize',14);
        %else
        %     title(strcat('V_{',num2str(dose_val),'}, ',...
        %         num2str(vol,4),' cc threshold, N_{Fx} = ',...
        %         num2str(num_fx)),'fontsize',14);
        %end
        
    end
