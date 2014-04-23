function [fig, pdx, cox_hr]=fPlotKaplanMeierCurve_DVx(cur_CGobj,volume_val,split_val,num_fx)
        screen_size=get(0,'ScreenSize');
        if ~exist('num_fx','var')
            num_fx=-1;
        end
        
        disp([num2str(volume_val),', ',num2str(split_val)]);
        
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
        if volume_val < 0 %max dose
            fvol_val = volume_val;
        else
            [~,fvol_val] = min(abs(CG.mBinsVol - volume_val));
        end
        f = cellfun(@(x) strcmpi('VDx',x),CG.mLogRank(:,1));
        
        
        dv=zeros(CG.mNumInGrp,1); % dose at volume v
        numstart=CG.mLogRankMinSize;
        pdx=1;
        for v=fvol_val:fvol_val
            % dose under v
            
            if v<0 %max dose
                v_bins = [CG.mGrp.mVolCum];
                d_bins = [CG.mGrp.mDoseBins_LQ];
                
                for j=1:size(v_bins,2) % for each patient
                    vol = v_bins(:,j);
                    dose = d_bins(:,j);
                    
                    vol(~vol)=nan;
                    nan_inds = find(isnan(vol));
                    if ~isempty(nan_inds)
                        min_ind = nan_inds(1)-1;
                    else
                        min_ind = length(vol);
                    end
                    %[~,min_ind] = min(vol);
                    dv(j) = dose(min_ind);
                end
                
            else
                
                dv(:)=0;
                for k=1:CG.mNumInGrp
                    dv(k) = CG.mGrp(k).fDoseAtVol( CG.mBinsVol(v) );
                end
                
            end
            
            g=find(dv>-1); % non-zeros volume cases
            g = intersect(g,find(flgfx'));
            % median volume
            if split_val<0
                dv_split = median(dv);
            else
                dv_split = split_val;
            end
            [~,fdose] = min(abs(CG.mBinsDose - dv_split));
            disp(['Median value of D_{',num2str(volume_val),'} = ',num2str(dv_split)]);
            
            
            % (di,vj)
            numend=length(g)-numstart;
            for d=fdose:fdose
                % check smaple size
                flg_volbelow1=dv(g)<=CG.mBinsDose(d); f=length(find(flg_volbelow1)); % group DVHs by (d,v)
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
                disp(['D_{',num2str(volume_val),'} split at ',num2str(dv_split),...
                    ' has HR: ',num2str(cox_hr)]);
                % compute survival curves and compare them
                sa=sa.fCalculateSurvivalCurve();
                sa=sa.fCombineSurvivalTime();
                sa=sa.fCompareSurvivalByLogrank();
                pdx=sa.mpValue;
            end
        end
        disp(['Log-Rank p-value of D_{',num2str(volume_val),...
                '} split at ',num2str(dv_split),...
                ' = ',num2str(pdx)]);
            
        [~,twoyr_high_val] = min(abs(sa.mSurvivalTimeSorted{1} - 24));
        [~,twoyr_low_val] = min(abs(sa.mSurvivalTimeSorted{2} - 24));
        
        disp(['Actuarial survival rate BELOW median at 2 yrs: ',...
            num2str(1-sa.mSurvivalCurve{1}(twoyr_high_val))]);
         disp(['Actuarial survival rate ABOVE median at 2 yrs: ',...
            num2str(1-sa.mSurvivalCurve{2}(twoyr_low_val))]);
        disp([]);
        fig=figure;  clf reset; hold on; % grid on;
        set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
        sa_times_below=sa.mSurvivalTimeSorted{1};
        sa_curve_below=1-sa.mSurvivalCurve{1};
        sa_times_below(sa_times_below<0)=0;%set negative onset to t=0
        
        sa_times_above=sa.mSurvivalTimeSorted{2};
        sa_curve_above=1-sa.mSurvivalCurve{2};
        sa_times_above(sa_times_above<0)=0;%set negative onset to t=0
        
        h_km(1)=stairs(sa_times_below./12,sa_curve_below,'LineWidth',2.5);
        plot(sa_times_below(sa.mCensorStatistics{1}(:,1))./12,...
            sa_curve_below(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',22,'LineWidth',1.1);
    
        h_km(2)=stairs(sa_times_above./12,sa_curve_above,'r','LineWidth',2.5);
        plot(sa_times_above(sa.mCensorStatistics{2}(:,1))./12,...
            sa_curve_above(sa.mCensorStatistics{2}(:,1)),'r+','MarkerSize',22,'LineWidth',1.1);
        
        str_pval2 = ['Log-Rank p-value = ',num2str(pdx,'%3.1e\n'),10,...
            'HR = ',num2str(cox_hr,'%3.1f')];
        %text(38,0.25,str_pval2,'FontSize',18,'Location',);
        if volume_val<0
%         lgnd=legend(h_km,...
%              strcat('D$_{\rm{max}}\leq',num2str(dv_split,'%3.1f'),'$~Gy$_{10}$'),...
%              strcat('D$_{\rm{max}} >',num2str(dv_split,'%3.1f'),'$~Gy$_{10}$'),...
%             'Location','Best');
        lgnd=legend(h_km,...
             strcat('D$_{\rm{max}}\leq',num2str(dv_split,'%3.1f'),'$~Gy'),...
             strcat('D$_{\rm{max}} >',num2str(dv_split,'%3.1f'),'$~Gy'),...
            'Location','Best');
            
        else
        lgnd=legend(h_km,...
            strcat('D$_{',num2str(volume_val,'%3.1f'),'\rm{cc}}\leq',num2str(dv_split,'%3.1f'),'$~Gy'),...
            strcat('D$_{',num2str(volume_val,'%3.1f'),'\rm{cc}} >',num2str(dv_split,'%3.1f'),'$~Gy'),...
            'Location','Best');
        end
        set(lgnd,'FontSize',22);
        set(lgnd,'Position',[0.6, 0.5, 0.18, 0.01]);
        h_lgnd=legend;
        set(h_lgnd,'interpreter','latex');
        
        textbp(str_pval2,'FontSize',22,'BackgroundColor','w','interpreter','latex');
        
        %tmp
        if num_fx>0
            ylim([0 0.51]);
            xlim([0 70]);
        end
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'FontSize',22);
        xlabel(['Years'],'fontsize',24);
        ylabel(['Probability of Complication'],'fontsize',24);


    end
