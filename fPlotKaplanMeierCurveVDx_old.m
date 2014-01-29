function fPlotKaplanMeierCurveVDx_old(cur_CGobj,dose_val,split_val)
        disp([num2str(dose_val),', ',num2str(split_val)]);
        
%         sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
%         
%         CG = cur_CGobj{m};
%         % survival/complication time
%         f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
%         f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
%         compdate = inf(CG.mNumInGrp,1);
%         lastfollowup = inf(CG.mNumInGrp,1);
%         compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
%         lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
%         compdate = min( lastfollowup, compdate );
%         flgcensor = [CG.mGrp.mFlgCensor]';
        
      
        [~,fdose_val] = min(abs(CGobj{m}.mBinsDose - dose_val));
        f = cellfun(@(x) strcmpi('VDx',x),CGobj{m}.mLogRank(:,1));
        pvx = CGobj{m}.mLogRank{f,2};        
        pvx = squeeze(pvx(fdose_val,:,:));
        f = pvx(:,6) < 2;
           
        vd(:)=0;
        for k=1:CG.mNumInGrp
            vd(k) = CG.mGrp(k).fVolAtDose( CG.mBinsDose(fdose_val) );
        end
        if split<0
            vol = median(vd);
        else
            vol = split;
        end
        
        [~,fvol] = min(abs(CGobj{m}.mBinsVol - vol));
        disp(['Median value of V_{',num2str(dose_val),'} = ',num2str(vol)]);
        disp(['Log-Rank p-value of V_{',num2str(dose_val),...
            '} split at ',num2str(vol),...
            ' = ',num2str(pvx(fvol,5))]);
        
        str_pval2 = {strcat('Log-Rank p-value = ',num2str(pvx(fvol,5),'%10.2e\n'))};
        
        sa = CGobj{m}.mLogRank{f,2};
        sa = sa{2}{fdose_val,fvol};
         
        fig=figure;  clf reset; hold on; % grid on;
        set(fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
        h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
        plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
        h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
        plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
            1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
        %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
        text(38,0.25,str_pval2,'FontSize',14);
        lgnd=legend(h_km,...
            strcat('V$_{',num2str(dose_val),'}\leq',num2str(vol,4),'$~cc'),...
            strcat('V$_{',num2str(dose_val),'}\geq',num2str(vol,4),'$~cc'),'Location','Best');
        set(lgnd,'FontSize',14);
        h_lgnd=legend;
        set(h_lgnd,'interpreter','latex');
    
        set(gca,'xminortick','on','yminortick','on');
        xlabel(['Months'],'fontsize',14);
        ylabel(['Probability of CW Pain'],'fontsize',14);
        title(strcat('V_{',num2str(dose_val),'}, ',num2str(vol),' cc threshold'),'fontsize',14);
        
        
    end