% Title: Functional brain connectivity predictors of prospective substance use initiation and their environmental correlates
% Contact: Omid Kardan omidk@med.umich.edu
% PLS analyses scripts for Study 1 results (Fig 1 and supp Figs S1 & S2)

% Requirements: 
% 1) Uses csv file generated from the R script
% "SUI_in_ABCD_compile.R" (keep Matlab custom function "read_df_abcd_sui.m" in the current directory) 

% 2) Uses the PLS scripts from https://www.rotman-baycrest.on.ca/index.php?section=345
% Download plscmd and plsgui and place in Pls_folder

% 3) Uses bluewhitered from Nathan Childress (2022). 
% (https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered),
% MATLAB Central File Exchange.

% 4) Uses export_fig from https://github.com/altmany/export_fig

%%
clear all
% read the csv files generated from the ABCD data using SUI_in_ABCD_compile.R
T_4 =read_df_abcd_sui('~/df_filled_demog_bcog_Y4.csv'); 
T_3 =read_df_abcd_sui('~/df_filled_demog_bcog_Y3.csv');
T_2 =read_df_abcd_sui('~/df_filled_demog_bcog_Y2.csv');
T_1 =read_df_abcd_sui('~/df_filled_demog_bcog_Y1.csv');
T_0 =read_df_abcd_sui('~/df_filled_demog_bcog_Y0.csv');

%% create a variable that includes the list of participants with passed qc and FD thresh rsFC data at both Y0 and Y2
rsFC_avail = zeros(11868,2);
for s = 1:11868
    for c = 1:2
        sub_id = char(T_4.src_subject_id(s));
        if c==1
            for rr=1:4
                ff = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_baselineYear1Arm1_rest_',num2str(rr),'_gordonsubc.mat'];
                if exist(ff)
                    rsFC_avail(s,1) = 1; continue %
                end
            end
        end
        
        if c==2
            for rr=1:4
                ff2 = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_2YearFollowUpYArm1_rest_',num2str(rr),'_gordonsubc.mat'];
                if exist(ff2)
                    rsFC_avail(s,2) = 1; continue%
                end
            end
        end
      
    end
end
with_2b_2024 = find(rsFC_avail(:,1)==1 & rsFC_avail(:,2)==1);

%%  %% Create matched control group for the prospSUI_brain participants
load('~/with_2b_2024.mat')

anyyear_su_ids = find(T_4.anytlfbSU >0 | T_3.anytlfbSU >0 | T_2.anytlfbSU >0 | T_1.anytlfbSU>0 | T_0.anytlfbSU >0 ...
 | T_4.mypiSU >0 | T_3.mypiSU >0 | T_2.mypiSU >0 | T_1.mypiSU >0 | T_0.mypiSU >0 );

year012_su_ids = find( T_2.anytlfbSU >0 | T_1.anytlfbSU>0 | T_0.anytlfbSU >0 ...
 | T_2.mypiSU >0 | T_1.mypiSU >0 | T_0.mypiSU >0 );

intersect(with_2b_2024,year012_su_ids); % subs in the sub-sample that were excluded due to SUI at Y2 or earlier 

prospSUI_ids = find(T_4.anytlfbSU >0 | T_4.mypiSU >0 | ...
    T_3.anytlfbSU >0 | T_3.mypiSU >0); 

y4_unknown_ids = find(isnan(T_4.interview_age)); % partiicpants with unknown Y4 SUI status to be excluded from possible control participants

prospSUI_brain = intersect(setdiff(prospSUI_ids,year012_su_ids),with_2b_2024); % make sure SUI group were negative in year1 and year 2
prospSUI_brain_n = length(prospSUI_brain);


T_2.raceth = T_2.White;
T_2.raceth(T_2.Black == 1) = 2;
T_2.raceth(T_2.Hispanic == 1) = 3;
T_2.raceth(T_2.Asian == 1) = 4;
T_2.raceth(T_2.Other == 1) = 5;
T_2.interview_age(isnan(T_2.interview_age)) = nanmean(T_2.interview_age);
T_2.agey = round(T_2.interview_age/6);
T_2.HighestEd(isnan(T_2.HighestEd)) = nanmean(T_2.HighestEd);
T_2.edu = round(T_2.HighestEd/5);

income_matches = cell(prospSUI_brain_n,1);
age_matches = cell(prospSUI_brain_n,1);
sex_matches = cell(prospSUI_brain_n,1);
raceth_matches = cell(prospSUI_brain_n,1);
edu_matches = cell(prospSUI_brain_n,1);

for k =1:prospSUI_brain_n
    
    if ~isnan(T_2.Income(prospSUI_brain(k)) )
        a = find(round(T_2.Income) == round(T_2.Income(prospSUI_brain(k))));
    else a = 1:11868;
    end
        income_matches{k} = intersect(with_2b_2024,a);
        
    if ~isnan(T_2.agey(prospSUI_brain(k)))
        a = find(T_2.agey == T_2.agey(prospSUI_brain(k)));
    else a = 1:11868
    end
        age_matches{k} = intersect(with_2b_2024,a);
    
    if ~isnan(T_2.Male_bin(prospSUI_brain(k)))
        a = find(T_2.Male_bin == T_2.Male_bin(prospSUI_brain(k)));
    else a = 1:11868;
    end
         sex_matches{k} = intersect(with_2b_2024,a);
         
    if ~isnan(T_2.raceth(prospSUI_brain(k)))
        a = find(T_2.raceth == T_2.raceth(prospSUI_brain(k)));
    else a = 1:11868;
    end
         raceth_matches{k} = intersect(with_2b_2024,a);
    
    if ~isnan(T_2.edu(prospSUI_brain(k)))
        a = find(T_2.edu == T_2.edu(prospSUI_brain(k)));
    else a = 1:11868;
    end
         edu_matches{k} = intersect(with_2b_2024,a);
            
end
cont_mj = cell(prospSUI_brain_n,1); cont_mj_full = cell(prospSUI_brain_n,1); cont_mj_bare = cell(prospSUI_brain_n,1);
for k = 1:prospSUI_brain_n
    gg = intersect( intersect(sex_matches{k},age_matches{k}), intersect(raceth_matches{k},edu_matches{k}));
    cont_mj_full{k} = setdiff(intersect( gg, income_matches{k} ),  union( anyyear_su_ids, y4_unknown_ids));
    cont_mj{k} = setdiff( gg,  union( anyyear_su_ids, y4_unknown_ids));
    gg2 = intersect( intersect(sex_matches{k},age_matches{k}), raceth_matches{k});
    cont_mj_bare{k} = setdiff( gg2,  union( anyyear_su_ids, y4_unknown_ids));
    
    if k == 17  % this one participant had no match so loosened the criteria for them to only age & sex & race
       gg3 = intersect(intersect(find(T_2.raceth == T_2.raceth(prospSUI_brain(k))),...
           find(T_2.agey - T_2.agey(prospSUI_brain(k)) < 1)), sex_matches{k});
       
      cont_mj_bare{k} = setdiff( intersect(with_2b_2024,gg3),  union( anyyear_su_ids, y4_unknown_ids));
    end
end

HCss = []; llls =[];
for ww = 1:10000
    HCs =NaN(prospSUI_brain_n,1);
    for k=1:prospSUI_brain_n
        if ~isempty(cont_mj_full{k})
            if length(cont_mj_full{k}) >1
                for jjk = 1: length(cont_mj_full{k})
                    dd = randperm(length(cont_mj_full{k}),1);
                    if ~ismember(cont_mj_full{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_full{k}(dd);
            
            
        elseif ~isempty(cont_mj{k})
            if length(cont_mj{k}) >1
                for jjk = 1: length(cont_mj{k})
                    dd = randperm(length(cont_mj{k}),1);
                    if ~ismember(cont_mj{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj{k}(dd);
            
            
        elseif ~isempty(cont_mj_bare{k})
            if length(cont_mj_bare{k}) >1
                for jjk = 1: length(cont_mj_bare{k})
                    dd = randperm(length(cont_mj_bare{k}),1);
                    if ~ismember(cont_mj_bare{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_bare{k}(dd);
        end
    end
    llls = [llls; length(unique(HCs))];
        
        HCss = [HCss  HCs];

end

 matches = HCss(:,1);
for j=1:10000
    [v, w] = unique( matches, 'stable' );
    duplicate_indices = setdiff( 1:numel(matches), w );
    if length(duplicate_indices)==0
        break
    end
    matches(duplicate_indices) = HCss(duplicate_indices,j);
end
        
matches_228 = unique(matches); % save these as your matched controls
SUI_233 = prospSUI_brain; % save these as your SUI group

%% Study 1 main PLS (Figure 1) includes all available fMRI runs per sub for all subs regardless of family

clear all

load('matches_228.mat')
load('SUI_233.mat')

% read the csv files generated from the ABCD data using SUI_in_ABCD_compile.R
T_4 =read_df_abcd_sui('~/df_filled_demog_bcog_Y4.csv'); 
T_3 =read_df_abcd_sui('~/df_filled_demog_bcog_Y3.csv');
T_2 =read_df_abcd_sui('~/df_filled_demog_bcog_Y2.csv');
T_1 =read_df_abcd_sui('~/df_filled_demog_bcog_Y1.csv');
T_0 =read_df_abcd_sui('~/df_filled_demog_bcog_Y0.csv');


addpath(genpath('~\Pls')); % folder containing the PLS scripts from https://www.rotman-baycrest.on.ca/index.php?section=345
rng('default');
groups{1} = [matches_228];
groups{2} = [SUI_233];
nTimes = 2;
nParcels = 333+54+31; badParcels=[]; goodps = setdiff(1:nParcels,badParcels); 

inds = find(tril(ones(length(goodps)),-1)==1);
datamat_lst = cell(length(groups),1); 

for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_4.src_subject_id(groups{g}(s)))
            
            vec =[]; vec2 =[];
            if c==1
                for rr=1:4
                    ff = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_baselineYear1Arm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff)
                        load(ff); vec1 = full_conn(inds)';  %
                        vec2 = [vec2; vec1];
                    end
                end
                
            end
            
            if c==2
                for rr=1:4
                    ff2 = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_2YearFollowUpYArm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff2)
                        load(ff2); vec1 = full_conn(inds)';
                        vec2 = [vec2; vec1];%
                    end
                end
            end
            if isempty(vec2)
                vec = NaN(1,55278);
            end
            if ~isempty(vec2)
                vec = nanmean(vec2,1);
            end
            datamat_lst{g} = [datamat_lst{g}; vec];
            s
            
        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nTimes;
option.method = 1; 
option.num_boot = 1000;
option.num_perm = 1000;
option.meancentering_type=[2]; % grand mean removed only (main analysis)
result = pls_analysis(datamat_lst, num_subj, num_cond, option);
result.datamat_lst = datamat_lst;
save('PLS_result_233any_228hc_rsFCs_y0_y2_multrun.mat','result');


% plotting the results
% Uses bluewhitered from Nathan Childress (2022) (https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered)
clear all
addpath(genpath('~\bluewhitered'))
% export_fig from https://github.com/altmany/export_fig
addpath(genpath('~\export_fig-master'));

load('PLS_result_233any_228hc_rsFCs_y0_y2_multrun.mat'); sucol = [.87 .77 .1]; ss = 'SU+';

ps = result.perm_result.sprob
norm_dist = (result.boot_result.distrib )/ ...
    (max(max(max(result.boot_result.distrib))) - min(min(min(result.boot_result.distrib))));
    
lv=1;
    
   gt_contrast = squeeze(norm_dist(:,lv,:))';
 gt_contrast2 = [gt_contrast(:,1), gt_contrast(:,3), gt_contrast(:,2), gt_contrast(:,4)];
 
brainPs = find(result.boot_result.compare_u(:,lv)>=3);
brainNs = find(result.boot_result.compare_u(:,lv)<=-3);
hcn = result.num_subj_lst(1); % length(groups{1})
sun = result.num_subj_lst(2); % length(groups{2})
beg = [1,hcn+1,2*hcn+1,2*hcn+sun+1];
en = [hcn, 2*hcn, 2*hcn+sun, 2*hcn+2*sun];
figure;
hold on
for k=1:4
if mod(k,2)==1; col = [.35 .65 0.35];  % color for columns 1 and 3 (HC)
else
    col = sucol;
end
   bar(k,mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(k,mean(gt_contrast2(:,k)),2*std(gt_contrast2(:,k)),'.k');
end
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'9-10 yrs old','9-10 yrs old','11-12 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS Weight'); 
hold off

% left panel in Figure 1
% export_fig -m5 -transparent rsFC_SU233multrun_abcd_pls_lv1_PLSweightbars.jpg

figure
subplot(1,2,1);
hold on
for k=1:4
if k<3 col = [.1 .2 .3];
else
    col = sucol;
end
   bar(k,mean(gt_contrast(:,k)),'FaceColor',col,'EdgeColor','k');
   errorbar(k,mean(gt_contrast(:,k)),2*std(gt_contrast(:,k)),'.k');
end
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'9-10 yrs old','11-12 yrs old','9-10 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',16);
ylabel('PLS Weight'); axis square
hold off  

    
    cbcovs = result.s.^2./sum(result.s.^2);
    title(['LV',num2str(lv),...
        '  \sigma_{XY} = ',num2str(round(cbcovs(lv),3)),...
        '  p = ',num2str(round(ps(lv),3))]);


subplot(1,2,2);
goodps = 1:418;
inds = find(tril(ones(length(goodps)),-1)==1);
WBcompare_u_1 = result.boot_result.compare_u(:,lv);
temp11=zeros(418,418);temp11(inds)=WBcompare_u_1;
temp1z = temp11+temp11'; 

temp1 = temp1z; 
temp1(abs(temp1z) <=3)=0;
temp1(temp1z >3)=1;
temp1(temp1z <-3)=-1;
imagesc(temp1,[-.4,.4]); hold on

colormap(bluewhitered), colorbar;

xnames = {'Auditory','Cingulo-Opercular', 'Cingulo-Parietal','Default', 'Dorsal Attn', 'Fronto-Parietal',...
    'Orbito & Temp Pole', 'Retro-Temporal', 'Soma-Motor-h','Soma-Motor-m' ,...
    ' ','Ventral Attn','Visual', 'Subcortical'  'Cerebellar'};
netsizes1 = [0 24    40     5    41    32    24    47     8    38     8     4    23  39 54];
netsizes = [24    40     5    41    32    24    47     8    38     8     4    23    39  54  31];
set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',16);
set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',16);

makans = cumsum(netsizes1);

xcords = []; ycords = []; cccs = [];
for jj=1:15
    for j=1:length(netsizes)
        ycord = [makans(j);  makans(j);   makans(j)+netsizes(j);  makans(j)+netsizes(j)];
        xcord = [makans(jj);  makans(jj)+netsizes(jj);     makans(jj)+netsizes(jj);      makans(jj)];
        ccc = mean(mean(temp1(makans(jj)+1:makans(jj)+netsizes(jj), makans(j)+1:makans(j)+netsizes(j))));
        xcords = [xcords  xcord];
        ycords = [ycords  ycord];
        cccs = [cccs; ccc];
        
    end
end
ff = fill(xcords,ycords,cccs); colormap(bluewhitered), colorbar;
for i=1:225
    ff(i,1).EdgeColor = [1,1,1];
end
fill([0;418;418],[0;0;418],'white');

makansb = cumsum(netsizes);
for j=1:length(netsizes)
line([makansb(j),makansb(j)],[makansb(j) ,418],'Color','black');
end
for j=1:length(netsizes)
line([0 ,makansb(j)],[makansb(j),makansb(j)],'Color','black');
end

line([0 ,418],[0 ,418],'Color','black');
ylabel('Contribution to the Latent Variable (prop |Z|>3)','Position',[470,214,1]); axis square
[r ~ ] = corr(result.usc(:,1),result.vsc(:,1))


% right panle in Figure 1:
% export_fig -m5 -transparent rsFC_SU233multrun_abcd_pls_lv1_PLSpropmap.jpg


%% Applying the brain saliences to the study 2 participants
clear all
% read the csv files generated from the ABCD data using SUI_in_ABCD_compile.R
T_4 =read_df_abcd_sui('~/df_filled_demog_bcog_Y4.csv'); 
T_3 =read_df_abcd_sui('~/df_filled_demog_bcog_Y3.csv');
T_2 =read_df_abcd_sui('~/df_filled_demog_bcog_Y2.csv');
T_1 =read_df_abcd_sui('~/df_filled_demog_bcog_Y1.csv');
T_0 =read_df_abcd_sui('~/df_filled_demog_bcog_Y0.csv');

load('PLS_result_233any_228hc_rsFCs_y0_y2_multrun.mat'); % can load the single-run or the no-family results here instead from supp sections of this script
U = result.u(:,1);
V = result.v(:,1);

anyyear_su_ids = find(T_4.anytlfbSU >0 | T_3.anytlfbSU >0 | T_2.anytlfbSU >0 | T_1.anytlfbSU>0 | T_0.anytlfbSU >0 ...
 | T_4.mypiSU >0 | T_3.mypiSU >0 | T_2.mypiSU >0 | T_1.mypiSU >0 | T_0.mypiSU >0 );


load('with_2b_2024.mat')
tab_fd1 = readtable('~\ABCD_stats_postFD.csv'); % frame displacement (FD) values for each run of each partiicpant is in ABCD_stats_postFD.csv 
tab_fd2 = readtable('~\ABCD_stats_postFD.csv');

nTimes = 2;
nParcels = 333+54+31; badParcels=[]; goodps = setdiff(1:nParcels,badParcels); 
inds = find(tril(ones(length(goodps)),-1)==1);

load('SUI_233.mat')
load('matches_228.mat')

Us =[];  subs =[];  V_diff_y0 =[]; V_diff_y2 =[];
for s = 1:11868
    sub_id = char(T_4.src_subject_id(s));
    
    if ismember(s,with_2b_2024)
        used_id = 0;
        if ismember(s,matches_228')
            used_id = 1;
        end
        if ismember(s,SUI_233')
            used_id = 2;
        end
        for c = 1:nTimes   
            
            if c==1
                vec11 =[]; pfd11 =[]; fd11 =[];
                for rr=1:4
                    ff = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_baselineYear1Arm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff)
                        load(ff); vec1 = U'*full_conn(inds); 
 pfd1 = tab_fd1.postmeanFD( string(tab_fd1.eventname) == 'baselineYear1Arm1' & string(tab_fd1.subid) == sub_id &  tab_fd1.run == rr);
 vec11 = [vec11; vec1];  pfd11 = [pfd11; pfd1];                       
 
 
                    end
                end
            end
            
            if c==2
                vec22 =[]; pfd22 =[]; fd22 =[];
                for rr=1:4
                    ff2 = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_2YearFollowUpYArm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff2)
                        load(ff2); vec2 = U'*full_conn(inds);
pfd2 = tab_fd2.postmeanFD( string(tab_fd2.eventname) == '2YearFollowUpYArm1' & string(tab_fd2.subid) == sub_id  & tab_fd2.run == rr);  
vec22 = [vec22; vec2];  pfd22 = [pfd22; pfd2];  

                    end
                end
            end
            
        end
        
    else
        vec11 = NaN; vec22 = NaN; fd11 = NaN; fd22 = NaN;  pfd11 = NaN; pfd22 = NaN; used_id = NaN;
    end
    Us = [Us; [nanmean(vec11,1)  nanmean(vec22,1)   nanmean(pfd11,1)  nanmean(pfd22,1)  used_id]];
    subs = [subs; string(sub_id)];
    V_diff_y0 = [V_diff_y0; V(3)-V(1)];
    V_diff_y2 = [V_diff_y2; V(4)-V(2)];
    
    s
end
any_year = zeros(11868,1);
any_year(anyyear_su_ids) = 1;
Tble = table(subs,Us(:,5),any_year,V_diff_y0,V_diff_y2,Us(:,1),Us(:,2),Us(:,3),Us(:,4),...
    'VariableNames',{'subid','used_ids','any_year_SU','V_diff_y0','V_diff_y2','U_y0','U_y2','pFD_y0','pFD_y2'});

% save the brain saliences (Us) and V_diffs as a csv for the study_2 regressions R script
writetable(Tble,'abcd_pls_rsFC_SU_multrun_scores_forstudy2.csv');


%% %%%%%%%%%%%%%%%%%% posthoc PLS with temporal component removed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all

load('matches_228.mat')
load('SUI_233.mat')

T_4 =read_df_abcd_sui('~/df_filled_demog_bcog_Y4.csv');

addpath(genpath('E:\Omid\nube\October\Pls'));
rng('default');
groups{1} = [matches_228];
groups{2} = [SUI_233];
nTimes = 2;
nParcels = 333+54+31; badParcels=[]; goodps = setdiff(1:nParcels,badParcels); 

inds = find(tril(ones(length(goodps)),-1)==1);
datamat_lst = cell(length(groups),1); 
% ####*****************####
for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_4.src_subject_id(groups{g}(s)))
            
            vec =[]; vec2 =[];
            if c==1
                for rr=1:4
                    ff = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_baselineYear1Arm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff)
                        load(ff); vec1 = full_conn(inds)';  %
                        vec2 = [vec2; vec1];
                    end
                end
                
            end
            
            if c==2
                for rr=1:4
                    ff2 = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_2YearFollowUpYArm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff2)
                        load(ff2); vec1 = full_conn(inds)';
                        vec2 = [vec2; vec1];%
                    end
                end
            end
            if isempty(vec2)
                vec = NaN(1,55278);
            end
            if ~isempty(vec2)
                vec = nanmean(vec2,1);
            end
            datamat_lst{g} = [datamat_lst{g}; vec];
            s
            
        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nTimes;
option.method = 1; 
option.num_boot = 1000;
option.num_perm = 1000;
option.meancentering_type=[1]; % 1 remove overall condition diffrences
result = pls_analysis(datamat_lst, num_subj, num_cond, option);

result.datamat_lst = datamat_lst;
save('PLS_result_233any_228hc_rsFCs_y0_y2_multrun_mc1.mat','result');


%% %%%%%%%%%%%%%%%%%%%%%%%% Sensitivity analysis PLS (Supp 1): equal one fMRI run per sub %%%%%%%%%%%%%%%%%
clear all
T_4 =read_df_abcd_sui('~/df_filled_demog_bcog_Y4.csv');

load('matches_228.mat')
load('SUI_233.mat')

addpath(genpath('~\Pls'));
rng('default');
groups{1} = [matches_228];
groups{2} = [SUI_233];
nRuns = 2;Runs=[1,2];
nParcels = 333+54+31; badParcels=[]; goodps = setdiff(1:nParcels,badParcels); 

inds = find(tril(ones(length(goodps)),-1)==1);
datamat_lst = cell(length(groups),1); 
% ####*****************####
for g = 1:numel(groups)
    for c = 1:nRuns
        for s = 1:numel(groups{g})         
          sub_id = char(T_4.src_subject_id(groups{g}(s)));  
          
          vec =[];
            if c==1  
                for rr=1:4
                 ff = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_baselineYear1Arm1_rest_',num2str(rr),'_gordonsubc.mat'];
                 if exist(ff)
                 load(ff); vec = full_conn(inds)'; continue % 
                 end
                end
            end
            
            if c==2  
                for rr=1:4
                ff2 = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_2YearFollowUpYArm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff2)
                load(ff2); vec = full_conn(inds)'; continue%
                    end
                end
            end 
            if isempty(vec)
           vec = NaN(1,55278);
            end
            datamat_lst{g} = [datamat_lst{g}; vec];
            s

        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nRuns;
option.method = 1; 
option.num_boot = 1000;
option.num_perm = 1000;
option.meancentering_type=[2]; 
result = pls_analysis(datamat_lst, num_subj, num_cond, option);

result.datamat_lst = datamat_lst;
save('PLS_result_233any_228hc_rsFCs_y0_y2_singlerun.mat','result');

% plotting the results
clear all
addpath(genpath('~\bluewhitered'))
load('PLS_result_233any_228hc_rsFCs_y0_y2_singlerun.mat'); sucol = [.87 .77 .1]; ss = 'SU+';

ps = result.perm_result.sprob
norm_dist = (result.boot_result.distrib )/ ...
    (max(max(max(result.boot_result.distrib))) - min(min(min(result.boot_result.distrib))));
    
lv=1;
    
   gt_contrast = squeeze(norm_dist(:,lv,:))';
 gt_contrast2 = [gt_contrast(:,1), gt_contrast(:,3), gt_contrast(:,2), gt_contrast(:,4)];
 
brainPs = find(result.boot_result.compare_u(:,lv)>=3);
brainNs = find(result.boot_result.compare_u(:,lv)<=-3);
hcn = result.num_subj_lst(1); % length(groups{1})
sun = result.num_subj_lst(2); % length(groups{2})
beg = [1,hcn+1,2*hcn+1,2*hcn+sun+1];
en = [hcn, 2*hcn, 2*hcn+sun, 2*hcn+2*sun];
figure;
hold on
for k=1:4
if mod(k,2)==1; col = [.35 .65 0.35];  % color for columns 1 and 3 (HC)
else
    col = sucol;
end
   bar(k,mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(k,mean(gt_contrast2(:,k)),2*std(gt_contrast2(:,k)),'.k');
end
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'9-10 yrs old','9-10 yrs old','11-12 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS Weight (au)'); 
hold off

% left panel in Figure S1
% export_fig -m5 -transparent rsFC_SU233singlerun_abcd_pls_lv1_PLSweightbars.jpg


figure
subplot(1,2,1);
hold on
for k=1:4
if k<3 col = [.1 .2 .3];
else
    col = sucol;
end
   bar(k,mean(gt_contrast(:,k)),'FaceColor',col,'EdgeColor','k');
   errorbar(k,mean(gt_contrast(:,k)),2*std(gt_contrast(:,k)),'.k');
end
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'9-10 yrs old','11-12 yrs old','9-10 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',16);
ylabel('PLS Weight (au)'); axis square
hold off  

    
    cbcovs = result.s.^2./sum(result.s.^2);
    title(['LV',num2str(lv),...
        '  \sigma_{XY} = ',num2str(round(cbcovs(lv),3)),...
        '  p = ',num2str(round(ps(lv),3))]);

subplot(1,2,2);
goodps = 1:418;
inds = find(tril(ones(length(goodps)),-1)==1);
WBcompare_u_1 = result.boot_result.compare_u(:,lv);
temp11=zeros(418,418);temp11(inds)=WBcompare_u_1;
temp1z = temp11+temp11'; 

temp1 = temp1z; 
temp1(abs(temp1z) <=3)=0;
temp1(temp1z >3)=1;
temp1(temp1z <-3)=-1;
imagesc(temp1,[-.2,.2]); hold on

colormap(bluewhitered), colorbar;

xnames = {'Auditory','Cingulo-Opercular', 'Cingulo-Parietal','Default', 'Dorsal Attn', 'Fronto-Parietal',...
    'Orbito & Temp Pole', 'Retro-Temporal', 'Soma-Motor-h','Soma-Motor-m' ,...
    ' ','Ventral Attn','Visual', 'Subcortical'  'Cerebellar'};
netsizes1 = [0 24    40     5    41    32    24    47     8    38     8     4    23  39 54];
netsizes = [24    40     5    41    32    24    47     8    38     8     4    23    39  54  31];
set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',16);
set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',16);

makans = cumsum(netsizes1);

xcords = []; ycords = []; cccs = [];
for jj=1:15
    for j=1:length(netsizes)
        ycord = [makans(j);  makans(j);            makans(j)+netsizes(j);  makans(j)+netsizes(j)];
        xcord = [makans(jj);  makans(jj)+netsizes(jj);     makans(jj)+netsizes(jj);      makans(jj)];
        ccc = mean(mean(temp1(makans(jj)+1:makans(jj)+netsizes(jj), makans(j)+1:makans(j)+netsizes(j))));
        xcords = [xcords  xcord];
        ycords = [ycords  ycord];
        cccs = [cccs; ccc];
        
    end
end
ff = fill(xcords,ycords,cccs); colormap(bluewhitered), colorbar;
for i=1:225
    ff(i,1).EdgeColor = [1,1,1];
end
fill([0;418;418],[0;0;418],'white');

makansb = cumsum(netsizes);
for j=1:length(netsizes)
line([makansb(j),makansb(j)],[makansb(j) ,418],'Color','black');
end
for j=1:length(netsizes)
line([0 ,makansb(j)],[makansb(j),makansb(j)],'Color','black');
end

line([0 ,418],[0 ,418],'Color','black');
ylabel('Contribution to the Latent Variable (prop |Z|>3)','Position',[470,214,1]); axis square
[r p ] = corr(result.usc(:,1),result.vsc(:,1))

% right panle in Figure S1:
% export_fig -m5 -transparent rsFC_SU233singlerun_abcd_pls_lv1_PLSpropmap.jpg


%%  PLS sensitivity analysis with one family member randomly kept (Supp 2)

clear all
% read the csv files generated from the ABCD data using SUI_in_ABCD_compile.R
T_4 =read_df_abcd_sui('~/df_filled_demog_bcog_Y4.csv'); 
T_3 =read_df_abcd_sui('~/df_filled_demog_bcog_Y3.csv');
T_2 =read_df_abcd_sui('~/df_filled_demog_bcog_Y2.csv');
T_1 =read_df_abcd_sui('~/df_filled_demog_bcog_Y1.csv');
T_0 =read_df_abcd_sui('~/df_filled_demog_bcog_Y0.csv');

% create new SUI and matched control groups with the no-family sample
famid = T_0.rel_family_id;
    [v, w] = unique( famid, 'stable' );
    famid_ex = setdiff( 1:numel(famid), w );

rsFC_avail = zeros(11868,2);
for s = 1:11868
    if ismember(s,famid_ex)
        continue;
    end
    for c = 1:2
        sub_id = char(T_4.src_subject_id(s));
        if c==1
            for rr=1:4
                ff = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_baselineYear1Arm1_rest_',num2str(rr),'_gordonsubc.mat'];
                if exist(ff)
                    rsFC_avail(s,1) = 1; continue %
                end
            end
        end
        
        if c==2
            for rr=1:4
                ff2 = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_2YearFollowUpYArm1_rest_',num2str(rr),'_gordonsubc.mat'];
                if exist(ff2)
                    rsFC_avail(s,2) = 1; continue%
                end
            end
        end
        
        if round(s/100) == s/100
            s
        end        
    end
end
with_2b_2024_nofam = find(rsFC_avail(:,1)==1 & rsFC_avail(:,2)==1);
%
load('with_2b_2024_nofam.mat')

anyyear_su_ids = find(T_4.anytlfbSU >0 | T_3.anytlfbSU >0 | T_2.anytlfbSU >0 | T_1.anytlfbSU>0 | T_0.anytlfbSU >0 ...
 | T_4.mypiSU >0 | T_3.mypiSU >0 | T_2.mypiSU >0 | T_1.mypiSU >0 | T_0.mypiSU >0 );

year012_su_ids = find( T_2.anytlfbSU >0 | T_1.anytlfbSU>0 | T_0.anytlfbSU >0 ...
 | T_2.mypiSU >0 | T_1.mypiSU >0 | T_0.mypiSU >0 );

prospSUI_ids = find(T_4.anytlfbSU >0 | T_4.mypiSU >0 | ...
    T_3.anytlfbSU >0 | T_3.mypiSU >0); mj_n = length(prospSUI_ids);

y4_unknown_ids = find(isnan(T_4.interview_age));

prospSUI_brain = intersect(setdiff(prospSUI_ids,year012_su_ids),with_2b_2024_nofam); % make sure SUI were negative in year1 and year 2
prospSUI_brain_n = length(prospSUI_brain);
%
T_2.raceth = T_2.White;
T_2.raceth(T_2.Black == 1) = 2;
T_2.raceth(T_2.Hispanic == 1) = 3;
T_2.raceth(T_2.Asian == 1) = 4;
T_2.raceth(T_2.Other == 1) = 5;
T_2.interview_age(isnan(T_2.interview_age)) = nanmean(T_2.interview_age);
T_2.agey = round(T_2.interview_age/6);
T_2.HighestEd(isnan(T_2.HighestEd)) = nanmean(T_2.HighestEd);
T_2.edu = round(T_2.HighestEd/5);

income_matches = cell(prospSUI_brain_n,1);
age_matches = cell(prospSUI_brain_n,1);
sex_matches = cell(prospSUI_brain_n,1);
raceth_matches = cell(prospSUI_brain_n,1);
edu_matches = cell(prospSUI_brain_n,1);

for k =1:prospSUI_brain_n
    
    if ~isnan(T_2.Income(prospSUI_brain(k)) )
        a = find(round(T_2.Income) == round(T_2.Income(prospSUI_brain(k))));
    else a = 1:11868;
    end
        income_matches{k} = intersect(with_2b_2024_nofam,a);
        
    if ~isnan(T_2.agey(prospSUI_brain(k)))
        a = find(T_2.agey == T_2.agey(prospSUI_brain(k)));
    else a = 1:11868
    end
            age_matches{k} = intersect(with_2b_2024_nofam,a);
    
    if ~isnan(T_2.Male_bin(prospSUI_brain(k)))
        a = find(T_2.Male_bin == T_2.Male_bin(prospSUI_brain(k)));
    else a = 1:11868;
    end
         sex_matches{k} = intersect(with_2b_2024_nofam,a);
         
    if ~isnan(T_2.raceth(prospSUI_brain(k)))
        a = find(T_2.raceth == T_2.raceth(prospSUI_brain(k)));
    else a = 1:11868;
    end
            raceth_matches{k} = intersect(with_2b_2024_nofam,a);
    
    if ~isnan(T_2.edu(prospSUI_brain(k)))
        a = find(T_2.edu == T_2.edu(prospSUI_brain(k)));
    else a = 1:11868;
    end
            edu_matches{k} = intersect(with_2b_2024_nofam,a);
            
end
cont_mj = cell(prospSUI_brain_n,1); cont_mj_full = cell(prospSUI_brain_n,1); cont_mj_bare = cell(prospSUI_brain_n,1);
for k = 1:prospSUI_brain_n
    gg = intersect( intersect(sex_matches{k},age_matches{k}), intersect(raceth_matches{k},edu_matches{k}));
    cont_mj_full{k} = setdiff(intersect( gg, income_matches{k} ),  union( anyyear_su_ids, y4_unknown_ids));
    cont_mj{k} = setdiff( gg,  union( anyyear_su_ids, y4_unknown_ids));
    gg2 = intersect( intersect(sex_matches{k},age_matches{k}), raceth_matches{k});
    cont_mj_bare{k} = setdiff( gg2,  union( anyyear_su_ids, y4_unknown_ids));
    
    if k == 17  % this partiicpant again had no matches so loosened the criteria
       gg3 = intersect(intersect(find(T_2.raceth == T_2.raceth(prospSUI_brain(k))),...
           find(T_2.agey - T_2.agey(prospSUI_brain(k)) < 1)), sex_matches{k});
       
      cont_mj_bare{k} = setdiff( intersect(with_2b_2024_nofam,gg3),  union( anyyear_su_ids, y4_unknown_ids));
    end
end

HCss = []; llls =[];
for ww = 1:10000
    HCs =NaN(prospSUI_brain_n,1);
    for k=1:prospSUI_brain_n
        if ~isempty(cont_mj_full{k})
            if length(cont_mj_full{k}) >1
                for jjk = 1: length(cont_mj_full{k})
                    dd = randperm(length(cont_mj_full{k}),1);
                    if ~ismember(cont_mj_full{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_full{k}(dd);
            
            
        elseif ~isempty(cont_mj{k})
            if length(cont_mj{k}) >1
                for jjk = 1: length(cont_mj{k})
                    dd = randperm(length(cont_mj{k}),1);
                    if ~ismember(cont_mj{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj{k}(dd);
            
            
        elseif ~isempty(cont_mj_bare{k})
            if length(cont_mj_bare{k}) >1
                for jjk = 1: length(cont_mj_bare{k})
                    dd = randperm(length(cont_mj_bare{k}),1);
                    if ~ismember(cont_mj_bare{k}(dd),HCs)
                        continue
                    end
                end
            else dd =1;
            end
            HCs(k) = cont_mj_bare{k}(dd);
        end
    end
    llls = [llls; length(unique(HCs))];
        
        HCss = [HCss  HCs];

end
 matches = HCss(:,1);
for j=1:10000
    [v, w] = unique( matches, 'stable' );
    duplicate_indices = setdiff( 1:numel(matches), w );
    if length(duplicate_indices)==0
        break
    end
    matches(duplicate_indices) = HCss(duplicate_indices,j);
end
       
matches_180_nofam = unique(matches);
SUI_187_nofam = prospSUI_brain;
% save the matches and SUI groups for the nofam version

% The PLS analysis with the new SUI and control groups from the no-family
% smaple
clear all
T_4 =read_df_abcd_sui('~/df_filled_demog_bcog_Y4.csv');

load('matches_180_nofam.mat')
load('SUI_187_nofam.mat')

addpath(genpath('~\Pls'));
rng('default');
groups{1} = [matches_180_nofam];
groups{2} = [SUI_187_nofam];
nTimes = 2;
nParcels = 333+54+31; badParcels=[]; goodps = setdiff(1:nParcels,badParcels); 

inds = find(tril(ones(length(goodps)),-1)==1);
datamat_lst = cell(length(groups),1); 

for g = 1:numel(groups)
    for c = 1:nTimes
        for s = 1:numel(groups{g})
            sub_id = char(T_4.src_subject_id(groups{g}(s)))
            
            vec =[]; vec2 =[];
            if c==1
                for rr=1:4
                    ff = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_baselineYear1Arm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff)
                        load(ff); vec1 = full_conn(inds)';  %
                        vec2 = [vec2; vec1];
                    end
                end
                
            end
            
            if c==2
                for rr=1:4
                    ff2 = ['rsFC_data_directory\NDARINV',sub_id(9:16),'_2YearFollowUpYArm1_rest_',num2str(rr),'_gordonsubc.mat'];
                    if exist(ff2)
                        load(ff2); vec1 = full_conn(inds)';
                        vec2 = [vec2; vec1];%
                    end
                end
            end
            if isempty(vec2)
                vec = NaN(1,55278);
            end
            if ~isempty(vec2)
                vec = nanmean(vec2,1);
            end
            datamat_lst{g} = [datamat_lst{g}; vec];
            s
            
        end
    end
end

num_subj = [length(groups{1})  length(groups{2})];
num_cond = nTimes;
option.method = 1; 
option.num_boot = 1000;
option.num_perm = 1000;
option.meancentering_type=[2]; 

result = pls_analysis(datamat_lst, num_subj, num_cond, option);
result.datamat_lst = datamat_lst;

save('PLS_result_187any_180hc_rsFCs_y0_y2_nofam.mat','result');

% plotting the results
clear all
addpath(genpath('~\bluewhitered'))
load('PLS_result_187any_180hc_rsFCs_y0_y2_nofam.mat'); sucol = [.87 .77 .1]; ss = 'SU+';

ps = result.perm_result.sprob
norm_dist = (result.boot_result.distrib )/ ...
    (max(max(max(result.boot_result.distrib))) - min(min(min(result.boot_result.distrib))));
    
lv=1;
    
   gt_contrast = squeeze(norm_dist(:,lv,:))';
 gt_contrast2 = [gt_contrast(:,1), gt_contrast(:,3), gt_contrast(:,2), gt_contrast(:,4)];
 
brainPs = find(result.boot_result.compare_u(:,lv)>=3);
brainNs = find(result.boot_result.compare_u(:,lv)<=-3);
hcn = result.num_subj_lst(1); % length(groups{1})
sun = result.num_subj_lst(2); % length(groups{2})
beg = [1,hcn+1,2*hcn+1,2*hcn+sun+1];
en = [hcn, 2*hcn, 2*hcn+sun, 2*hcn+2*sun];
figure;
hold on
for k=1:4
if mod(k,2)==1; col = [.35 .65 0.35];  % color for columns 1 and 3 (HC)
else
    col = sucol;
end
   bar(k,mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(k,mean(gt_contrast2(:,k)),2*std(gt_contrast2(:,k)),'.k');
end
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'9-10 yrs old','9-10 yrs old','11-12 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS Weight (au)'); 
hold off

% left panel in Figure S2
% export_fig -m5 -transparent rsFC_SU187nofam_abcd_pls_lv1_PLSweightbars.jpg


figure
subplot(1,2,1);
hold on
for k=1:4
if k<3 col = [.1 .2 .3];
else
    col = sucol;
end
   bar(k,mean(gt_contrast(:,k)),'FaceColor',col,'EdgeColor','k');
   errorbar(k,mean(gt_contrast(:,k)),2*std(gt_contrast(:,k)),'.k');
end
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'9-10 yrs old','11-12 yrs old','9-10 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',16);
ylabel('PLS Weight (au)'); axis square
hold off  

    
    cbcovs = result.s.^2./sum(result.s.^2);
    title(['LV',num2str(lv),...
        '  \sigma_{XY} = ',num2str(round(cbcovs(lv),3)),...
        '  p = ',num2str(round(ps(lv),3))]);


subplot(1,2,2);
goodps = 1:418;
inds = find(tril(ones(length(goodps)),-1)==1);
WBcompare_u_1 = result.boot_result.compare_u(:,lv);
temp11=zeros(418,418);temp11(inds)=WBcompare_u_1;
temp1z = temp11+temp11'; 

temp1 = temp1z; 
temp1(abs(temp1z) <=3)=0;
temp1(temp1z >3)=1;
temp1(temp1z <-3)=-1;
imagesc(temp1,[-.2,.2]); hold on

colormap(bluewhitered), colorbar;

xnames = {'Auditory','Cingulo-Opercular', 'Cingulo-Parietal','Default', 'Dorsal Attn', 'Fronto-Parietal',...
    'Orbito & Temp Pole', 'Retro-Temporal', 'Soma-Motor-h','Soma-Motor-m' ,...
    ' ','Ventral Attn','Visual', 'Subcortical'  'Cerebellar'};
netsizes1 = [0 24    40     5    41    32    24    47     8    38     8     4    23  39 54];
netsizes = [24    40     5    41    32    24    47     8    38     8     4    23    39  54  31];
set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',16);
set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',16);

makans = cumsum(netsizes1);

xcords = []; ycords = []; cccs = [];
for jj=1:15
    for j=1:length(netsizes)
        ycord = [makans(j);  makans(j);            makans(j)+netsizes(j);  makans(j)+netsizes(j)];
        xcord = [makans(jj);  makans(jj)+netsizes(jj);     makans(jj)+netsizes(jj);      makans(jj)];
        ccc = mean(mean(temp1(makans(jj)+1:makans(jj)+netsizes(jj), makans(j)+1:makans(j)+netsizes(j))));
        xcords = [xcords  xcord];
        ycords = [ycords  ycord];
        cccs = [cccs; ccc];
        
    end
end
ff = fill(xcords,ycords,cccs); colormap(bluewhitered), colorbar;
for i=1:225
    ff(i,1).EdgeColor = [1,1,1];
end
fill([0;418;418],[0;0;418],'white');

makansb = cumsum(netsizes);
for j=1:length(netsizes)
line([makansb(j),makansb(j)],[makansb(j) ,418],'Color','black');
end
for j=1:length(netsizes)
line([0 ,makansb(j)],[makansb(j),makansb(j)],'Color','black');
end

line([0 ,418],[0 ,418],'Color','black');
ylabel('Contribution to the Latent Variable (prop |Z|>3)','Position',[470,214,1]); axis square
[r p ] = corr(result.usc(:,1),result.vsc(:,1))

% right panle in Figure S2:
% export_fig -m5 -transparent rsFC_SU187nofam_abcd_pls_lv1_PLSpropmap.jpg


