% Title: Functional brain connectivity predictors of prospective substance use initiation and their environmental correlates
% Contact: Omid Kardan omidk@med.umich.edu
% Revised for Biological Psychiatry: Cognitive Neuroscience and
% Neuroimaging

% PLS analyses scripts for Study 1 results (Fig 2 and supp Figs S1 , S2 , and S3)

% Requirements: 

% Requires running the SUI_abcd_PLS_main_and_supp.m
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
addpath(genpath('~\bluewhitered'))

load('~/PLS_result_233any_228hc_rsFCs_y0_y2_multrun.mat'); sucol = [.87 .77 .1]; ss = 'SU+';

ps = result.perm_result.sprob
load('~/with_2b_2024.mat')
load('~/SUI_233.mat')
load('~/matches_228.mat')
tab_fds = readtable('~\abcd_pls_rsFC_SUscores_and_subs_2024_233_multrun.csv');
site_ids = readtable('~/Sites.csv');
dems_table0 = read_demy0('~\df_filled_mh_cog_Y0.csv');
dems_table2 = read_demy2('~\df_filled_mh_cog_Y2.csv');

Sites = zeros(11868,22);
for K=1:22
    Sites((site_ids.site_num == K),K) = 1;
end

lv=1;   %~#

adjy0cont = zeros(922,1);
adjy0cont(1:228) = 1;

adjy2cont = zeros(922,1);
adjy2cont(229:456) = 1;

adjy0sui = zeros(922,1);
adjy0sui(457:689) = 1;

adjy2sui = zeros(922,1);
adjy2sui(690:922) = 1;

dd = tab_fds.used_ids;
covariates_y0cont = [tab_fds.pFD_y0(dd == 1)   Sites(dd==1,:)  dems_table0.White(dd == 1) dems_table0.Black(dd == 1) dems_table0.Hispanic(dd == 1) dems_table0.HighestEd(dd == 1) dems_table0.Income(dd == 1) dems_table0.Male_bin(dd == 1)  dems_table0.interview_age(dd == 1)];
covariates_y2cont = [tab_fds.pFD_y2(dd == 1)  Sites(dd==1,:) dems_table2.White(dd == 1) dems_table2.Black(dd == 1) dems_table2.Hispanic(dd == 1) dems_table2.HighestEd(dd == 1) dems_table2.Income(dd == 1) dems_table2.Male_bin(dd == 1)  dems_table2.interview_age(dd == 1)];
covariates_y0sui = [tab_fds.pFD_y0(dd == 2)   Sites(dd==2,:)  dems_table0.White(dd == 2) dems_table0.Black(dd == 2) dems_table0.Hispanic(dd == 2) dems_table0.HighestEd(dd == 2) dems_table0.Income(dd == 2) dems_table0.Male_bin(dd == 2)  dems_table0.interview_age(dd == 2)];
covariates_y2sui = [tab_fds.pFD_y2(dd == 2)  Sites(dd==2,:) dems_table2.White(dd == 2) dems_table2.Black(dd == 2) dems_table2.Hispanic(dd == 2) dems_table2.HighestEd(dd == 2) dems_table2.Income(dd == 2) dems_table2.Male_bin(dd == 2)  dems_table2.interview_age(dd == 2)];

covariates = [covariates_y0cont; covariates_y2cont; covariates_y0sui; covariates_y2sui];

[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates, 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates, 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0sui], covariates, 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2sui], covariates, 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[];
for jj=1:1000
   rand_inds1 = randi([1,228],228,1); rand_inds11 = randi([1,228],228,1); 
   rand_inds2 = randi([1,233],233,1); rand_inds22 = randi([1,233],233,1);
   inds = [rand_inds1; 228+ rand_inds11; 456+ rand_inds2; 456+233+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0sui(inds)], covariates(inds), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2sui(inds)] , covariates(inds), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];

end



    
 gt_contrast2 = [a1s, a2s, a3s, a4s];
 
brainPs = find(result.boot_result.compare_u(:,lv)>=3);
brainNs = find(result.boot_result.compare_u(:,lv)<=-3);
hcn = result.num_subj_lst(1); % length(groups{1})
sun = result.num_subj_lst(2); % length(groups{2})

figure;
hold on
xss = [1,1.8,3.2,4];
for k=1:4
if k<3; col = [.35 .65 0.35];  % color for columns 1 and 2 (HC)
else
    col = sucol;
end
   mm = mean(gt_contrast2(:,k));
   negs = mm - min(gt_contrast2(:,k));
   pos = max(gt_contrast2(:,k)) - mm;
   bar(xss(k),mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(xss(k),mean(gt_contrast2(:,k)),negs,pos,'.k');
end
set(gca,'XTick',[1,1.8,3.2,4],'XTickLabels',{'9-10 yrs old','11-12 yrs old','9-10 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS loading (partial corr)'); ylim([-.5,0.5]);
hold off


figure
subplot(1,2,1);
hold on
for k=1:4
if k<3 col = [.1 .2 .3];
else
    col = sucol;
end
   bar(k,mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor','k');
   errorbar(k,mean(gt_contrast2(:,k)),2*std(gt_contrast2(:,k)),'.k');
end
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'9-10 yrs old','11-12 yrs old','9-10 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',16);
ylabel('PLS loading (partial corr)'); axis square
hold off  

    
    cbcovs = result.s.^2./sum(result.s.^2);
    title(['LV',num2str(lv),...
        '  \sigma_{XY} = ',num2str(round(cbcovs(lv),3)),...
        '  p = ',num2str(round(ps(lv),3))]);


subplot(1,2,2);
goodps = 1:418;
inds = find(tril(ones(length(goodps)),-1)==1);
WMcompare_u_1 = result.boot_result.compare_u(:,lv);
temp11=zeros(418,418);temp11(inds)=WMcompare_u_1;
temp1z = temp11+temp11'; 
temp1 = temp1z; 
temp1(abs(temp1z) <=3)=0;
temp1(temp1z >3)=1;
temp1(temp1z <-3)=-1;
imagesc(temp1,[-.2,.2]); hold on   %#############

colormap(bluewhitered), colorbar;
load('E:\Omid\abcd\UMich_connectomeDev/GAnetnames.mat')
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
ylabel('Contribution to the Latent Variable (mean Z)','Position',[470,214,1]); axis square
[r p ] = partialcorr([result.usc(:,1),result.vsc(:,1)],covariates, 'Rows', 'pairwise')
%%
addpath(genpath('~\export_fig-master'));
export_fig -m5 -transparent rsFC_SU233multrun_abcd_pls_lv1_PLSweightbars_partial.jpg
%export_fig -m5 -transparent rsFC_SU233multrun_abcd_pls_lv1_PLSpropmap.jpg

%% supp 1 (one run per participant)

% plotting the results
clear all
addpath(genpath('~\bluewhitered'))

load('~/PLS_result_233any_228hc_rsFCs_y0_y2.mat'); sucol = [.87 .77 .1]; ss = 'SU+';

ps = result.perm_result.sprob
load('~/with_2b_2024.mat')
load('~/SUI_233.mat')
load('~/matches_228.mat')
tab_fds = readtable('~\abcd_pls_rsFC_SUscores_and_subs_2024_233.csv');
site_ids = readtable('~/Sites.csv');
dems_table0 = read_demy0('~\df_filled_mh_cog_Y0.csv');
dems_table2 = read_demy2('~\df_filled_mh_cog_Y2.csv');

Sites = zeros(11868,22);
for K=1:22
    Sites((site_ids.site_num == K),K) = 1;
end

lv=1;   %~#

adjy0cont = zeros(922,1);
adjy0cont(1:228) = 1;

adjy2cont = zeros(922,1);
adjy2cont(229:456) = 1;

adjy0sui = zeros(922,1);
adjy0sui(457:689) = 1;

adjy2sui = zeros(922,1);
adjy2sui(690:922) = 1;

dd = tab_fds.used_ids;
covariates_y0cont = [tab_fds.pFD_y0(dd == 1)   Sites(dd==1,:)  dems_table0.White(dd == 1) dems_table0.Black(dd == 1) dems_table0.Hispanic(dd == 1) dems_table0.HighestEd(dd == 1) dems_table0.Income(dd == 1) dems_table0.Male_bin(dd == 1)  dems_table0.interview_age(dd == 1)];
covariates_y2cont = [tab_fds.pFD_y2(dd == 1)  Sites(dd==1,:) dems_table2.White(dd == 1) dems_table2.Black(dd == 1) dems_table2.Hispanic(dd == 1) dems_table2.HighestEd(dd == 1) dems_table2.Income(dd == 1) dems_table2.Male_bin(dd == 1)  dems_table2.interview_age(dd == 1)];
covariates_y0sui = [tab_fds.pFD_y0(dd == 2)   Sites(dd==2,:)  dems_table0.White(dd == 2) dems_table0.Black(dd == 2) dems_table0.Hispanic(dd == 2) dems_table0.HighestEd(dd == 2) dems_table0.Income(dd == 2) dems_table0.Male_bin(dd == 2)  dems_table0.interview_age(dd == 2)];
covariates_y2sui = [tab_fds.pFD_y2(dd == 2)  Sites(dd==2,:) dems_table2.White(dd == 2) dems_table2.Black(dd == 2) dems_table2.Hispanic(dd == 2) dems_table2.HighestEd(dd == 2) dems_table2.Income(dd == 2) dems_table2.Male_bin(dd == 2)  dems_table2.interview_age(dd == 2)];

covariates = [covariates_y0cont; covariates_y2cont; covariates_y0sui; covariates_y2sui];

[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates, 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates, 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0sui], covariates, 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2sui], covariates, 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[];
for jj=1:1000
   rand_inds1 = randi([1,228],228,1); rand_inds11 = randi([1,228],228,1); 
   rand_inds2 = randi([1,233],233,1); rand_inds22 = randi([1,233],233,1);
   inds = [rand_inds1; 228+ rand_inds11; 456+ rand_inds2; 456+233+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0sui(inds)], covariates(inds), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2sui(inds)] , covariates(inds), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];

end



    
 gt_contrast2 = [a1s, a2s, a3s, a4s];
 
brainPs = find(result.boot_result.compare_u(:,lv)>=3);
brainNs = find(result.boot_result.compare_u(:,lv)<=-3);
hcn = result.num_subj_lst(1); % length(groups{1})
sun = result.num_subj_lst(2); % length(groups{2})

figure;
hold on
xss = [1,1.8,3.2,4];
for k=1:4
if k<3; col = [.35 .65 0.35];  % color for columns 1 and 2 (HC)
else
    col = sucol;
end
   mm = mean(gt_contrast2(:,k));
   negs = mm - min(gt_contrast2(:,k));
   pos = max(gt_contrast2(:,k)) - mm;
   bar(xss(k),mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(xss(k),mean(gt_contrast2(:,k)),negs,pos,'.k');
end
set(gca,'XTick',[1,1.8,3.2,4],'XTickLabels',{'9-10 yrs old','11-12 yrs old','9-10 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS loading (partial corr)'); ylim([-.5,0.5]);
hold off
[r p ] = partialcorr([result.usc(:,1),result.vsc(:,1)],covariates, 'Rows', 'pairwise')
%%
addpath(genpath('~\export_fig-master'));
export_fig -m5 -transparent rsFC_SU233onerun_abcd_pls_lv1_PLSweightbars_partial.jpg
%% supp 2 (one member per family)

% plotting the results
clear all
addpath(genpath('~\bluewhitered'))

load('~/PLS_result_187any_180hc_rsFCs_y0_y2_nofam.mat'); sucol = [.87 .77 .1]; ss = 'SU+';

ps = result.perm_result.sprob
load('~/with_2b_2024_nofam.mat')
load('~/SUI_187_nofam.mat')
load('~/matches_180_nofam.mat')
tab_fds = readtable('~\abcd_pls_rsFC_SUscores_and_subs_2024_187_nofam.csv');
site_ids = readtable('~/Sites.csv');
dems_table0 = read_demy0('~\df_filled_mh_cog_Y0.csv');
dems_table2 = read_demy2('~/df_filled_mh_cog_Y2.csv');

Sites = zeros(11868,22);
for K=1:22
    Sites((site_ids.site_num == K),K) = 1;
end

lv=1;   %~#

adjy0cont = zeros(734,1);
adjy0cont(1:180) = 1;

adjy2cont = zeros(734,1);
adjy2cont(181:360) = 1;

adjy0sui = zeros(734,1);
adjy0sui(361:547) = 1;

adjy2sui = zeros(734,1);
adjy2sui(548:734) = 1;

dd = tab_fds.used_ids;
covariates_y0cont = [tab_fds.pFD_y0(dd == 1)   Sites(dd==1,:)  dems_table0.White(dd == 1) dems_table0.Black(dd == 1) dems_table0.Hispanic(dd == 1) dems_table0.HighestEd(dd == 1) dems_table0.Income(dd == 1) dems_table0.Male_bin(dd == 1)  dems_table0.interview_age(dd == 1)];
covariates_y2cont = [tab_fds.pFD_y2(dd == 1)  Sites(dd==1,:) dems_table2.White(dd == 1) dems_table2.Black(dd == 1) dems_table2.Hispanic(dd == 1) dems_table2.HighestEd(dd == 1) dems_table2.Income(dd == 1) dems_table2.Male_bin(dd == 1)  dems_table2.interview_age(dd == 1)];
covariates_y0sui = [tab_fds.pFD_y0(dd == 2)   Sites(dd==2,:)  dems_table0.White(dd == 2) dems_table0.Black(dd == 2) dems_table0.Hispanic(dd == 2) dems_table0.HighestEd(dd == 2) dems_table0.Income(dd == 2) dems_table0.Male_bin(dd == 2)  dems_table0.interview_age(dd == 2)];
covariates_y2sui = [tab_fds.pFD_y2(dd == 2)  Sites(dd==2,:) dems_table2.White(dd == 2) dems_table2.Black(dd == 2) dems_table2.Hispanic(dd == 2) dems_table2.HighestEd(dd == 2) dems_table2.Income(dd == 2) dems_table2.Male_bin(dd == 2)  dems_table2.interview_age(dd == 2)];

covariates = [covariates_y0cont; covariates_y2cont; covariates_y0sui; covariates_y2sui];

[a1 ] = partialcorr([result.usc(:,lv),adjy0cont],covariates, 'Rows', 'pairwise');
[a2 ] = partialcorr([result.usc(:,lv),adjy2cont], covariates, 'Rows', 'pairwise');
[a3 ] = partialcorr([result.usc(:,lv),adjy0sui], covariates, 'Rows', 'pairwise');
[a4 ] = partialcorr([result.usc(:,lv),adjy2sui], covariates, 'Rows', 'pairwise');

a1s =[]; a2s =[]; a3s =[]; a4s =[];
for jj=1:1000
   rand_inds1 = randi([1,180],187,1); rand_inds11 = randi([1,180],180,1); 
   rand_inds2 = randi([1,187],187,1); rand_inds22 = randi([1,187],187,1);
   inds = [rand_inds1; 180+ rand_inds11; 360+ rand_inds2; 360+187+ rand_inds22];
   a11 = partialcorr([result.usc(inds,lv),adjy0cont(inds)], covariates(inds), 'Rows', 'pairwise');
   a1s = [a1s; a11(1,2)];
   a22 = partialcorr([result.usc(inds,lv),adjy2cont(inds)], covariates(inds), 'Rows', 'pairwise');
a2s = [a2s; a22(1,2)];
a33 = partialcorr([result.usc(inds,lv),adjy0sui(inds)], covariates(inds), 'Rows', 'pairwise');
a3s = [a3s; a33(1,2)];
a44 = partialcorr([result.usc(inds,lv),adjy2sui(inds)] , covariates(inds), 'Rows', 'pairwise');
a4s = [a4s; a44(1,2)];

end

    
gt_contrast2 = [a1s, a2s, a3s, a4s];
 
brainPs = find(result.boot_result.compare_u(:,lv)>=3);
brainNs = find(result.boot_result.compare_u(:,lv)<=-3);
hcn = result.num_subj_lst(1); % length(groups{1})
sun = result.num_subj_lst(2); % length(groups{2})

figure;
hold on
xss = [1,1.8,3.2,4];
for k=1:4
if k<3; col = [.35 .65 0.35];  % color for columns 1 and 2 (HC)
else
    col = sucol;
end
   mm = mean(gt_contrast2(:,k));
   negs = mm - min(gt_contrast2(:,k));
   pos = max(gt_contrast2(:,k)) - mm;
   bar(xss(k),mean(gt_contrast2(:,k)),'FaceColor',col,'EdgeColor',col);
   errorbar(xss(k),mean(gt_contrast2(:,k)),negs,pos,'.k');
end
set(gca,'XTick',[1,1.8,3.2,4],'XTickLabels',{'9-10 yrs old','11-12 yrs old','9-10 yrs old','11-12 yrs old'},...
    'XTickLabelRotation',45,'FontSize',22);
ylabel('PLS loading (partial corr)'); ylim([-.5,0.5]);
hold off
[r p ] = partialcorr([result.usc(:,1),result.vsc(:,1)],covariates, 'Rows', 'pairwise')
%%
addpath(genpath('~\export_fig-master'));
export_fig -m5 -transparent rsFC_SU187nofam_abcd_pls_lv1_PLSweightbars_partial.jpg
%%
