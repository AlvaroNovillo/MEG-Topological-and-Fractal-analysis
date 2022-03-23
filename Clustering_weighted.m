clear all
%Generate the data
load("..\plv_78ROIs\sexy_cn_312s_plv_78_rois_alpha.mat")
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array
%%

boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix

%Node number 
%Computation of the clustering coefficient and spl. 
for i=1:length(find(boys))
    g_boys = graph(boy_rois(:,:,i),'upper','omitselfloops');
    ad_boys = adjacency(g_boys,'weighted');
    C_boys(i) = mean(clustering_coef_wu(ad_boys));
    D_boys = distance_wei(ad_boys);
    D_boys(D_boys==Inf)=0;
    dist_boys(i)= mean(D_boys(D_boys>0));
end
for i=1:length(find(girls))
    g_girls= graph(girl_rois(:,:,i),'upper','omitselfloops');
    ad_girls = adjacency(g_girls,'weighted');
    C_girls(i) = mean(clustering_coef_wu(ad_girls));
    D_girls = distance_wei(ad_girls);
    D_girls(D_girls==Inf)=0;
    dist_girls(i)= mean(D_girls(D_girls>0));
end

%There was a subject whose values where to high
dist_girls(dist_girls == max(dist_girls)) = nan;
C_girls(C_girls == max(C_girls)) = nan;

%Save the values for further use 
save('C_boys.mat','C_boys');
save('C_girls.mat','C_girls');
save('dist_boys.mat','dist_boys');
save('dist_girls.mat','dist_girls');
%%
% Alternatively, you can load the data from GitHub
load('data_cluster_length_weighted\C_boys.mat');
load("data_cluster_length_weighted\C_girls.mat");
load("data_cluster_length_weighted\dist_boys.mat");
load("data_cluster_length_weighted\dist_girls.mat");
%%
%eTIV and age/sex classification 
[age_boys,eTIV_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls] = neuro_sex_load(girls,sample);

%%
paint_topo(C_boys,dist_boys,C_girls,dist_girls, age_boys,eTIV_boys,age_girls, eTIV_girls)
%% 
% 

function sex_ROIs =load_ROIs(sex, fcmatrix)
%%%
%Load the ROI's matrix and format them
%Arguments: sex -> Logical array
%fcmatrix: 4D array (Database of MEG)
%%%
    %Load
    sex_ROIs = zeros(78,78,length(find(sex)));
    sex_ROIs = fcmatrix(:,:,find(sex));
    %Format the matrix (Remove the diagonal)
    for i = 1:length(find(sex))
        sex_ROIs(:,:,i) = sex_ROIs(:,:,i) - diag(diag(sex_ROIs(:,:,i)));
    end   
end

function paint_topo(C_boys,dist_boys,C_girls,dist_girls,age_boys,eTIV_boys,age_girls,eTIV_girls)
%%%
%Sequence of plots and stadistics
%%%
    figure();
    %C vs dist 
    hold on;
    plot(C_boys,dist_boys,'o')
    plot(C_girls,dist_girls,'s')
    legend('Boys', 'Girls')
    xlabel('Clustering Coeff')
    ylabel('Shortest path length')
    hold off;
    
    figure();
    %Histogram of C %Add fit distribution
    hold on;
    h1 = histogram(C_boys);
    h2 = histogram(C_girls);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    legend('Boys', 'Girls')
    xlabel('Clustering Coeff')
    hold off;
    sprintf('The mean of the distribution for C is %.3f for boys and %.3f for girls', mean(full(C_boys)), mean(full(C_girls),'omitnan'))
    sprintf('The std of the distribution for C is %.3f for boys and %.3f for girls', std(full(C_boys)), std(full(C_girls),'omitnan'))
    
    figure();
    %Histogram of d
    hold on;
    h1 = histogram(dist_boys);
    h2 = histogram(dist_girls);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    legend('Boys', 'Girls')
    xlabel('Shortest path length')
    sprintf('The mean of the distribution for dist is %.3f for boys and %.3f for girls', mean(full(dist_boys)), mean(full(dist_girls),'omitnan'))
    sprintf('The std of the distribution for dist is %.3f for boys and %.3f for girls', std(full(dist_boys)), std(full(dist_girls),'omitnan'))
    hold off
    
    %Clustering
    figure();
    plot(age_boys,C_boys,'o')
    hold on;
    plot(age_girls,C_girls,'s')
    xlabel('Age')
    ylabel('Clustering')
    title('Clustering Coeff vs age')
    legend('Boys', 'Girls')
    hold off;
    
    
    figure();
    plot(eTIV_boys,C_boys,'o')
    hold on;
    plot(eTIV_girls,C_girls,'s')
    xlabel('eTIV')
    ylabel('Clustering')
    title('Clustering Coeff vs eTIV')
    legend('Boys', 'Girls')
    hold off;
    
    %Distance
    
    figure();
    plot(age_boys,dist_boys,'o')
    hold on;
    plot(age_girls,dist_girls,'s')
    xlabel('Age')
    ylabel('Shortest path length')
    title('Shortest path length vs age')
    legend('Boys', 'Girls')
    hold off;
    
    
    figure();
    plot(eTIV_boys,dist_boys,'o')
    hold on;
    plot(eTIV_girls,dist_girls,'s')
    xlabel('eTIV')
    ylabel('Shortest path length')
    title('Shortest path length vs eTIV')
    legend('Boys', 'Girls')
    hold off;
end

function [age_sex,eTIV_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
end