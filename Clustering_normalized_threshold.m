%Calculation of the Clustering and Shortest Path Lenght normalizing the
%matrix using a threshold (5% of the highest weighted links)
clear all
%Generate the data
load("..\plv_78ROIs\sexy_cn_312s_plv_78_rois_alpha.mat")
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array

boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix

%%
for i=1:length(find(boys))
     %Boys
    g_boys = graph(threshold_norm(boy_rois(:,:,i)),'upper','omitselfloops');
    ad_boys = adjacency(g_boys);
    C_boys(i) = mean(clustering_coef_bu(ad_boys));
    D_boys = distances(g_boys,'Method',"unweighted");
    D_boys(D_boys==Inf)=0;
    dist_boys(i)= mean(D_boys(D_boys>0));
end

for i=1:length(find(girls))
     %Boys
    g_girls = graph(threshold_norm(girl_rois(:,:,i)),'upper','omitselfloops');
    ad_girls = adjacency(g_girls);
    C_girls(i) = mean(clustering_coef_bu(ad_girls));
    D_girls = distances(g_girls,'Method',"unweighted");
    D_girls(D_girls==Inf)=0;
    dist_girls(i)= mean(D_girls(D_girls>0));
end
%Transform to arrays 
C_boys = full(C_boys);
C_girls = full(C_girls);
dist_boys = full(dist_boys);
dist_girls = full(dist_girls);

%Save the values for further use 
save('C_boys.mat','C_boys');
save('C_girls.mat','C_girls');
save('dist_boys.mat','dist_boys');
save('dist_girls.mat','dist_girls');
%%
% Alternatively, you can load the data from GitHub
load('data_cluster_length_normalized_threshold\C_boys.mat');
load("data_cluster_length_normalized_threshold\C_girls.mat");
load("data_cluster_length_normalized_threshold\dist_boys.mat");
load("data_cluster_length_normalized_threshold\dist_girls.mat");
%%
%eTIV and age/sex classification 
[age_boys,eTIV_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls] = neuro_sex_load(girls,sample);

%%
paint_topo(C_boys,dist_boys,C_girls,dist_girls, age_boys,eTIV_boys,age_girls, eTIV_girls)
%%
function sex_ROIs =load_ROIs(sex, fcmatrix)
%%%
%Load the ROI's matrix and format them
%Arguments: sex -> Logical array
%fcmatrix: 4D array (Database of MEG)
%%%
    %Load
    sex_ROIs = zeros(78,78,length(find(sex)));
    sex_ROIs = fcmatrix(:,:,find(sex));
    %Format the matrix (Remove the diagonal and linear normalization
    for i = 1:length(find(sex))
        sex_ROIs(:,:,i) = sex_ROIs(:,:,i) - diag(diag(sex_ROIs(:,:,i)));
    end   
end

function thresh_plv = threshold_norm(plvmatrix)
%%%
%Creates N instances of the input PLV matrix taking the weights as
%probabilites
%%%

%Check if the graph is connected
%  g = graph(threshold_norm(boy_rois(:,:,1)),'upper');
%  bins = conncomp(g)
%  isConnected = all(bins == 1)

%It is also connected with only the 10% of the most relevant values.

    threshold = prctile(plvmatrix,80);
    thresh_plv = plvmatrix >= threshold;
end

function C=clustering_coef_bu(G)

%CLUSTERING_COEF_BU     Clustering coefficient
%
%   C = clustering_coef_bu(A);
%
%   The clustering coefficient is the fraction of triangles around a node
%   (equiv. the fraction of node?s neighbors that neighbors of each other).
%
%   Input:      A,      binary undirected connection matrix
%
%   Output:     C,      clustering coefficient vector
%
%   Reference: Watts and Strogatz (1998) Nature 393:440-442.
%
%
%   Mika Rubinov, UNSW, 2007-2010
n=length(G);
C=zeros(n,1);
for u=1:n
    V=find(G(u,:));
    k=length(V);
    if k>=2                %degree must be at least 2
        S=G(V,V);
        C(u)=sum(S(:))/(k^2-k);
    end
end
end

function paint_topo(C_boys,dist_boys,C_girls,dist_girls,age_boys,eTIV_boys,age_girls,eTIV_girls)
%%%
%Sequence of plots and stadistics
%%%

    %C vs dist 
    figure();
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