clear all
load("..\plv_78ROIs\sexy_cn_312s_plv_78_rois_alpha.mat")
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array
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
%Alternatively, you can load the data from GitHub
% load('data_cluster_length\C_boys.mat');
% load('data_cluster_length\C_girls.mat');
% load("data_cluster_length\dist_boys.mat");
% load("data_cluster_length\dist_girls.mat");
%%
paint_topo(C_boys,dist_boys,C_girls,dist_girls)
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

function paint_topo(C_boys,dist_boys,C_girls,dist_girls)
%%%
%Plot C vs d for boys and girls 
%%%
    figure();
    hold on;
    plot(C_boys,dist_boys,'o')
    plot(C_girls,dist_girls,'s')
    legend('Boys', 'Girls')
    xlabel('Clustering Coeff')
    ylabel('Shortest path length')
    hold off;
    
    figure();
    hold on;
    h1 = histogram(C_boys);
    h2 = histogram(C_girls);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    legend('Boys', 'Girls')
    xlabel('Clustering Coeff')
    hold off
    
    figure();
    hold on;
    h1 = histogram(dist_boys);
    h2 = histogram(dist_girls);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    legend('Boys', 'Girls')
    xlabel('Shortest path length')
    hold off
end