%Calculation of the Clustering and Shortest Path Lenght normalizing the
%matrix using weights as probabilities.
clear all
%Generate the data
load("..\plv_78ROIs\sexy_cn_312s_plv_78_rois_alpha.mat")
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array

boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix
%%
%MAIN
N = 10; %Number of instances
rng(123); %random seed
%Boys
for i=1:length(find(boys))
    boys_norm = binom_norm(boy_rois(:,:,i),N);
    for j = 1:N
        g_boys = graph(boys_norm(:,:,j),'upper','omitselfloops');
        N_boys = numnodes(g_boys);
        ad_boys = adjacency(g_boys);
        C_boys(i) = mean(clustering_coef_bu(ad_boys));
        D_boys = distances(g_boys,'Method',"unweighted");
        D_boys(D_boys==Inf)=0;
        dist_boys(i)= mean(D_boys(D_boys>0))/(N_boys-1);
    end
    errC_boys(i) = std(C_boys)./sqrt(length(C_boys));
    C_boys_mean(i) = mean(C_boys);
    
    errdist_boys(i) = std(dist_boys)./sqrt(length(dist_boys));
    dist_boys_mean(i) = mean(dist_boys);    
end
C_boys_mean = full(C_boys_mean);
errC_boys = full(errC_boys);

%Girls
for i=1:length(find(girls))
    girls_norm = binom_norm(girl_rois(:,:,i),N);
    for j =1:N
        g_girls = graph(girls_norm(:,:,j),'upper','omitselfloops');
        N_girls = numnodes(g_girls);
        ad_girls = adjacency(g_girls);
        C_girls(i) = mean(clustering_coef_bu(ad_girls));
        D_girls = distances(g_girls,'Method',"unweighted");
        D_girls(D_girls==Inf)=0;
        dist_girls(i)= mean(D_girls(D_girls>0))/(N_girls-1);
    end
    errC_girls(i) = std(C_girls)./sqrt(length(C_girls));
    C_girls_mean(i) = mean(C_girls);
    
    errdist_girls(i) = std(dist_girls)./sqrt(length(dist_girls));
    dist_girls_mean(i) = mean(dist_girls);    
end
C_girls_mean = full(C_girls_mean);
errC_girls = full(errC_girls);
%%
%Save the values for further use
save('C_boys.mat','C_boys_mean');
save('errC_boys.mat','errC_boys');
save('dist_boys.mat','dist_boys_mean');
save('errdist_boys.mat','errdist_boys');

save('C_girls.mat','C_girls_mean');
save('errC_girls.mat','errC_girls');
save('dist_girls.mat','dist_girls_mean');
save('errdist_girls.mat','errdist_girls');
%%
% Alternatively, you can load the data from GitHub
load("data_cluster_length_normalized_binom\C_boys.mat")
load("data_cluster_length_normalized_binom\C_girls.mat")
load("data_cluster_length_normalized_binom\dist_boys.mat")
load("data_cluster_length_normalized_binom\dist_girls.mat")
load("data_cluster_length_normalized_binom\errC_boys.mat")
load("data_cluster_length_normalized_binom\errC_girls.mat")
load("data_cluster_length_normalized_binom\errdist_boys.mat")
load("data_cluster_length_normalized_binom\errdist_girls.mat")
%%
%eTIV and age/sex classification 
[age_boys,eTIV_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls] = neuro_sex_load(girls,sample);

%%
paint_topo(C_boys_mean,dist_boys_mean,C_girls_mean,dist_girls_mean, age_boys,eTIV_boys,age_girls, eTIV_girls)
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
    %Format the matrix (Remove the diagonal)
    for i = 1:length(find(sex))
        sex_ROIs(:,:,i) = sex_ROIs(:,:,i) - diag(diag(sex_ROIs(:,:,i)));
    end   
end

function norm_plv = binom_norm(plvmatrix,N)
%%%
%Creates N instances of the input PLV matrix taking the weights as
%probabilites
%%%
    norm_plv = zeros(78,78,N);

    for k=1:N
        norm_plv(:,:,k) = binornd(1,plvmatrix);
    end
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