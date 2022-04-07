%Weigthed clustering coefficient C(p_ij)
clear all
%Generate the data
load("..\plv_nodes\sexy_cn_312s_plv_1202nodes_alpha.mat")
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array
k = 1;
%%
boy_nodes = load_nodes(boys, fcmatrix); %Boys ROIs matrix
girl_nodes = load_nodes(girls, fcmatrix); %Girls ROIs matrix
N = 10;
%%
%Computation of the clustering coefficient and spl. 
for i=1:length(find(boys))
    g_boys = graph(boy_nodes(:,:,i),'upper','omitselfloops');
    C_boys(i) = mean(clustering_coef_wu(boy_nodes(:,:,i)));
    D_boys = distances(graph(1./boy_nodes(:,:,i)),'Method','positive');
    D_boys(D_boys==Inf)=0;
    dist_boys(i)= mean(D_boys(D_boys>0));
    deg_boys(i) = mean(degree(g_boys));
    for k = 1:N
        %For each network, we produce a random network with tbe same
        %properties
        n = g_boys.numnodes;
        e = g_boys.numedges;
        T=table(randi(n,e,2),'Var',{'EndNodes'});
        G=graph(T,'omitselfloops');
        %Proof that the graphs are connected.
        if length(unique(conncomp(G)))== 1
            k = k;
        else
            k = k - 1;
        end
        ad_rand_b = adjacency(simplify(G),'weighted');
        C_rand_b(k) = mean(clustering_coef_wu(ad_rand_b));
        D_rand_b = distances(graph(1./ad_rand_b),'Method','positive');
        D_rand_b(D_rand_b==Inf)=0;
        dist_rand_b(k)= mean(D_rand_b(D_rand_b>0));
        deg_rand_b(k) = mean(degree(G));
    end 
    Cmean_rand_b(i) = mean(C_rand_b);
    distmean_rand_b(i) = mean(dist_rand_b);
    degmean_rand_b(i) = mean(deg_rand_b);
    
end
for i=1:length(find(girls))
    g_girls= graph(girl_nodes(:,:,i),'upper','omitselfloops');
    C_girls(i) = mean(clustering_coef_wu(girl_nodes(:,:,i)));
    D_girls =distances(graph(1./girl_nodes(:,:,i)),'Method','positive');
    D_girls(D_girls==Inf)=0;
    dist_girls(i)= mean(D_girls(D_girls>0));
    deg_girls(i) = mean(degree(g_girls));
    for k = 1:N
        %For each network, we produce a random network with tbe same
        %properties
        n = g_girls.numnodes;
        e = g_girls.numedges;
        T=table(randi(n,e,2),'Var',{'EndNodes'});
        G=graph(T,'omitselfloops');
        %Proof that the graphs are connected.
        if length(unique(conncomp(G)))== 1
            k = k;
        else
            k = k - 1;
        end
        ad_rand_g = adjacency(simplify(G),'weighted');
        C_rand_g(k) = mean(clustering_coef_wu(ad_rand_g));
        D_rand_g = distances(graph(1./ad_rand_g),'Method','positive');
        D_rand_g(D_rand_g==Inf)=0;
        dist_rand_g(k)= mean(D_rand_g(D_rand_g>0));
        deg_rand_g(k) = mean(degree(G));
    end 
    Cmean_rand_g(i) = mean(C_rand_g);
    distmean_rand_g(i) = mean(dist_rand_g);
    degmean_rand_g(i) = mean(deg_rand_g);
end
%Transform to arrays 
C_boys = full(C_boys);
C_girls = full(C_girls);
dist_boys = full(dist_boys);
dist_girls = full(dist_girls);

Cmean_rand_b = full(Cmean_rand_b);
distmean_rand_b = full(distmean_rand_b);


Cmean_rand_g = full(Cmean_rand_g);
distmean_rand_g = full(distmean_rand_g);
%%

C_boys = C_boys./Cmean_rand_b;
C_girls = C_girls./Cmean_rand_g;
dist_boys = dist_boys./distmean_rand_b;
dist_girls = dist_girls./distmean_rand_g;
deg_boys = deg_boys./degmean_rand_b;
deg_girls = deg_girls./degmean_rand_g;

% Save the values for further use 
save('C_boys.mat','C_boys');
save('C_girls.mat','C_girls');
save('dist_boys.mat','dist_boys');
save('dist_girls.mat','dist_girls');
save('deg_boys.mat',"deg_boys");
save('deg_girls.mat',"deg_girls");
%%
% Alternatively, you can load the data from GitHub
load('data_cluster_nodes\C_boys.mat');
load("data_cluster_nodes\C_girls.mat");
load("data_cluster_nodes\dist_boys.mat");
load("data_cluster_nodes\dist_girls.mat");
load("data_cluster_nodes\deg_boys.mat");
load("data_cluster_nodes\deg_girls.mat");
%%
%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);

%%
paint_topo(C_boys,dist_boys,C_girls,dist_girls, age_boys,eTIV_boys,age_girls, eTIV_girls,deg_boys,deg_girls)

%% 
% 

function sex_nodes =load_nodes(sex, fcmatrix)
%%%
%Load the ROI's matrix and format them
%Arguments: sex -> Logical array
%fcmatrix: 4D array (Database of MEG)
%%%
    %Load
    sex_nodes = zeros(78,78,length(find(sex)));
    sex_nodes = fcmatrix(:,:,find(sex));
    %Format the matrix (Remove the diagonal)
    for i = 1:length(find(sex))
        sex_nodes(:,:,i) = sex_nodes(:,:,i) - diag(diag(sex_nodes(:,:,i)));
         %Linear normalization of the weights
        sex_nodes(:,:,i) = (sex_nodes(:,:,i)- min(min(sex_nodes(:,:,i))))/(max(max(sex_nodes(:,:,i))) - min(min(sex_nodes(:,:,i))));
    end   
end


function paint_topo(C_boys,dist_boys,C_girls,dist_girls,age_boys,eTIV_boys,age_girls,eTIV_girls,deg_boys, deg_girls)
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
%     sprintf('The mean of the distribution for C is %.3f for boys and %.3f for girls', mean(full(C_boys)), mean(full(C_girls),'omitnan'))
%     sprintf('The std of the distribution for C is %.3f for boys and %.3f for girls', std(full(C_boys)), std(full(C_girls),'omitnan'))
    
    figure();
    %Histogram of d
    hold on;
    h1 = histogram(dist_boys);
    h2 = histogram(dist_girls);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    legend('Boys', 'Girls')
    xlabel('Shortest path length')
%     sprintf('The mean of the distribution for dist is %.3f for boys and %.3f for girls', mean(full(dist_boys),'omitnan'), mean(full(dist_girls),'omitnan'))
%     sprintf('The std of the distribution for dist is %.3f for boys and %.3f for girls', std(full(dist_boys), 'omitnan'), std(full(dist_girls),'omitnan'))
    hold off
    
    %Mean degree
    figure();
    hold on;
    h1 = histogram(deg_boys);
    h2 = histogram(deg_girls);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    legend('Boys', 'Girls')
    xlabel('Mean degree')
    hold off;
    
    
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
    
    %Degree
    figure();
    plot(age_boys,deg_boys,'o')
    hold on;
    plot(age_girls,deg_girls,'s')
    xlabel('Age')
    ylabel('Mean degree')
    title('Mean degree vs age')
    legend('Boys', 'Girls')
    hold off;
    
    
    figure();
    plot(eTIV_boys,deg_boys,'o')
    hold on;
    plot(eTIV_girls,deg_girls,'s')
    xlabel('eTIV')
    ylabel('Mean degree')
    title('Mean degree vs eTIV')
    legend('Boys', 'Girls')
    hold off;
    
end

function [age_sex,eTIV_sex,study_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
    study_sex = sample.neuro_vals(find(sex),3);
end