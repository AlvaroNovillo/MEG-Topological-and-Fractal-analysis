%Weigthed clustering coefficient C(w_ij)
clear all
%Generate the data
band = ["alpha"];
dum = sprintf('../plv_78ROIs/sexy_cn_312s_plv_78_rois_%s.mat',band);
files=dir(dum);
filename=horzcat(files.folder,'\',files.name);

load(filename); 
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array

boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix
k = 1;

%%
N = 100;
%Computation of the clustering coefficient and spl. 
for i=1:length(find(boys))
    g_boys = graph(boy_rois(:,:,i),'upper','omitselfloops');
    ad_boys = adjacency(g_boys,'weighted');
    C_boys(i) = mean(clustering_coef_wu(ad_boys));
    D_boys = distances(graph(1./ad_boys),'Method','positive');
    D_boys(D_boys==Inf)=0;
    dist_boys(i)= mean(D_boys(D_boys>0));
    deg_boys(i) = mean(degree(g_boys));
    for k = 1:N
        %For each network, we produce a random network with the same
        %properties
        n = g_boys.numnodes;
        e = g_boys.numedges;
        T=table(randi(n,e,2),abs(random('Normal',mean(g_boys.Edges.Weight), std(g_boys.Edges.Weight),1,length(g_boys.Edges.Weight))'),'VariableNames',{'EndNodes','Weights'});
        G=simplify(graph(T,'omitselfloops')); %Para eliminar enlaces repetidos
        %Proof that the graphs are connected.
        if length(unique(conncomp(G)))== 1
            k = k;
        else
            k = k - 1;
        end
        ad_rand_b = adjacency(G,G.Edges.Weights);
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
    g_girls= graph(girl_rois(:,:,i),'upper','omitselfloops');
    ad_girls = adjacency(g_girls,'weighted');
    C_girls(i) = mean(clustering_coef_wu(ad_girls));
    D_girls =distances(graph(1./ad_girls),'Method','positive');
    D_girls(D_girls==Inf)=0;
    dist_girls(i)= mean(D_girls(D_girls>0));
    deg_girls(i) = mean(degree(g_girls));
    for k = 1:N
        %For each network, we produce a random network with tbe same
        %properties
        n = g_girls.numnodes;
        e = g_girls.numedges;
        T=table(randi(n,e,2),abs(random('Normal',mean(g_girls.Edges.Weight), std(g_girls.Edges.Weight),1,length(g_girls.Edges.Weight))'),'VariableNames',{'EndNodes','Weights'});
        G=simplify(graph(T,'omitselfloops'));
        %Proof that the graphs are connected.
        if length(unique(conncomp(G)))== 1
            k = k;
        else
            k = k - 1;
        end
        ad_rand_g = adjacency(G,G.Edges.Weights);
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
%Raw data representation 
figure();
%C vs dist 
hold on;
plot(C_boys,dist_boys,'o')
plot(C_girls,dist_girls,'s')
legend('Man', 'Woman')
xlabel('Clustering Coeff')
ylabel('Shortest path length')
title('Raw data')
hold off;
%%

C_boys = C_boys./Cmean_rand_b;
C_girls = C_girls./Cmean_rand_g;
dist_boys = dist_boys./distmean_rand_b;
dist_girls = dist_girls./distmean_rand_g;
deg_boys = deg_boys./degmean_rand_b;
deg_girls = deg_girls./degmean_rand_g;
% % 
% % Save the values for further use 
% save('C_boys.mat','C_boys');
% save('C_girls.mat','C_girls');
% save('dist_boys.mat','dist_boys');
% save('dist_girls.mat','dist_girls');
% save('deg_boys.mat',"deg_boys");
% save('deg_girls.mat',"deg_girls");
%%
%Alternatively, you can load the data from GitHub
%alpha band 
load('data_cluster_length_weighted\alpha\C_boys.mat');
load('data_cluster_length_weighted\alpha\C_girls.mat');
load('data_cluster_length_weighted\alpha\dist_boys.mat');
load('data_cluster_length_weighted\alpha\dist_girls.mat');
load('data_cluster_length_weighted\alpha\deg_boys.mat');
load('data_cluster_length_weighted\alpha\deg_girls.mat');
%beta band
% load('data_cluster_length_weighted\beta\C_boys.mat');
% load("data_cluster_length_weighted\beta\C_girls.mat");
% load("data_cluster_length_weighted\beta\dist_boys.mat");
% load("data_cluster_length_weighted\beta\dist_girls.mat");
% load("data_cluster_length_weighted\beta\deg_boys.mat");
% load("data_cluster_length_weighted\beta\deg_girls.mat");

%%
%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);

%%
paint_topo(C_boys,dist_boys,C_girls,dist_girls, age_boys,eTIV_boys,age_girls, eTIV_girls,deg_boys,deg_girls,study_boys,study_girls)

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
         %Linear normalization of the weights
        sex_ROIs(:,:,i) = (sex_ROIs(:,:,i)- min(min(sex_ROIs(:,:,i))))/(max(max(sex_ROIs(:,:,i))) - min(min(sex_ROIs(:,:,i))));
    end   
end


function paint_topo(C_boys,dist_boys,C_girls,dist_girls,age_boys,eTIV_boys,age_girls,eTIV_girls,deg_boys, deg_girls,study_boys,study_girls)
%%%
%Sequence of plots and stadistics
%%%
    figure();
    %C vs dist 
    hold on;
    plot(C_boys,dist_boys,'o')
    plot(C_girls,dist_girls,'s')
    legend('Man', 'Woman')
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
    legend('Man', 'Woman')
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
    legend('Man', 'Woman')
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
    legend('Man', 'Woman')
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
    legend('Man', 'Woman')
    hold off;
    
    
    figure();
    plot(eTIV_boys,C_boys,'o')
    hold on;
    plot(eTIV_girls,C_girls,'s')
    xlabel('eTIV')
    ylabel('Clustering')
    title('Clustering Coeff vs eTIV')
    legend('Man', 'Woman')
    hold off;
    
    %Distance
    
    figure();
    plot(age_boys,dist_boys,'o')
    hold on;
    plot(age_girls,dist_girls,'s')
    xlabel('Age')
    ylabel('Shortest path length')
    title('Shortest path length vs age')
    legend('Man', 'Woman')
    hold off;
    
    
    figure();
    plot(eTIV_boys,dist_boys,'o')
    hold on;
    plot(eTIV_girls,dist_girls,'s')
    xlabel('eTIV')
    ylabel('Shortest path length')
    title('Shortest path length vs eTIV')
    legend('Man', 'Woman')
    hold off;
    
    %Degree
    figure();
    plot(age_boys,deg_boys,'o')
    hold on;
    plot(age_girls,deg_girls,'s')
    xlabel('Age')
    ylabel('Mean degree')
    title('Mean degree vs age')
    legend('Man', 'Woman')
    hold off;
    
    
    figure();
    plot(eTIV_boys,deg_boys,'o')
    hold on;
    plot(eTIV_girls,deg_girls,'s')
    xlabel('eTIV')
    ylabel('Mean degree')
    title('Mean degree vs eTIV')
    legend('Man', 'Woman')
    hold off;
    
    %Box-Whiskers plots 
    
    %C 
    C_box = [C_boys';C_girls'];
    g = [ones(size(C_boys')); 2*ones(size(C_girls'))]; %Boys and girls class.
    binEdges = 1:3;
    bins = {'Man', 'Woman'};
    groupSex = discretize(g,binEdges,'categorical',bins);
    
    figure();
    hold on;
    boxchart(groupSex,C_box,'GroupByColor',groupSex)
    ylabel('Clustering ')
    xlabel('gender')
    legend('Location',"best")
    title('Clustering Coefficient')
    hold off;

    %dist
    dist_box = [dist_boys';dist_girls'];
    g = [ones(size(dist_boys')); 2*ones(size(dist_girls'))]; %Boys and girls class.
    binEdges = 1:3;
    bins = {'Man', 'Woman'};
    groupSex = discretize(g,binEdges,'categorical',bins);
    
    figure();
    hold on;
    boxchart(groupSex,dist_box,'GroupByColor',groupSex)
    ylabel('Shortest path length')
    xlabel('gender')
    legend('Location',"best")
    title('Shortest path length')
    hold off;
    
    %degree
    k_box = [deg_boys';deg_girls'];
    g = [ones(size(deg_boys')); 2*ones(size(deg_girls'))]; %Boys and girls class.
    binEdges = 1:3;
    bins = {'Man', 'Woman'};
    groupSex = discretize(g,binEdges,'categorical',bins);
    
    figure();
    hold on;
    boxchart(groupSex,k_box,'GroupByColor',groupSex)
    ylabel('K')
    xlabel('gender')
    legend('Location',"best")
    title('K')
    hold off;
    
    
    %C per age
    figure();
    ages = [age_boys;age_girls];
    binEdges = 50:10:90;
    bins = {'50s','60s','70s','80s'};
    groupAge = discretize(ages,binEdges,'categorical',bins);
    boxchart(groupAge,C_box,'GroupByColor',groupSex)
    title('Clustering vs age')
    legend('Location',"best")
    ylabel('Clustering')
    xlabel('age')
    
    %dist per age
    figure();
    boxchart(groupAge,dist_box,'GroupByColor',groupSex)
    title('Shortest path length vs age')
    legend('Location',"best")
    ylabel('Shortest path length')
    xlabel('age')
    
    %K per age 
    figure();
    boxchart(groupAge,k_box,'GroupByColor',groupSex)
    title('K vs age')
    legend('Location',"best")
    ylabel('K')
    xlabel('age')
    
    
    %C per study_year

    figure();
    study_period = [study_boys;study_girls];
    study_period(isnan(study_period)) = 0;
    binEdges = 0:5:35;
    bins = {'-5','5','10' ,'15','20','25','+30'};
    groupStudy = discretize(study_period,binEdges,'categorical',bins);
    boxchart(groupStudy,C_box,'GroupByColor',groupSex)
    title('Clustering vs Study period')
    legend('Location',"best")
    ylabel('Clustering')
    xlabel('Study period')
    
    %dist per study_year
    figure();
    groupStudy = discretize(study_period,binEdges,'categorical',bins);
    boxchart(groupStudy,dist_box,'GroupByColor',groupSex)
    title('Shortest path length vs Study period')
    legend('Location',"best")
    ylabel('Shortest path length')
    xlabel('Study period')
    
    %K per study_year
    figure();
    groupStudy = discretize(study_period,binEdges,'categorical',bins);
    boxchart(groupStudy,k_box,'GroupByColor',groupSex)
    title('K  vs Study period')
    legend('Location',"best")
    ylabel('K')
    xlabel('Study period')
    
    %Years of study and age vs C
    figure();
    binEdges = 0:10:30;
    bins = {'0','10','20'};
    groupStudy = discretize(study_period,binEdges,'categorical',bins);
    boxchart(groupStudy,C_box,'GroupByColor',groupAge)
    title('Clustering vs Study period per age')
    legend('Location',"best")
    ylabel('Clustering')
    xlabel('Study period')
    
    %Years of study and age vs dist
    figure();
    boxchart(groupStudy,dist_box,'GroupByColor',groupAge)
    title('Shortest Path Length vs Study period per age')
    legend('Location',"best")
    ylabel('Shortest path length')
    xlabel('Study period')
end

function [age_sex,eTIV_sex,study_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
    study_sex = sample.neuro_vals(find(sex),3);
end