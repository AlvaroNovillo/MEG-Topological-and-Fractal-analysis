%Intra - ROI study
clear all;
%Generate the data
band = ["alpha"];
dum = sprintf('../plv_nodes/sexy_cn_312s_plv_1202nodes_%s.mat',band);
files=dir(dum);
filename=horzcat(files.folder,'\',files.name);
load(filename); %DESCARGAR ARCHIVO BETA
load("..\fcinfo.mat") %Map nodes 2 ROIs
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array
% 
%%
%Variables needed
ROI_names = fcinfo.rois_names;
nodes2rois = fcinfo.nodes1202_in_rois;
nodes = fcinfo.nodes;
ROIs = fcinfo.rois_in;
boy_nodes = load_nodes(boys, fcmatrix); %Boys nodes matrix
girl_nodes = load_nodes(girls, fcmatrix); %Girls nodes matrix

names = {'ROIs_names', 'Clustering', 'C_error','Shortest_Path_Length', 'Dist_error','Full_C','Full_dist'};
c = cell(length(names),1);
results_boys = cell2struct(c,names);
results_girls = cell2struct(c,names);


%%
%Boys
for i = 1:length(ROIs)
    for j = 1:length(find(boys))
        g_b_ROI = graph(boy_nodes(nodes2rois{i},nodes2rois{i},j),'upper','omitselfloops');
        C_boys(j) = mean(clustering_coef_wu(boy_nodes(nodes2rois{i},nodes2rois{i},j)));
        D_boys = distances(graph(1./boy_nodes(nodes2rois{i},nodes2rois{i},j)),'Method','positive');
        D_boys(D_boys==Inf)=0;
        dist_boys(j)= mean(D_boys(D_boys>0));
    end
    results_boys.ROIs_names(i) = ROI_names(i);
    results_boys.Clustering(i) = mean(C_boys);
    results_boys.Shortest_Path_Length(i) = mean(dist_boys);
    results_boys.Full_C(:,i) = C_boys;
    results_boys.Full_dist(:,i) = dist_boys;
    
    %Error
    results_boys.C_error(i) = std(C_boys)./sqrt(length(C_boys));
    results_boys.Dist_error(i) = std(dist_boys)./sqrt(length(dist_boys));
end


%Girls
for i = 1:length(ROIs)
    for j = 1:length(find(girls))
        g_g_ROI = graph(girl_nodes(nodes2rois{i},nodes2rois{i},j),'upper','omitselfloops');
        C_girls(j) = mean(clustering_coef_wu(girl_nodes(nodes2rois{i},nodes2rois{i},j)));
        D_girls = distances(graph(1./girl_nodes(nodes2rois{i},nodes2rois{i},j)),'Method','positive');
        D_girls(D_girls==Inf)=0;
        dist_girls(j)= mean(D_girls(D_girls>0));
    end
    results_girls.ROIs_names(i) = ROI_names(i);
    results_girls.Clustering(i) = mean(C_girls);
    results_girls.Shortest_Path_Length(i) = mean(dist_girls);
    results_girls.Full_C(:,i) = C_girls;
    results_girls.Full_dist(:,i) = dist_girls;
    
    %Error
    results_girls.C_error(i) = std(C_girls)./sqrt(length(C_girls));
    results_girls.Dist_error(i) = std(dist_girls)./sqrt(length(dist_girls));
end
%%
dum = sprintf('../plv_78ROIs/sexy_cn_312s_plv_78_rois_%s.mat',band);
files=dir(dum);
filename=horzcat(files.folder,'\',files.name);
load(filename) %ROIs
%Connectivity diagonal
boys_diag = load_connectivity(boys, fcmatrix);
girls_diag = load_connectivity(girls, fcmatrix);

%Add to the result cell
results_boys.Connect = boys_diag.Connect'; results_boys.ConnErr = boys_diag.Error';
results_girls.Connect = girls_diag.Connect'; results_girls.ConnErr = girls_diag.Error';
%%
%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);


%Add to the result cell
results_boys.age_boys = age_boys'; results_boys.eTIV_boys = eTIV_boys'; results_boys.study_boys = study_boys';
results_girls.age_girls = age_girls'; results_girls.eTIV_girls = eTIV_girls'; results_girls.study_girls = study_girls';
%%
%Save the results 

save(sprintf('Intra_ROI_results/%s/boys.mat',band),"results_boys")
save(sprintf('Intra_ROI_results/%s/girls.mat',band),"results_girls")
%%
%Load the results 
clear all; 

band = ["alpha"];

load(sprintf('Intra_ROI_results/%s/boys.mat',band))
load(sprintf('Intra_ROI_results/%s/girls.mat',band))
%%
paint_topo(results_boys.Clustering, results_boys.Shortest_Path_Length,results_girls.Clustering,results_girls.Shortest_Path_Length,results_boys.ROIs_names,results_boys.C_error,results_boys.Dist_error,results_girls.C_error,results_girls.Dist_error,results_boys.Connect, results_girls.Connect)
%%
select_paint(3,results_boys.ROIs_names{3},results_boys.Full_C,results_boys.Full_dist,results_girls.Full_C,results_girls.Full_dist,results_boys.age_boys,results_boys.eTIV_boys,results_girls.age_girls,results_girls.eTIV_girls,results_boys.study_boys,results_girls.study_girls)
%%
function [age_sex,eTIV_sex,study_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
    study_sex = sample.neuro_vals(find(sex),3);
end

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


function sex_diag_avg = load_connectivity(sex, fcmatrix)
%%%
%Load the ROI's matrix and extract the diagonal elements
%Arguments: sex -> Logical array
%fcmatrix: 4D array (Database of MEG)
%%%

    names = {'Connect', 'Error'};
    c = cell(length(names),1);
    sex_diag_avg= cell2struct(c,names);

    %Load
    sex_nodes = zeros(78,78,length(find(sex)));
    sex_nodes = fcmatrix(:,:,find(sex));
    %Extract the diagonal
    for i = 1:length(find(sex))
        sex_diag(:,i) = diag(sex_nodes(:,:,i));
        
    end   
    
    sex_diag_avg.Connect = mean(sex_diag,2);
    sex_diag_avg.Error = std(sex_diag,0,2)./sqrt(78);
    
end

function paint_topo(C_boys,dist_boys,C_girls,dist_girls,ROI_tag,C_b_err,dist_b_err,C_g_err,dist_g_err,conn_boys,conn_girls)
%%%
%Sequence of plots and stadistics
%%%
    figure();
    %C vs dist 
    hold on;
    errorbar(C_boys,dist_boys,C_b_err,'o')
    errorbar(C_girls,dist_girls,dist_g_err,'s')
    legend('Man', 'Woman')
    xlabel('Clustering Coeff')
    ylabel('Shortest path length')
    hold off;
    
    
    %C and dist per ROI
    figure();
    hold on;
    tag = categorical(unique(ROI_tag));
    errorbar(tag,C_boys,C_b_err, '-o')
    errorbar(tag,dist_boys,dist_b_err,'-o')
    
    errorbar(tag,C_girls,C_g_err, '-s')
    errorbar(tag,dist_girls,dist_g_err,'-s')
    
    legend('Man Clustering', 'Man Distance' ,'Woman Clustering', 'Woman Distance');
    
    xlabel('ROIs names')
    ylabel('C and Shortest path')
    
    
    %Connectivity vs C
    figure();
    hold on;
    plot(C_boys,conn_boys,'o')
    plot(C_girls,conn_girls,'s')
    legend('Man', 'Woman')
    xlabel('Clustering Coeff')
    ylabel('Connectivity')
    hold off;
    
    
    %Connectivity vs dist
    figure();
    hold on;
    plot(dist_boys,conn_boys,'o')
    plot(dist_girls,conn_girls,'s')
    legend('Man', 'Woman')
    xlabel('Shortest path length')
    ylabel('Connectivity')
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
end

function select_paint(ROI,ROI_name,C_boys,dist_boys,C_girls,dist_girls,age_boys,eTIV_boys,age_girls,eTIV_girls,study_boys,study_girls)

    %C 
    C_box = [double(C_boys(:,ROI))',double(C_girls(:,ROI))'];
    g = [ones(size(C_boys(:,ROI)')), 2*ones(size(C_girls(:,ROI)'))]; %Boys and girls class.
    binEdges = 1:3;
    bins = {'Man', 'Woman'};
    groupSex = discretize(g,binEdges,'categorical',bins);
    
    
    %dist
    dist_box = [dist_boys(:,ROI)',dist_girls(:,ROI)'];
    g = [ones(size(dist_boys(:,ROI)')), 2*ones(size(dist_girls(:,ROI)'))]; %Boys and girls class.
    binEdges = 1:3;
    bins = {'Man', 'Woman'};
    groupSex = discretize(g,binEdges,'categorical',bins);
    
    %C per age
    figure();
    ages = [age_boys,age_girls];
    binEdges = 50:10:90;
    bins = {'50s','60s','70s','80s'};
    groupAge = discretize(ages,binEdges,'categorical',bins);
    boxchart(groupAge,C_box,'GroupByColor',groupSex)
    title('Clustering vs age')
    legend('Location',"best")
    ylabel(sprintf('%s Clustering',ROI_name))
    xlabel('age')
    
    %dist per age
    figure();
    boxchart(groupAge,dist_box,'GroupByColor',groupSex)
    title('Shortest path length vs age')
    legend('Location',"best")
    ylabel(sprintf('%s Shortest Path Length',ROI_name))
    xlabel('age')
    
     %C per study_year

    figure();
    study_period = [study_boys,study_girls];
    study_period(isnan(study_period)) = 0;
    binEdges = 0:5:35;
    bins = {'-5','5','10' ,'15','20','25','+30'};
    groupStudy = discretize(study_period,binEdges,'categorical',bins);
    boxchart(groupStudy,C_box,'GroupByColor',groupSex)
    title('Clustering vs Study period')
    legend('Location',"best")
    ylabel(sprintf('%s Clustering',ROI_name))
    xlabel('Study period')
    
    %dist per study_year
    figure();
    groupStudy = discretize(study_period,binEdges,'categorical',bins);
    boxchart(groupStudy,dist_box,'GroupByColor',groupSex)
    title('Shortest path length vs Study period')
    legend('Location',"best")
    ylabel(sprintf('%s Shortest Path Length',ROI_name))
    xlabel('Study period')
    
    %Years of study and age vs C
    figure();
    binEdges = 0:10:30;
    bins = {'0','10','20'};
    groupStudy = discretize(study_period,binEdges,'categorical',bins);
    boxchart(groupStudy,C_box,'GroupByColor',groupAge)
    title('Clustering vs Study period per age')
    legend('Location',"best")
    ylabel(sprintf('%s Clustering',ROI_name))
    xlabel('Study period')
    
    %Years of study and age vs dist
    figure();
    boxchart(groupStudy,dist_box,'GroupByColor',groupAge)
    title('Shortest Path Length vs Study period per age')
    legend('Location',"best")
    ylabel(sprintf('%s Shortest Path Length',ROI_name))
    xlabel('Study period')

end