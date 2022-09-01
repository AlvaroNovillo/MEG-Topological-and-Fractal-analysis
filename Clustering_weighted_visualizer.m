%% VISUALIZATION OF RESULTS

% %Run from here if you have already downloaded/generated the results
% 
%Generate the data
clear all;
band = ["alpha"]; %You can try 'alpha' or 'beta'
gband='\alpha ';%You can try '\alpha' or '\beta';

type='rois';
% type='215nodes';
% type='1202nodes';
%% 78X78 REPRESENTATION

if(strcmp(type,'rois')==1)
dum = sprintf('../plv_78ROIs/sexy_cn_312s_plv_78_%s_%s.mat',type,band);
%% 215X215 REPRESENTATION
% 

elseif (strcmp(type,'215nodes')==1)
dum = sprintf('../plv_215/sexy_cn_312s_plv_%s_%s.mat',type,band);
%% 1202 x 1202 REPRESENTATION

elseif(strcmp(type,'1202nodes')==1)
dum = sprintf("../plv_nodes/sexy_cn_312s_plv_%s_%s.mat",type,band);
end
%%
%Load the selected representation
files=dir(dum);
filename=horzcat(files.folder,'/',files.name);

load(filename);

load(sprintf('%s/%s/boys.mat',type,band));
load(sprintf('%s/%s/girls.mat',type,band));

%Boys
C_boys = results_boys.Raw_C;
dist_boys = results_boys.Raw_L; 
C_b = results_boys.C;
L_b = results_boys.L; 
%Girls
C_girls = results_girls.Raw_C;
dist_girls = results_girls.Raw_L; 
C_g = results_girls.C;
L_g = results_girls.L; 

figure();
%C vs dist 
subplot(1,2,1),hold on;
plot(C_boys,dist_boys,'o')
plot(C_girls,dist_girls,'s')
%axis([0 0.5 3 9])
legend('Man', 'Woman')
xlabel('C')
ylabel('L')
tit=sprintf('Raw data %s-%s',gband,type);
title(tit)
axis square
hold off;
% 
subplot(1,2,2),hold on
plot(C_b,L_b,'o')
plot(C_g,L_g,'s')
%axis([0.95 1.05 1. 2])

legend('Man', 'Woman','Location','southeast')
xlabel('C/C_{rand}')
ylabel('L/L_{rand}')

tit=sprintf('Norm data %s-%s',gband,type);
title(tit)
axis square
hold off;

% filename=sprintf('L_vs_C_%s_%s_surrog%d.pdf',band,type,N);
% print('-dpdf',filename,'-bestfit')
toc
%% 

%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array


% eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);


paint_topo(C_b,L_b,C_g,L_g, age_boys,eTIV_boys,age_girls,eTIV_girls,study_boys,study_girls)

%%
function sex_ROIs =load_ROIs(sex, fcmatrix)
%Load the ROI's matrix and format them
%Arguments: sex -> Logical array
%fcmatrix: 4D array (Database of MEG)
    %Load
    nROIs=size(fcmatrix,1);%Network size
    nsubjects=sum(sex);%Number of subjects of each sex type
    sex_ROIs = zeros(nROIs,nROIs,nsubjects);%esta linea se puede borrar
    sex_ROIs = fcmatrix(:,:,sex);
    ind=~logical(eye(nROIs));%indexes of the non diagonal elements     
    id=logical(eye(nROIs));%indexes of the diagonal elements
    %Format the matrix (Remove the diagonal)
    for i = 1:nsubjects
        
         %Linear normalization of the weights
        dum=sex_ROIs(:,:,i);        
        minval=min(dum(ind));%min value among the non diagonal elements
        maxval=max(dum(ind));%max value among the non diagonal elements
        meanval=mean(dum(ind));%mean value of all weights 

        %sex_ROIs(:,:,i) = (sex_ROIs(:,:,i)- min(min(sex_ROIs(:,:,i))))/(max(max(sex_ROIs(:,:,i))) - min(min(sex_ROIs(:,:,i))));
        dum = (dum-minval)./(maxval-minval);    
        %dum=dum./meanval; %testing another normalization 
        dum(id)=0;%setting diagonal elements to 0 -Note that in the ROI representation the element (40,40) is a NaN
        %sex_ROIs(:,:,i) = sex_ROIs(:,:,i) - diag(diag(sex_ROIs(:,:,i)));
        sex_ROIs(:,:,i)=dum;
    end   
end


function paint_topo(C_boys,dist_boys,C_girls,dist_girls,age_boys,eTIV_boys,age_girls,eTIV_girls,study_boys,study_girls)
%Sequence of plots and stadistics
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
    
   
    %% Box-Whiskers plots 
    
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

function [str] = strengths_und(CIJ)
%STRENGTHS_UND        Strength
%
%   str = strengths_und(CIJ);
%
%   Node strength is the sum of weights of links connected to the node.
%
%   Input:      CIJ,    undirected weighted connection matrix
%
%   Output:     str,    node strength
%
%
%   Olaf Sporns, Indiana University, 2002/2006/2008

% compute strengths
str = sum(CIJ);        % strength
end

function [age_sex,eTIV_sex,study_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
    study_sex = sample.neuro_vals(find(sex),3);
end