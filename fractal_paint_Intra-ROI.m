%Intra-ROI correlation dimension.
clear all
%Generate the data

band = ["alpha"];
dum = sprintf('../plv_nodes/sexy_cn_312s_plv_1202nodes_%s.mat',band);
files=dir(dum);
filename=horzcat(files.folder,'\',files.name);

load(filename);

%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array


load(sprintf("Cm_Intra-ROI/%s/Cm_boys.mat",band));
load(sprintf("Cm_Intra-ROI/%s/Cm_girls.mat",band));
%%
%Previous variables
ROIs = results_boys.ROIs_names;
r=exp(0:0.13:5);
%%

%Fractal dimension and error of each ROI per m
for i = 1:length(ROIs)
    for m = 1:5 %m = 7, 8, 9, 10, 11
        Cm = results_boys.Cm(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_boys(m,i),delta_boys(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
end


for i = 1:length(ROIs)
    for m = 1:5 %m = 7, 8, 9, 10, 11
        Cm = results_girls.Cm(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_girls(m,i),delta_girls(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
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
%Visualization of a man result
a = 2; %ROI 
[r_int,int,A] = arrayfun(@(x) new_filter(r,smoothdata(results_boys.Cm(x,:,a),'lowess',6)),1:5,'UniformOutput',false);
figure();
hold on;
arrayfun(@(x) paint(r,smoothdata(results_boys.Cm(x,:,a),'lowess',6),int{x}),1:5)
legend('m=7','','m=8','','m=9','','m=10','', 'm=11','Fit interval','Location','Best');
title('Correlation sum of a man')
hold off;
%Visualization of the algorithim
figure();
plot(A{5})
hold on;
plot(int{5},A{5}(int{5}),'r-o')
legend('Cm differences (m=11)','Fit interval')
hold off;
%Computation of the fractal dimension 
[frac_dim,delta] = arrayfun(@(x) fractalfit(r_int{x},smoothdata(results_boys.Cm(x,:,a),'lowess',6),int{x}),1:5,'UniformOutput',false)
%\beta vs m
figure();
m = 7:11;
hold on;
arrayfun(@(x) errorbar(m(x),frac_dim{x},delta{x},'ko'), 1:5)
hold off;
xlabel('m');
ylabel('\beta'); 
legend('\beta');
hold off;
%%
%Visualization of a woman result
[r_int,int,A] = arrayfun(@(x) new_filter(r,smoothdata(results_girls.Cm(x,:,a),'lowess',6)),1:5,'UniformOutput',false);
figure();
hold on;
arrayfun(@(x) paint(r,smoothdata(results_girls.Cm(x,:,a),'lowess',6),int{x}),1:5)
legend('m=7','','m=8','','m=9','','m=10','', 'm=11','Fit interval','Location','Best');
title('Correlation sum of a woman')
hold off;
%Visualization of the algorithim
figure();
plot(A{5})
hold on;
plot(int{5},A{5}(int{5}),'r-o')
legend('Cm differences (m=11)','Fit interval')
hold off;
%Computation of the fractal dimension 
[frac_dim,delta] = arrayfun(@(x) fractalfit(r_int{x},smoothdata(results_girls.Cm(x,:,a),'lowess',6),int{x}),1:5,'UniformOutput',false)
%\beta vs m
figure();
m = 7:11;
hold on;
arrayfun(@(x) errorbar(m(x),frac_dim{x},delta{x},'ko'), 1:5)
hold off;
xlabel('m');
ylabel('\beta'); 
legend('\beta');
hold off;
%%
%Figures 

%frac_dim
frac_box = [frac_dim_boys(5,:)';frac_dim_girls(5,:)'];
g = [ones(size(frac_dim_boys(5,:)')); 2*ones(size(frac_dim_girls(5,:)'))]; %Boys and girls class.
binEdges = 1:3;
bins = {'Man','Woman'};
groupSex = discretize(g,binEdges,'categorical',bins);

figure();
hold on;
boxchart(groupSex,frac_box,'GroupByColor',groupSex)
ylabel('\beta')
xlabel('gender')
legend('Location',"best")
hold off;
%%
% Correlation dimension vs ROIs
figure();
hold on;
arrayfun(@(x) errorbar(1:78,frac_dim_boys(x,:),delta_boys(x,:),'-o'), 5)
arrayfun(@(x) errorbar(1:78,frac_dim_girls(x,:),delta_girls(x,:),'-s'), 5)
xlabel('ROIs')
ylabel('\beta')
legend('Man','Woman')
hold off;
%Connectivity vs Correlation dimension 
figure();
hold on;
arrayfun(@(x) errorbar(frac_dim_boys(x,:),results_boys.Connect,results_boys.ConnErr,'o'), 5)
arrayfun(@(x) errorbar(frac_dim_girls(x,:),results_girls.Connect,results_girls.ConnErr,'s'), 5)
ylabel('Connectivity')
xlabel('\beta')
legend('Man','Woman')
hold off;
%%
%Lets see if there is any relation betweeen topology and geometry
topo_boys = load(sprintf('../Clustering and Shortest Path Length/Intra_ROI_results/%s/boys.mat',band));
topo_girls = load(sprintf('../Clustering and Shortest Path Length/Intra_ROI_results/%s/girls.mat',band));
%%
%Clustering vs Corr dim

figure();
hold on;
arrayfun(@(x) errorbar(frac_dim_boys(x,:),topo_boys.results_boys.Clustering,topo_boys.results_boys.C_error,'o'), 5)
arrayfun(@(x) errorbar(frac_dim_girls(x,:),topo_girls.results_girls.Clustering,topo_girls.results_girls.C_error,'s'), 5)
ylabel('Clustering')
xlabel('\beta')
legend('Man','Woman')
hold off;

%dist vs Corr dim

figure();
hold on;
arrayfun(@(x) errorbar(frac_dim_boys(x,:),topo_boys.results_boys.Shortest_Path_Length,topo_boys.results_boys.Dist_error,'o'), 5)
arrayfun(@(x) errorbar(frac_dim_girls(x,:),topo_girls.results_girls.Shortest_Path_Length,topo_girls.results_girls.Dist_error,'s'), 5)
ylabel('Shortest Path Length')
xlabel('\beta')
legend('Man','Woman')
hold off;
%%
%Box plots 

%Clustering vs Corr dim

%C 
C_box = [topo_boys.results_boys.Clustering';topo_girls.results_girls.Clustering'];


figure();
frac = [frac_dim_boys(5,:),frac_dim_girls(5,:)]';
binEdges = 2:7;
bins = {'2','3','4','5','6'};
groupFrac = discretize(frac,binEdges,'categorical',bins);
boxchart(groupFrac,C_box,'GroupByColor',groupSex)
title('Clustering vs Correlation Dimension')
legend('Location',"best")
ylabel('Clustering')
xlabel('\beta')

%Shortest path vs Corr Dim

%dist
dist_box = [topo_boys.results_boys.Shortest_Path_Length';topo_girls.results_girls.Shortest_Path_Length'];
figure();
boxchart(groupFrac,dist_box,'GroupByColor',groupSex)
title('Shortest Path Length vs Correlation Dimension')
legend('Location',"best")
ylabel('Shortest Path Length')
xlabel('\beta')


%Connectivity vs Corr Dim

conn_box = [results_boys.Connect';results_girls.Connect'];
figure();
boxchart(groupFrac,conn_box,'GroupByColor',groupSex)
title('Connectivity vs Correlation Dimension')
legend('Location',"best")
ylabel('Connectivity')
xlabel('\beta')

%%
function paint(r,Cm,int)
    %Representation of Cm(r) vs r
    loglog(r,Cm(1,1:end),'-o','MarkerSize',4)
    loglog(r(int), Cm(1,int), 'b-o','MarkerSize',4)
    xlabel('r');
    ylabel('C_m(r)');
end
function [r_int,int,A] = new_filter(r,Cm)
    A = diff(Cm);
    index = (A>= max(A)/4); %Threshold
    gt=find(index~=0);
    lower = min(gt);
    upper = max(gt);
    int = lower:upper; %Interval where the fit is performed
    r_int = log(r(int)); %Differential section of r where the fit is performed 

end

function [frac_dim,delta] = fractalfit(r_int,Cm,int)
    %Computation of the fit
    [P_5,S] = polyfit(r_int,log(Cm(1,int)),1);
    frac_dim = P_5(1);
    %Estimation of the standard error of the slope
    Sy = sqrt(sum((log(Cm(1,int)) - P_5(1)*r_int - P_5(2)).^2)/(length(r_int)- 2));
    delta = Sy*sqrt(length(r_int)/(length(r_int)*sum(r_int.^2)-sum(r_int)^2));
    
    %Linear fit statistics
    Rscore = 1 - (S.normr/norm(log(Cm(1,int)) - mean(log(Cm(1,int)))))^2;
%     sprintf('The Correlation Dimension is %.3f +/- %.3f, with R^2 of %.3f',P_5(1),delta,Rscore)


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