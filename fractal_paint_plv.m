% %Run from here if you have already downloaded/generated the results
% 
%Generate the data
clear all;
band = ["alpha"]; %You can try 'alpha' or 'beta'
gband='\alpha ';%You can try '\alpha' or '\beta';

% type='rois'; NOT IMPLEMENTED YET!
type='215nodes';
% type='1202nodes'; NOT IMPLEMENTED YET!
%% 78X78 REPRESENTATION

if(strcmp(type,'rois')==1)
dum = sprintf('../plv_78ROIs/sexy_cn_312s_plv_78_%s_%s.mat',type,band);
%% 215X215 REPRESENTATION

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


%%
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array

%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);

Cm_boys = results_boys.Cm;
Cm_girls = results_girls.Cm;

Cm_null_boys = results_boys.Cm_rand;
Cm_null_girls = results_girls.Cm_rand;

r=exp(0:0.07:5.2);

% Fractal dimension and error per m
for i = 1:length(find(boys))
    for m = 1:5 %m = 7, 8, 9, 10, 11
        Cm = Cm_boys(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_boys(m,i),delta_boys(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
end

for i = 1:length(find(boys))
    for m = 1:5 %m = 7, 8, 9, 10, 11
        Cm = Cm_null_boys(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_null_boys(m,i),delta_null_boys(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
end

for i = 1:length(find(girls))
    for m = 1:5 %m = 7, 8, 9, 10, 11
        Cm = Cm_girls(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_girls(m,i),delta_girls(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
end

for i = 1:length(find(girls))
    for m = 1:5 %m = 7, 8, 9, 10, 11
        Cm = Cm_null_girls(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_null_girls(m,i),delta_null_girls(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
end
%%
%%% Null model 
frac_dim_boys = frac_dim_boys./frac_dim_null_boys ;
frac_dim_girls = frac_dim_girls./frac_dim_null_girls;
delta_boys = delta_boys./delta_null_boys;
delta_girls = delta_girls./delta_null_girls;
%%
%Visualization of a man result
a = 106; %Patient 
[r_int,int,A] = arrayfun(@(x) new_filter(r,smoothdata(Cm_boys(x,:,a),'lowess',6)),1:5,'UniformOutput',false);
figure();
hold on;
arrayfun(@(x) paint(r,smoothdata(smoothdata(Cm_boys(x,:,a),'lowess',6),'lowess',6),int{x}),1:5)
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
[frac_dim,delta] = arrayfun(@(x) fractalfit(r_int{x},smoothdata(Cm_boys(x,:,a),'lowess',6),int{x}),1:5,'UniformOutput',false)
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
a = 1;
[r_int,int,A] = arrayfun(@(x) new_filter(r,smoothdata(Cm_girls(x,:,a),'lowess',6)),1:5,'UniformOutput',false);
figure();
hold on;
arrayfun(@(x) paint(r,smoothdata(Cm_girls(x,:,a),'lowess',6),int{x}),1:5)
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
[frac_dim,delta] = arrayfun(@(x) fractalfit(r_int{x},smoothdata(Cm_girls(x,:,a),'lowess',6),int{x}),1:5,'UniformOutput',false)
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


%frac_dim per age
figure();
ages = [age_boys;age_girls];
binEdges = 50:10:90;
bins = {'50s','60s','70s','80s'};
groupAge = discretize(ages,binEdges,'categorical',bins);
boxchart(groupAge,frac_box,'GroupByColor',groupSex)
title('Corr. dimension vs age')
legend('Location',"best")
ylabel('\beta')
xlabel('age')

%frac_dim per study_year

figure();
study_period = [study_boys;study_girls];
study_period(isnan(study_period)) = 0;
binEdges = 0:5:35;
bins = {'-5','5','10' ,'15','20','25','+30'};
groupStudy = discretize(study_period,binEdges,'categorical',bins);
boxchart(groupStudy,frac_box,'GroupByColor',groupSex)
title('Corr. dimension vs Study period')
legend('Location',"best")
ylabel('\beta')
xlabel('Study period')

%frac_dim vs eTIV

figure();
plot(eTIV_boys,frac_dim_boys(5,:),'o')
hold on;
plot(eTIV_girls,frac_dim_girls(5,:),'s')
title('Corr. dimension vs eTIV')
legend('Man','Woman','Location',"best")
ylabel('\beta')
xlabel('eTIV')
hold off;

%%
function [age_sex,eTIV_sex,study_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
    study_sex = sample.neuro_vals(find(sex),3);
end
function paint(r,Cm,int)
    %Representation of Cm(r) vs r
    loglog(r,Cm(1,1:end),'-o','MarkerSize',4)
    loglog(r(int), Cm(1,int), 'b-o','MarkerSize',4)
    xlabel('r');
    ylabel('C_m(r)');
end
function [r_int,int,A] = new_filter(r,Cm)
    A = diff(Cm);
    index = (A>= max(A)/2); %Threshold
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