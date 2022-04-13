%Correlation dimension plvmatrix
clear all
%Generate the data
load("..\plv_78ROIs\sexy_cn_312s_plv_78_rois_alpha.mat")
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array

%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);

Cm_boys = load("Cm\Cm_boys.mat");
Cm_girls = load("Cm\Cm_girls.mat");

r=exp(-1:0.05:5);

%Fractal dimension and error per m
for i = 1:length(find(boys))
    for m = 1:5 %m = 4, 5, 6, 7, 8
        Cm = Cm_boys.Cm_boys(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_boys(m,i),delta_boys(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
end
for i = 1:length(find(girls))
    for m = 1:5 %m = 4, 5, 6, 7, 8
        Cm = Cm_girls.Cm_girls(m,:,i);
        %Filter the data
        Cm = smoothdata(Cm,'lowess',6); 
        [r_int,int] = new_filter(r,Cm);
        [frac_dim_girls(m,i),delta_girls(m,i)] = fractalfit(r_int,Cm,int); %Computes the corr.dimension as a function of m.
    end
end
%%
%Visualization of a man result
[r_int,int,A] = arrayfun(@(x) new_filter(r,Cm_boys.Cm_boys(x,:,1)),1:5,'UniformOutput',false);
figure();
hold on;
arrayfun(@(x) paint(r,Cm_boys.Cm_boys(x,:,1),int{x}),1:5)
legend('m=4','','m=5','','m=6','','m=7','', 'm=8','Fit interval','Location','Best');
title('Correlation sum of a man')
hold off;
%Visualization of the algorithim
figure();
plot(A{5})
hold on;
plot(int{5},A{5}(int{5}),'r-o')
legend('Cm differences (m=8)','Fit interval')
hold off;
%Computation of the fractal dimension 
[frac_dim,delta] = arrayfun(@(x) fractalfit(r_int{x},Cm_boys.Cm_boys(x,:,1),int{x}),1:5,'UniformOutput',false)
%\beta vs m
figure();
m = 4:8;
hold on;
arrayfun(@(x) errorbar(m(x),frac_dim{x},delta{x},'ko'), 1:5)
hold off;
xlabel('m');
ylabel('\beta'); 
legend('\beta');
hold off;
%%
%Visualization of a woman result
[r_int,int,A] = arrayfun(@(x) new_filter(r,Cm_girls.Cm_girls(x,:,1)),1:5,'UniformOutput',false);
figure();
hold on;
arrayfun(@(x) paint(r,Cm_girls.Cm_girls(x,:,1),int{x}),1:5)
legend('m=4','','m=5','','m=6','','m=7','', 'm=8','Fit interval','Location','Best');
title('Correlation sum of a woman')
hold off;
%Visualization of the algorithim
figure();
plot(A{5})
hold on;
plot(int{5},A{5}(int{5}),'r-o')
legend('Cm differences (m=8)','Fit interval')
hold off;
%Computation of the fractal dimension 
[frac_dim,delta] = arrayfun(@(x) fractalfit(r_int{x},Cm_girls.Cm_girls(x,:,1),int{x}),1:5,'UniformOutput',false)
%\beta vs m
figure();
m = 4:8;
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
bins = {'Boy','Girl'};
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
legend('Boys','Girls','Location',"best")
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
%     delta = sqrt(1/(length(r_int)-2)*(sum((log(Cm(1,int)) - mean(log(Cm(1,int)))).^2)))/sqrt(sum((r_int - mean(r_int)).^2));
    Sy = sqrt(sum((log(Cm(1,int)) - P_5(1)*r_int - P_5(2)).^2)/(length(r_int - 2)));
    delta = Sy*sqrt(length(r_int)/(length(r_int)*sum(r_int.^2)-sum(r_int)^2));
    
    %Linear fit statistics
    Rscore = 1 - (S.normr/norm(log(Cm(1,int)) - mean(log(Cm(1,int)))))^2;
%     sprintf('The Correlation Dimension is %.3f +/- %.3f, with R^2 of %.3f',P_5(1),delta,Rscore)


end