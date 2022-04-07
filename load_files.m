%Filter betweeen boys and Girls
load("..\plv_78ROIs\sexy_cn_312s_plv_78_rois_alpha.mat")
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array
boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix

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
         %Linear normalization of the weights
        sex_ROIs(:,:,i) = (sex_ROIs(:,:,i)- min(min(sex_ROIs(:,:,i))))/(max(max(sex_ROIs(:,:,i))) - min(min(sex_ROIs(:,:,i))));
    end   
end

%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);
