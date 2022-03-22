%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array
boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix

%Ahora, calcular el clustering y el shortest path con 
%Las funciones de BTC
%%
function sex_ROIs =load_ROIs(sex, fcmatrix)
%%%
    %Load the ROI's matrix and format them
    %Arguments: sex -> Logical array
    %fcmatrix: 4D array (Database of MEG)
%%%
    %Load
    load("plv_78ROIs\sexy_cn_312s_plv_78_rois_alpha.mat")
    sex_ROIs = zeros(78,78,length(find(sex)));
    sex_ROIs = fcmatrix(:,:,find(sex));
    %Format the matrix (Remove the diagonal)
    for i = 1:length(find(sex))
        sex_ROIs(:,:,i) = sex_ROIs(:,:,i) - diag(diag(sex_ROIs(:,:,i)));
    end   
end