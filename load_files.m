    %CHOOSE THE REPRESENTATION YOU WANT TO VISUALIZE BY UNCOMMENTING THE
    %CODE LINES
    %THEN, SIMPLY RUN THE CODE.
clear all;


%Generate the data
band = ["alpha"]; %You can try 'alpha' or 'beta'


type='rois';
% type='215nodes';
% type='1202nodes';




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

%Load the selected representation
files=dir(dum);
filename=horzcat(files.folder,'/',files.name);

load(filename);

%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array

%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);

%PLV matrix load
boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix

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

function [age_sex,eTIV_sex,study_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
    study_sex = sample.neuro_vals(find(sex),3);
end


