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


%%
%Variables needed

ROI_names = fcinfo.rois_names;
nodes2rois = fcinfo.nodes1202_in_rois;
nodes = fcinfo.nodes;
ROIs = fcinfo.rois_in;
boy_nodes = load_nodes(boys, fcmatrix); %Boys nodes matrix
girl_nodes = load_nodes(girls, fcmatrix); %Girls nodes matrix

names = {'ROIs_names', 'Cm'};
c = cell(length(names),1);
results_boys = cell2struct(c,names);
results_girls = cell2struct(c,names);

results_boys.ROIs_names = ROI_names;
results_girls.ROIs_names = ROI_names;
r=exp(0:0.13:5);
Cm_boys = zeros(5,length(r),length(find(boys)));
Cm_girls = zeros(5,length(r),length(find(girls)));
%%
%Boys
tic
for a=1:length(ROIs)
    X = fcinfo.nodespos_mm(nodes2rois{a},1);
    Y = fcinfo.nodespos_mm(nodes2rois{a},2);
    Z = fcinfo.nodespos_mm(nodes2rois{a},3);
    for i=1:length(find(boys))
        g_boys = graph(boy_nodes(nodes2rois{a},nodes2rois{a},i),'upper','omitselfloops');
        ad_boys = adjacency(g_boys,double(g_boys.Edges.Weight));
        
    %             Main cluster
        [bin,binsize] = conncomp(g_boys);
        idx = binsize(bin) == max(binsize);
        GCC_boys = subgraph(g_boys, idx);
        X = X(find(idx));
        Y = Y(find(idx));
        Z = Z(find(idx));
        adGCC_boys = adjacency(GCC_boys,double(GCC_boys.Edges.Weight));
    %     %Random Walker
        n_0 = 0;
        try
        n_0 = initial_node(adGCC_boys,n_0); %Aleatory initial node
        end
        N = 7*numnodes(GCC_boys); %Length of the RW
        try
        walk = rand_walk(adGCC_boys,N,n_0); %Random Walker initiation
        end
        try
        Cm_boys(:,:,i) = EmbDim(11,N,walk,X,Y,Z); %Computes the Corr. Sum 
        end
        %Compute the mean
    end
    results_boys.Cm(:,:,a) = mean(Cm_boys,3);
    
end
toc
save(sprintf('Cm_Intra-ROI/%s/Cm_boys.mat',band),"results_boys");

%Girls
tic
for a=1:length(ROIs)
    X = fcinfo.nodespos_mm(nodes2rois{a},1);
    Y = fcinfo.nodespos_mm(nodes2rois{a},2);
    Z = fcinfo.nodespos_mm(nodes2rois{a},3);
    for i=1:length(find(girls))
        g_girls = graph(girl_nodes(nodes2rois{a},nodes2rois{a},i),'upper','omitselfloops');
        ad_girls = adjacency(g_girls,double(g_girls.Edges.Weight));
        
    %             Main cluster
        [bin,binsize] = conncomp(g_girls);
        idx = binsize(bin) == max(binsize);
        GCC_girls = subgraph(g_girls, idx);
        X = X(find(idx));
        Y = Y(find(idx));
        Z = Z(find(idx));
        adGCC_girls = adjacency(GCC_girls,double(GCC_girls.Edges.Weight));
    %     %Random Walker
        n_0 = 0;
        try
        n_0 = initial_node(adGCC_girls,n_0); %Aleatory initial node
        end
        N = 7*numnodes(GCC_girls); %Length of the RW
        try
        walk = rand_walk(adGCC_girls,N,n_0); %Random Walker initiation
        end
        try
        Cm_girls(:,:,i) = EmbDim(11,N,walk,X,Y,Z); %Computes the Corr. Sum 
        end
        %Compute the mean
    end
    results_girls.Cm(:,:,a) = mean(Cm_girls,3);
    
end
toc
save(sprintf('Cm_Intra-ROI/%s/Cm_girls.mat',band),"results_girls");

%%
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

function n_0 = initial_node(AdGCC,n_0)
    [row_0,~] = find(AdGCC);
    n_0 = row_0(randi(length(row_0))); %Aleatory initial node
    
end
  
function Cm = EmbDim(Emb,N,walk,X,Y,Z)
%Algorithm to compute the Correlation Sum
    for m = 7:Emb
        j = 1;
        V = zeros(N-m,m);
        for i = 1:N-m
            V(i,:) = walk(i:i+m-1); %Serie temporal
        end
        for r=exp(0:0.13:5)
            drawnow;
            d = 0;
            heaviside = 0;
            for i = 1:N-m
                xi = X(V(i,:))';   yi = Y(V(i,:))'; zi = Z(V(i,:))';
                Xj = X(V(i+1:N-m,:)); Yj = Y(V(i+1:N-m,:)); Zj = Z(V(i+1:N-m,:));
                if i == N-m-1
                    Xj = Xj'; Yj = Yj'; Zj = Zj';
                end
                Xi = repmat(xi,size(Xj,1),1); Yi = repmat(yi,size(Xj,1),1);Zi = repmat(zi,size(Xj,1),1);
                A = abs(Xi-Xj); B = abs(Yi-Yj);D = abs(Zi-Zj);
                C = [A,B,D];
                aux = r-max(C');
                heaviside = heaviside+sum(aux >= 0);
            end
            Cm(m-6,j) = 2*heaviside/((N-m)*(N-m+1));
            j = j + 1;
        end
    end
end


function walk = rand_walk(R,N,n)
% R --> matriz de adyacencia
% N --> numero de elementos del random walker
% n --> nodo inicial

    connectivity = false;
    walk = zeros(1,N);
    
    while connectivity == false
    R_prob = binornd(1,R);
    R_prob = R_prob - diag(diag(R_prob));
    G = graph(R_prob,'omitselfloops','upper');
    
    %Proof that the graph is connected.
        if length(unique(conncomp(G)))== 1
            connectivity = true;
        else
            connectivity = false;
        end
    end
    %Random walker
    for i = 1:N
        neigh = find(R_prob(:,n));
        walk(i) = neigh(randi(length(neigh)));
        n = walk(i);
    end
end