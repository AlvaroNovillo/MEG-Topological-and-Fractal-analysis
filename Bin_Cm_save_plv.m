%Correlation dimension plvmatrix

%In this method, we create M=10 binary instances for each patient
%The corr. sum of each patient will be the mean of those M instances.

%Then, the fractality of each patient is normalized by N = 5 random weighted graphs
%with the same links and nodes, and the same distribution of weights.

clear all
%Generate the data

band = ["alpha"];
dum = sprintf('../plv_78ROIs/sexy_cn_312s_plv_78_rois_%s.mat',band);
files=dir(dum);
filename=horzcat(files.folder,'\',files.name);

load(filename) 
%Filter betweeen boys and Girls
boys = sample.neuro_vals(:,2) == 1; %Boys logical array
girls = sample.neuro_vals(:,2) == 2; %Girls logical array


boy_rois = load_ROIs(boys, fcmatrix); %Boys ROIs matrix
girl_rois = load_ROIs(girls, fcmatrix); %Girls ROIs matrix

%eTIV and age/sex classification 
[age_boys,eTIV_boys,study_boys] = neuro_sex_load(boys,sample);
[age_girls, eTIV_girls,study_girls] = neuro_sex_load(girls,sample);


%ROIs' position
X = aal78rois.rois_pos(:,1);
Y = aal78rois.rois_pos(:,2);
Z = aal78rois.rois_pos(:,3);

r=exp(0:0.05:5.2);
%%
%Boys
M = 5; %Binary instances 
tic
for i=1:length(find(boys))
    g_boys = graph(boy_rois(:,:,i),'upper','omitselfloops');
    ad_boys = adjacency(g_boys,'weighted');
    
%             Main cluster
    [bin,binsize] = conncomp(g_boys);
    idx = binsize(bin) == max(binsize);
    GCC_boys = subgraph(g_boys, idx);
    X = X(find(idx));
    Y = Y(find(idx));
    Z = Z(find(idx));
    adGCC_boys = adjacency(GCC_boys,'weighted');
%     %Random Walker
    for k = 1:M
        n_0 = 0;
        n_0 = initial_node(adGCC_boys,n_0); %Aleatory initial node
        N = 7*numnodes(GCC_boys); %Length of the RW
        walk = rand_walk(adGCC_boys,N,n_0); %Random Walker initiation
        Cm(:,:,k) = EmbDim(11,N,walk,X,Y,Z); %Computes the Corr. Sum 
    end
    %Compute the mean
    for j = 1:5
        Cm_boys(j,:,i) = arrayfun(@(x) mean(Cm(j,x,:)),1:length(r));
    end
end 
toc
save(sprintf('Cm_new/%s/Cm_boys.mat',band),"Cm_boys");

%Girls 
tic
for i=1:length(find(girls))
    g_girls = graph(girl_rois(:,:,i),'upper','omitselfloops');
    ad_girls = adjacency(g_girls,'weighted');
    
%             Main cluster
    [bin,binsize] = conncomp(g_girls);
    idx = binsize(bin) == max(binsize);
    GCC_girls = subgraph(g_girls, idx);
    X = X(find(idx));
    Y = Y(find(idx));
    Z = Z(find(idx));
    adGCC_girls = adjacency(GCC_girls,'weighted'); 
%     %Random Walker
    for k = 1:M
        n_0 = 0;
        n_0 = initial_node(adGCC_girls,n_0); %Aleatory initial node
        N = 7*numnodes(GCC_girls); %Length of the RW
        walk = rand_walk(adGCC_girls,N,n_0); %Random Walker initiation
        Cm(:,:,k) = EmbDim(11,N,walk,X,Y,Z); %Computes the Corr. Sum        
    end
    for j = 1:5
        Cm_girls(j,:,i) = arrayfun(@(x) mean(Cm(j,x,:)),1:length(r));
    end
end
toc

save(sprintf('Cm_new/%s/Cm_girls.mat',band),"Cm_girls");
%% Null model computation

%Boys
M = 5; %Binary instances 
tic
for i=1:length(find(boys))
    g_boys = graph(boy_rois(:,:,i),'upper','omitselfloops');
    ad_boys = adjacency(g_boys,'weighted');
    
%             Main cluster
    [bin,binsize] = conncomp(g_boys);
    idx = binsize(bin) == max(binsize);
    GCC_boys = subgraph(g_boys, idx);
    X = X(find(idx));
    Y = Y(find(idx));
    Z = Z(find(idx));
    adGCC_boys = adjacency(GCC_boys,'weighted');
%     %Random Walker
    for k = 1:M
        %For each network, we produce a random network with the same
        %properties
        n = GCC_boys.numnodes;
        e = GCC_boys.numedges;
        T=table(randi(n,e,2),abs(random('Normal',mean(GCC_boys.Edges.Weight), std(GCC_boys.Edges.Weight),1,length(GCC_boys.Edges.Weight))'),'VariableNames',{'EndNodes','Weights'});
        G=simplify(graph(T,'omitselfloops'));
        
        %Proof that the graphs are connected.
        if length(unique(conncomp(G)))== 1
            k = k;
        else
            k = k - 1;
        end
        ad_rand_b = adjacency(G,G.Edges.Weights);
        n_0 = 0;
        n_0 = initial_node(ad_rand_b,n_0); %Aleatory initial node
        N = 7*numnodes(G); %Length of the RW
        walk = rand_walk(ad_rand_b,N,n_0); %Random Walker initiation
        Cm(:,:,k) = EmbDim(11,N,walk,X,Y,Z); %Computes the Corr. Sum 
    end
    %Compute the mean
    for j = 1:5
        Cm_rand_b(j,:,i) = arrayfun(@(x) mean(Cm(j,x,:)),1:length(r));
    end
end 
toc
save(sprintf('Cm_new/%s/Cm_rand_b.mat',band),"Cm_rand_b");

%Girls 
tic
for i=137:length(find(girls))
    g_girls = graph(girl_rois(:,:,i),'upper','omitselfloops');
    ad_girls = adjacency(g_girls,'weighted');
    
%             Main cluster
    [bin,binsize] = conncomp(g_girls);
    idx = binsize(bin) == max(binsize);
    GCC_girls = subgraph(g_girls, idx);
    X = X(find(idx));
    Y = Y(find(idx));
    Z = Z(find(idx));
    adGCC_girls = adjacency(GCC_girls,'weighted'); 
%     %Random Walker
    for k = 1:M
        %For each network, we produce a random network with tbe same
        %properties
        n = GCC_girls.numnodes;
        e = GCC_girls.numedges;
        T=table(randi(n,e,2),abs(random('Normal',mean(GCC_girls.Edges.Weight), std(GCC_girls.Edges.Weight),1,length(GCC_girls.Edges.Weight))'),'VariableNames',{'EndNodes','Weights'});
        G=simplify(graph(T,'omitselfloops'));
        %Proof that the graphs are connected.
        if length(unique(conncomp(G)))== 1
            k = k;
        else
            k = k - 1;
        end
        ad_rand_g = adjacency(G,G.Edges.Weights);
        n_0 = 0;
        n_0 = initial_node(ad_rand_g,n_0); %Aleatory initial node
        N = 7*numnodes(G); %Length of the RW
        walk = rand_walk(ad_rand_g,N,n_0); %Random Walker initiation
        Cm(:,:,k) = EmbDim(11,N,walk,X,Y,Z); %Computes the Corr. Sum        
    end
    for j = 1:5
        Cm_rand_g(j,:,i) = arrayfun(@(x) mean(Cm(j,x,:)),1:length(r));
    end
end
toc

save(sprintf('Cm_new/%s/Cm_rand_g.mat',band),"Cm_rand_g");

%%
% %Plot Cm vs r
% %Plot Cm(r) vs r
% [r_int,int,A] = arrayfun(@(x) new_filter(r,Cm_boys(x,:)),1:5,'UniformOutput',false);
% figure();
% hold on;
% arrayfun(@(x) paint(r,Cm_boys(x,:),int{x}),1:5)
% legend('m=4','','m=5','','m=6','','m=7','', 'm=8','Fit interval','Location','Best');
% hold off;
% %Visualization of the algorithim
% figure();
% plot(A{5})
% hold on;
% plot(int{5},A{5}(int{5}),'r-o')
% legend('Cm differences (m=8)','Fit interval')
% hold off;
% %Computation of the fractal dimension 
% [frac_dim,delta] = arrayfun(@(x) fractalfit(r_int{x},Cm(x,:),int{x}),1:5,'UniformOutput',false)
% %\beta vs m
% figure();
% m = 7:11;
% hold on;
% arrayfun(@(x) errorbar(m(x),frac_dim{x},delta{x},'ko'), 1:5)
% hold off;
% xlabel('m');
% ylabel('\beta'); 
% legend('\beta');
% hold off;
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


function [age_sex,eTIV_sex,study_sex] = neuro_sex_load(sex,sample)
    age_sex = sample.neuro_vals(find(sex),1);
    eTIV_sex = sample.vols_vals(find(sex),1);
    study_sex = sample.neuro_vals(find(sex),3);
end

function n_0 = initial_node(AdGCC,n_0)
    AdGCC = binornd(1,AdGCC);
    AdGCC = AdGCC - diag(diag(AdGCC));
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
        for r=exp(0:0.05:5.2)
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
        neigh = find (R_prob(:,n));
        walk(i) = neigh(randi(length(neigh)));
        n = walk(i);
    end
end

function paint(r,Cm,int)
    %Representation of Cm(r) vs r
    loglog(r,Cm(1,1:end),'-o','MarkerSize',4)
    loglog(r(int), Cm(1,int), 'b-o','MarkerSize',4)
    xlabel('r');
    ylabel('C_m(r)');
    title('Correlation sum');
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