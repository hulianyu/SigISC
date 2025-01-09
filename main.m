clearvars
addpath([cd '/']);
addpath([cd '/Sequence2BinaryData']);
addpath([cd '/Evaluation']);
filename = char('activity','aslbu','auslan2','context','epitope','gene',...
    'news','pioneer','question','reuters','robot','skating','unix','webkb');
text_data_name = char('news', 'reuters', 'webkb');
IS = size(filename,1);
Results = zeros(IS,3);
eTime = zeros(IS,1);
numPi = zeros(IS,2);
Depths = zeros(IS,2);
for I = 1:IS
    load([strtrim(filename(I,:)), '_label', '.mat']); % Label
    GT = double(Label)';
    %% Sequential Pattern Discovery under Multiple Constraints: Binary_data
    load([strtrim(filename(I,:)), '_binary', '.mat']); % Binary_data
    X = double(Binary_data);
    %% Splitting
    [N,M] = size(X);
    objsID = 1:N;
    X = [objsID' X];
    pi_Node = zeros(N+1,1);
    h = 0;
    tic
    [Node, pi_Node] = SigISC(X,M,pi_Node,h);
    eTime(I,1) = toc;
    pi = pi_Node(1:end-1);
    %% Clustering quality
    results = ClusteringMeasure(GT, pi);
    Results(I,1) = results(3); % Purity
    Results(I,2) = results(2); % NMI
    Results(I,3) = results(7); % F-score
    disp(Results);
    %% Depths
    numPi(I,1) = length(unique(GT));
    numPi(I,2) = length(unique(pi));
    Depths(I,1) = treeDepth(Node);
    Depths(I,2) = averageLeafDepth(Node);
    drawTree_seq(Node, filename(I,:));
end