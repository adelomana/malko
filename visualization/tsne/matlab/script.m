%
% This script runs a tSNE on the expression data
% 

% 0. user defined variables
fileName='/Volumes/omics4tb/alomana/projects/mscni/data/single.cell.data.xlsx';
figureFile='figure.tsne.matlab.pdf';

% 1. read data
fileName='/Volumes/omics4tb/alomana/projects/mscni/data/single.cell.data.xlsx';
A=xlsread(fileName);
E=A(2:end,2:end);
labels=A(2:end,1);

% 2. analysis
% labels, dim, initDim, perplexity
%Y = tsne(E,[],2,50,10); 

Y = tsne(E); 
% 3. plot figure
gscatter(Y(:,1),Y(:,2),labels,'brgm');
%print(figureFile,'-dpdf');

%scatter3(Y(:,1),Y(:,2),Y(:,3),15,labels,'filled')

