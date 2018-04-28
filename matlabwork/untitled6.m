load fisheriris

%rng default % for reproducibility
%Y = tsne(meas);
%gscatter(Y(:,1),Y(:,2),species)

% D = textscan(fileID,formatSpec,'HeaderLines',2,'Delimiter','\t')

fileName='/Users/alomana/Google Drive File Stream/My Drive/projects/mscni/data/expression/single cell data for Serdar and Adrian.txt';
f=fopen(fileName,'r');

format='%f %f %f';
D = textscan(f,format,'HeaderLines',1,'Delimiter','\t');
fclose(f);