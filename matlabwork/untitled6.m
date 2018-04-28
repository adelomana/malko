load fisheriris

%rng default % for reproducibility
%Y = tsne(meas);
%gscatter(Y(:,1),Y(:,2),species)

fileName='/Users/alomana/Google Drive File Stream/My Drive/projects/mscni/data/expression/test.txt';
f=fopen(fileName,'r');
T = readtable(f);
fclose(f);