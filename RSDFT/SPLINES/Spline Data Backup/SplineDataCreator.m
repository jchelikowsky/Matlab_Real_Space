Al;
Ar;
B;
Be;
C;
Cl;
F;
H;
He;
Li;
Mg;
N;
Na;
Ne;
O;
P;
S;
Si;
%call script to load data into a variable

save('splineData.mat','dataAl','dataAr','dataB','dataBe','dataC','dataCl');% no append on first save
% because file must be created first
save('splineData.mat','dataF','dataH','dataHe','dataLi','dataMg','dataN','-append');
save('splineData.mat','dataNa','dataNe','dataO','dataP','dataS','dataSi','-append');
%call save with the file splineData.mat and the name of the variable,
%don't forget the -append flag
