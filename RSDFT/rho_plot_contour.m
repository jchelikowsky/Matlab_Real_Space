clear all
wfnid = fopen('./wfn.dat')


[nx]=fread(wfnid, 1, 'uint32')
[domain.radius]=fread(wfnid, 1, 'double')
[domain.h]=fread(wfnid, 1, 'double')

n2=nx;

[pot_length]= fread(wfnid, 1, 'uint32')


[pot]=fread(wfnid, pot_length, 'double');
x=[1:n2];
y=[1:n2];
z=[1:n2];


[rho_length]=fread(wfnid, 1, 'uint32')



[rho]=fread(wfnid, rho_length, 'double');

V=reshape(rho,n2,n2,n2);

rho_plane(:,:)=(V(:,:,n2/2+1)+V(:,:,n2/2-1))*0.5;

rho_plot=rot90(rho_plane);

clf('reset')


% FigHandle = figure;
%  set(FigHandle, 'Position', [100, 100, 500, 500]);
hold on;
contour(rho_plot,8,'k')
hold on;
plot(42.35,36.5,'o','MarkerFaceColor','k','Markersize',15)
hold on;
plot(30.65,36.5,'o','MarkerFaceColor','k','Markersize',15)






fclose(wfnid);