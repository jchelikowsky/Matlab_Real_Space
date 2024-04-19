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
x=(x-n2*0.5-0.5)*domain.h;
y=(y-n2*0.5-0.5)*domain.h;
z=(z-n2*0.5-0.5)*domain.h;
[rho_length]=fread(wfnid, 1, 'uint32')

[rho]=fread(wfnid, rho_length, 'double');

V=reshape(rho,n2,n2,n2);

rho_plane(:,:)=(V(:,:,n2/2+1)+V(:,:,n2/2))*0.5;

rho_plot=rot90(rho_plane);

clf('reset')

% FigHandle = figure;
%  set(FigHandle, 'Position', [100, 100, 500, 500]);
hold on;
contour(x,y,rho_plot,8,'k')

[w_length]=fread(wfnid, 1, 'uint32' );
[nev]=fread(wfnid, 1, 'uint32');
for i = 1:nev
    [W]=fread(wfnid, w_length, 'double');
end

[N_types]=fread(wfnid, 1, 'uint32');
for i=1:N_types
    [na]=fread(wfnid, 1, 'uint32');
    for j=1:na
        xyz=fread(wfnid,3,"double");
        if (xyz(3) < domain.h * 0.5 && xyz(3) > -domain.h * 0.5)
            hold on;
            plot(xyz(1),xyz(2),'o','MarkerFaceColor','k','Markersize',4)
        end
    end
end

fclose(wfnid);