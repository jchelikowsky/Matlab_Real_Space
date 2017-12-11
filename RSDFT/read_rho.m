clear all
wfnid = fopen('./wfn.dat')


[little_big_test]=fread(wfnid, 1, 'uint32')
[domain.radius]=fread(wfnid, 1, 'double')
[domain.h]=fread(wfnid, 1, 'double')

[pot_length]= fread(wfnid, 1, 'uint32')


[pot]=fread(wfnid, pot_length, 'double');
x=[1:52];
y=[1:52];
z=[1:52];


[rho_length]=fread(wfnid, 1, 'uint32')



[rho]=fread(wfnid, rho_length, 'double');

V=reshape(rho,52,52,52);

xm=max(rho)
fv=isosurface(x,y,z,V,xm*0.05);

p = patch(fv);
isonormals(x,y,z,V,p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1]);
view(3); axis tight
camlight 
lighting gouraud







fclose(wfnid);