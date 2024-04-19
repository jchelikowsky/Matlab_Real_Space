clear all
wfnid = fopen('./wfn.dat')


[little_big_test]=fread(wfnid, 1, 'uint32')
[domain.radius]=fread(wfnid, 1, 'double')
[domain.h]=fread(wfnid, 1, 'double')

[pot_length]= fread(wfnid, 1, 'uint32')


[pot]=fread(wfnid, pot_length, 'double');
pot_length_cr = round(pot_length^(1/3))
x=[1:pot_length_cr];
y=[1:pot_length_cr];
z=[1:pot_length_cr];


[rho_length]=fread(wfnid, 1, 'uint32')



[rho]=fread(wfnid, rho_length, 'double');

V=reshape(rho,pot_length_cr,pot_length_cr,pot_length_cr);

xm=max(rho)
% C60 used fv=isosurface(x,y,z,V,xm*0.05);
fv=isosurface(x,y,z,V,xm*0.5);
p = patch(fv);
isonormals(x,y,z,V,p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1]);
view(3); axis tight
camlight 
lighting gouraud


fclose(wfnid);
