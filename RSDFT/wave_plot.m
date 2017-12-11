clear all;
clear all;
wfnid = fopen('./wfn.dat')


[nx]=fread(wfnid, 1, 'uint32')
[radius]=fread(wfnid, 1, 'double')
[h]=fread(wfnid, 1, 'double')

n2=nx;

[pot_length]= fread(wfnid, 1, 'uint32')


[pot]=fread(wfnid, pot_length, 'double');
x=[1:n2];
y=[1:n2];
z=[1:n2];


[rho_length]=fread(wfnid, 1, 'uint32')



[rho]=fread(wfnid, rho_length, 'double');



[w_length]= fread(wfnid,1, 'uint32' );
[nev]= fread(wfnid, 1, 'uint32');
display('number of states')
display(nev)
prompt = 'What state would you like to plot? ';
nstate = input(prompt)
prompt = 'How many atoms in the plot? ';
natms = input(prompt)
for natm=1:natms
    promptx='x-coordinate in a.u.';
    prompty='y-coordinate in a.u.';
    xatom(natm)=input(promptx);
    yatom(natm)=input(prompty);
end
for i = 1:nev
[W(:,i)]=fread(wfnid,w_length, 'double');
end

W_state=W(:,nstate);

%% check norm
w_norm=W_state' *  W_state
Wave=reshape(W_state,n2,n2,n2);
n_2=n2/2

for i=1:n2
    for j=1:n2
Wave_plot(j,i)=0.5*(Wave(i,j,n_2)+Wave(i,j,n_2+1));
  
    end
end
wmax2=max(max(Wave_plot))
wmin2=min(min(Wave_plot))
test=abs(Wave_plot);
wmax=max(max(test))

for i=1:n2
    for j=1:n2
if abs(Wave_plot(i,j)) < wmax/100
    Wave_plot(i,j) =0.;
end
    end
end

ncont=10
v=-wmax*(1:ncont)/ncont
v2=wmax*(1:ncont)/ncont
%contour(Wave_plot,10)
clf('reset')

sign=wmax2*wmin2
if sign > 0
        contour(Wave_plot,v2,'-k');
        hold on;
        contour(Wave_plot,v,'-k');
end
if sign < 0
     contour(Wave_plot,v2,'-k');
        hold on;
        contour(Wave_plot,v,'--k');
end

%% put atoms in plot
for natm=1:natms
    xx=xatom(natm)/h+n2/2+0.5
    yy=yatom(natm)/h+n2/2+0.5
hold on;
plot(xx,yy,'o','MarkerFaceColor','k','Markersize',15)
hold on;
end


%hold on;
%plot(42.35,36.5,'o','MarkerFaceColor','k','Markersize',15)





fclose(wfnid);