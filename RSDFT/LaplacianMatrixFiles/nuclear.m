function [enuc]=nuclear(Domain, Atoms)

% Find ion-ion core repulsion
 

indx = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_types = length(Atoms);
for at_typ=1:N_types
    typ = Atoms(at_typ).typ;
    elem=importdata('elements_new.csv');
    for i=1:length(elem.data)
        if strcmp(typ,elem.textdata{i})
            index=i;
            break
        end
    end
    xyz=Atoms(at_typ).coord;
    %retrieve corresponding information from elements.csv
    Z=elem.data(index,2);
    natoms = size(xyz,1);
%%%%
    for at = 1: natoms
%% scan all points for each atom (not optimal but OK)--copied from
%% ppot.m with a small change of formula to find potential
       indx = indx+1; 
        xx(indx)=xyz(at,3);
        yy(indx)=xyz(at,2);
        zz(indx)=xyz(at,1);
		zt(indx)=Z;		 
    end
end

%%%
%%%  Nuclear repulsion term --  For now assume homopolar molecule
%%%

enuc=0;
for at=1:indx
    xx1=xx(at);
    yy1=yy(at);
    zz1=zz(at);
    Z1=zt(at);
    
    atomIndex2=1:indx;
    xx2=xx(atomIndex2);
    yy2=yy(atomIndex2);
    zz2=zz(atomIndex2);
    Z2=zt(atomIndex2);
    d2=(xx1-xx2).^2+(yy1-yy2).^2+(zz1-zz2).^2;
    Z2=Z2(d2~=0);%these 2 lines remove the 0 elements to avoid a division by 0 on line 54
    d2=d2(d2~=0);
    R2=sqrt(d2);
    enuc=enuc+sum(Z1*Z2./R2);
end

%%%  
%%%  Fix double counting, but convert to Ry--no factor of two 
%%%
%%%

