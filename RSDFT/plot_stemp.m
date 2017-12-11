
XREF=EM(120)
for j=1:150
   X1(j)=EQ(j)
%X1(j)=EM(j)-XREF
    Y1(j)=1.0 
end
stem(X1,Y1)