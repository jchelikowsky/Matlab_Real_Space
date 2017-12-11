function elemExists=elementExists(potentialElementName,elementData)
elementNames=cell(size(elementData.data,1),1);

for i=1:size(elementData.data,1)
    elementNames(i)=elementData.textdata(i);
end
elemExists=1;


if (any(strcmp(potentialElementName,elementNames))==0)
       elemExists=0;
end   
   

return;

