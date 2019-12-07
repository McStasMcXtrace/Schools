function Vol=loadvol(datadir,maxslices)

M=[];
for j=0:maxslices
    M=[M iData([datadir '/qa_vs_qb_' num2str(j) '.dat'])];
end

Vol=cat(M,3);
slice(log(Vol));
