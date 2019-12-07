function Vol=loadvol(datadir)

M=[];
for j=0:23
    M=[M iData([datadir '/qa_vs_qb_' num2str(j) '.dat'])];
end

Vol=cat(M,3);
slice(log(Vol));
