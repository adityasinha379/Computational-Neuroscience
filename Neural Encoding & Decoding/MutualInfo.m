function MI = MutualInfo(confusion)
MI=0;
confusion=confusion/size(confusion,1);
for i=1:size(confusion,1)
    for j=1:size(confusion,2)
        if(confusion(i,j)~=0)
            MI=MI+confusion(i,j)*log2(confusion(i,j)*size(confusion,1)/sum(confusion(:,j)));          %confusion matrix has entries of p(y/x)
        end
    end
end
end