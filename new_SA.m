function [value] = new_SA(Clean,Seg,k,th)
[m,n] = size(Clean);
Clean_index = unique(Clean);
new_Clean = zeros(m,n);
new_Seg = zeros(m,n);
for i=1:k
    new_Clean(Clean == Clean_index(i))= i;
    if i == 1
        new_Seg(Seg<th(i)) = i;
    elseif i == k
        new_Seg(Seg>=th(i-1)) = i;
    else
        index1 = Seg >= th(i-1);
        index2 = Seg < th(i);
        new_Seg(logical(index1.*index2)) = i;
    end
end
SA_index = new_Seg == new_Clean;
value = sum(SA_index(:))/(m*n);
