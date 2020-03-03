function assignment=randomize_equal_cols(assignment,HKRR_energy_b,uEn)
%deal with identical columns, i.e. effectively identical RRs

for i=1:size(uEn,1)%loop over unique cols - i.e. rows of uEn
    cols=find(ismember(HKRR_energy_b',uEn(i,:),'rows'));
    %indices of the rows in HKRR_energy_b' that are equal to uEn(i,:)
    %ismember(HKRR_energy_b',uEn(i,:),'rows') returns a vector with 1 if the row is equal to uEn(i,:) and 0 if the row is different from uEn(i,:)
    if size(cols,1)>1
        %the column is repeated. Make a random permutation of these cols in assignment
        %first find which rows these columns are assigned to
        rows=find(ismember(assignment,cols));
        %then allow for a random permutation of these rows
        p=randperm(size(rows,1));
        assignment(rows)=assignment(rows(p));
    end
end

end

