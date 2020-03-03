function assignment=randomize_equal_rows(assignment,HKRR_energy_b,uEn)
%deal with identical rows, i.e. effectively identical HKs

for i=1:size(uEn,1)%loop over unique rows
    rows=find(ismember(HKRR_energy_b,uEn(i,:),'rows'));
    %indices of the rows in HKRR_energy_b that are equal to uEn(i,:)
    %ismember(HKRR_energy_b,uEn(i,:),'rows') returns a vector with 1 if the row is equal to uEn(i,:) and 0 if the row is different from uEn(i,:)
    if size(rows,1)>1
        %the row is repeated. Make a random permutation of these rows in assignment
        p=randperm(size(rows,1));
        assignment(rows)=assignment(rows(p));
    end
end

end

