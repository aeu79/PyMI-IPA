function Results = Predict_pairs(encoded_focus_alignment, PMIs, LengthA, table_count_species)
%makes pairing predictions on encoded_focus_alignment using PMI scores

[N, alignment_width] = size(encoded_focus_alignment);
L=alignment_width-2; % last 2 columns contain species index and initial sequence index

%initialize the Results array, used for saving data
Results=zeros(N-2,5);
%col 1: species
%col 2: HK index in initial alignment
%col 3: RR index in initial alignment
%col 4: score of pairing
%col 5: gap 

%total pair counter
totcount = 0; 


%loop over species
for i=1:size(table_count_species,1)
    
    test_seqs = encoded_focus_alignment(table_count_species(i,2):table_count_species(i,3),:);
    NSeqs = table_count_species(i,3)-table_count_species(i,2)+1;
    species_id=table_count_species(i,1);

    %now compute the PMI-based pairing score of all the HK-RR pairs within the species corresponding to i
    Pairing_scores = Compute_pairing_scores(test_seqs,NSeqs,PMIs, LengthA, L);
    if NSeqs==1
        assignment=1;
        Pairing_scores_b=Pairing_scores-min(min(Pairing_scores)); %ensure that all elements are >=0
    elseif isequal(min(Pairing_scores(:)),max(Pairing_scores(:)))
        assignment=randperm(NSeqs); %avoid spurious positive results
        Pairing_scores_b=Pairing_scores-min(min(Pairing_scores));
    else %use the Hungarian algorithm
        Pairing_scores_b=Pairing_scores-min(min(Pairing_scores));
        [assignment, score] = assignmentoptimal(Pairing_scores_b);
        %deal with identical rows, i.e. effectively identical HKs
        uEn = unique(Pairing_scores_b, 'rows');
        if size(uEn,1) < size(Pairing_scores_b,1)
            assignment=randomize_equal_rows(assignment,Pairing_scores_b,uEn);
        end
        %deal with identical cols, i.e. effectively identical RRs
        uEn = unique(Pairing_scores_b', 'rows'); %transpose to deal with columns
        if size(uEn,1) < size(Pairing_scores_b,2) %uEn is transposed
            assignment=randomize_equal_cols(assignment,Pairing_scores_b,uEn);
        end
    end
    
    bigval=1e3*abs(max(Pairing_scores_b(:)));
    
    for j=1:NSeqs
        totcount=totcount+1;
        Results(totcount,1)=species_id;
        Results(totcount,2)=test_seqs(j,L+2);  %initial index of HK sequence (HK: line)
        Results(totcount,3)=test_seqs(assignment(j),L+2); %initial index of RR sequence (RR: col)
        Results(totcount,4)=Pairing_scores(j,assignment(j)); %absolute score of the pairing
        if NSeqs==1
            Results(totcount,5)=abs(Pairing_scores); %no real gap... consider that absolute energy is gap
        elseif isequal(min(Pairing_scores(:)),max(Pairing_scores(:)))            
            Results(totcount,5)=0; %no gap for this assignment
        else
            %calculate gap for this assignment
            Pairing_scores_mod=Pairing_scores_b;
            Pairing_scores_mod(j,assignment(j))=bigval;
            [~, score_mod] = assignmentoptimal(Pairing_scores_mod);
            Results(totcount,5)=score_mod-score;
        end
    end
    
end

end
