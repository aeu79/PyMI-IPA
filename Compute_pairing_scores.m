function Pairing_scores = Compute_pairing_scores(test_seqs,NSeqs,PMIs, LengthA, L)
%calculate PMI-based pairing scores between all pairs of HKs and RRs in test_seqs
%Line: HK; Col: RR.

Pairing_scores = zeros(NSeqs,NSeqs);

for i = 1:NSeqs %to choose the HK
    for j = 1:NSeqs %to choose the RR
        for a = 1:LengthA %sites in HK
            for b = LengthA+1:L %sites in RR
                aa1=test_seqs(i,a); %aa in HK i at site a
                aa2=test_seqs(j,b); %aa in RR j at site b
                Pairing_scores(i,j) = Pairing_scores(i,j) + PMIs(a,b,aa1,aa2);
            end
        end
    end
end

end
