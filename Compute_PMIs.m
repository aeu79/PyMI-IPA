function [PMIs, Meff] = Compute_PMIs(encoded_focus_alignment, pseudocount_weight, theta)

q=21;

%compute the empirical frequencies from the sequence data
[Meff, Pij_true, Pi_true, alignment_width] = count_alignment(encoded_focus_alignment, theta, q);

%include pseudocounts
[Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q);

%compute PMIs
PMIs = Get_PMIs(Pij, Pi, alignment_width, q); 

end



%%auxiliary function definitions

function [Meff, Pij_true, Pi_true, alignment_width] = count_alignment(encoded_focus_alignment, theta, q)

%calculate Meff
[alignment_height,alignment_width] = size(encoded_focus_alignment);
W = ones(1, alignment_height);
if(theta > 0.0) % whether you weight or not
    W = (1./(1+sum(squareform(pdist(encoded_focus_alignment, 'hamm')<theta))));
end
Meff=sum(W);

%compute the frequencies
Pij_true = zeros(alignment_width, alignment_width, q, q); %a 4-dim matrix (equivalent to a q^2*alignment_width^2 square matrix)
Pi_true = zeros(alignment_width, q);

%single-site frequencies
for j=1:alignment_height %sequence index in the sequence alignment
	for i=1:alignment_width %site index in the sequence alignment
		Pi_true(i, encoded_focus_alignment(j, i)) = Pi_true(i, encoded_focus_alignment(j, i)) + W(j); %increment the proba to have this residue at position i by the weight W(j) of the sequence j considered
	end
end
Pi_true = Pi_true/Meff; %normalization

%two-site frequencies
for l=1:alignment_height %sequence index
	for i=1:alignment_width-1 %site index
		for j=i+1:alignment_width %site index
			Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) + W(l);
			Pij_true(j, i, encoded_focus_alignment(l, j), encoded_focus_alignment(l, i)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j));
		end
	end
end
Pij_true = Pij_true/Meff;

%fill the diagonal of Pij_true
scra = eye(q, q); %has to be the same amino acid
for i=1:alignment_width
	for alpha=1:q
		for beta=1:q
			Pij_true(i, i, alpha, beta) = Pi_true(i, alpha) * scra(alpha, beta);
		end
	end
end

end


function [Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q)
%add pseudocounts to deal with some finite size effects

Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(alignment_width, alignment_width, q, q);
Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(alignment_width, q);

%correct things on the diagonal
scra = eye(q);
for i=1:alignment_width
	for alpha = 2:q
		for beta = 2:q
			Pij(i, i, alpha, beta) = (1.-pseudocount_weight)*Pij_true(i, i, alpha, beta) + pseudocount_weight/q*scra(alpha, beta);
		end
	end
end

end



function C = Get_PMIs(Pij, Pi, alignment_width, q)
%Compute PMIs from the matrices of frequencies. 
%Here, the 1st aa type(=gap) is ignored (slight improvement)

C=zeros(alignment_width,alignment_width,q,q);
for i=1:alignment_width
	for j=1:alignment_width
		for alpha=2:q
			for beta=2:q
                C(i,j,alpha,beta) = log(Pij(i, j, alpha, beta)/(Pi(i, alpha)*Pi(j, beta)));
			end
		end
	end
end

end
