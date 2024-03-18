function [A_PARA, B_PARA] = BEST_UNBIASED_ESTIMATOR(S, n_vector)
	K = length(n_vector);
	A_PARA = zeros(K);
	B_PARA = zeros(K);
	for k = 1 : K
		for kp = 1 : K
			row_index = [(sum(n_vector(1 : (k - 1))) + 1) : sum(n_vector(1 : k))];
			col_index = [(sum(n_vector(1 : (kp - 1))) + 1) : sum(n_vector(1 : kp))];
			SUB_matrix = S(row_index, col_index);
			if kp == k
				SUB_on_diag = mean(diag(SUB_matrix));
				SUB_off_diag = (sum(SUB_matrix, 'all') - sum(diag(SUB_matrix))) / (n_vector(k) * n_vector(kp) - n_vector(k));
				A_PARA(k, kp) = SUB_on_diag - SUB_off_diag;
				B_PARA(k, kp) = SUB_off_diag;
			else
				B_PARA(k, kp) = sum(SUB_matrix, 'all') / (n_vector(k) * n_vector(kp));
			end
		end
    end
    B_PARA = (B_PARA + B_PARA') / 2;
end

