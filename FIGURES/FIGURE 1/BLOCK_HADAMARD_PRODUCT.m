function N_Matrix = BLOCK_HADAMARD_PRODUCT(A, B, n_vec) 
	K = length(n_vec);
    COL_temp = [];
	for k = 1 : K
		ROW_temp = [];
		for kp = 1 : K
			if kp == k
				ROW_temp = [ROW_temp, A(k, k) * eye(n_vec(k)) + B(k, k) * ones(n_vec(k), n_vec(k))];
            else
                ROW_temp = [ROW_temp, B(k, kp) * ones(n_vec(k), n_vec(kp))];
			end
		end
		COL_temp = [COL_temp; ROW_temp];
	end
	N_Matrix = COL_temp;
end