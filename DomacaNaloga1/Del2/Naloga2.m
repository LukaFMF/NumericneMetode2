% za tocke p_0,p_1,...,p_m,p_(m+1) (m + 2 tock) generiramo
% novo zaporedje, ki ga sestavlja 2^k*m + 2 tock
P = [
	4 2; 
	1 4; 
	-2 4; 
	-4 -1;
	-2 -3;
	2 -2
];

P_2 = chaikGen(P,2);
razdalijaMedTockama = norm(P_2(1,:) - P_2(18,:),2)

function P_k = chaikGen(P,k)
	m = size(P,1) - 2;

	P_p = P;
	for l = 1:k
		P_k = zeros(2^l*m + 2,2);

		for j = 0:(2^(l - 1)*m)
			P_k(2*j + 1,:) = 3/4*P_p(j + 1,:) + 1/4*P_p(j + 2,:); 
			P_k(2*j + 2,:) = 1/4*P_p(j + 1,:) + 3/4*P_p(j + 2,:);
		end
		P_p = P_k;
	end
	P_k = P_p;
end