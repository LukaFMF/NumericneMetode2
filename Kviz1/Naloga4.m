% za tocke p_0,p_1,...,p_m,p_(m+1) (m + 2 tock) generiramo
% novo zaporedje, ki ga sestavlja 2^k*m + 2 tock
d = 82/67;
P = [
	5 2; 
	1 4; 
	d d; 
	-2 4;
	-4 -2;
	-2*d, -3*d;
	2 -2
];

P_1 = chaikGen(P,1);
razdalijaMedTockama = norm(P_1(11,:) - P_1(12,:),2)

P_2 = chaikGen(P,2);
dolzLomljenke = dolzinaLomljenke(P_2)

function dolz = dolzinaLomljenke(P)
	numPts = size(P,1);

	dolz = 0;
	for i = 1:(numPts - 1)
		dolz = dolz + norm(P(i,:) - P(i + 1,:),2);
	end
end

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