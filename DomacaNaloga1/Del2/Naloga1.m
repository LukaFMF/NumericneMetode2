% Imamo tocke b_i; i = 0,1,...,n (Bezierjeva krivulja stopnje n)
% Z zvisanjem stopnje zelimo dobiti tocke c_i; i = 0,1,...,n, n + 1 (Bezierjeva krivulja stopnje n + 1)
%    x  , y 
B = [
	-5/3,0;
	-16/3,4/3;
	-2,16/3;
	0,18/5;
	8/3,2
];

C = zvisajStopnjo(B);

vsotaAbscis = sum(C(:,1))
[najRazdalija,najInx] = max(vecnorm(C'))

function C = zvisajStopnjo(B)
	numPts = size(B,1);

	C = [B(1,:)];

	for i = 2:numPts
		fct = (i-1)/numPts;
		C = [C;fct*B(i-1,:) + (1 - fct)*B(i,:)];
	end
	C = [C;B(end,:)];
	size(B)
	size(C)
end