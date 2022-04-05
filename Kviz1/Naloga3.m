% Imamo tocke b_i; i = 0,1,...,n (Bezierjeva krivulja stopnje n)
% Z zvisanjem stopnje zelimo dobiti tocke c_i; i = 0,1,...,n, n + 1 (Bezierjeva krivulja stopnje n + 1)
%    x  , y
c = 106/97; 
B = [
	-2,1/2;
	-16/3,c;
	-2,16/3;
	0,18/5
];

C = zvisajStopnjo(B);
D = zvisajStopnjo(C);

vsotaOrdinat = sum(C(:,2))
d_1 = D(2,:); % d_0 = D(1,:)
oddaljenostOdTocke = norm(d_1 - [-4,4],2)

function C = zvisajStopnjo(B)
	numPts = size(B,1);

	C = [B(1,:)];

	for i = 2:numPts
		fct = (i-1)/numPts;
		C = [C;fct*B(i-1,:) + (1 - fct)*B(i,:)];
	end
	C = [C;B(end,:)];
end