format long;
A = 18/53;
f = @(x) exp(A*x);
E = [0,1/2,1];
int = [0,1];

[a,b] = remesLin(f,E,int,1);

function [a,b] = remesLin(f,E,int,k)

	dim = size(E,2);

	for i = 1:k
		% matrika sistema
		col1 = (-ones(dim,1)).^((0:(dim - 1))');
		col2 = ones(dim,1);
		col3 = E';
		A = [col1,col2,col3];

		% desna stran
		b = f(E)';

		% za resitev matrike velikosti visine 3
		koef = A\b;
		m = koef(1);
		b = koef(2);
		a = koef(3);

		% aproksimacijski polinom
		p = @(x) a*x + b;

		residual = @(x) f(x) - p(x);

		% pripravimo na fminbnd
		resAbs = @(x) -abs(residual(x));

		y = fminbnd(resAbs,int(1),int(end));

		inx = -1;
		for j = 1:(dim - 1)
			if E(j) <= y && y <= E(j + 1)
				inx = j;
				break;
			end
		end

		% vsota koeficientov po prvem koraku
		razlikaKoefPoPrvem = b - a

		lowerR = residual(E(inx));
		%upperR = residual(E(inx + 1));
		residualY = residual(y);

		% ce sta predznaka enaka, bo produkt pozitiven
		if lowerR*residualY > 0
			E(j) = y;
		else
			E(j + 1) = y;
		end
	end
	napakaPoPrvem = abs(f(y) - p(y))
end