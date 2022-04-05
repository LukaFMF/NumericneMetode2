% % Na vajah smo na tablo naredili en korak Remesovega postopka za funkcijo
% % f(x) = sin(x) na intervalu [0,pi] z začetno množico točk E1 = [0,pi/3,(2*pi)/3,pi]
% % Drugi korak postopka izvedemo v Matlabu.
% f = @(x) sin(x);

% % Določimo množico E2:
% % E2 = [0,pi/3,2.7327,pi];
% E2 = [0,pi/3,2*pi/3,pi];

% % Določimo matriko sistema enačb
% A = [1 1 0 0;
%      -1 1 pi/3 pi*pi/9;
%      1 1 2.7327 2.7327^2;
%      -1 1 pi pi^2];
 
% % Določimo desno stran sistema enačb
% b = [f(E2(1));f(E2(2));f(E2(3));f(E2(4))];

% % Rešimo sistem enačb in dobimo vektor koeficientov: [m2;c2;b2;a2]
% koef = A\b;
% m2 = koef(1);
% c2 = koef(2);
% b2 = koef(3);
% a2 = koef(4);

% % Določimo aproksimacijski polinom
% p2 = @(x) a2*x.^2 + b2*x + c2;

% % Določimo residual in ga narišemo
% r2 = @(x) f(x)-p2(x);
% x = linspace(0,pi,1001);
% figure
% plot(x,r2(x))

% % Vzamemo absolutno vrednosti residuala
% r2_abs = @(x) abs(f(x)-p2(x));

% % Absolutno vrednost še negiramo, da bomo lahko uporabili fminbnd
% r2_abs2 = @(x) -r2_abs(x);

% % Narišemo
% figure
% hold on
% plot(x,r2_abs(x));
% plot(x,r2_abs2(x));
% hold off

% % Poiščemo y2, pri katerem r2 doseže maksimalno absolutno vrednost
% y2 = fminbnd(r2_abs2,0,pi)

% % Preverimo, za koliko se še razlikujeta absolutni vredsnoti r2(y2) in m2
% abs(r2(y2))-abs(m2)
format long;
f = @(x) log(x);
E = [1,3/2,2];
int = [1,2];

[a,b] = remesLin(f,E,int,2);

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
		if i == 1
			vstKoefPoPrvem = a + b
		end

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
	napakaPoDrugem = abs(m)
end

