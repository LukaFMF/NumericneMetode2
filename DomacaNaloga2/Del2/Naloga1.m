format long;
f = @(x,y) y + 15*exp(x).*cos(15*x);
y0 = 0;
intX = [0 1];
fTocna = @(x) exp(x).*sin(15*x);

RK4_a = [0 1/2 1/2 1];
RK4_b = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
RK4_g = [1/6 2/6 2/6 1/6];
fprintf("RK4:\n");
for r = 0:4
	fprintf("r = %d:",r);
	h = .1*2^(-r);
	x = intX(1):h:intX(2);
	y = rungeKutta(x,f,y0,RK4_a,RK4_b,RK4_g);

	maxGlobalnaNapaka = max(abs(fTocna(x) - y))
end

heun_a = [0 1];
heun_b = [0 0; 1 0];
heun_g = [1/2 1/2];
fprintf("Huen:\n");
for r = 0:4
	fprintf("r = %d:",r);
	h = .1*2^(-r);
	x = intX(1):h:intX(2);
	y = rungeKutta(x,f,y0,heun_a,heun_b,heun_g);

	maxGlobalnaNapaka = max(abs(fTocna(x) - y))
end

function y = rungeKutta(x,f,y0,alpha,beta,gamma)
	% Opis:
	%  Funkcija rungeKutta2 vrne numericno resitev navadne diferencialne enacbe
	%  y' = f(x,y) pri pogoju y(x(1)) = y0, ki je izracunana z metodo podana z 
	% eksplicitno Butcherjevo shemo (alpha,beta,gamma)
	%
	% Definicija:
	%  y = rungeKutta2(x,f,y0,alpha,beta,gamma)
	%
	% Vhodni podatki:
	%  x     vrstica delilnih tock,
	%  f     funkcija f v obliki @(x,y) f(x,y),
	%  y0    zacetna vrednost resitve.
	%  alpha vektor koeficientov alpha_i 
	%  beta  strogo spodnje trikotna matrika koeficientov beta_ij 
	%  gamma vektor koeficientov gamma_i 
	%
	% Izhodni podatek:
	%  y    vrstica numericnih priblizkov za vrednosti tocne resitve v delilnih
	%       tockah.

	dimMetode = size(alpha,2);
	stTock = size(x,2);

	y = zeros(1,stTock);
	y(1) = y0;
	for r = 2:stTock
		h = x(r) - x(r-1);
		k = zeros(1,dimMetode);
		for i = 1:dimMetode
			sumK = dot(beta(i,:),k);
			k(i) = f(x(r-1) + h*alpha(i),y(r-1) + h*sumK);
		end

		gammaSumK = dot(gamma,k);
		y(r) = y(r-1) + h*gammaSumK;
	end
end