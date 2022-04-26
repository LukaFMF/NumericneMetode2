format long;

b = 642/41;
f = @(x,y) b*exp(x).*sin(x) + y;
y0 = 1 - b;
intX = [0 5];
fTocna = @(x) exp(x) - b*exp(x).*cos(x);
h = .01;
xi = 0:h:5;

heun_a = [0 1];
heun_b = [0 0; 1 0];
heun_g = [1/2 1/2];
yi = rungeKutta(xi,f,y0,heun_a,heun_b,heun_g);
y500Napaka = abs(fTocna(xi(end)) - yi(end))

najmanjsaRazdalija = min(vecnorm([xi; yi] - repmat([4;0],1,501))) % norma na vsakem vektorju do tocke (4,0)


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