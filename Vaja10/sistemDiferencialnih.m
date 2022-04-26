format long;

x = 0:.1:.9;
y2Prime = @(x,y,yPrime) yPrime*y.^2 - y;
y0 = [1, 0]; % y(0) = 1, y'(0) = 0

RK3_a = [0 1/2 1];
RK3_b = [0 0 0; 1/2 0 0; -1 2 0];
RK3_g = [1/6 2/3 1/6];
fs = cell(1,2);
fs{1} = @(x,Y) Y(2);
fs{2} = @(x,Y) y2Prime(x,Y(1),Y(2));
yRK3 = rungeKuttaSistem(x,fs,y0,RK3_a,RK3_b,RK3_g)'


F = @(x,Y) [Y(2);y2Prime(x,Y(1),Y(2))];
sol = ode45(F,[0,1],[1,0]); % 2. argument poda interval [0,1]

Y = deval(sol,x)' % prvi stolpec - vrednsti y, drugi stolpec - vrednosti y' 

function y = rungeKuttaSistem(x,fs,y0,alpha,beta,gamma)
	% Opis:
	%  Funkcija rungeKutta vrne numericno resitev navadne diferencialne enacbe
	%  y' = f(x,y) pri pogoju y(x(1)) = y0, ki je izracunana z metodo podana z 
	% eksplicitno Butcherjevo shemo (alpha,beta,gamma)
	%
	% Definicija:
	%  y = rungeKutta(x,f,y0,alpha,beta,gamma)
	%
	% Vhodni podatki:
	%  x     vrstica delilnih tock,
	%  fs    celicno polje funcij, ki definirajo sistem oblike y_i' = fs{i}(x,y_1,y_2,...,y_n)
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
	stEnacb = size(fs,2);

	y = [y0'];
	for r = 2:stTock
		h = x(r) - x(r-1);
		k = zeros(stEnacb,dimMetode);
		for i = 1:dimMetode
			Y = zeros(1,stEnacb);
			for j = 1:stEnacb
				sumK = dot(beta(i,:),k(j,:));
				Y(j) = y(j,r-1) + h*sumK;
			end

			for j = 1:stEnacb
				% st parametrov funckij fs je odvisno od stevila enacb
				% 2 enacbe = 3 argumenti (x,[y,y'])
				% 3 enacbe = 4 argumenti (x,[y,y',y''])
				% naracunamo jih v prejsnji zanki
				k(j,i) = fs{j}(x(r-1) + h*alpha(i),Y);
			end
		end
		yPrib = zeros(stEnacb,1);
		for i = 1:stEnacb
			gammaSumK = dot(gamma,k(i,:));
			yPrib(i) = y(i,r-1) + h*gammaSumK;
		end
		y = [y yPrib];
	end
end


