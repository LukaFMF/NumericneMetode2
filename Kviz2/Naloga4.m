format long;

RK4_a = [0 1/2 1/2 1];
RK4_b = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
RK4_g = [1/6 2/6 2/6 1/6];

alpha = 182/17;
beta = 3;
gamma = 8;
delta = 2;
y0 = [2,12];

fs = cell(1,2);
fs{1} = @(x,Y) alpha*Y(1) - beta*Y(1)*Y(2);
fs{2} = @(x,Y) delta*Y(1)*Y(2) - gamma*Y(2);

h = .1;
t = h*(0:100);

% prvi stolpec - pleni, drugi stolpec - plenilci
yRK4 = rungeKuttaSistem(t,fs,y0,RK4_a,RK4_b,RK4_g)'

stZajcevPriT50 = yRK4(51,1)

najmanjsaRazlikaVPopulaciji = min(abs(yRK4(:,1) - yRK4(:,2)))

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