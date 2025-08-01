clear ; clc; close all;set_figure_background()

lambda = 1.55*10^(-6);
textColor = [0.9 0.9 0.9]; 

[Beta,N_effective,kx_valid,gamma_valid,W] = Set_WG_param(1.75,1.65,10,lambda,[]);
fprintf("\n\n")
[E_MM,x] = Set_ElectricField(Beta,N_effective,kx_valid,gamma_valid,W);
for i=1:length(E_MM(:,1))
    fprintf("normalized Ey: %f \n",max(abs(E_MM(i,:))));
end
fprintf("\n\n")
[Beta2,N_effective2,kx_valid2,gamma_valid2,~] = Set_WG_param(1.75,1.65,1,[],0.5*10^-6);
fprintf("\n\n")
[E_SM,x2]= Set_ElectricField(Beta2,N_effective2,kx_valid2,gamma_valid2,0.5*10^(-6));

[psi,z] = overlap_integeral(E_MM,E_SM,x,Beta);
H = calc_thick(2,1.5,1,1.75,lambda,1);
t = calc_thick(2,1.5,1,1.65,lambda,1);
Rib_width = calc_thick(1.75,1.65,1.65,N_effective,lambda,(1.75/1.65)^2);
fprintf("\n\n")
fprintf("The thickness in Rib %f (um) \n",H*10^6);
fprintf("The thickness in slab %f (um)\n",t*10^6);
fprintf("The Width of Rib %f (um)\n",Rib_width*10^6);

%% plotting for overlap integeral :
 % at L we have triple image getting it by Trial and error
% also we can get approximate same value by noting from graph that self image
% around 220um which will equal to 3*L
L=2*pi/(Beta(1)- Beta(3)); % Beta(3) as Beta(2) is odd and won't propagate

idx_Z0 = 1;
figure(),plot(x*10^6,abs(psi(idx_Z0,:)));
xlabel("X(um)")
ylabel("|E(x,z)|")
title(sprintf("Input field:Z = 0 = %0.1f(um)",0));

idx_Z1 = find(z>=L/3,1);
figure(),plot(x*10^6,abs(psi(idx_Z1,:)));
xlabel("X(um)")
ylabel("|E(x,z)|")
title(sprintf("triple image:Z = Lpi = %0.1f(um)",L/3*10^6));

idx_Z2 = find(z>=0.5*L,1);
figure(),plot(x*10^6,abs(psi(idx_Z2,:)));
title(sprintf("double image:Z = 1.5Lpi = %0.1f(um)",0.5*L*10^6));
xlabel("X(um)")
ylabel("|E(x,z)|")

idx_Z3 = find(z>=L,1);
figure(),plot(x*10^6,abs(psi(idx_Z3,:)));
title(sprintf("Z = 3Lpi = %0.1f(um)",L*10^6));
xlabel("X(um)")
ylabel("|E(x,z)|")



function [Beta,N_effective,kx_valid,gamma_valid,W] = Set_WG_param(n1,n2,m,lambda,W)
    if nargin < 4 || isempty(lambda)
    lambda = 1.55*10^(-6); %defaul value ;
    end
    Ko = 2*pi/lambda;
    NA = (n1^2-n2^2)^0.5;
    V_cutoff = m*pi/2;
    if nargin < 5 || isempty(W)
         W = floor(((lambda/(2*pi))*V_cutoff/NA)*10^7)/10^7;
    end

fprintf("core Width is(2W): %0.4fu \n",W*2*10^6); 
V = (2*pi/lambda)*W*(n1^2-n2^2)^0.5;
step = 10^(-5);
Kx_W =[0:step:5.1]*pi; % for ten modes Kx will be from 0 to 5*pi  
gamma_V = (-(Kx_W).^2 + (V).^2 ).^0.5;
% we define Kx_W to 5.1*pi then for Kx_W bigger than 5*pi we will have g as imgainary so
% we will cut the img part
[~,index_img] =find(imag(gamma_V));
gamma_V(index_img) = [];
Kx_W(index_img) = [];
%% define dispersion relation:
gamma_W_even = tan(Kx_W).*Kx_W;
gamma_W_odd = -cot(Kx_W).*Kx_W;
gamma_W_even(gamma_W_even <0) = 0;
gamma_W_odd(gamma_W_odd <0) = 0;
gamma_W = gamma_W_odd + gamma_W_even;

%% find intersection between Vnumber curve and disperstion relation:
diff = gamma_V - gamma_W;
intersection_index =zeros(10,1);
ii = 0;
for i = 2:length(diff)-1
if diff(i) > 0 && diff(i+1) < 0
    ii = ii+1;
    intersection_index(ii) = i;
end
end
intersection_index(intersection_index == 0) =[]; 
%% define Beta & Neffective:
Beta = ((n1*2*pi/lambda)^2 - (Kx_W(intersection_index)./W).^2).^0.5;
N_effective = Beta /Ko;
kx_valid = Kx_W(intersection_index)./W;
gamma_valid = gamma_W(intersection_index)./W;
figure(),plot(Kx_W/pi,gamma_V,"LineWidth", 1.5,"color",'r'),hold on;
plot(Kx_W/pi,gamma_W,"lineWidth",1.5,"color",'b'), ylim([0,max(gamma_V)+1]);hold on
plot(Kx_W(intersection_index)/pi, gamma_V(intersection_index), 'ko', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor',[0, 0.75, 1], ... 
    'LineWidth', 1); hold off; grid on;
xlabel("Kx*W")
ylabel("gamma*W")
title("dispersion relation")
xt = get(gca, 'XTick');              
xticklabels(arrayfun(@(x) [num2str(x) '\pi'], xt, 'UniformOutput', false));
for i = 1:length(N_effective)
fprintf("effective index of mode %d = %0.4f\n",i-1 ,N_effective(i))
end
end



function [E_y,x] = Set_ElectricField(Beta,N_effective,kx_valid,gamma_valid,W)
%% set x and E steps :
x_step = 10^(-8);
W = ceil(W*10^8)/10^8;
X_boundary = 13*10^-6;
x = -X_boundary:x_step:X_boundary ;
core_idx1 = find(x>=-W,1);
core_idx2 = find(x>=W,1);
x_region_1 =x(1:core_idx1-1); 
x_core = x(core_idx1:core_idx2);
x_region_2 = x(core_idx2+1:end);
E_region_1 =zeros(length(N_effective),length(x_region_1));
E_core = zeros(length(N_effective),length(x_core));
E_region_2 =zeros(length(N_effective),length(x_region_2));
E_y =zeros(length(N_effective),length(x));
%% electric field :
even_num = 0;
odd_num = 0;
evenFigure = figure('Tag', 'even');
oddFigure = figure('Tag', 'odd');
len = length(N_effective);
Row = ceil(len/4);
col = abs(floor(len/2) - Row) ;
for i=1:len
if mod(i-1,2) ==0 %%even mode 
 % electric field region 1:
    even_num =even_num+1; 
    E_region_1(i,:) = exp(1*gamma_valid(i).*(x_region_1 + W)).*cos(kx_valid(i).*W);
    E_core(i,:) = cos(kx_valid(i).*x_core);
    E_region_2(i,:) = exp(-1*gamma_valid(i).*(x_region_2 - W)).*cos(kx_valid(i).*W);
    E_y(i,:) = [E_region_1(i,:), E_core(i,:),E_region_2(i,:)];
    Am =trapz(x, abs(E_y(i,:)).^2);
    Am = 1/Am.^(0.5);
    E_y(i,:) = Am.*E_y(i,:);
    figure(evenFigure),subplot(Row,col,even_num),plot(x*10^6,E_y(i,:),"linewidth",1.5),title(sprintf("TE%d",i-1)),grid on,hold on;
    xlabel("X(um)"),ylabel("E(V/m)")
    max_amp = max(abs(E_y(i, :)))*1.1; % Scale impulse to max E_y_core
    stem([W, W,-W,-W]*10^6, [max_amp, -max_amp,max_amp, -max_amp], 'Color', [0.0 0.9 0.9], 'LineWidth', 1.2, 'Marker', 'none');
    ylim([-max_amp,max_amp]);
    hold off;
  %  legend_label_1(i)=sprintf("Mode number : %d",i-1);
elseif(mod(i-1,2) == 1 )

    odd_num =odd_num+1; 
    E_region_1(i,:) = -1*exp(1*gamma_valid(i).*(x_region_1+W)).*sin(kx_valid(i).*W);
    E_core(i,:) = sin(kx_valid(i).*x_core);
    E_region_2(i,:) = exp(-1*gamma_valid(i).*(x_region_2-W)).*sin(kx_valid(i).*W);
    E_y(i,:) =[E_region_1(i,:), E_core(i,:),E_region_2(i,:)];
    Am =trapz(x, abs(E_y(i,:)).^2);
    Am = 1/Am.^(0.5);
    E_y(i,:) = Am.*E_y(i,:);
    figure(oddFigure),subplot(3,2,odd_num),plot(x*10^6,E_y(i,:),"linewidth",1.5),hold on,title(sprintf("TE%d",i-1));grid on;
    xlabel("X(um)"),ylabel("E(V/m)")
    max_amp = max(abs(E_y(i, :)))*1.1; % Scale impulse to max E_y_core
    stem([W, W,-W,-W]*10^6, [max_amp, -max_amp,max_amp, -max_amp], 'Color', [0 0.9 0.9], 'LineWidth', 1.2, 'Marker', 'none');
    ylim([-max_amp,max_amp]);
    hold off;
end
end
end

function [psi,z]= overlap_integeral(E_MM,E_SM,x,Beta) 
E_MM_conj = conj(E_MM);
N = length(E_MM(:,1));
Cv =zeros(N,1);
z_step = 5*10^(-8);
z = 0:z_step:5*10^-4;
psi = zeros(length(z),length(x));
for i = 1:N
 Cv(i,1) = trapz(x,E_SM.*E_MM_conj(i,:));   
 fprintf("OverLap Integeral of TE%d : %f \n",i-1,Cv(i));
end
for i = 1:length(z)
    for ii=1:N
    psi(i,:)= Cv(ii,1)*E_MM(ii,:)*exp(-1j*Beta(ii)*z(i))+ psi(i,:);
    end
end
figure(),
imagesc(z.*10^6,x.*10^6,abs(psi).'),colorbar;
title("E(x,y)");
xlabel("z(um)");ylabel("x(um)")
end


function [H] = calc_thick(n1,n2,n3,N_effective,lambda,eta)
% eta will equal to one in case TE 
% in case TM eta equal n1^2 /n2^2
Neff = N_effective(1);
b = (Neff^2 -n2^2)/(n1^2 - n2^2);
a = (n2^2 - n3^2)/(n1^2 - n2^2);
V = (atan(eta*sqrt(b/(1-b))) + atan(eta*sqrt((b+a)/(1-b))))/sqrt(1-b);
H = (lambda/(2*pi))*V/sqrt(n1^2- n2^2);
end
function set_figure_background()
set(0, 'DefaultFigureColor', [0.1, 0.1, 0.2]);
set(0, 'DefaultFigureColor', [0.1 0.1 0.2]); 
set(0, 'DefaultAxesColor', [0.2 0.2 0.3]); 
textColor = [0.9 0.9 0.9]; 
set(0, 'DefaultAxesXColor', textColor);     
set(0, 'DefaultAxesYColor', textColor);          
set(0, 'DefaultTextColor', textColor);       
set(0, 'DefaultLineColor', textColor);
set(0, 'DefaultAxesGridColor', textColor*0.7);
set(0, 'DefaultAxesMinorGridColor', textColor*0.35);
end
