clc
clear all
global Vt % Global constant

% Constants
k = 1.3806503e-23;   %Boltzmann [J/K]
q = 1.60217646e-19;  %Electron charge [C]

% Experimental values and data
V=  [0, 0.0057, 0.0646, 0.1185, 0.1678, 0.2132, 0.2545, 0.2924, 0.3269, 0.3585, 0.3873, 0.4137, 0.4373, 0.4507, 0.459, 0.4784, 0.496, 0.5119, 0.5265, 0.5398, 0.5521, 0.5633, 0.572692511]; 
I_exp= [0.7605, 0.7605, 0.76, 0.759, 0.757, 0.757, 0.7555, 0.754, 0.7505, 0.7465, 0.7385, 0.728, 0.7065, 0.6894, 0.6755, 0.632, 0.573, 0.499, 0.413, 0.3165, 0.212, 0.1035, 0];

T = 33 + 273.15;     %Nominal operating temperature [K]
n = 1;               %Number of cells in the array
Vt = n*k*T/q;        %Thermal Voltage   


% Minimize Least squares
[umin,fval]=fminsearch(@(u)RECT(u,V,I_exp),[1,1e-8,1,10,1]);

% Results: parameters of equivalent circuit

Ipv=umin(1)
I0=umin(2)
Rs=umin(3)
Rsh=umin(4) 
a=umin(5)

% plot results

 I_modelo = zeros(size(V,2),1)';
 for i=1:size(V,2)
    I_modelo(i) = Panel_Current(umin,V(i));
end
 
 figure (1)
 plot (V,I_exp, 'o');
 hold on
 plot (V,I_modelo,'k');
