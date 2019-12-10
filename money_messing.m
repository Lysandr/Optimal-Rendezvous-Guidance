%% money messing
clc; clear all; close all;

base_sal = 95000;
stock  = 167000;

sp500 = 1.1;
spx = 1.2;

sp_inv = 68680.97 + 12000;
spx_inv = stock;

years_T = 5;
years = 1:years_T;
years_plus = 1:(years_T +1);
for k = years
    sp_inv(k+1) = sp_inv(k)*(sp500);
    spx_inv(k+1) = spx_inv(k)*(spx);
end 



figure
format longg
plot(years_plus, sp_inv); grid on; hold on;
plot(years_plus, spx_inv);
ax = gca;
ax.XRuler.Exponent = 0;





























