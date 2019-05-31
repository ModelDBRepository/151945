function P = NaKpumprate(c_Ko, c_Nai);
% The pump rate for the Na/K pump

K_mNai=10.0;% K+ threshold concentration for Na/K-pump (mol/m^3);
K_mKo=1.5; %Na+ threshold concentration for Na/K-pump (mol/m^3);
Pmax = 1.115e-6; %Maximum pumprate, mol/(m^2*s);
P = Pmax*(c_Ko./(c_Ko+K_mKo)).*c_Nai.^1.5./(c_Nai.^1.5+K_mNai.^1.5); %(mol/(m^2*s))

