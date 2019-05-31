function f_Kirfactor  = Kirf(c_Ko,c_Ki,c_Ko0,c_Ki0, v_m, v_m0, psifac)
% Returns the Kir-factor for the K+ conductance

e_K=psifac.*log(c_Ko./c_Ki); % K+ reversal potential (V)
e_Ko0 = psifac*log(c_Ko0./c_Ki0); % K+ reversal potential during rest (V)

% converting to mV
delta_v = 1000*(v_m - e_K);
v_m_mil = 1000*v_m;
e_Ko0_mil = 1000*e_Ko0;

fakt1=(1+exp(18.5/42.4))./(1+exp((delta_v + 18.5)/42.4));
fakt2=(1+exp(-(118.6+e_Ko0_mil)/44.1))./(1+exp(-(118.6+v_m_mil)/44.1));

f_Kirfactor = sqrt(c_Ko./c_Ko0).*fakt1.*fakt2;


