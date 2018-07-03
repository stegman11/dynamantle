% Mass balance type approach to enrichment of Iron in melt... 
%
%

% McDonough & Sun 95 'Pyrolitic' BSE composition table 5
%     Fe = 6.26 wt%
%     Mg = 22.8 wt%
%     Si = 10.65 wt%
%     Al = 2.35 wt%
%     Ca = 2.53 wt%
%      O = 100% - rest  ?? , = 55.41 wt%
amtO=100-6.26-22.8-10.65-2.53-2.35
%molar masses, appx g/mol:
mm_Mg = 24.305;
mm_Si = 28.085;
mm_O = 15.999;
mm_Fe = 55.845;
mm_Ca = 40.08;

%convert wt-% to mol-%
M_mantle = 4e24 % kg
M_mantle = 4e27 % g

m_Ca = .0253 * M_mantle
mol_Ca = m_Ca / mm_Ca

m_Fe = .0626 * M_mantle
mol_Fe = m_Fe / mm_Fe

m_Mg = .2280 * M_mantle
mol_Mg = m_Mg / mm_Mg

m_Si = .1065 * M_mantle
mol_Si = m_Si / mm_Si

m_O = .5541 * M_mantle
mol_O = m_O / mm_O

mol_tot = mol_Fe + mol_Mg + mol_Si + mol_O + mol_Ca
molpct_Fe = mol_Fe/mol_tot
molpct_O = mol_O/mol_tot
molpct_Si = mol_Si/mol_tot
molpct_Mg = mol_Mg/mol_tot
molpct_Ca = mol_Ca/mol_tot

F = [0:.01:1];
k_D = 0.3; %LZ 4/27 - need k_D to be element specific, 0.3 for Iron
k_Dca = 0.2; %LZ 4/27 for Ca in Mg Peroviskite
cL_Fe = molpct_Fe ./ ( F.*(1-k_D) + k_D );

cL_Ca = molpct_Ca ./ ( F.*(1-k_Dca) + k_Dca );
plot(F,cL_Fe)
%LZ 4/27 is 1.5 is entropy fusion of the "solvent" (i.e. the mantle).
%cL_UTH = 100e-9 ./ ( F.*(1-k_D) + k_D );
%plot(F,cL_UTH)
%  T_m  - T_m^o = k_B * T_m / delS_A^o (Cx^s - Cx^l)
T_mfe = 5000 ./ (1 - (1/1.5).*(cL_Fe*k_D - cL_Fe));
%figure
plot(F,T_mfe)
T_mca = 5000 ./ (1 - (1/1.5).*(cL_Ca*k_Dca - cL_Ca));
figure
plot(F,T_mca)
%T_m = 5000 ./ (1 - (1/1.5).*(cL_UTH*0 - cL_UTH));
%figure
%plot(F,T_m)

%try adding the effects of Fe and Ca together
delTmfe = 5000-T_mfe;
delTmca = 5000-T_mca;
delTm = delTmfe + delTmca;

figure
plot(F,delTmfe,'r'),hold on
plot(F,delTmca,'b'),hold on
plot(F,delTm,'k'),hold off
%mass_pv = Mg + Si + 3*O % = 100.3870 g/mol

% Mantle composition approximation = Perovskite =  (Mg,Fe)SiO3

%Murakami et al Nature 2012, "up to 93% of lower mantle perovskite"
% Though really the Mg/Si ratio is unknown, and I'm not sure that is
% representative of WHOLE mantle, but rather, 'lower' mantle, where the
% abstract doesn't put a number on that.   Upper Mantle Olivine (Mg2SiO4?) 

%Say Partition Coeff (Kd) is 0.4.  Kd = c_s/c_l (Jackson et al 2013)

%------
% calculating c_o, original concentration of iron in otherwise mgsi03
% mantle, this is a gross approximation, not a detailed composition
%------







% Labrosse et al 2007 and ZigSteg2013 used delS = 300 J / (kg K)
%Melting Point depression from Alfe et al
%  T_m  - T_m^o = k_B * T_m / delS_A^o (Cx^s - Cx^l)
%   T_m = melting Temperature of liquid
%   T_m^o = Melting temperature of PURE SOLVENT (e.g. Perovs.)
%   k_B = 1.3806488 × 10-23 (m2 kg)/(s2 K);  Boltzman Constant
%      note (m2 kg) / (s2 K) = J/K
%   delS_A^o == S_A^oL - S_A^oS == Entropy of fusion of pure solvent =
%   J/(mol K), Joule = (kg m2)/s2  * in Alfe et al it is per atom.
%  Cx^s = concentration in solid (%)
%  Cx^l = concentration in liquid (%)
%  concentrations in units of ?
%   k_B/delS_A^o = units of mol?
%  so K - K = mol_p * K * (mol/mol - mol/mol)?

% E_fus = 1.5Nk Stixrude 2009 - what is Nk, N = num atoms, k =
% boltzman?
% 
%E_Fus = 300  % J/(kg K) Labrosse et al 2007
%E_fus_atom = 300


%mass_pv = Mg + Si + 3*O % = 100.3870 g/mol
% What is N for 1 kg?
%Av = 6.022e23;
%kb = 1.38006488e-23;
%convertS = 300 * (1/1000) * (100.3870) * (1/Av)
%test = convertS / (5*kb)
%ninkilo = Av * 1000 / (100/5)
%onenk = 1.0*kb*ninkilo
%onefivenk = 1.5*kb*ninkilo

%T_m = 4000 % K
%delS_kb = 1.5 % N = 1, kb cancels with that above




