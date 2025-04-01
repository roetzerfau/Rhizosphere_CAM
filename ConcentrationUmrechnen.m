%2.425500 / C  0.026    1.055000e+03 
C = 0.0539; n = round(45 * 5.9340); %g/cm³ 
% C = 1; n = 1.055000e+03 ;
C = 10^-4; n = 26415;
% C = 0.3168; n = 655;
% C = 333.8 * 1e-6; n = 1;
porosity = 0.411;
soil_particleNNZ = 250 * 250 * (1-porosity);
soilParticleDensity = 2.36;

soilfactor = ( soil_particleNNZ * soilParticleDensity/1000 )
%massdrySoil = 2.65 * 10^-12* 500* 500 * (1- 0.45) %gsoil
massdrySoil = soilParticleDensity * 10^-12* 500* 500 * (1- porosity) %gsoil
%eine zelle: 2 * 2 * 1 = 4 (micron)^3; 
% 10^12 / 4  = 2.5*10^11
celltocm = 2.5*10^11; %1/10^12 * 4 (
massabs1cell = C / celltocm * 1000 %(g/cm^3) -> (mg C/1cell) = (mgC/4micron^3)
masspersoil1cell = massabs1cell/massdrySoil %(mg C) -> (mgC/gsoil)

massabsdomain = massabs1cell * n %(mg C)
masspersoildomain = masspersoil1cell * n %(mgC/gsoil)

%masspersoildomain = masspersoildomain * 1000 %(microgC/gsoil)

%% Anders rum %micron g/gsoil
masspersoildomain = 158/1000 %(micron g/gsoil) -> (mgC/gsoil)
massabsdomain = masspersoildomain  * massdrySoil %(mgC/gsoil) -> (mgC)
massabs1cell = massabsdomain/1;%(250^2*0.41)
C = massabs1cell *  celltocm  / 1000 %(mgC/cell) -> (gC/cm^3) von domain auf cm
%C = 0.0539 * 2.7041e+03; n = 1; 
C = C/(250^2 * 0.411)


%%Fläche
C = 3.6;%(micro g C/(h cm²)
risemassabs1cell = C / 10^8 * 4 / 10^6 %(micro g C/(h cm²) -> =(g C/(h cellwand) durch die Wand (4 microm^2) kommt so viel C rein
%cm2 to micrometer2 -> 1e+8
%jeder micron^3 würfel in cm^3 bekommt eigenen exudate
C = risemassabs1cell *  celltocm %(g C/(h cell) -> (g C/(h cm³)
C_total = C * 24 /soilfactor