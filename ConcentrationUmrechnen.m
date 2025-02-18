%2.425500 / C  0.026    1.055000e+03 
C = 0.0539; n = round(45 * 5.9340); %g/cm³ 
% C = 1; n = 1.055000e+03 ;
C = 10^-4; n = 26415;
% C = 0.3168; n = 655;
% C = 333.8 * 1e-6; n = 1;


massdrySoil = 2.65 * 10^-12* 500* 500 * (1- 0.45) %gsoil
%eine zelle: 2 * 2 * 1 (micron)^3; 
% 10^12 / 4  = 2.5*10^11
cmtocell = 2.5*10^11; %1/10^12 * 4
massabs1cell = C / cmtocell * 1000 %(g/cm^3) -> (mg C/1cell) = (mgC/4micron^3)
masspersoil1cell = massabs1cell/massdrySoil %(mg C) -> (mgC/gsoil)

massabsdomain = massabs1cell * n %(mg C)
masspersoildomain = masspersoil1cell * n %(mgC/gsoil)

%masspersoildomain = masspersoildomain * 1000 %(microgC/gsoil)

%% Anders rum %micron g/gsoil
masspersoildomain = 158/1000; %(micron g/gsoil) -> (mgC/gsoil)
massabsdomain = masspersoildomain  * massdrySoil %(mgC/gsoil) -> (mgC)
massabs1cell = massabsdomain/1%(250^2*0.41)
C = massabs1cell *  cmtocell  / 1000 %(mgC) -> (gC/cm^3) von domain auf cm
%C = 0.0539 * 2.7041e+03; n = 1; 


%%Fläche
C = 3.6;%(micro g C/(h cm²)
risemassabs1cell = C / 10^8 * 4 / 10^6 %(micro g C/(h cm²) -> =(g C/(h cell) durch die Wand kommt so viel C rein

C = risemassabs1cell*  cmtocell 