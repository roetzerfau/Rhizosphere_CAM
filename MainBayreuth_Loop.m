
% parfor i = 1:4
% if(i == 1)
%    Main_Bayreuth3(true,10,"move_mucilageCN10_10x5thdayShoot_rangeDiff_noCBdepend_por45",false);
% end
% if(i == 2)
% Main_Bayreuth3(true,100,"move_mucilageCN100_10x5thdayShoot_rangeDiff_noCBdepend_por45",false);
% end
% if(i == 3)
% Main_Bayreuth3(true,10,"move_mucilageCN10_10x5thdayShoot_rangeDiff_por45",true);
% end
% if(i == 4)
% Main_Bayreuth3(true,100, "move_mucilageCN100_10x5thdayShoot_rangeDiff_por45",true);
% end
% end
%Main_Bayreuth3(true,100, "move_mucilageCN100_10x5thdayShoot_rangeDiff_por45",true);
%Main_Bayreuth3(true,100,"move_mucilageCN100_10x5thdayShoot_rangeDiff_noCBdepend_por45",false);
%Main_Bayreuth3(false,100,"Nomove_test",false);

%Main_Bayreuth3(true,100, "move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN",true, true);
%Main_Bayreuth3(true,100, "move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_nostartMN",true, false);

%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN1",true, 1);
%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0",true, 0);
%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0_5",true, 0.5);


%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0_5_domainSPP",true, 0.5);

%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_PaperReady_leak1percent",true, 0.5);
%Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_10percentPOM",true, 0.5);
%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_PaperReady_2nd",true, 0.5);
%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_30_1",true, 0.5);
%Main_Bayreuth3(false,100, "test",true, 0.5);
%Main_Bayreuth3(false,100, "Nomove_mucilageCN100_PaperReady_zeroRM_RG_2nd",true, 0.5);

%Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_2POM_leak_0_0001",true, 0.5, 2, 0.0001);
%Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_2POM_leak_0_0042",true, 0.5, 2, 0.0042);
%Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_2POM_leak_0_01",true, 0.5, 2, 0.01);
exudatesDays = [200:5:249];
%Main_Bayreuth3(true,50,[200], "paperReady_test_19_02_POM10_zuwenigexudate",true, 0.5, 10, 0.0001);
%Main_Bayreuth3(true,50,[200], "paperReady_test_19_02_POM10",true, 0.5, 10, 0.0001);
%Main_Bayreuth3(false,50,[200], "paperReady_test_19_02_POM10_movefalse",true, 0.5, 10, 0.0001);
%Main_Bayreuth3(true,50,exudatesDays, "paperReady_test_19_02_POM10_exudationDays",true, 0.5, 10, 0.0001);
%Main_Bayreuth3(true,50,[200], "paperReady_test_19_02_POM10_duengen_0",true, 0.5, 10, -10^-5);
%Main_Bayreuth3(false,50,[200], "paperReady_test_18_03_POM10_noMove",true, 0.5, 10, 0.0001);

%Main_Bayreuth3(false,50,[200], "paperReady_referenceSetting_without_particle_Movement",true, 0.5, 10, 0.0001);

%Main_Bayreuth3(true,50,[200], "paperReady_test_18_03_POM10_kPom_0_1",true, 0.5, 10, 0.0001);
%Main_Bayreuth3(true,50,[200], "paperReady_test_18_03_POM10_kPom_0_5",true, 0.5, 10, 0.0001);
%Main_Bayreuth3(true,50,[200], "paperReady_test_18_03_POM10_kPom_0_9",true, 0.5, 10, 0.0001);


%Main_Bayreuth3(false,50,[200], "paperReady_test_18_03_POM10_noMove_oldMB",true, 0.5, 10, 0.0001);
%Main_Bayreuth3(true,50,[200], "paperReady_test_18_03_POM10",true, 0.5, 10, 0.0001);


parfor i = 1:3
if(i == 1)
Main_Bayreuth3(false,40,[200], "paperReady_referenceSetting_CN40_without_particle_Movement",true, 0.5, 10, 0.0001);
end
if(i == 2)
Main_Bayreuth3(true,40,[200], "paperReady_CN40_referenceSetting",true, 0.5, 10, 0.0001);
end
if(i == 3)
Main_Bayreuth3(true,40,exudatesDays, "paperReady_CN40_referenceSetting_exudationDays",true, 0.5, 10, 0.0001);
end
% if(i == 4)
% Main_Bayreuth3(true,10,[200], "paperReady_referenceSetting_CN10",true, 0.5, 10, 0.0001);
% end
end



% parfor i = 1:4
% if(i == 1)
% Main_Bayreuth3(true,50,[200], "paperReady_referenceSetting_leaked_0",true, 0.5, 10, 0);
% end
% if(i == 2)
% Main_Bayreuth3(true,50,[200], "paperReady_referenceSetting_duengen",true, 0.5, 10, -10^-6);
% end
% if(i == 3)
% Main_Bayreuth3(true,100,[200], "paperReady_referenceSetting_CN100",true, 0.5, 10, 0.0001);
% end
% if(i == 4)
% Main_Bayreuth3(true,50,[200], "paperReady_referenceSetting_leaked_0_01",true, 0.5, 10, 0.01);
% end
% end
% parfor i = 1:3
% if(i == 1)
% Main_Bayreuth3(true,50,[200], "paperReady_referenceSetting_POM2",true, 0.5, 2, 0.0001);
% end
% if(i == 2)
% Main_Bayreuth3(true,50,[200], "paperReady_referenceSetting_duengen",true, 0.5, 10, -10^-6/24);
% end
% if(i == 3)
% Main_Bayreuth3(true,50,[200:50:499], "paperReady_referenceSetting_50dayspulstimesteps",true, 0.5, 10, 0.0001);
% end
% end
%   Main_Bayreuth3(true,50,[200], "paperReady_referenceSetting_again_again",true, 0.5, 10, 0.0001);
% parfor i = 1:4
% if(i == 1)
% Main_Bayreuth3(true,40,[200], "paperReady_referenceSetting_CN40",true, 0.5, 10, 0.0001);
% end
% if(i == 2)
% Main_Bayreuth3(true,30,[200], "paperReady_referenceSetting_CN30",true, 0.5, 10, 0.0001);
% end
% if(i == 3)
% Main_Bayreuth3(true,20,[200], "paperReady_referenceSetting_CN20",true, 0.5, 10, 0.0001);
% end
% if(i == 4)
% Main_Bayreuth3(true,22,[200], "paperReady_referenceSetting_CN22",true, 0.5, 10, 0.0001);
% end
% end
% parfor i = 1:3
% if(i == 1)
% Main_Bayreuth3(true,75,[200], "paperReady_referenceSetting_CN75",true, 0.5, 10, 0.0001);
% end
% if(i == 2)
% Main_Bayreuth3(true,50,[], "paperReady_referenceSetting_noInput",true, 0.5, 1000, 0.0001);
% end
% if(i == 3)
%     Main_Bayreuth3(true,50,[200], "paperReady_referenceSetting_again",true, 0.5, 10, 0.0001);
% end
% end
%Main_Bayreuth3(true,80,[200], "paperReady_referenceSetting_CN80",true, 0.5, 10, 0.0001);
% parfor i = 1:3
% if(i == 1)
% Main_Bayreuth3(true,10,exudatesDays, "paperReady_referenceSetting_exudatesDays_CN10",true, 0.5, 10, 0.0001);
% end
% if(i == 2)
% Main_Bayreuth3(true,51,[200], "paperReady_referenceSetting_CN51",true, 0.5, 10, 0.0001);
% end
% if(i == 3)
%     Main_Bayreuth3(true,49,[200], "paperReady_referenceSetting_CN49",true, 0.5, 10, 0.0001);
% end
% end
%Main_Bayreuth3(true,50,[200], "paperReady_test_19_02_POM10_leak_0_01",true, 0.5, 10, 0.01);
%Main_Bayreuth3(true,50,[200], "paperReady_test_19_02_POM10_leak_0",true, 0.5, 10, 0);

%Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_Vmax_11_02",true, 0.5, 2, 0.0001);

%Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM10_leakMineral",true, 0.5, 10, 0.01);

% Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM10_leak0001",true, 0.5, 10, 0.0001);
% Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM10_leak0042",true, 0.5, 10, 0.0042);
% Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM10_leak01",true, 0.5, 10, 0.01);
% Main_Bayreuth3(false,100, "Nomove_mucilageCN100_PaperReady_1puls",true, 0.5, 10, 0.0001);

% parfor i = 1:4
% if(i == 1)
%     Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM10_leak0001",true, 0.5, 10, 0.0001);
% end
% if(i == 2)
% Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM10_leak0042",true, 0.5, 10, 0.0042);
% end
% if(i == 3)
% Main_Bayreuth3(true,100, "move_mucilageCN100_PaperReady_1puls_POM10_leak01",true, 0.5, 10, 0.01);
% end
% if(i == 4)
% Main_Bayreuth3(false,100, "Nomove_mucilageCN100_PaperReady_1puls",true, 0.5, 10, 0.0001);
% end
% end