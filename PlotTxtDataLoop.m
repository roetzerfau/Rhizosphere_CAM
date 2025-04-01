%PlotTxtdata2("Nomove_mucilageCN100_PaperReady_zeroRM_RG_2nd","Nomove_mucilageCN100_PaperReady_zeroRM_RG_2nd")
%PlotTxtdata2("test", "test")
%PlotTxtdata2("Nomove_mucilageCN100_PaperReady", "Nomove_mucilageCN100_PaperReady")
%PlotTxtdata2("Nomove_mucilageCN100_PaperReady_spreadNoDivide", "Nomove_mucilageCN100_PaperReady_spreadNoDivide")
%PlotTxtdata2("Nomove_mucilageCN100_PaperReady_spreadNoDivide_leak_noMaint", "Nomove_mucilageCN100_PaperReady_spreadNoDivide_leak_noMaint")
%move_mucilageCN100_PaperReady_1puls

%PlotTxtdata2("Nomove_mucilageCN100_PaperReady","Nomove_mucilageCN100_PaperReady")
%PlotTxtdata2("Nomove_mucilageCN100_PaperReady_leak1percent","Nomove_mucilageCN100_PaperReady_leak1percent")


% files = ["paperReady_referenceSetting", "paperReady_referenceSetting_duengen","paperReady_referenceSetting_leaked_0","paperReady_referenceSetting_leaked_0_01"];
% names_appendix = [" CN50"," duengen"," leaked0", " leaked001"];
% PlotTxtdataCompare(files, "compareLeak", names_appendix)
% 
% files = ["paperReady_referenceSetting", "paperReady_referenceSetting_POM2"];
% names_appendix = [" 10"," 2"];
% PlotTxtdataCompare(files, "comparePOM", names_appendix)
% 
% files = ["paperReady_referenceSetting_CN40", "paperReady_referenceSetting_CN10","paperReady_referenceSetting_50dayspulstimesteps","paperReady_referenceSetting_exudationDays"];
% names_appendix = [" CN40"," CN10"," 50dayspuls", " exudationDays"];
% PlotTxtdataCompare(files, "compareExudateTiming2", names_appendix)
% 
files = ["paperReady_referenceSetting_CN40", "paperReady_referenceSetting_CN10","paperReady_referenceSetting_exudationDays"];
names_appendix = [" CN40"," CN10", " exudationDays"];
PlotTxtdataCompare(files, "compareExudateTiming", names_appendix)
% % 
% 
% % 
% files = ["paperReady_referenceSetting", "paperReady_referenceSetting_without_particle_Movement"];
% names_appendix = [" withMovement"," withoutMovement"];
% PlotTxtdataCompare(files, "compareParticleMovement", names_appendix)

files = ["paperReady_referenceSetting_CN10","paperReady_referenceSetting_CN20","paperReady_referenceSetting_CN40", "paperReady_referenceSetting_CN100"];
names_appendix = [ " CN10"," CN20", " CN40" ," CN100"];
PlotTxtdataCompare(files, "compareExudateQuality", names_appendix)

files = ["paperReady_referenceSetting_CN10","paperReady_referenceSetting_CN20",...
    "paperReady_referenceSetting_CN40", "paperReady_referenceSetting_CN100", ...
    "paperReady_referenceSetting_exudationDays"];%, "paperReady_referenceSetting_CN100"
names_appendix = [" CN10", " CN20", " CN40", " CN100", " exudationDays"];%," CN100"
PlotTxtdataCompare(files, "compareAll", names_appendix)

% 
% files = ["paperReady_referenceSetting_CN10","paperReady_referenceSetting", "paperReady_referenceSetting_CN100"];
% names_appendix = [" reference"," CN10", " CN100"];
% PlotTxtdataCompare(files, "compareAll", names_appendix)

PlotTxtdata2("paperReady_referenceSetting", "paperReady_referenceSetting")
PlotTxtdata2("paperReady_referenceSetting_CN40", "paperReady_referenceSetting_CN40")
PlotTxtdata2("paperReady_referenceSetting_without_particle_Movement", "paperReady_referenceSetting_without_particle_Movement")
% PlotTxtdata2("paperReady_referenceSetting_exudationDays", "paperReady_referenceSetting_exudationDays")
% PlotTxtdata2("paperReady_referenceSetting_CN10", "paperReady_referenceSetting_CN10")
% 
% PlotTxtdata2("paperReady_referenceSetting_leaked_0", "paperReady_referenceSetting_leaked_0")
% PlotTxtdata2("paperReady_referenceSetting_duengen", "paperReady_referenceSetting_duengen")
% PlotTxtdata2("paperReady_referenceSetting_CN100", "paperReady_referenceSetting_CN100")
PlotTxtdata2("paperReady_referenceSetting_leaked_0_01", "paperReady_referenceSetting_leaked_0_01")
PlotTxtdata2("paperReady_referenceSetting_POM2", "paperReady_referenceSetting_POM2")



PlotTxtdata2("paperReady_test_18_03_POM10", "paperReady_test_18_03_POM10")
PlotTxtdata2("paperReady_test_18_03_POM10_alterMBfactor", "paperReady_test_18_03_POM10_alterMBfactor")
PlotTxtdata2("paperReady_test_18_03_POM10_kPom_0_1", "paperReady_test_18_03_POM10_kPom_0_1")
PlotTxtdata2("paperReady_test_18_03_POM10_kPom_0_5", "paperReady_test_18_03_POM10_kPom_0_5")
PlotTxtdata2("paperReady_test_18_03_POM10_kPom_0_9", "paperReady_test_18_03_POM10_kPom_0_9")

PlotTxtdata2("paperReady_test_19_02_POM10_CN100", "paperReady_test_19_02_POM10_CN100")
PlotTxtdata2("paperReady_test_19_02_POM10_CN10", "paperReady_test_19_02_POM10_CN10")
PlotTxtdata2("paperReady_test_19_02_POM10_CN10_exudationDays", "paperReady_test_19_02_POM10_CN10_exudationDays")
PlotTxtdata2("paperReady_test_19_02_POM10_leak_0","paperReady_test_19_02_POM10_leak_0")
PlotTxtdata2("paperReady_test_19_02_POM10_duengen_0","paperReady_test_19_02_POM10_duengen_0")


PlotTxtdata2("paperReady_test_19_02_POM10_leak_0_01","paperReady_test_19_02_POM10_leak_0_01")


PlotTxtdata2("paperReady_test_19_02_POM10","paperReady_test_19_02_POM10")
PlotTxtdata2("paperReady_test_19_02_POM10_movefalse","paperReady_test_19_02_POM10_movefalse")
PlotTxtdata2("paperReady_test_19_02_POM10_exudationDays","paperReady_test_19_02_POM10_exudationDays")





PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_Vmax_11_02","move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_Vmax_11_02")
PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_10_02","move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_10_02")
% PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_oldPOM_10_02","move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_oldPOM_10_02")
% PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_oldBimass_10_02","move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_oldBimass_10_02")
% PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_oldK_CLiquid_10_02","move_mucilageCN100_PaperReady_1puls_POM2_testValues_newbase_oldK_CLiquid_10_02")

PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_K_Cliquid_less","move_mucilageCN100_PaperReady_1puls_POM2_testValues_K_Cliquid_less")
PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_K_Cliquid_less_POM_less","move_mucilageCN100_PaperReady_1puls_POM2_testValues_K_Cliquid_less_POM_less")
PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_basis","move_mucilageCN100_PaperReady_1puls_POM2_testValues_basis")
PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues_biomassbefore","move_mucilageCN100_PaperReady_1puls_POM2_testValues_biomassbefore")

PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM2_testValues","move_mucilageCN100_PaperReady_1puls_POM2_testValues")
PlotTxtdata2("move_mucilageCN100_PaperReady_2POM_leak_0_0042","move_mucilageCN100_PaperReady_2POM_leak_0_0042")
PlotTxtdata2("move_mucilageCN100_PaperReady_2POM_leak_0_0001","move_mucilageCN100_PaperReady_2POM_leak_0_0001")
PlotTxtdata2("move_mucilageCN100_PaperReady_2POM_leak_0_01","move_mucilageCN100_PaperReady_2POM_leak_0_01")

PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM10_leak0001", "move_mucilageCN100_PaperReady_1puls_POM10_leak0001")
PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM10_leak0042", "move_mucilageCN100_PaperReady_1puls_POM10_leak0042")
PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM10_leak01", "move_mucilageCN100_PaperReady_1puls_POM10_leak01")
PlotTxtdata2("Nomove_mucilageCN100_PaperReady_1puls", "Nomove_mucilageCN100_PaperReady_1puls")
PlotTxtdata2("move_mucilageCN100_PaperReady_1puls_POM10_leakMineral", "move_mucilageCN100_PaperReady_1puls_POM10_leakMineral")

PlotTxtdata2("move_mucilageCN100_PaperReady_1puls", "move_mucilageCN100_PaperReady_1puls")
PlotTxtdata2("move_mucilageCN100_PaperReady_2POM_leak_0_0001", "move_mucilageCN100_PaperReady_2POM_leak_0_0001")
PlotTxtdata2("move_mucilageCN100_PaperReady_2POM_leak_0_0042", "move_mucilageCN100_PaperReady_2POM_leak_0_0042")
PlotTxtdata2("move_mucilageCN100_PaperReady_2POM_leak_0_01", "move_mucilageCN100_PaperReady_2POM_leak_0_01")


PlotTxtdata2("move_mucilageCN100_PaperReady_2nd", "move_mucilageCN100_PaperReady_2nd")

PlotTxtdata2("Nomove_mucilageCN100_30_1", "Nomove_mucilageCN100_30_1")
PlotTxtdata2("move_mucilageCN100_PaperReady_2nd", "move_mucilageCN100_PaperReady_2nd")
PlotTxtdata2("Nomove_mucilageCN100_PaperReady_leak1percent", "Nomove_mucilageCN100_PaperReady_leak1percent")
PlotTxtdata2("Nomove_mucilageCN100_PaperReady_spreadNoDivide_leak", "Nomove_mucilageCN100_PaperReady_spreadNoDivide_leak")
PlotTxtdata2("move_mucilageCN100_PaperReady_10percentPOM", "move_mucilageCN100_PaperReady_10percentPOM")

% PlotTxtdata2("Nomove_mucilageCN100_PaperReady_spreadNoDivide_leak", "Nomove_mucilageCN100_PaperReady_spreadNoDivide_leak")
%PlotTxtdata2("move_mucilageCN100_PaperReady", "move_mucilageCN100_PaperReady")
%PlotTxtdata2("move_mucilageCN100_PaperReady_zeroRM_RG", "move_mucilageCN100_PaperReady_zeroRM_RG")


PlotTxtdata2("Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN1","Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN1")
PlotTxtdata2("Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0","Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0")
PlotTxtdata2("Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0_5","Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0_5")
PlotTxtdata2("Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0_5_domainSPP","Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN0_5_domainSPP")




PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_nostartMN","move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_nostartMN")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN","move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_factor50_startMN")


PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por45","move_mucilageCN100_10x5thdayShoot_rangeDiff_por45")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_noCBdepend_por45","move_mucilageCN100_10x5thdayShoot_rangeDiff_noCBdepend_por45")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_Kx10","move_mucilageCN100_10x5thdayShoot_rangeDiff_por45_Kx10")



PlotTxtdata2("Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45", "Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45")
PlotTxtdata2("Nomove_mucilageCN10_10x5thdayShoot_rangeDiff_por45", "Nomove_mucilageCN10_10x5thdayShoot_rangeDiff_por45")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por45", "move_mucilageCN100_10x5thdayShoot_rangeDiff_por45")
PlotTxtdata2("move_mucilageCN10_10x5thdayShoot_rangeDiff_por45", "move_mucilageCN10_10x5thdayShoot_rangeDiff_por45")



PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por45", "move_mucilageCN100_10x5thdayShoot_rangeDiff_por45")

%PlotTxtdata2("Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por40_02", "Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por40_02")
%PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por40_02", "move_mucilageCN100_10x5thdayShoot_rangeDiff_por40_02")

%PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff", "move_mucilageCN100_10x5thdayShoot_rangeDiff")


PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff", "move_mucilageCN100_10x5thdayShoot_rangeDiff_por50")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_homoDiff", "move_mucilageCN100_10x5thdayShoot_homoDiff_por50")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por45", "move_mucilageCN100_10x5thdayShoot_rangeDiff_por45")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_homoDiff_por45", "move_mucilageCN100_10x5thdayShoot_homoDiff_por45")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_rangeDiff_por50", "move_mucilageCN100_10x5thdayShoot_rangeDiff_por50")
PlotTxtdata2("move_mucilageCN100_10x5thdayShoot_homoDiff_por50", "move_mucilageCN100_10x5thdayShoot_homoDiff_por50")

PlotTxtdata2("test_C_NRoot10_Negativereturn", "test_C_NRoot10_Negativereturn")
PlotTxtdata2("test_C_NRoot100_Negativereturn", "test_C_NRoot100_Negativereturn")
PlotTxtdata2("test_C_NRoot100_Negativeclamp", "test_C_NRoot100_Negativeclamp")


PlotTxtdata2("noMove_mucilageCN100_10x5thdayShoot_factor_plus_one", " default parameter, CN100_10x5thdayShoot factor_plus_one")
PlotTxtdata2("noMove_mucilageCN100_10x5thdayShoot_factor_10"," default parameter, CN100_10x5thdayShoot factor_10")





PlotTxtdata2("noMove_mucilageCN1_2yearCycle_2nd", " default parameter, no Movement of particles, mucilage input with CN 1 for 20 days, 2 year cycle ")
PlotTxtdata2("noMove_2years"," default parameter, no Movement of particles ")
PlotTxtdata2("noMove_mucilageCN01_2yearCycle", "default parameter, no Movement of particles, mucilage input with CN 0.1 for 20 days, 2 year cycle ")
% PlotTxtdata2("noMove_mucilageCN1_2yearCycle", "2000 days, default parameter, no Movement of particles, mucilage input with CN 1 for 20 days, 2 year cycle ")
PlotTxtdata2("noMove_mucilageC_2yearCycle", "default parameter, no Movement of particles, mucilage input with CN 100 for 20 days, 2 year cycle ")
PlotTxtdata2("noMove_HighN_2yearCycle", "no Movement of particles, N Mucilage input (amount of N as the default carbon amount, CN 0.01) for 20 days, 2 year cycle ")
PlotTxtdata2("Move_normal", "default parameter, Movement of particles ")
PlotTxtdata2("noMove_noPOMDecay", "default parameter, no Movement of particles,no POM decay")






















% PlotTxtdata2("Move_mucilageC_2yearCycle", "2000 days, default parameter, Movement of particles, mucilage input for 20 days, 2 year cycle ")
% PlotTxtdata2("noMove_mucilageC_after100steps", "1000 days default parameter, no Movement of particles, mucilage input after 100 days for 10 days")
% PlotTxtdata2("noMove_mucilageC_after100stepsCN01", "1000 days default parameter, no Movement of particles, mucilage input after 100 days for 10 days with high N amount")

% PlotTxtdata2("noMove_mucilageC_after100stepsCN01")
% PlotTxtdata2("noMove_mucilageC_after100steps")
% PlotTxtdata2("noMove_C_N_S1000")
% PlotTxtdata2("noMove_C_N_S01")
% PlotTxtdata2("DOC0001")
% PlotTxtdata2("noMove_noMNDecay")
% PlotTxtdata2("noMove_noPOMDecay")
% PlotTxtdata2("noMove_noPOMDecay_startMucilage001")
% PlotTxtdata2("noMove_startMucilage001")