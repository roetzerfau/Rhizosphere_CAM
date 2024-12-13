
parfor i = 1:4
if(i == 1)
    Main_Bayreuth3(false,10,"Nomove_mucilageCN10_10x5thdayShoot_rangeDiff_por45");
end
if(i == 2)
Main_Bayreuth3(false,100,"Nomove_mucilageCN100_10x5thdayShoot_rangeDiff_por45");
end
if(i == 3)
Main_Bayreuth3(true,10,"move_mucilageCN10_10x5thdayShoot_rangeDiff_por45");
end
if(i == 4)
Main_Bayreuth3(true,100,"move_mucilageCN100_10x5thdayShoot_rangeDiff_por45");
end
end