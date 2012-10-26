#!/bin/csh
rm fetchdata_mat.log
unsetenv DISPLAY
foreach ista (`cat stacount`)
	foreach itime (`cat timecount`);
		echo "fetch_data("$ista"," $itime")" > tempmatscript.m
		cat tempmatscript.m
		nohup matlab < tempmatscript.m >> fetchdata_mat.log
	end
end

