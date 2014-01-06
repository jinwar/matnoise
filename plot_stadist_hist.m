clear;

load stainfo_BHZ.mat

stlas = [stainfo.lat];
stlos = [stainfo.lon];

sta_pair_num = 0;

for ista = 1:length(stlas)
	for jsta = ista+1:length(stlas)
		sta_pair_num = sta_pair_num+1;
		sta_dist(sta_pair_num) = deg2km(distance(stlas(ista),stlos(ista),stlas(jsta),stlos(jsta)));
	end
end

figure(68)
clf
hist(sta_dist,30);
