function f = T_BB(height, percent_cloud_cover)
	if height >= 6000
		f = 204.4;
	else
		f = 214.4 - 0.20*percent_cloud_cover;
	end
end