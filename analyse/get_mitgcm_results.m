function ss = get_mitgcm_results(gamma_T, ensembles,members,timeslices,oceanout_folder)
% 02/06/22 
%
% Return a solution structure of the mitgcm results.
%
ss = struct;
secs_per_year = 365*24*60^2;
density_ice = 918.0;


for ig = 1:length(gamma_T)
    for ie = 1:length(ensembles)
        for im = 1:length(members) %loop over ensemble members
            for it = 1:length(timeslices)
		%filename
	   	run_no = strcat(num2str(gamma_T(ig), '%02.f'), num2str(ensembles(ie)), num2str(members(im), '%02.f'), num2str(timeslices(it), '%03.f'));
		fname = strcat(oceanout_folder, 'ATTRmit_', run_no, '/run/state2D.nc');
                %if simulation exists, store melt rate at timeslices
		if exist(fname, 'file')
                    melt = ncread(fname, 'SHIfwFlx', [1, 1, 1], [Inf, Inf, Inf]);
                    melt = -melt * secs_per_year / density_ice;
                    melt = squeeze(melt(:,:,end));
                    szm  = size(melt);

		    %put the melt on the same grid as wavi. Final 40 grid
                    %points are ocean only (restoring boundary at 340km,
                    %ice front at 300km)
                    m = zeros(300,50);
                    nf = szm - 40; %40 ocean ocean only grid points
                    m(end - nf + 1:end,:) = melt(1:nf, :);                
                    ss(ig,ie, im, it).m = m;
		end %end if simulation exists
	    end %end loop over timeslices
	end %end loop over ensemble members
    end %end loop over ensembles
end %end loop over gamma value
