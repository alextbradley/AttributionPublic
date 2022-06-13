function ss = get_wavi_results(gamma_T,ensembles,members,timeslices, iceout_folder,dx,dy)
% 02/06/22
%
% Pull the wavi results. Return as a solution strucutre with size
% length(gamma_T) x length(ensembles) x length(members) x length(timeslices)
density_ice = 918.0;
ss = struct;
for ig = 1:length(gamma_T)
    for ie = 1:length(ensembles)
        for im = 1:length(members)
            %generate run number
            run_no = strcat(num2str(gamma_T(ig), '%02.f'), num2str(ensembles(ie)), num2str(members(im), '%02.f'));

            %if run exists, load it and loop over timeslices
            fname = strcat(iceout_folder, 'ATTR_', run_no, '/run/outfile.nc');

            if exist(fname, 'file') %if file exists, load info
                h = ncread(fname, 'h');
                t = ncread(fname, 'TIME');
                gr = ncread(fname ,'grfrac');
                m = ncread(fname, 'm');
                x = ncread(fname, 'x');
                y = ncread(fname, 'y');   
                b = ncread(fname, 'b');   %bathymetry
                haf = h - (1028.0/918.0)*(-b);
                haf0 = haf(:,:,1); %haf at first time
                vaf0 = sum(sum(haf0(haf0 > 0))) * dx *dy; %vaf at first time


                for it = 1:length(timeslices)
                    ss(ig, ie, im, it).run_no = run_no;

                    [~, idx] = min(abs(timeslices(it) - t)); %get index of wavi outfile closest to specified timeslice

                    ss(ig, ie, im, it).h = squeeze(h(:,:,idx)); %ice thickness at this timeslice 
                    ss(ig, ie, im, it).m = squeeze(m(:,:,idx)); %melt at this timeslice 
                    ss(ig, ie, im, it).gr = squeeze(gr(:,:,idx)); %grounded fraction at this timeslice
			
                    %compute the SLR
                    hafT = haf(:,:,idx); %height aboove floatatioon at this timeslice
                    vafT = sum(sum(hafT(hafT > 0)))*dx *dy; %and corresponding vaf
                    dmaf = (vaf0 - vafT)*density_ice;   %change in mass above floatation
                    slr = dmaf/1e9/361.8;               %associated sea level rise
		    ss(ig,ie,im,it).slr = slr;
                end %end loop over timeslice 
	    else %no solution exists for this parameter set
		    ss(ig,ie,im,it).gr = nan;
		    ss(ig,ie,im,it).m = nan;
		    ss(ig,ie,im,it).slr = nan;
	   	    ss(ig,ie,im,it).h = nan; 
            end %end if statement on solution existing
        end %end loop over ensemble members
    end %end loop over ensembles
end %end loop over gamma_T



