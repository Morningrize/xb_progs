#! /bin/octave --no-gui
#This tiny script puts the fit in an excitation energy vs. counts
#spectrum (and saves it). Strength conversion will happen elsewhere

args = argv();

disp( ['Loading ', args{1}] );
load( args{1} );

nrgs = zeros( numel( fin ), 1 );
for ii=1:numel( fin )
    if ~isempty( strfind( fin{ii}, 'twoplus' ) )
        nrgs(ii) = 4;
    else
        ee = sscanf( fin{ii}, 'sims/r%d/events_%fMeV_r%d.kb.dp.xb' );
        nrgs(ii) = ee(2);
    end
end

[nrgs, sorting] = sort( nrgs );
spc_pees = spc_pees( sorting )(:);
spc_errs = spc_errs( sorting )(:);

%NOTE: this is where I can include efficiency, specificity and the errors in the discourse.
%      not done yet, but soon to be done. Duh.
counts = 1e5*spc_pees;
c_errs = 1e5*spc_errs;

data = xb_load_data( '/home/gatto/PhD_local/xb_data/xb/tin-132/lead-1134/9xx_AR/tin132-to132/v05/r9xx_cuts_narrow.xb' );
%NOTE: this is the only place where the Tpats will have some sort of importance. I need to give it a brughing over
%      and settle on a bunch of triggers I actually want to use for event scaling (I already don't use them
%      for selection... I just don't want to be screwed up by them).
scalers = xb_tpat2scaler( data, 0x20 );

%I could ignore the scalers here, but I'm still not sure. Anyhow, they are saved
%and I can divide the RP and XS by them.
rp = xb_counts2rp( scalers, counts );%/scalers.a;
rp_errs = xb_counts2rp( scalers, c_errs );%/scalers.a;

xs = xb_rp2xs( rp, 1.132, 'lead' );
xs_errs = xb_rp2xs( rp_errs, 1.132, 'lead' );

xb_save_spc( [args{1},'_exnrg_counts'], { counts }, { nrgs }, { c_errs } );
xb_save_spc( [args{1},'_exnrg_xs'], { xs }, { nrgs }, { xs_errs } );

%this could be made better by actually using the incoming beta distribution
%instead of just the central value. But what can you do.
[vpn_e1, vpn_e2] = xb_virtualphotons( 1e3*nrgs, xb_beam2beta( 520, 132, 50 ), 50, 132, 82, 207.2 );
phxs = xs./vpn_e1(:).*nrgs;
phxs_errs = xs_errs./vpn_e1(:).*nrgs;
strength = xb_xs2strength( 1e3*nrgs, ones( size( nrgs ) ), phxs, 1 );
strength_errs = xb_xs2strength( 1e3*nrgs, ones( size( nrgs ) ), phxs_errs, 1 );

xb_save_spc( [args{1},'_exnrg_phxs'], { phxs }, { nrgs }, { phxs_errs } );
xb_save_spc( [args{1},'_exnrg_strength'], { strength }, { nrgs }, { strength_errs } );

alpha_D = 8*pi/9*sum( strength./nrgs.^2 );
alpha_D_err = 8*pi/9*sum( sum( strength_errs.^2 )./nrgs.^2 );
disp( ['alpha_D = ',num2str( alpha_D ),'+/-',num2str( alpha_D_err )] );

clear data ee;
save( [args{1},'_procd'] );
