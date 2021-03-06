%this thing is supposed to drive a fit. Manually. The interface shall be the same as "fitter"
%The experience: very different.
%
% [pees, pees_errs, chisq] = manual_fitter( spc_pees, spc_model, h_data, ...
%                                           extrermes, binZ, minopts )
%
%Also parameters and return values are the same.

function [pees, pee_errs, chisq] = manual_fitter( spc_pees, spc_model, h_data, ...
                                                  extremes, binZ, minopts )
    if ~exist( 'minopts', 'var' )
        minopts = { 'lr', 1e-2, 'z', 1e-9, 'M', 1e5 };
    end
    preweight = 1e-4;

    pees = spc_pees;
    %realistically, this is only going to happen between 0 and 1
    pcage = [ones( 1, numel( pees ) );zeros( 1, numel( pees ) )];
    
    warning off;

    go_on = 1;
    while go_on
        user_says = input( 'mf> ', 's' );
        [cmd, opts] = mf_parser( user_says );
        L_model = __make_model( h_data, spc_model, extremes, preweight );
        switch( cmd )
            case { 'h', 'help', 'HALP!' }
                __display_halp();
            case { 'd', 'draw' }
                %redraw the plot
                __draw( binZ, spc_model( preweight*pees ), h_data, extremes,  'redraw' );
            case { 'show' }
                disp( ['Jval           ',num2str( L_model( pees ))] );
                disp( ['Chi2           ',num2str( __chisq( h_data, spc_model, preweight*pees, ...
                                                           extremes ))] );
                disp( ['Pees           ',num2str( preweight*pees )] );
                disp( ['Preweight      ',num2str( preweight )] );
                Jcov = xb_covariance( L_model, pees  );
                disp( 'Pees standard errors: ' );
                disp( num2str( abs( sqrt( diag( Jcov )' )/sqrt( extremes(2) - extremes(1) ) ) ) );
            case { 'p', 'parameters' }
                oldpees = pees;
                if isempty( opts ); disp( pees );
                else pees = __set_pees( oldpees, opts ); end
                if numel( pees ) ~= numel( oldpees )
                    disp( ['The correct number of pees is ',num2str( numel( oldpees ) )] );
                    pees = oldpees;
                end
            case { 'pw', 'preweight' }
                if isempty( opts ); disp preweight;
                else preweight = str2num( opts{1} ); end
            case { 'c', 'constraint' }
                if isempty( opts ); disp( pcage );
                else try
                    pcage = [opts{2}*ones( 1, numel( pees ) );opts{1}*ones( 1, numel( pees ) )];
                catch
                    warning( 'Inconsistent constraint format!' );
                end end
            case { 'r', 'run' }
                [pees, jval, rc] = xb_constrained_gradient_descent( L_model, pees, pcage, minopts );
            case { 'anr', 'animate-run' }
                max_iter = minopts{6};
                minopts{6} = ceil( max_iter/120 );
                pees_last = pees;
                for ii=1:ceil(max_iter/120):max_iter
                    disp( ['Iterations to ',num2str( ii+ceil(max_iter/120) )] );
                    [pees_last, jval, rc] = xb_constrained_gradient_descent( L_model, pees_last, pcage, minopts );
                    disp( ['Current jval :', num2str( jval )] );
                    disp( ['Current rc   :', num2str( rc )] );
                    if ~isempty( opts ) && strcmp( opts{1}, 'sticky' )
                        __draw( binZ, spc_model( preweight*pees_last ) , h_data, extremes );
                    else
                        __draw( binZ, spc_model( preweight*pees_last ), h_data, extremes, 'redraw' );
                    end
                    sleep( 0.5 );
                    if rc == 0; break; end
                end
                isok = input( 'accept?> ', 's' );
                switch( isok )
                    case { 'y', 'Y', 'yes', 'Yes' }
                        pees = pees_last;
                end
                minopts{6} = max_iter;
            case { 'ama', 'animate-manual' }
                pees_last = __set_pees( opts(1) );
                pees_matrix = linspace( pees, pees_last, 120 );
                for ii=1:120
                    disp( [ 'Parameters: ',num2str( pees_matrix(:,ii)' ) ] );
                    jval = L_model( pees_matrix(:,ii) );
                    disp( ['Jval: ',num2str( jval )] );
                    if ~isempty( opts ) && strcmp( opts{end}, 'sticky' );
                        __draw( binZ, spc_model( preweight*pees_matrix(:,ii) ), h_data, extremes );
                    else
                        __draw( binZ, spc_model( preweight*pees_matrix(:,ii) ) , h_data, extremes, 'redraw' );
                    end
                    sleep( 0.5 )
                end
                isok = input( 'accept?> ', 's' );
                switch( isok )
                    case { 'y', 'Y', 'yes', 'Yes' }
                        pees = pees_last;
                end
            case { 'o', 'minopts' }
                switch( opts{1} )
                    case 'show'
                        disp( minopts );
                    case 'lr'
                        minopts(2) = str2num( opts{2} );
                    case 'zr'
                        minopts(4) = str2num( opts{2} );
                    case 'mi'
                        minopts(6) = str2num( opts{2} );
                    otherwise
                        disp( [opts{1},' is not a valid minimizer option'] );
                    end
            case { 'x', 'extremes' }
                if isempty( opts )
                    disp( extremes )
                else
                    extremes(1) = str2num( opts{1} );
                    extremes(2) = str2num( opts{2} );
                end
            case { 'ok', 'OK' }
                go_on = 0;
	    case { 'die!' }
		exit;
            otherwise
                warning( 'Unknown command. Read the HALP!' );
        end
    end
    if exist( 'pees_las', 'var' );
        pees = pees_last;
    end
    J_cov = xb_covariance( L_model, pees );
    pee_errs = preweight*abs( sqrt( diag( J_cov )' )/sqrt( extremes(2) - extremes(1) ) );
    pees *= preweight;
    chisq = __chisq( h_data, spc_model, pees, extremes );

    warning on;
end

%utility to set parameters
function pees = __set_pees( pees, opts )
    if numel( opts ) == 2
        pees(str2num( opts{1} )) = str2num( opts{2} );
    elseif numel( opts ) == 1
        eval( ['pees = ',opts{1},';'] );
    end
end

%do the chisquared
function chisq = __chisq( hdata, hmodel, pees, extremes )
    hmodel = hmodel( pees );
    chisq = xb_goodness_of_fit( hdata(extremes(1):extremes(2)), ...
	                        hmodel(extremes(1):extremes(2)), ...
			        numel( pees ) );
end

%utility to draw
function __draw( binZ, hmodel, hdata, extremes, mm )
    if exist( 'mm' ) && strcmp( mm, 'redraw' ); hold off; end
    
    stairs( binZ, hdata, 'linewidth', 3, 'k--' );
    hold on;
    stairs( binZ, round( hmodel ), 'linewidth', 2 );
    plot( [binZ(extremes(1)), binZ(extremes(1))], [1, max( hdata )], 'linewidth', 3,'r--' );
    plot( [binZ(extremes(2)), binZ(extremes(2))], [1, max( hdata )], 'linewidth', 3,'r--' );
    ylabel( ['#/',num2str(binZ(2)-binZ(1)),'KeV'] );
    xlabel( 'KeV' );
    set( gca, 'yscale', 'log', 'fontsize', 24, 'linewidth', 2 ); grid on;
end

%utility to make the model
function model = __make_model( hdata, spc_model, extremes, preweight )
    L_hdata = log( max( hdata, 1 ) );
    L_spc_model = @( p ) log( max( spc_model( preweight*p ), 1 ) );
    
    if ~exist( 'extremes' )
        model = @( p ) sum( norm( L_hdata - L_spc_model( p )) )/numel( h_data );
    else
        model = @( p ) sum( norm( L_hdata(extremes(1):extremes(2)) - ...
                                  L_spc_model( p )(extremes(1):extremes(2)) ) )/ ...
                                  numel( extremes(1):extremes(2) );
    end
end

%diplay the help
function __display_halp()
    disp( 'h|help|HALP!       --> displays this help.' );
    disp( 'd|draw             --> refresh the canvas if any and draw the current status.' );
    disp( 'show               --> display some statistics.' );
    disp( 'r|run              --> run the minimizer.' );
    disp( 'anr|animate-run    --> run the minimizer and show an animation of it.' );
    disp( 'ama|animate-manual --> set the parameters manually and animate it.' );
    disp( 'p|parameters [index p|p-vec]    --> set the parameters.' );
    disp( 'pw|preweight [number]           --> get/set the preweight.' );
    disp( 'c|constraint [lower upper]      --> get/set the constraint.' );
    disp( 'x|extremes [e1 e2] --> get/set extremes.' );
    disp( 'o|minopts [lr zr mi]    --> minimizer options' );
    disp( 'ok|OK              --> return.' );
end
