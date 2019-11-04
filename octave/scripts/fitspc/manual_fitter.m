%this thing is supposed to drive a fit. Manually. The interface shall be the same as "fitter"
%The experience: very different.
%
% [pees, pees_errs] = manual_fitter( spc_pees, spc_model, h_data, ...
%                                    extrermes, binZ, minopts )
%
%Also parameters and return values are the same.

function [pees, pee_errs] = manual_fitter( spc_pees, spc_model, h_data, ...
                                           extremes, binZ, minopts )
    if ~exist( 'minopts', 'var' )
        minopts = { 'lr', 1e-2, 'z', 1e-9, 'M', 1e5 };
    end
    
    pees = spc_pees;
    L_model = __make_model( h_data, spc_model, extremes );
    
    warning off;

    go_on = 1;
    while go_on
        user_says = input( 'mf> ', 's' );
        [cmd, opts] = mf_parser( user_says );
        
        switch( cmd )
            case { 'h', 'help', 'HALP!' }
                __display_halp();
            case { 'd', 'draw' }
                %redraw the plot
                __draw( binZ, spc_model( pees ), h_data, extremes,  'redraw' );
            case { 'show' }
                L_model = __make_model( h_data, spc_model,  extremes );
                disp( ['Jval           ',num2str( L_model( pees ))] );
                disp( ['Pees           ',num2str( pees )] );
                Jcov = xb_covariance( L_model, pees  );
                disp( 'Pees STD errors ' );
                disp( num2str( sqrt( diag( Jcov ) ) ) );
            case { 'p', 'parameters' }
                pees = __set_pees( opts );
            case { 'r', 'run' }
                L_model = @( p ) log( max( spc_model( p ), 1 ) );
                [pees, jval, rc] = xb_gradient_descent( L_model, pees, minopts );
            case { 'anr', 'animate-run' }
                max_iter = minopts{6};
                minopts{6} = round( max_iter/120 );
                pees_last = pees;
                for ii=1:round(max_iter/120):max_iter
                    disp( ['Iterations to ',num2str( ii+10 )] );
                    L_model = __make_model( h_data, spc_model,  extremes );
                    [pees_last, jval, rc] = xb_gradient_descent( L_model, pees_last, minopts );
                    disp( ['Current jval :', num2str( jval )] );
                    if ~isempty( opts ) && strcmp( opts{1}, 'sticky' )
                        __draw( binZ, spc_model( pees_last ) , h_data, extremes );
                    else
                        __draw( binZ, spc_model( pees_last ), h_data, extremes, 'redraw' );
                    end
                    sleep( 0.5 );
                end
                isok = input( 'accept?> ', 's' );
                switch( isok )
                    case { 'y', 'Y', 'yes', 'Yes' }
                        pees = pees_last;
                end
            case { 'ama', 'animate-manual' }
                pees_last = __set_pees( opts(1) );
                pees_matrix = linspace( pees, pees_last, 120 );
                for ii=1:120
                    disp( [ 'Parameters: ',num2str( pees_matrix(:,ii)' ) ] );
                    L_model = __make_model( h_data, spc_model,  extremes );
                    jval = L_model( pees_matrix(:,ii) );
                    disp( ['Jval: ',num2str( jval )] );
                    if ~isempty( opts ) && strcmp( opts{end}, 'sticky' );
                        __draw( binZ, spc_model( pees_matrix(:,ii) ), h_data, extremes );
                    else
                        __draw( binZ, spc_model( pees_matrix(:,ii) ) , h_data, extremes, 'redraw' );
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
    pee_errs = sqrt( diag( J_cov ) );

    warning on;
end

%utility to set parameters
function pees = __set_pees( opts )
    if numel( opts ) == 2
        pees(str2num( opts{1} )) = str2num( opts{2} );
    elseif numel( opts ) == 1
        eval( ['pees = ',opts{1},';'] );
    end
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
function model = __make_model( hdata, spc_model, extremes )
    L_hdata = log( max( hdata, 1 ) );
    L_spc_model = @( p ) log( max( spc_model( p ), 1 ) );
    
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
    disp( 'x|extremes [e1 e2] --> view/set extremes.' );
    disp( 'o|minopts [lr zr mi]    --> minimizer options' );
    disp( 'ok|OK              --> return.' );
end
