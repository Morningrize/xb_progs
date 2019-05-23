%this thing is supposed to drive a fit. Manually. The interface shall be the same as "fitter"
%The experience: very different.
%
% [pees, pees_errs] = manual_fitter( spc_pees, spc_model, h_data, ...
%                                    extrermes, binZ, offset, minopts )
%
%Also parameters and return values are the same.

function function [pees, pee_errs] = fitter( spc_pees, spc_model, h_data, ...
                                             extremes, binZ, offset, minopts )
    if ~exist( 'offset', 'var' )
        offset = zeros( size( binZ ) );
    end
    if ~exist( 'minopts', 'var' )
        minopts = { 'lr', 1e-2, 'z', 1e-9, 'M', 1e5 };
    end
    
    pees = spc_pees;
    
    go_on = 1;
    while go_on
        user_says = input( 'mf> ', 'S' );
        [cmd, opts] = mf_parser( user_says );
        
        switch( cmd )
            case { 'h', 'help', 'HALP!' }
                __display_halp();
            case { 'd', 'draw' }
                %redraw the plot
                
                __draw( binZ, spc_model, h_data, extremes, offset, 'redraw' );
            case { 'show' }
                %display some relevant things.
            case { 'p', 'parameters' }
                pees = __set_pees( pees, opts );
            case { 'r', 'run' }
                L_spc_model = @( p ) log( max( spc_model( p ) + offset, 1 ) );
                [pees, jval, rc] = xb_gradient_descent( model, pees, minopts );
            case { 'anr', 'animate-run' }
                max_iter = minopts{6};
                minopts{6} = 10;
                pees_last = pees;
                for ii=1:10:max_iter
                    disp( ['Iterations to ',num2str( ii+10 )] );
                    model = __make_model( spc_model, offset, extremes );
                    [pees_last, jval, rc] = xb_gradient_descent( model, pees_last, minopt );
                    disp( ['Current jval :', num2str( jval )] );
                    if opts{1} == 'sticky'
                        __draw( binZ, model( pees_last ), h_data, extremes, offset );
                    else
                        __draw( binZ, model( pees_last ), h_data, extremes, offset, 'redraw' );
                    end
                    sleep( 0.5 );
                end
                isok = input( 'accept? >' );
                switch( isok )
                    case { 'y', 'Y', 'yes', 'Yes' }
                        pees = pees_last;
                end
            case { 'ama', 'animate-manual' }
                pees_last = __set_pees( pees, opts );
                pees_matrix = linspace( pees, pees_last, 120 );
                for ii=1:120
                    disp( [ 'Parameters: ',num2str( pees_matrix(:,ii)' ) ] );
                    model = __make_model( spc_model, offset, extremes );
                    jval = model( pees_matrix(:,ii) );
                    if opts{end} == 'sticky'
                        __draw( binZ, model( pees_matrix(:,ii) ), h_data, extremes, offset );
                    else
                        __draw( binZ, model( pees_matrix(:,ii) ), h_data, extremes, offset, 'redraw' );
                    end
                    sleep( 0.5 )
                end
                isok = input( 'accept? >' );
                switch( isok )
                    case { 'y', 'Y', 'yes', 'Yes' }
                        pees = pees_last;
                end
            case { 'x', 'extremes' }
                extremes(1) = str2num( opts{1} );
                extremes(2) = str2num( opts{2} );
            case { 'ok', 'OK' }
                go_on = 0;
            otherwise
                warning( 'Unknown command. Read the HALP!' );
        end
    pees = pees_last;
    J_cov = xb_covariance( model, pees );
    pee_errs = sqrt( diag( J_cov ) );
end

%utility to set parameters
function pees = __set_pees( opts )
    if numel( opts ) == 2
        pees(str2num( opts{1} )) = str2num( opts{2} );
    elseif numel( opts ) == 1
        eval( ['pees = ',opts{1}] );
    end
end

%utility to draw
function __draw( hmodel, binZ, hdata, extremes, offset, mm )
    %...
end

%utility to make the model
function model = __make_model( spc_model, offset, extremes )
    %...
end

%diplay the help
function __display_halp()
    %...
end
