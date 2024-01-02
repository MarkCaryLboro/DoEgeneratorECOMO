classdef SobolSequence < DoEgeneratorECOMO
    %----------------------------------------------------------------------
    % A class to generate an experimental design, using a sobol sequence
    % for the ECOMO model identification approach using Bayesian 
    % Optimisation.
    %----------------------------------------------------------------------
    properties ( SetAccess = protected )
        Leap        (1,1) double    = max( primes( 49 ) )
        Skip        (1,1) double    = max( primes( 7301 ) )
    end % protected methods

    methods
        function obj = SobolSequence( Factors, Con )
            %--------------------------------------------------------------
            % Class constructor
            %
            % obj = SobolSequence( Factors, Con );
            %
            % Input Arguments:
            %
            % Factors --> (struct) Multidimensional Structure defining factor 
            %                      properties with fields:
            %
            %             Name   - (string) Name of factor
            %             Units  - (string) Factor units
            %             Fixed  - (logical) True if fixed factor. False
            %                      if distributed factor.
            %             Lo     - (double) Low natural limit for factor
            %             Hi     - (double) High natural limit for factor
            %             Sz     - (int64) Size of corresponding lookup
            %                      table (only for distributed parameters).
            %                      Set to [ 1, 1 ] if the factor is fixed.
            %             Type   - (string) Set to "Parameter" to denote am 
            %                      identifable coefficient of "Boundary" to
            %                      denote a boundary condition.
            %             Spline - (struct) configuration for a distributed
            %                      parameter, with fields:
            %
            %                      X   - (string) Input factor name(s)
            %                      M   - (int8) Spline 
            %                      K   - (cell)
            %                      Xlo - (double) Low limit(s) for 
            %                            x-factor(s) range
            %                      Xhi - (double) High limit(s) for 
            %                            x-factor(s) range
            % Con   --> (struct) Structure defining constraint properties
            %                    with fields
            %
            %   Name        - (string) Name of factor
            %   derivative  - set to 0,1 or 2 {0} to specify the spline
            %                 derivative to which the constraint applies.
            %   type        - set to '==','>=' or '<='
            %   value       - constraint bound value
            %   x           - x-ordinates at which constraints apply.
            %                 Leave empty to specify all training
            %                 x-ordinates.
            %
            %--------------------------------------------------------------
            arguments
                Factors (1,:) struct  = struct.empty  
                Con     (1,:) struct  = struct.empty;
            end
            S = warning;
            if isempty( Factors )
                warning( "off" );
            end
            if nargin > 0 && all( obj.fieldCheck( Factors ) )
                obj = obj.addFactor( Factors, Con );
            else
                warning( "Missing fields in input structure" );
            end
            warning( S );
        end % SobolSequence
    end % Constructor method

    methods 
        function obj = generate( obj, N, Opts )
            %--------------------------------------------------------------
            % Generate a design
            %
            % obj = obj.generate( N );
            %
            % Input Arguments:
            %
            % N     --> (int64) Number of design points.
            %
            % Optional Arguments:
            %
            % Leap      --> First sample in Sobol sequence. Set this to a
            %               big prime number.
            % Skip      --> Sample rate in sequence. Set to a prime number
            % Scramble  --> Set this to true if scrambled sequence is
            %               required (recommended)
            %--------------------------------------------------------------
            arguments
                obj           (1,1)   SobolSequence
                N             (1,1)   int64        { mustBePositive( N ) } = 101
                Opts.Leap     (1,1)   double       = obj.Leap
                Opts.Skip     (1,1)   double       = obj.Skip
                Opts.Scramble (1,1)   logical      = false
            end
            obj = obj.clearDesign();
            %--------------------------------------------------------------
            % Determine number of parameters
            %--------------------------------------------------------------
            D = table2cell( obj.DesignInfo );
            D = cellfun( @max, D );
            D = max( D( : ) );                                              % Nan is ignored.
            %--------------------------------------------------------------
            % Now create the sobol set
            %--------------------------------------------------------------
            obj.Leap = Opts.Leap;
            obj.Skip = Opts.Skip;
            P = sobolset( D + 1, "Leap", obj.Leap, "Skip", obj.Skip,...
                          "PointOrder", "standard" );
            if Opts.Scramble
                %----------------------------------------------------------
                % Scramble the design
                %----------------------------------------------------------
                P = scramble( P, 'MatousekAffineOwen' );
                obj.Scramble = true;
            end 
            %--------------------------------------------------------------
            % Set initial size of experiment
            %--------------------------------------------------------------
            obj.InitialSize = N;
            %--------------------------------------------------------------
            % Apply constraints if defined
            %--------------------------------------------------------------
            if obj.Constrained
                %----------------------------------------------------------
                % Identify feasible combinations for the distributed
                % factors
                %----------------------------------------------------------
                Des = obj.applyConstraints( N, P );                         % Retain only feasible combinations 
                obj.NumPoints_ = size( Des, 1 );
            else
                obj.NumPoints_ = N;
                Des = net( P, double( obj.NumPoints ) );
            end
            Des = Des( :, 2:end );
            obj.NumColDes_ = D;
            %--------------------------------------------------------------
            % Decode the design
            %--------------------------------------------------------------
            Des = obj.decodeDesign( Des );
            obj.Design = Des;
            obj.ExportReady = true;
        end % generate
    end % concrete abstract method signatures

    methods 
        function V = evalSplineConstraint( obj, Name, S, Sz, N )
            %--------------------------------------------------------------
            % Return a constrained design of the appropriate size for a 
            % single distributed factor. All constraints applied to a
            % factor will be satisfied.
            %
            % V =  obj.evalSplineConstraint( S, Sz, N );
            %
            % Input Arguments:
            %
            % Name  --> (string) name of constrained factor
            % S     --> (struct) spline constraint definition structure
            % Sz    --> (int64) Desired design size
            % N     --> (int64) number of points to evaluate spline
            %           constraint at.
            %--------------------------------------------------------------
            arguments
                obj     (1,1) SobolSequence { mustBeNonempty( obj ) }
                Name    (1,1) string        { mustBeNonempty( Name ) }
                S       (:,:) struct        { mustBeNonempty( S ) }
                Sz      (1,1) int64         { mustBeNonempty( Sz ) }
                N       (1,1) int64         { mustBeNonempty( N ) } = 101
            end
            NumCon = max( size( S ) );                                      % Number of constraints
            %--------------------------------------------------------------
            % Fetch the B-spline object
            %--------------------------------------------------------------
            B = obj.Bspline{ Name, "Object" };
            if iscell( B )
                B = B{ : };
            end
            %--------------------------------------------------------------
            % Create the evaluation vector (spline input)
            %--------------------------------------------------------------
            X = linspace( B.a, B.b, N ).';
            %--------------------------------------------------------------
            % Determine the number of parameters & establish pointers
            %--------------------------------------------------------------
            D = B.nb + B.k;
            Cidx = 1:B.nb;
            Kidx = ( B.nb + 1 ):( B.nb + B.k ); 
            V = zeros( Sz, D );
            Finish = 0;
            %--------------------------------------------------------------
            %  Generate a constrained design
            %--------------------------------------------------------------
            W = waitbar( 0, "Selecting Feasible Points");
            while ( Finish < Sz )
                %----------------------------------------------------------
                % Create a seperate Sobol set for the constrained variable
                %----------------------------------------------------------
                P = sobolset( D, "Leap", obj.Leap, "Skip", obj.Skip );
                if obj.Scramble
                    % Apply a scramble if specified
                    P = scramble( P, 'MatousekAffineOwen' );
                end
                %----------------------------------------------------------
                % Generate a design & decode the parameters
                %----------------------------------------------------------
                Des = net( P, Sz );                                         % Coded on interval [ 0,1 ]
                Coeff = obj.decodeSplineCoeff( Name, Des( :,Cidx) );
                Knot = B.decode( Des( :, Kidx ) );
                %----------------------------------------------------------
                % Initialise feasible design point pointer
                %----------------------------------------------------------
                Ptr = false( Sz, 1 );
                for Idx = 1:Sz
                    %------------------------------------------------------
                    % Evaluate the spline constraint for each point in
                    % the design
                    %------------------------------------------------------
                    B.alpha = Coeff( Idx, : );
                    B.n = Knot( Idx, : );
                    Feasible = false( N, NumCon );
                    for II = 1:NumCon
                        Dx = S( II ).derivative;
                        Type = S( II ).type;
                        Val = S( II ).value;
                        C = B.calcDerivative( X, Dx );
                        switch Type
                            case { ">=", "=>" }
                                Feasible( :, II ) = ( C >= Val );
                            case { "<=", "=<"}
                                Feasible( :, II ) = ( C <= Val );
                            case "=="
                                Feasible( :, II ) = ( C == Val );
                            otherwise
                                Msg = 'Unrecognised constraint type "%s"';
                                error( Msg, Type );
                        end
                    end
                    Feasible = all( Feasible, 2 );
                    if all( Feasible )
                        %--------------------------------------------------
                        % Add a feasible point to the design
                        %--------------------------------------------------
                        Ptr( Idx ) = true;
                    end
                end %/Idx
                Start = Finish + 1;
                Finish = Start + sum( Ptr ) - 1;
                V ( Start:Finish, ( 1:D ) ) = Des( Ptr, : );
                Msg = sprintf( "Selected %4.0f out of %4.0f feasible points", Finish, Sz);
                waitbar( double( Finish ) / double( Sz ), W, Msg );
            end %/while
            %--------------------------------------------------------------
            % Clip the design to the right size
            %--------------------------------------------------------------
            V = V( 1:Sz, : );
            delete( W );
        end % evalSplineConstraint
    end % ordinary methods
end % SobolSequence