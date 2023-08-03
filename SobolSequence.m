classdef SobolSequence < DoEgeneratorECOMO
    %----------------------------------------------------------------------
    % A class to generate an experimental design, using a sobol sequence
    % for the ECOMO model identification approach using Bayesian 
    % Optimisation.
    %----------------------------------------------------------------------
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
                Factors (1,:) struct  { mustBeNonempty( Factors ) }  
                Con     (1,:) struct  = struct.empty;
            end
            if nargin > 0 && all( obj.fieldCheck( Factors ) )
                obj = obj.addFactor( Factors, Con );
            else
                warning( "Missing fields in input structure" );
            end
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
            % N     --> (int64) Number of design loints.
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
                N             (1,1)   int64        { mustBePositive( N ) } = 101;
                Opts.Leap     (1,1)   double       = max( primes( 7301 ) );
                Opts.Skip     (1,1)   double       = max( primes( 49 ) );
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
            % Now create the design
            %--------------------------------------------------------------
            P = sobolset( D, "Leap", Opts.Leap, "Skip", Opts.Skip );
            if Opts.Scramble
                %----------------------------------------------------------
                % Scramble the design
                %----------------------------------------------------------
                P = scramble( P, 'MatousekAffineOwen' );
                obj.Scramble = true;
            end
            %--------------------------------------------------------------
            % Apply constraints if defined
            %--------------------------------------------------------------
            if obj.Constrained
                %----------------------------------------------------------
                % Identify feasible combinations for the distributed
                % factors
                %----------------------------------------------------------
                Des = net( P, N );                                          % Coded on interval [ 0,1 ]
                Des = obj.applyConstraints( N, Des );                       % Retain only feasible combinations 
                obj.NumPoints_ = size( Des, 1 );
            else
                obj.NumPoints_ = N;
                Des = net( P, obj.NumPoints );
            end
            obj.NumColDes_ = D;
            %--------------------------------------------------------------
            % Decode the design
            %--------------------------------------------------------------
            Des = obj.decodeDesign( Des );
            obj.Design = Des;
            obj.ExportReady = true;
        end % generate
    end % concrete abstract method signatures
end % SobolSequence