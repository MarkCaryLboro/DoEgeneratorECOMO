classdef SobolSequence < DoEgeneratorECOMO
    %----------------------------------------------------------------------
    % A class to generate an experimental design, using a sobol sequence
    % for the ECOMO model identification approach using Bayesian 
    % Optimisation.
    %----------------------------------------------------------------------
    methods
        function obj = SobolSequence( Factors, Opts )
            %--------------------------------------------------------------
            % Class constructor
            %
            % obj = SobolSequence( Factors );
            %
            % Input Arguments:
            %
            % Factors   --> (struct) Structure defining factor properties
            %               with fields:
            %
            %               Name  - (string) Name of factor
            %               Units - (string) Factor units
            %               Fixed - (logical) True if fixed factor. False
            %                       if distributed factor.
            %               Lo    - (double) Low natural limit for factor
            %               Hi    - (double) High natural limit for factor
            %
            % Optional Arguments:
            % 
            % Name   --> (string), may be either "M" for spline order or
            %            "K" for number of knots.
            % Value -->  (double) If a scaler, then all distributed
            %            parameters will be assigned the same order or 
            %            number of knots. If a vector, with elements equal 
            %            to the number of distributed factors, then each 
            %            distributed factor is assigned a unique order or 
            %            number of knots respectively.
            %
            % Examples: obj = SobolSequence( Factors, "M", 4, "K", 2 )
            %           obj = SobolSequence( Factors, "M", [3, 4], "K", 2 )
            %           obj = SobolSequence( Factors, "K", [2, 3], "M", 4 )
            %--------------------------------------------------------------
            arguments
                Factors (1,:) struct
                Opts.M  (:,1) int8      = 4
                Opts.K  (:,1) int8      = 2
            end
            if nargin > 0 && all( obj.fieldCheck( Factors ) )
                obj = obj.addFactor( Factors, "K", Opts.K, "M", Opts.M );
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
            % "
            %--------------------------------------------------------------
            arguments
                obj           (1,1)   SobolSequence
                N             (1,1)   int64        { mustBePositive( N ) } = 101;
                Opts.Leap     (1,1)   double       = max( primes( 7301 ) );
                Opts.Skip     (1,1)   double       = max( primes( 49 ) );
                Opts.Con      (1,1)   logical      = false
                Opts.Scramble (1,1) logical        = false
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
            if Opts.Con
                %----------------------------------------------------------
                % Identify feasible combinations for the distributed
                % factors
                %----------------------------------------------------------
                Des = net( P, N );                                          % Coded on interval [ 0,1 ]
                obj = obj.applyConstraints( N, Des );                       % Retain only feasible combinations 
                obj.NumPoints_ = size( obj.Design, 1 );
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