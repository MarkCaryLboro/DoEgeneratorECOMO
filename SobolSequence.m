classdef SobolSequence < DoEgeneratorECOMO

    
    methods
        function obj = SobolSequence( Factors )
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
            %--------------------------------------------------------------
            if nargin > 0 && all( obj.fieldCheck( Factors ) )
                obj = obj.addFactor( Factors );
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
                obj         (1,1)   SobolSequence
                N           (1,1)   int64          { mustBePositive( N ) } = 101;
                Opts.Leap   (1,1)   double         = max( primes( 7301 ) );
                Opts.Skip   (1,1)   double         = max( primes( 49 ) );
                Opts.Con    (1,1)   logical        = false
            end
            obj.NumPoints_ = N;
            D = sum( obj.Bspline.NumPar ) + double( obj.NumFixed );
            P = sobolset( D, "Leap", Opts.Leap, "Skip", Opts.Skip );
            Des = net( P, obj.NumPoints );
            if Opts.Con
                %----------------------------------------------------------
                % Identify feasible combinations for the distributed
                % factors
                %----------------------------------------------------------
                Des = obj.applyConstraints( Des );
            end
            obj.Design = Des;
        end % generate
        
        function obj = scrambleDesign( obj, varargin )
        end % scrambleDesign
    end % concrete abstract method signatures
    
end % SobolSequence