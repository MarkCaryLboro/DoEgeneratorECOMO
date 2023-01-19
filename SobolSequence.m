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
        function obj = generate( obj, Opts )
            %--------------------------------------------------------------
            % Generate a design
            %
            % obj = obj.generate( "Name1", "value1", ..., "Name#", value# )
            %--------------------------------------------------------------
            arguments
                obj         (1,1)   SobolSequence
                Opts.Size   (1,1)   double         = 101
                Opts.Leap   (1,1)   double         = 7301
                Opts.Skip   (1,1)   double         = 49
            end
            D = height( obj.Factors );
            P = sobolset( D, "Leap", Opts.Leap, "Skip", Opts.Skip );
            Des = net( P, Opts.Size );

        end % generate
        
        function obj = scrambleDesign( obj, varargin )
        end % scrambleDesign
    end % concrete abstract method signatures
    
end % SobolSequence