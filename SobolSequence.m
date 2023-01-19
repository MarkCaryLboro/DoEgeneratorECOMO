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
        function obj = generate( obj, varargin )
        end % generate
        
        function obj = scrambleDesign( obj, varargin )
        end % scrambleDesign
    end % concrete abstract method signatures
    
end % SobolSequence