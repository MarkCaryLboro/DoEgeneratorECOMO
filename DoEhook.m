classdef DoEhook < handle
    %----------------------------------------------------------------------
    % A class to generate configuration data required to run the ECOMO
    % model. The class is a listener for the SobolSequence class, which
    % generates a DoE for the fixed and distributed parameters identified
    % from data in the model.
    %----------------------------------------------------------------------

    properties ( Access = protected)
        Lh (1,1)                                                            % Listener handle
        N  (1,1)    int64                                                   % Number of lookup table cells
    end % Protected properties

    methods
        function obj = DoEhook( Src )
            %--------------------------------------------------------------
            % Define the DoEhook to listen for and process the 
            % DESIGN_AVAILABLE event
            %
            % obj = DoEhook( Src );
            %
            % Input Arguments:
            %
            % Src --> Event source object.
            %--------------------------------------------------------------
            arguments
                Src (1,1) SobolSequence { mustBeNonempty( Src ) }
            end
            obj.Lh = addlistener( Src, "DESIGN_AVAILABLE",...
                        @( SrcObj, Evnt )obj.eventCb( SrcObj, Evnt ) );
        end % DoEhook
    end % Constructor method

    methods
        function eventCb( obj, Src, E )
            arguments
                obj   (1,1) DoEhook { mustBeNonempty( obj )} 
                Src   (1,1)         { mustBeNonempty( Src ) }
                E     (1,1)         { mustBeNonempty( E ) }
            end
            Src.Factors
        end % eventCB
    end % ordinary methods
end % DoEhook