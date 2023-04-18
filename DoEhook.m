classdef DoEhook < handle
    %----------------------------------------------------------------------
    % A class to generate configuration data required to run the ECOMO
    % model. The class is a listener for the SobolSequence class, which
    % generates a DoE for the fixed and distributed parameters identified
    % from data in the model.
    %----------------------------------------------------------------------
    
    events
        RUN_EXPERIMENT                                                      % Run the full experiment
    end

    properties ( SetAccess = protected)
        Lh          (1,1)                                                   % Listener handle for DESIGN_AVAILABLE event
        Uh          (1,1)                                                   % Listener handle for UPDATE event
        ConfigFile  (1,1) string                                            % ECOMO model configuration file
        ParTable    (:,:) table                                             % Parameter table in engineering units
        TubeLength  (1,1) double                                            % Length of the tube [mm]
        TubeIntDia  (1,1) double                                            % Clean inner diameter of tube [mm]
        NumFactors  (1,1) int64                                             % Number of factors
        NumFixed    (1,1) int64                                             % Number of fixed factors
        NumDist     (1,1) int64                                             % Number of B-spline factors
        DistIdx     (1,:) logical                                           % Logical index to distributed parameters
        Design      (:,:) double                                            % Durrent experimental design in engineering units
        Type        (1,:) string                                            % String array of variable types (Parameter or Boundary)
    end % Protected properties

    properties ( SetAccess = protected, Dependent = true )
        NumPoints   (1,1) int64                                             % Number of points in the design
    end % dependent properties

    methods
        function obj = DoEhook()
            %--------------------------------------------------------------
            % Define the DoEhook to listen for and process the 
            % DESIGN_AVAILABLE and UPDATE events
            %
            % obj = DoEhook();
            %--------------------------------------------------------------

            %--------------------------------------------------------------
            % Capture the configuration file for the ECOMO model
            %--------------------------------------------------------------
            [ File, Path ] = uigetfile(".m","Select ECOMO model configuration file", "MultiSelect","off");
            obj.ConfigFile = fullfile( Path, File );
            Ok = true;
            if ( File == 0 )
                Ok = false;
            end
            assert( Ok, "Must select ECOMO model configuration file!" );
        end % DoEhook
    end % Constructor method

    methods
        function obj = addDesignAvailableListener( obj, Src )
            %--------------------------------------------------------------
            % Add a listener for the DESIGN_AVAILABLE event broadcast by a
            % SobolSequence object.
            %
            % obj = obj.addUpdateListener( Src )
            %
            % Input Arguments:
            %
            % Src --> SobolSequence object
            %--------------------------------------------------------------
            arguments
                obj (1,1) DoEhook       
                Src (1,1) SobolSequence { mustBeNonempty( Src ) }
            end
            %--------------------------------------------------------------
            % Define the listener
            %--------------------------------------------------------------
            obj.Lh = addlistener( Src, "DESIGN_AVAILABLE",...
                        @( SrcObj, Evnt )obj.eventCb( SrcObj, Evnt ) );
            obj.Design = Src.Design;
        end % addDesignAvailableListener

        function runSimulation( obj )
            %--------------------------------------------------------------
            % Trigger the RUN_EXPERIMENT event to execute the simulation
            %
            % obj.runSimulation();
            %--------------------------------------------------------------
            arguments
                obj (1,1) DoEhook { mustBeNonempty( obj ) }
            end
            notify( obj, 'RUN_EXPERIMENT' );
        end % runSimulation

        function obj = addProcessNewQueryListener( obj, Src )
            %--------------------------------------------------------------
            % Add a listener for the PROCESS_NEW_QUERY event broadcast by
            % object
            %
            % obj = obj.addProcessNewQueryListener( Src )
            %
            % Input Arguments:
            %
            % Src --> ecomoInterface object
            %--------------------------------------------------------------
            arguments
                obj (1,1)   DoEhook 
                Src (1,1)   ecomoInterface  { mustBeNonempty( Src ) }
            end
            %--------------------------------------------------------------
            % Define the listener
            %--------------------------------------------------------------
            obj.Uh = addlistener( Src, "PROCESS_NEW_QUERY",...
                     @( Src, Evnt )obj.eventCbProcessNewQuery( Src,...
                                                                  Evnt ) );
        end % addProcessNewQueryListener

        function obj = eventCbProcessNewQuery( obj, Src, E )
            %--------------------------------------------------------------
            % PROCESS_NEW_QUERY event listener
            %
            % This function creates a table listing the parameters in the
            % experiment. Distributed parameters are converted to lookup
            % tables.
            %--------------------------------------------------------------     
            arguments
                obj   (1,1) DoEhook { mustBeNonempty( obj )}                % DoEhook object
                Src   (1,1)         { mustBeNonempty( Src ) }               % SobolSequence object
                E     (1,1)         { mustBeNonempty( E ) }                 % EventData object 
            end
            %--------------------------------------------------------------
            % Event check
            %--------------------------------------------------------------
            Ename = string( E.EventName );
            Ok = contains( "PROCESS_NEW_QUERY", Ename );
            assert( Ok, 'Not processing the %s event supplied', Ename );
            obj = obj.augmentParTable( Src );
            notify( obj, 'RUN_EXPERIMENT' );
        end % eventCbProcessNewQuery

        function obj = eventCb( obj, Src, E )
            %--------------------------------------------------------------
            % DESIGN_AVAILABLE event listener
            %
            % This function creates a table listing the parameters in the
            % experiment. Distributed parameters are converted to lookup
            % tables.
            %--------------------------------------------------------------
            arguments
                obj   (1,1) DoEhook { mustBeNonempty( obj ) }               % DoEhook object
                Src   (1,1)         { mustBeNonempty( Src ) }               % SobolSequence object
                E     (1,1)         { mustBeNonempty( E ) }                 % EventData object 
            end
            %--------------------------------------------------------------
            % Event check
            %--------------------------------------------------------------
            Ename = string( E.EventName );
            Ok = contains( "DESIGN_AVAILABLE", Ename );
            assert( Ok, 'Not processing the %s event supplied', Ename );
            %--------------------------------------------------------------
            % Set geometric parameters
            %--------------------------------------------------------------
            obj.TubeIntDia = Src.TubeIntDia / 1000;                         % convert to [m]
            obj.TubeLength = Src.TubeLength / 1000;                         % Convert to [m]
            %--------------------------------------------------------------
            % Set DoE Parameters
            %--------------------------------------------------------------
            obj.NumFactors = Src.NumFactors;                                % Number of factors
            obj.NumFixed = Src.NumFixed;                                    % Number of fixed factors
            obj.NumDist = Src.NumDist;                                      % Number of B-spline factors
            obj.DistIdx = Src.DistIdx;                                      % Logical index to distributed parameters
            %--------------------------------------------------------------
            % Create the table of physical parameter values
            %--------------------------------------------------------------
            obj = obj.createParTable( Src );
        end % eventCB

        function obj = setSimulated( obj, N, State)
            %--------------------------------------------------------------
            % Set the Nth row of the "Simulated" column in the ParTable
            % property to the desired State
            %
            % obj = obj.setSimulated( N, State );
            %
            % Input Arguments:
            %
            % N         --> (int64) Row number(s) to set
            % State     --> (logical) State to set {false}
            %--------------------------------------------------------------
            arguments
                obj
                N     (1,:)  int64      { mustBeNonempty( N ) }
                State (1,:)  logical    = false
            end
            N = double( N );
            obj.ParTable.Simulated( N ) = State;
        end % setSimulated
    end % ordinary methods

    methods
        function N = get.NumPoints( obj )
            N = int64( size( obj.Design, 1 ) );
        end
    end % Get/Set methods

    methods ( Access = protected )
        function obj = augmentParTable( obj, Src )
            %--------------------------------------------------------------
            % Augment the parameter table with the new query
            %
            % obj = obj.augmentParTable( Src );
            %
            % Input Arguments
            %
            % Src --> Event source object (ecomoInterface object).
            %--------------------------------------------------------------
            S = obj.Lh.Source{ : };
            S = S.addDesignPoint( Src.B.Xnext );   
            obj.Design = S.Design;
            Simulated = [ obj.ParTable.Simulated; false ];
            obj = obj.createParTable( S );
            obj.ParTable.Simulated = Simulated;
        end

        function obj = createParTable( obj, Src )
            %--------------------------------------------------------------
            % Creates the parameter table for running the experiment
            %
            % obj = obj.createParTable( Src )
            %
            % Input Arguments
            %
            % Src --> Event source object.
            %--------------------------------------------------------------
            Didx = Src.DistIdx;
            Fnames = string( Src.Factors.Name );
            Npts = Src.NumPoints;
            VarTypes = obj.createVarTypes( Didx );
            T = table( 'Size', [ Npts, numel( VarTypes ) ],...
                'VariableTypes', VarTypes );
            T.Properties.VariableNames = [ Fnames; "Simulated" ];
            obj.Type = obj.getParameterTypes( Src ); 
            for R = 1:Npts
                %----------------------------------------------------------
                % Fill out the table a row at a time
                %----------------------------------------------------------
                for Q = 1:Src.NumFactors
                    if Didx( Q )
                        %--------------------------------------------------
                        % Distributed factor. Calculate lookup table
                        %--------------------------------------------------
                        LookUp = obj.makeLookUp( Src, Fnames( Q ), R );
                        T( R, Q ) = { LookUp }; 
                    else
                        %--------------------------------------------------
                        % Fixed parameter
                        %--------------------------------------------------
                        Col = Src.DesignInfo{ Fnames( Q ), "Coefficients" };
                        if iscell( Col )
                            Col = Col{ : };
                        end
                        Sz = Src.Factors.Sz( Q,: );
                        if iscell( Sz )
                            Sz = cell2mat( Sz );
                        end
                        D = { Src.Design( R, Col ) };
                        if iscell( D )
                            D = cell2mat( D );
                        end
                        D = reshape( D, Sz );
                        T( R, Q ) = array2table( { D } );
                    end
                end % Q
            end % R
            obj.ParTable = T;
        end % createParTable
    end % protected methods

    methods ( Access = private, Static = true )
        function V = getParameterTypes( Src )
            %--------------------------------------------------------------
            % Create a string array detailing whether each design factor is
            % either an identifiable parameter or a boundary condition
            %
            % V = obj.getParameterTypes( Src )
            %
            % Input arguments:
            %
            % Src       --> Event source object.
            %--------------------------------------------------------------
            V = Src.Factors.Type;
            V = reshape( V, 1, numel( V ) );
        end % getParameterTypes

        function VarTypes = createVarTypes( Didx )
            %--------------------------------------------------------------
            % Create a cell array of variable types
            %
            % VarTypes = obj.createVarTypes( Didx );
            %
            % Input Arguments:
            %
            % Didx  --> (Logical) true if distributed factor
            %--------------------------------------------------------------
            N = numel( Didx );
            VarTypes = cell( 1, N );
            for Q = 1:N
                VarTypes{ Q } = 'cell';
            end % Q
            %--------------------------------------------------------------
            % Add the simulated column
            %--------------------------------------------------------------
            VarTypes{ end + 1 } = 'logical';
        end % createVarTypes

        function LookUp = makeLookUp( Src, Name, RunNumber )
            %--------------------------------------------------------------
            % Evaluate the B-spline and calculate the lookup table values
            %
            % LookUp = obj.makeLookUp( Src, Name );
            %
            % Input Arguments
            %
            % Src       --> Event source object.
            % Name      --> Name of parameter to process
            % RunNumber --> Current experimental design run
            %--------------------------------------------------------------
            Info = Src.DesignInfo( Name, : );
            Idx = contains( Src.Factors.Name, Name );
            Sz = Src.Factors.Sz( Idx, : );
            LookUp = zeros( Sz );
            LookUp( 1,: ) = linspace( 0, Src.TubeLength, max( Sz ) );       % Tube axial dimension to evaluate spline
            %--------------------------------------------------------------
            % Point to the columns in the design table defining the spline
            % parameters
            %--------------------------------------------------------------
            Coeff = Info.Coefficients;
            Knot = Info.Knots;
            if iscell( Coeff )
                Coeff = Info.Coefficients{ : };
            end
            if iscell( Knot )
                Knot = Info.Knots{ : };
            end
            %--------------------------------------------------------------
            % Capture the requested spline parameter values
            %--------------------------------------------------------------
            Coeff = Src.Design( RunNumber, Coeff );
            Knot = Src.Design( RunNumber, Knot );
            %--------------------------------------------------------------
            % Calculate the response
            %--------------------------------------------------------------
            Y = Src.evalSpline( LookUp( 1,: ), Name, Coeff, Knot );
            LookUp( 2,: ) = reshape( Y, 1, numel( Y ) );
            LookUp( 1,: ) = 0.001 * LookUp( 1,: );                          % Convert to [m] for simulation
            %--------------------------------------------------------------
            % Transpose if a Boundary Condition
            %--------------------------------------------------------------
%             if contains( Src.Factors.Type( Idx ), "Boundary" )
% %                 LookUp = LookUp.';
%             end
             LookUp = LookUp.';
        end
    end % private methods
end % DoEhook