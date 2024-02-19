classdef DoEhook < handle
    %----------------------------------------------------------------------
    % A class to generate configuration data required to run the ECOMO
    % model. 
    %----------------------------------------------------------------------
    properties ( SetAccess = protected)
        ParTable    (:,:) table                                             % Parameter table in engineering units
        TubeLength  (1,1) double                                            % Length of the tube [mm]
        TubeIntDia  (1,1) double                                            % Clean inner diameter of tube [mm]
        NumFactors  (1,1) int64                                             % Number of factors
        NumFixed    (1,1) int64                                             % Number of fixed factors
        NumDist     (1,1) int64                                             % Number of B-spline factors
        DistIdx     (1,:) logical                                           % Logical index to distributed parameters
        Type        (1,:) string                                            % String array of variable types (Parameter or Boundary)
        DesObj      (1,1)                                                   % Experimental design object
    end % Protected properties

    properties ( SetAccess = protected, Dependent = true )
        Design      (:,:)                                                   % Durrent experimental design in engineering units
        NumPoints   (1,1)                                                   % Number of points in the design
    end % dependent properties

    methods
        function obj = transDesign( obj, Src )
            %--------------------------------------------------------------
            % Translate the design into the corresponding ECOMO model
            % simulation format. This function creates a table listing the 
            % parameters in the experiment. Distributed parameters are 
            % converted to lookup tables.
            %--------------------------------------------------------------
            arguments
                obj   (1,1) DoEhook             { mustBeNonempty( obj ) }   % DoEhook object
                Src   (1,1) SobolSequence       { mustBeNonempty( Src ) }   % SobolSequence object
            end
            %--------------------------------------------------------------
            % Save the SobolSequence object
            %--------------------------------------------------------------
            obj.DesObj = Src;
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
        end % transDesign

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
                obj   (1,1) DoEhook             { mustBeNonempty( obj ) }   % DoEhook object
                N     (1,:)  int64      { mustBeNonempty( N ) }
                State (1,:)  logical    = false
            end
            N = double( N );
            obj.ParTable.Simulated( N ) = State;
        end % setSimulated

        function obj = augmentParTable( obj, Data )
            %--------------------------------------------------------------
            % Augment the parameter table with the new query
            %
            % obj = obj.augmentParTable( Data );
            %
            % Input Arguments
            %
            % Data --> DoE points to be translated to ECOMO model
            %          parameters
            %--------------------------------------------------------------
            arguments
                obj   (1,1) DoEhook             { mustBeNonempty( obj ) }   % DoEhook object
                Data  (:,:) double              { mustBeNonempty( Data ) }  % Data points to add
            end
            %--------------------------------------------------------------
            % Check the data has the appropriate number of columns
            %--------------------------------------------------------------
            N = size( obj.Design, 2 );
            Ok = ( size( Data, 2 ) == N );
            assert( Ok, "Data must be numeric and have %3.0f columns", N );
            %--------------------------------------------------------------
            % Add the data point
            %--------------------------------------------------------------
            S = obj.DesObj.addDesignPoint( Data );
            Simulated = [ obj.ParTable.Simulated; false ];
            obj = obj.createParTable( S );
            obj.ParTable.Simulated = Simulated;
        end % augmentParTable
    end % ordinary methods

    methods
        function N = get.NumPoints( obj )
            N = obj.DesObj.NumPoints;
        end

        function D = get.Design( obj )
            D = obj.DesObj.Design;
        end
    end % Get/Set methods

    methods ( Access = protected )
        function obj = overwritePartable( obj, X )
            %--------------------------------------------------------------
            % Overwrite the last entry in the parameter table
            %
            % obj = obj.overwritePartable( X );
            %
            % Input Arguments:
            %
            % X --> (double) Replacement parameter values
            %--------------------------------------------------------------
            S = obj.DesObj;
            S = S.overwriteDesignPoint( X );
            obj.Design = S.Design;
            Simulated = [ obj.ParTable.Simulated( 1:( end-1 ) ); false ];
            obj = obj.createParTable( S );
            obj.ParTable.Simulated = Simulated;
        end % overwritePartable

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
            Fnames = string( Src.Factors.Properties.RowNames );
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
                        LookUp = obj.makeLookUp( Fnames( Q ), R );
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

        function LookUp = makeLookUp( obj, Name, RunNumber )
            %--------------------------------------------------------------
            % Evaluate the B-spline and calculate the lookup table values
            %
            % LookUp = obj.makeLookUp( Src, Name, RunNumber );
            %
            % Input Arguments
            %
            % Name      --> Name of parameter to process
            % RunNumber --> Current experimental design run
            %--------------------------------------------------------------
            Src = obj.DesObj;
            Idx = matches( Src.Bspline.Properties.RowNames, Name );
            Tensor = Src.Bspline.Tensor( Idx );
            if Tensor
                %----------------------------------------------------------
                % Two-dimensional tensor product B-spline
                %----------------------------------------------------------
                LookUp = obj.tensorProdBsplineMakeLookUp( Name, RunNumber );
            else               
                %----------------------------------------------------------
                % One-dimensional B-spline
                %----------------------------------------------------------
                LookUp = obj.bSplineMakeLookUp( Name, RunNumber );
            end
        end % makeLookUp

        function LookUp = bSplineMakeLookUp( obj, Name, RunNumber )
            %--------------------------------------------------------------
            % Generate a 1-d lookup
            %
            % LookUp = obj.bSplineMakeLookUp( Name, RunNumber );
            %
            % Input Arguments:
            %
            % Name      --> (string) Variable name
            % RunNumber --> ( double) Current experimental design run
            %--------------------------------------------------------------
            Src = obj.DesObj;
            Info = Src.DesignInfo( Name, : );
            Idx = matches( Src.Factors.Properties.RowNames, Name );
            Sz = Src.Factors.Sz( Idx, : );
            %--------------------------------------------------------------
            % Retrieve Low and High input limits & define lookup table
            % input vector
            %--------------------------------------------------------------
            B = Src.Bspline{ Name, "Object" };
            if iscell( B )
                B = B{ : };
            end
            Lo = B.a;
            Hi = B.b;
            LookUp = zeros( Sz );
            LookUp( 1,: ) = linspace( Lo, Hi, max( Sz ) );                  % inputs at which to evaluate spline
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
            InputVar = Src.Bspline{ Name, "Xname" };
            Idx = matches( InputVar,"" );
            InputVar = InputVar( ~Idx );
            if matches( InputVar, "x", 'IgnoreCase', true)
                %----------------------------------------------------------
                % Convert to [m] for simulation
                %----------------------------------------------------------
                LookUp( 1,: ) = 0.001 * LookUp( 1,: );
            end
            LookUp = LookUp.';
        end % bSplineMakeLookUp

        function LookUp = tensorProdBsplineMakeLookUp( obj, Name, RunNumber )
            %--------------------------------------------------------------
            % Generate a 2-d lookup based on the tensor product spline
            %
            % LookUp = obj.tensorProdBsplineMakeLookUp( Name );
            %
            % Input Arguments:
            %
            % Name      --> (string) Variable name
            % RunNumber --> ( double) Current experimental design run
            %--------------------------------------------------------------
            Src = obj.DesObj;
            Info = Src.DesignInfo( Name, : );
            Idx = matches( Src.Factors.Properties.RowNames, Name );
            Sz = Src.Factors.Sz( Idx, : );
            T = Src.Bspline{ Name, "Object" };
            if iscell( T )
                T = T{ : };
            end
            Lo = T.A;
            Hi = T.B;
            %--------------------------------------------------------------
            % Create data points to evaluate the spline at
            %--------------------------------------------------------------
            X = linspace( Lo( 1 ), Hi( 1 ), Sz( 1 ) );
            Y = linspace( Lo( 2 ), Hi( 2 ), Sz( 2 ) );
            [ X, Y ] = meshgrid( X, Y );
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
            Knot = T.convertKnotSequences( Knot );
            %--------------------------------------------------------------
            % Calculate the response
            %--------------------------------------------------------------
            T = T.setAlpha( Coeff );
            T = T.setKnotSequences( Knot );
            Z = T.eval( [ X( :) Y( : ) ] );
            LookUp = [ X( : ), Y( : ), Z ];
        end % tensorProdBsplineMakeLookUp
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
    end % private & static methods
end % DoEhook