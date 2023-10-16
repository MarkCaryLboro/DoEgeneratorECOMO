classdef DoEgeneratorECOMO < handle
    % A class to handle mixed fixed and distributed parameter DoE
    % generation for the ECOMO model

    events
        DESIGN_AVAILABLE
    end % events

    properties ( Constant = true, Access = protected)
        Expected    string = [ "Name", "Units", "Fixed", "Lo", "Hi", "Sz", "Type", "Spline" ];
        ExpectCon = [ "name", "derivative", "type", "value", "x" ];
    end % constant properties

    properties ( Access = protected )
        NumPoints_  (1,1)    double                                         % Number of points in the design
        NumColDes_  (1,1)    double                                         % Number of colums in the design matrix
    end % protected properties

    properties ( SetAccess = protected )
        Design      (:,:)    double                                         % Design array
        DesignInfo  (:,:)    table                                          % Table of pointers to make it easy to populate the design table
        Bspline     (:,:)    table                                          % Table of bSplineTools objects (one row for each distributed parameter)
        Factors     (:,:)    table                                          % Factor details and type
        InitialSize (1,1)    double                                         % Initial size of experiment
        TubeLength  (1,1)    double  = 185.00                               % Length of the tube [mm]
        TubeIntDia  (1,1)    double  = 4.5                                  % Clean inner diameter of tube [mm]
        Constrained (1,1)    logical = false                                % True if design is constrained
        Scramble    (1,1)    logical   = false
    end % SetAccess protected

    properties ( SetAccess = protected )
        ExportReady (1,1)    logical = false                                % Set to true if experimental design is ready for export.
    end % private properties

    properties ( SetAccess = protected, Dependent = true )
        NumPoints           int64                                           % Number of points in the design
        NumFactors          int64                                           % Number of factors
        NumFixed            int64                                           % Number of fixed factors
        NumDist             int64                                           % Number of B-spline factors
        DistIdx             logical                                         % Logical index to distributed parameters
    end % Accessible dependent properties

    properties ( Access = private, Dependent = true )
    end % private dependent properties

    methods ( Abstract = true )
        obj = generate( obj, varargin )
    end % abstract method signatures

    methods
        function obj = addDesignPoint( obj, Data )
            %--------------------------------------------------------------
            % Add a data point to the current design
            %
            % obj = obj.addDesignPoint( Data );
            %
            % Input Arguments:
            %
            % Data  --> (double) row vector of new design data in
            %           engineering units.
            %--------------------------------------------------------------
            arguments
                obj     (1,1)               { mustBeNonempty( obj ) }
                Data    (1,:) double        
            end
            %--------------------------------------------------------------
            % Check the data format
            %--------------------------------------------------------------
            Ok = obj.checkDataFormat( Data, obj.NumColDes_ );
            assert( Ok, "Data must be numeric and have %3.0f columns",...
                            obj.NumColDes_ );
            obj.Design = [ obj.Design; Data ];
            obj.NumPoints_ = size( obj.Design, 1 );
        end % addDesignPoint

        function X = decode( obj, Xc, Name )
            %--------------------------------------------------------------
            % Decode from the interval [0,1] --> [a,b] for a fixed factor
            %
            % X = obj.decode( Xc, Name );
            %
            % Input Arguments:
            %
            % Xc    --> (double) Vector of coded data
            % Name  --> (string) Name of fixed parameter
            %--------------------------------------------------------------
            arguments
                obj     (1,1)        { mustBeNonempty( obj ) }
                Xc      (:,:) double { mustBeGreaterThanOrEqual( Xc, 0 ),...
                                       mustBeLessThanOrEqual( Xc, 1 ) }
                Name    (1,1) string { mustBeNonempty( Name ) }
            end
            %--------------------------------------------------------------
            % Find an exact match
            %--------------------------------------------------------------
            Idx = matches( obj.Factors.Properties.RowNames, Name );
            Ok = obj.Factors{ Idx, "Fixed" };
            assert( Ok, 'Factor "%s" is not a fixed parameter', Name);
            [ A, B ] = obj.getLimits( Idx );
            X = ( B - A ) .* Xc + A;
        end % decode

        function Xc = code( obj, X, Name )
            %--------------------------------------------------------------
            % Code from the interval [ a, b] --> [ 0, 1 ] for a fixed
            % factor.
            % 
            % Xc = obj.code( X, Name );
            %
            % Input Arguments:
            %
            % X     --> (double) Vector of uncoded data  
            % Name  --> (string) Name of fixed parameter
            %--------------------------------------------------------------
            arguments
                obj     (1,1)        { mustBeNonempty( obj ) }
                X      (:,1) double  { mustBeNonempty( X ) }
                Name    (1,1) string { mustBeNonempty( Name ) }
            end
            Idx = contains( obj.Factors.Name, Name );
            Ok = obj.Factors{ Idx, "Fixed" };
            assert( Ok, 'Factor "%s" is not a fixed parameter', Name);
            [ A, B ] = obj.getLimits( Idx );
            Xc = ( X - A ) ./ ( B - A );
        end % code

        function export( obj )
            %--------------------------------------------------------------
            % Export the design to the ECOMO model configuration class
            %
            % obj = export();
            %--------------------------------------------------------------
            arguments
                obj  (1,1)          { mustBeNonempty( obj ) }
            end
            if obj.ExportReady
                notify( obj, 'DESIGN_AVAILABLE' );
            else
                warning( 'Design not generated!' );
            end
        end % export
        
        function D = applyConstraints( obj, Sz, P )
            %--------------------------------------------------------------
            % Apply derivative constraints to the distributed parameter.
            % Spline evaluations outside the defined range are removed
            % from the experiment.
            % 
            % D = obj.applyConstraints( Sz, P );
            %
            % Input Arguments
            %
            % Sz  --> Desired size of design
            % P   --> Sobol set object
            %--------------------------------------------------------------
            arguments
                obj (1,1)                   { mustBeNonempty( obj ) }
                Sz  (1,1)  double
                P   (:,:)  sobolset                           
            end
            %--------------------------------------------------------------
            % Generate an unconstrained design to begin with
            %--------------------------------------------------------------
            D = net( P, Sz );
            %--------------------------------------------------------------
            % Point to the constrained parameters
            %--------------------------------------------------------------
            ConPtr = obj.Bspline.Constrained;
            %==============================================================
            % Apply the constraints
            %==============================================================
            ConFactors = obj.Bspline.Properties.RowNames( ConPtr );         % Retieve factor names with constraints
            ConFactors = string( ConFactors );                                                           
            for Q = 1:numel( ConFactors )
                %----------------------------------------------------------
                % Retrieve the constraint structure
                %----------------------------------------------------------
                Name = ConFactors( Q );
                S = obj.Bspline{ Name, "Constraint" };
                if iscell( S )
                    S = S{ : };
                end
                V = obj.evalSplineConstraint( Name, S, Sz, 101 );
                %----------------------------------------------------------
                % Add constrained points to the design
                %----------------------------------------------------------
                [ Cidx, Kidx ] = obj.getParameterPointers( Name );
                Start = min( Cidx );
                Finish = max( Kidx );
                D( :, Start:Finish ) = V;
            end % /Q
        end % applyConstraints

        function Y = evalSpline( obj, X, Name, Coeff, Knot )
            %--------------------------------------------------------------
            % Evaluate the B-spline at the coordinates specified
            %
            % Y = obj.evalSpline( Name, Coeff, Knots );
            %
            % Input Arguments:
            %
            % X     --> (double) Vector of axial dimensions at which to 
            %                    evaluate the B-spline  
            % Name  --> (string) Name of distributed parameter
            % Coeff --> (double) Basis function coefficients
            %--------------------------------------------------------------
            arguments
                obj   (1,1)          { mustBeNonempty( obj ) }
                X     (:,1)  double  { mustBeNonempty( X ) }
                Name  (1,1)  string  { mustBeNonempty( Name ) }
                Coeff (:,1)  double  { mustBeNonempty( Coeff ) }
                Knot  (:,1)  double  { mustBeNonempty( Knot ) }
            end
            %--------------------------------------------------------------
            % Check name of distributed parameter is valid
            %--------------------------------------------------------------
            Ok = matches( Name, obj.Factors.Properties.RowNames );
            assert( Ok, 'Parameter "%s" not defined', Name );
            %--------------------------------------------------------------
            % Check that parameter is of type "distributed"
            %--------------------------------------------------------------
            Idx = matches( obj.Factors.Properties.RowNames, Name );
            Ok = ~obj.Factors{ Idx, "Fixed" };
            assert( Ok, 'Parameter "%s" cannot be of type "Fixed"', Name );
            %--------------------------------------------------------------
            % Capture the relevant B-spline object
            %--------------------------------------------------------------
            B = obj.Bspline{ Name, "Object"};
            %--------------------------------------------------------------
            % Clip input range to the interval [B.a, B.b]
            %--------------------------------------------------------------
            X( X < B.a ) = B.a;
            X( X > B.b ) = B.b;
            %--------------------------------------------------------------
            % Check coefficient vector dimensionality is correct & assign
            %--------------------------------------------------------------
            Ok = ( numel( Coeff ) == B.nb );
            assert( Ok, '"%s" spline must have %2.0f coefficients', ...
                                                            Name, B.nb);
            B.alpha = Coeff;
            %--------------------------------------------------------------
            % Check knot vector dimensionality is correct & assign
            %--------------------------------------------------------------
            Ok = ( numel( Knot ) == B.k );
            assert( Ok, '"%s" spline must have %2.0f knots', ...
                                                            Name, B.k);
            B.n = Knot;
            %--------------------------------------------------------------
            % Calculate the response variables 
            %--------------------------------------------------------------
            Y = B.eval( X );
        end % evalSpline

        function obj = setPipeGeometry( obj, D, L )
            %--------------------------------------------------------------
            % Define the pipe geometric properties: the clean internal
            % diameter and the tube lenghth. Both dimensions are assumed to
            % be in [mm].
            %--------------------------------------------------------------
            arguments
                obj (1,1)            { mustBeNonempty( obj ) }
                D   (1,1)   double   { mustBePositive( D ) } = 4.50
                L   (1,1)   double   { mustBePositive( L ) } = 185.00
            end
            obj.TubeLength = L;
            obj.TubeIntDia = D;
            obj.clearDesign();
        end % setPipeGeometry

        function obj = clearDesign( obj )
            %--------------------------------------------------------------
            % Set the design array to empty and reset Scramble property to
            % false
            %
            % obj = obj.clearDesign();
            %--------------------------------------------------------------
            obj.Design = [];
            obj.Scramble = false;
        end % clearDesign

        function obj = addFactor( obj, S, C )
            %--------------------------------------------------------------
            % Add a factor to the factor table property. If factor is
            % already defined, details are overwritten with the new
            % information supplied.
            %
            % obj = obj.addFactor( S, C );
            %
            % Input Arguments:
            %
            % S      --> (struct) Multidimensional Structure defining factor 
            %                     properties with fields:
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
            %
            % Note each dimension of S must define a different factor.
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   
                S       (1,:)   struct   { mustBeNonempty( S ) }
                C       (1,:)   struct = struct.empty  
            end
            Q = max( size( S ) );
            %--------------------------------------------------------------
            % Check for unique factors
            %--------------------------------------------------------------
            Name = [ S.Name ];
            Ok = ( numel( unique( Name ) ) == Q );
            assert( Ok, "Factor names may not be repeated!");
            for N = 1:Q
                S( N ).Name = replace( S( N ).Name, " ", "" );              % Deblank channel names
            end
            %--------------------------------------------------------------
            % Clear the design
            %--------------------------------------------------------------
            obj = obj.clearDesign();
            %--------------------------------------------------------------
            % Check all necessary fields are present & parse
            %--------------------------------------------------------------
            obj.Factors = obj.parseFactor( S );                       
            %--------------------------------------------------------------
            % Create B-spline array for distributed factors
            %--------------------------------------------------------------
            obj = obj.createBsplineTable( C );
            %--------------------------------------------------------------
            % Generate the design information linking columns of the design
            % matrix to factor values
            %--------------------------------------------------------------
            obj = obj.genDesignInfo();
        end % addFactor

        function Out = decodeDesign( obj, Des )
            %--------------------------------------------------------------
            % Decode the design and return the result in engineering units
            %
            % Des = obj.decodeDesign( Des );
            %
            % Input Arguments:
            %
            % Des --> (double) Coded design on the interval [0,1]
            %--------------------------------------------------------------
            Out = zeros( size( Des ) );
            for Q = 1:obj.NumFactors
                Name = obj.Factors.Properties.RowNames( Q );
                %----------------------------------------------------------
                % Decode the columns by factor
                %----------------------------------------------------------
                Col = obj.DesignInfo{ Name, "Coefficients" };
                if iscell( Col )
                    Col = Col{ : };
                end
                if obj.DistIdx( Q )
                    %------------------------------------------------------
                    % Distributed factor. Process coefficients and then
                    % knots
                    %------------------------------------------------------
                    Coeffc = Des( :, Col );
                    Coeff = obj.decodeSplineCoeff( Name, Coeffc );
                    Out( :, Col ) = Coeff;
                    %------------------------------------------------------
                    % Now decode the knots
                    %------------------------------------------------------
                    Col = obj.DesignInfo{ Name, "Knots" };
                    if iscell( Col )
                        Col = Col{ : };
                    end
                    Kc = Des( :, Col );
                    B = obj.Bspline{ Name, "Object" };
                    K = B.decode( Kc );
                    Out( :, Col ) = K;
                else
                    %------------------------------------------------------
                    % Fixed factor
                    %------------------------------------------------------
                    Out( :, Col ) = obj.decode( Des( :, Col ), Name );
                end
            end % Q
        end % decodeDesign
    end % ordinary methods

    methods ( Access = protected )
        function Idx = intializeConstraintPtr( obj, N )
            %--------------------------------------------------------------
            % Return a logical array with columns corresponding to
            % constraints and rows to design points.
            %
            % Idx = obj.intializeConstraintPtr( N );
            %
            % Input Arguments:
            %
            % N --> (double) number of design points
            %--------------------------------------------------------------
            C = obj.Bspline.Constraint;
            %--------------------------------------------------------------
            % Determine the number of constraints
            %--------------------------------------------------------------
            NumCon = 0;
            for Q = 1:max( size( C ) )
                S = C{ Q };
                if ~ isempty( S )
                    NumCon = NumCon + max( size( S ) );
                end
            end
            %--------------------------------------------------------------
            % Define the logical output array
            %--------------------------------------------------------------
            Idx = false( N, NumCon );
        end % intializeConstraintPtr

        function obj = genDesignInfo( obj )
            %--------------------------------------------------------------
            % Create a table decoding the design table information
            %
            % obj = obj.genDesignInfo();
            %--------------------------------------------------------------
            N = obj.NumFactors;
            D = cell( N, 2 );
            Finish = 0;
            for Q = 1:N
                %----------------------------------------------------------
                % Parse the column information
                %----------------------------------------------------------
                if obj.Factors.Fixed( Q )
                    %------------------------------------------------------
                    % Parse a fixed factor
                    %------------------------------------------------------
                    Start = Finish + 1;
                    Sz = obj.Factors.Sz( Q, : );
                    Finish = prod( Sz ) + Start - 1;
                    Coeff = ( Start:Finish );
                    Knots = nan;
                    D( Q,: ) = { Coeff, Knots };
                else
                    %------------------------------------------------------
                    % Parse the distributed factor
                    %------------------------------------------------------
                    Name = string( obj.Bspline.Properties.RowNames{ Q } );
                    [ D( Q,: ), Finish ] = obj.parseDistributed( Name,...
                                                         Finish );
                end
            end
            %--------------------------------------------------------------
            % Update the design info table property with the pointers
            %--------------------------------------------------------------
            D = cell2table( D );
            D.Properties.VariableNames = [ "Coefficients", "Knots" ];
            D.Properties.RowNames = obj.Factors.Properties.RowNames;
            obj.DesignInfo = D;
        end % genDesignInfo

        function obj = createBsplineTable( obj, C )
            %--------------------------------------------------------------
            % Create a B-spline representation for each distributed
            % parameter
            %
            % obj = obj.createBsplineArray( C );
            %
            % Input Arguments:
            %
            % C     --> (struct) Structure defining constraint properties
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
            %--------------------------------------------------------------
            arguments
                obj (1,1)           { mustBeNonempty( obj ) }
                C   (1,:) struct = struct.empty
            end
            %--------------------------------------------------------------
            % Set the constrained flag to the appropriate state
            %--------------------------------------------------------------
            obj.Constrained = ~isempty( C );
            %--------------------------------------------------------------
            % Fetch names of distributed factors
            %--------------------------------------------------------------
            RowNames = obj.Factors.Properties.RowNames( obj.DistIdx );      % Names of distributed factors
            RowNames = string( RowNames );
            N = obj.NumDist;
            %--------------------------------------------------------------
            % Initialise list of input variables
            %--------------------------------------------------------------
            Xname = strings( N, 1 );
            %--------------------------------------------------------------
            % Parse the spline data and create the necessary Bspline
            % objects.
            %--------------------------------------------------------------
            if ( N > 0 )
                B = bSplineTools.empty( N, 0 );                             % Create array of Bsplines
                for Q = 1:N
                    %------------------------------------------------------
                    % Retrieve spline parameter structure
                    %------------------------------------------------------
                    S = table2struct( obj.Factors( RowNames( Q ),...
                                                               "Spline") );
                    S = S.Spline;
                    %------------------------------------------------------
                    % Form list of input variable names
                    %------------------------------------------------------
                    Xname( Q, 1 ) = S.X;
                    %------------------------------------------------------
                    % Parse one-dimensional spline
                    %------------------------------------------------------
                    B( Q ) = obj.makeOneDimensionalSpline( S );
                end
                NumBasis = double( [B.nb].' );
                NumKnots = double( [B.k].' );
                NumPar = NumBasis + NumKnots;
                B = table( B( : ), NumBasis, NumKnots, NumPar, Xname );
                B.Properties.RowNames = RowNames;
                B.Properties.VariableNames = [ "Object", "NumBasis",...
                                           "NumKnots", "NumPar", "Xname"];
                %----------------------------------------------------------
                % Parse derivative constraints if supplied
                %----------------------------------------------------------
                T = obj.parseBsplineConstraints( C );
                obj.Bspline = horzcat( B, T );
            end
        end % createBsplineTable

        function T = parseBsplineConstraints( obj, C )
            %--------------------------------------------------------------
            % Process any constraints to the Bspline table property
            %
            % T = obj.parseBsplineConstraints( C );
            %
            % Input Arguments:
            %
            % C     --> (struct) Structure defining constraint properties
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
            %--------------------------------------------------------------
            if isempty( C )
                C = obj.getDefaultConstraints();
            end
            Ok = all( obj.conFieldCheck( C ) );                         % Check necessary fields are present
            assert( Ok, "Missing information for factor %s", C.name);   % Throw an error if there is missing information
            %----------------------------------------------------------
            % Parse the constraint information. Allow for multiple
            % constraints being applicable to a given factor.
            %----------------------------------------------------------
            Names = string( obj.Factors.Properties.RowNames );
            Names = Names( obj.DistIdx );
            T = table( 'Size',  [ obj.NumDist, 2 ], ...
                'VariableTypes', {'cell', 'logical'} );
            T.Properties.VariableNames = [ "Constraint", "Constrained" ];
            T.Properties.RowNames = Names;
            Cnames = string( { C(:).name } );
            for Q = 1:numel( Names )
                %----------------------------------------------------------
                % Add the constraint structure to the table
                %----------------------------------------------------------
                Idx = matches( Cnames, Names( Q ) );
                Con = ~isempty( C( Idx ) );
                T( Names( Q ), : ) = cell2table( { { C( Idx ) }, Con });
            end %/Q
        end % parseBsplineConstraints

        function T = parseFactor( obj, S )
            %--------------------------------------------------------------
            % A method to parse the factor definition structure
            %
            % T = obj.parseFactor( S );
            %
            % Input Arguments:
            %
            % S     --> (struct) Structure defining factor properties
            %                    with fields:
            %
            %            Name  - (string) Name of factor
            %            Units - (string) Factor units
            %            Fixed - (logical) True if fixed factor. False
            %                    if distributed factor.
            %            Lo    - (double) Low natural limit for factor
            %            Hi    - (double) High natural limit for factor            
            %--------------------------------------------------------------
            Ok = all( obj.fieldCheck( S ) );                                % Check necessary fields are present
            assert( Ok, "Missing information for factor %s", S.Name);       % Throw an error if there is missing information
            if ( max( size( S ) ) == 1 )
                % Scalar structure
                T = struct2table( S, "AsArray", true );
            else
                T = struct2table( S );
            end
            T.Properties.RowNames = string( T.Name );
            Idx = ~contains( T.Properties.VariableNames, "Name" );
            T = T( :, obj.Expected( Idx ) );
        end % parseFactor

        function Ok = conFieldCheck( obj, C )
            %--------------------------------------------------------------
            % Output logical value to indicate expected field is present.
            %
            % Ok = obj.conFieldCheck( C );
            %
            % C     --> (struct) Structure defining factor properties
            %                    with expected fields:
            %
            %   Name  - (string) Name of factor
            %   derivative  --> set to 0,1 or 2 {0} to specify the spline
            %                   derivative to which the constraint applies.
            %   type          --> set to '==','>=' or '<='
            %   value         --> constraint bound value
            %   x             --> x-ordinates at which constraints apply.
            %                     Leave empty to specify all training
            %                     x-ordinates.
            %--------------------------------------------------------------
            Ok = false( size( obj.ExpectCon ) );
            F = fieldnames( C );
            N = numel( obj.ExpectCon );
            for Q = 1:N
                Ok( Q ) = contains( obj.ExpectCon( Q ), F, "IgnoreCase",...
                                                            true );
            end                     
        end % conFieldCheck

        function Ok = fieldCheck( obj, S )
            %--------------------------------------------------------------
            % Output logical value to indicate expected field is present.
            %
            % Ok = obj.fieldCheck( S );
            %
            % S     --> (struct) Structure defining factor properties
            %                    with expected fields:
            %
            %            Name  - (string) Name of factor
            %            Units - (string) Factor units
            %            Fixed - (logical) True if fixed factor. False
            %                    if distributed factor.
            %            Lo    - (double) Low natural limit for factor
            %            Hi    - (double) High natural limit for factor
            %            Type  - (string) Dentoes "Parameter" or "Boundary"
            %                    condition
            %--------------------------------------------------------------
            N = numel( obj.Expected );
            Ok = false( size( obj.Expected ) );
            F = fieldnames( S );
            for Q = 1:N
                Ok( Q ) = contains( obj.Expected( Q ), F, "IgnoreCase",...
                    true );
            end           
        end % fieldCheck

        function C = getDefaultConstraints( obj )
            %--------------------------------------------------------------
            % Return constraints structure with defined fields all empty,
            % except the name
            %
            % Default field names are:
            %
            %   name        - (string) Name of factor
            %   derivative  - set to 0,1 or 2 {0} to specify the spline
            %                 derivative to which the constraint applies.
            %   type        - set to '==','>=' or '<='
            %   value       - constraint bound value
            %   x           - x-ordinates at which constraints apply.
            %                 Leave empty to specify all training
            %                 x-ordinates.
            %--------------------------------------------------------------
            Names = string( obj.Factors.Properties.RowNames( obj.DistIdx ) );
            for Q = numel( Names ):-1:1
                %----------------------------------------------------------
                % Define the default constraint structure
                %----------------------------------------------------------
                C( Q ).name = Names( Q );
                C( Q ).derivative = [];
                C( Q ).type = [];
                C( Q ).value = [];
                C( Q ).x = [];
            end
        end % getDefaultConstraints

        function [ Coeff ] = decodeSplineCoeff( obj, Name, Coeffc )
            %--------------------------------------------------------------
            % Decode the spline coefficients for the cited distributed
            % parameter
            %
            % Coeff = decodeSplineCoeff( obj, Name, Coeffc )
            %
            % Input Arguments:
            %
            % Name   --> (string) Name of spline factor
            % Coeffc --> (double) Coded coefficients [0,1]
            %--------------------------------------------------------------
            Idx = matches( obj.Factors.Properties.RowNames, Name );
            A = obj.Factors{ Idx, "Lo" };
            if iscell( A )
                A = A{ : };
            end            
            B = obj.Factors{ Idx, "Hi" };
            if iscell( B )
                B = B{ : };
            end
            Coeff =  ( B - A ) .* Coeffc +  A; 
        end % decodeSplineCoeff

        function [ Cidx, Kidx ] = getParameterPointers( obj, Name )
            %--------------------------------------------------------------
            % Return pointers to the coefficient and knot locations in the
            % design matrix. If a fixed parameters Knots is a NaN.
            %
            % [ Coeff, Knots ] = obj.getParameterPointers( Name );
            %
            % Input Arguments:
            % 
            % Name --> (string) Name of parameter
            %
            % Output Arguments:
            %
            % Cidx   --> (double) columns of the design matrix containing
            %            the parameter coefficients
            % Kidx   --> (double) columns corresponding to the knots
            %            associated with a distributed parameter. Will be
            %            NaN if a fixed parameter.
            %--------------------------------------------------------------
            Kidx = obj.DesignInfo{ Name, "Knots" };                         % retrieve the knots
            if iscell( Kidx )
                Kidx = Kidx{ : };
            end
            Cidx = obj.DesignInfo{ Name, "Coefficients" };                  % retrieve the coefficients
            if iscell( Cidx )
                Cidx = Cidx{ : };
            end
        end % getParameterPointers
    end % protected methods

    methods ( Access = private )
        function [ Out, Finish ] = parseDistributed( obj, Name, Finish )
            %--------------------------------------------------------------
            % Output a cell array of dimension {1,2}. First cell contains
            % the columns for the spline coefficients. The second cell
            % contains the columns pertaiing to the knots.
            %
            % [ Out, Finish ] = obj.parseDistributed( Name, Finish );
            %
            % Input Arguments:
            %
            % Name   - (string) Name of variable
            % Finish - (double) Last column pointer value known
            %
            % Output Arguments:
            %
            % Out    - (cell) array of output columns
            % Finish - (doulbe) Updated column pointer
            %--------------------------------------------------------------
            Str = [ "NumBasis", "NumKnots" ];
            Out = cell( 1,2 );
            for Q = 1:2
                Start = Finish + 1;
                Finish = Start + obj.Bspline{ Name, Str( Q ) } - 1;
                Out{ Q } = Start:Finish;
            end
        end % parseDistributed

        function [ A, B ] = getLimits( obj, Idx )
            %--------------------------------------------------------------
            % Return the Hi & Lo limits for a factor
            %
            % [ A, B ] = obj.getLimits( Idx );
            %
            % Input Arguments:
            %
            % Idx   --> (logical) Pointer to limit information
            %
            % Output Arguments:
            %
            % A     --> (double) Low limits
            % B     --> (double) Hi limits
            %--------------------------------------------------------------
            A = obj.Factors.Lo( Idx );
            if iscell( A )
                A = A{ : };
            end
            B = obj.Factors.Hi( Idx );
            if iscell( B )
                B = B{ : };
            end     
            %--------------------------------------------------------------
            % Handle the matrix case by reshaping the A and B matrices
            %--------------------------------------------------------------
            A = reshape( A, 1, numel( A ) );
            B = reshape( B, 1, numel( B ) );            
        end % getLimits
    end % private methods

    methods ( Static = true, Access = protected )
        function Ok = checkDataFormat( X, C )
            %----------------------------------------------------------------------
            % Check to see if the number of data columns is correct and all data is
            % numeric.
            %
            % Ok = obj.checkDataFormat( X, C )
            %
            % Input Arguments:
            %
            % X     --> Data to augment the design with
            % C     --> Number of expected columns
            %----------------------------------------------------------------------
            Ok = true;
            Ok = Ok & all( isnumeric( X ) );
            Ok = Ok & ( numel( X ) == C );
        end % checkDataFormat

        function B = makeOneDimensionalSpline( S )
            %--------------------------------------------------------------
            % Return a correctly configured one-dimensional bSplineTools 
            % object
            %
            % B = obj.makeOneDimensionalSpline( S );
            %
            % Input Arguments:
            %
            % S     --> (struct) B-spline configuration structure with
            %                    fields:
            %
            %                    X   - (string) Input factor name(s)
            %                    M   - (int8) Spline 
            %                    K   - (cell)
            %                    Xlo - (double) Low limit(s) for 
            %                          x-factor(s) range
            %                    Xhi - (double) High limit(s) for 
            %                          x-factor(s) range
            %--------------------------------------------------------------
            DK = linspace( S.Xlo, S.Xhi, S.K + 2 ).';
            DK = DK( 2:end-1 );
            B = bSplineTools( S.M - 1 , DK, S.Xlo, S.Xhi );
        end % makeOneDimensionalSpline
    end % static and protected methods

    methods
        function N = get.NumPoints( obj )
            N = int64( obj.NumPoints_ );                                    % Number of points in the design
        end

        function N = get.NumFactors( obj )                             
            N = height( obj.Factors );                                      % Number of factors
        end

        function N = get.NumFixed( obj )                                           
            N = sum( obj.Factors{ :, "Fixed" } );                           % Number of fixed factors
        end

        function N = get.NumDist( obj )
            N = sum( ~obj.Factors{ :, "Fixed" } );                          % Number of B-spline factors
        end

        function Idx = get.DistIdx( obj )
            Idx = ~obj.Factors{ :, "Fixed" };                               % Point to the distributed parameters
        end
    end % set/get methods

end % DoEgeneratorECOMO