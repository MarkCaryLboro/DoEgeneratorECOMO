classdef DoEgeneratorECOMO < handle
    % A class to handle mixed fixed and distributed parameter DoE
    % generation for the ECOMO model

    events
        DESIGN_AVAILABLE
    end % events

    properties ( Constant = true, Access = protected)
        Expected    string = [ "Name", "Units", "Fixed", "Lo", "Hi", "Sz", "Type" ];
    end % constant properties

    properties ( Access = protected )
        NumPoints_  (1,1)    double                                         % Number of points in the design
        NumColDes_  (1,1)    double                                         % Number of colums in the design matrix
        BOptLh      (1,1)                                                   % UPDATE event listener handle
    end % protected properties

    properties ( SetAccess = protected )
        Design      (:,:)    double                                         % Design array
        DesignInfo  (:,:)    table                                          % Table of pointers to make it easy to populate the design table
        Bspline     (:,4)    table                                          % Table of bSplineTools objects (one row for each distributed parameter)
        Factors     (:,:)    table                                          % Factor details and type
        Scramble    (1,1)    logical = false                                % Set to true to apply scramble to design
        TubeLength  (1,1)    double  = 185.00                               % Length of the tube [mm]
        TubeIntDia  (1,1)    double  = 4.5                                  % Clean inner diameter of tube [mm]
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
            Idx = startsWith( obj.Factors.Name, Name ) & ...
                  endsWith( obj.Factors.Name, Name );
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

        function obj = applyConstraints( obj, Sz, Des )
            %--------------------------------------------------------------
            % Apply interval constraints to the distributed parameters.
            % Spline evaluations outside the supplied range are removed
            % from the experiment.
            % 
            % obj = obj.applyConstraints( Des )
            %
            % Input Arguments
            %
            % Des --> Current design set
            % Sz  --> Desired size of design
            %--------------------------------------------------------------
            arguments
                obj (1,1)         { mustBeNonempty( obj ) }
                Sz  (1,1)  double
                Des (:,:)  double                           = obj.Design
            end
            N = 1:obj.NumFactors;
            X = linspace( 0, obj.TubeLength, 1001 );                        % Define axial dimension
            Idx = true( size( Des, 1 ), obj.NumFactors );
            for Q = N
                %----------------------------------------------------------
                % Evaluate the splines
                %
                % 1. Extract coefficients and decode
                % 2. Extract knots and decode
                % 3. Evaluate spline
                %----------------------------------------------------------
                if ~obj.Factors{ Q, "Fixed" }
                    %------------------------------------------------------
                    % Distributed factor
                    %------------------------------------------------------
                    Name = obj.Factors{ Q, "Name" };
                    try
                        Kidx = obj.DesignInfo{ Q, "Knots" }{:};
                        Cidx = obj.DesignInfo{ Q, "Coefficients" }{:};
                    catch
                        Kidx = obj.DesignInfo{ Q, "Knots" };
                        Cidx = obj.DesignInfo{ Q, "Coefficients" };
                    end
                    B = obj.Bspline{ Name, "Object" };
                    K = B.decode( Des( :,Kidx ) );
                    
                    C = obj.decodeSplineCoeff( Name, Des( :, Cidx ) );
                    Lo = obj.Factors.Lo( Q );
                    if iscell( Lo )
                        Lo = Lo{ : };
                    end
                    Hi = obj.Factors.Hi( Q );
                    if iscell( Hi )
                        Hi = Hi{ : };
                    end                    
                    Ok = false( size( Des, 1 ), 1 );
                    for R = 1:size( C, 1 )
                        %--------------------------------------------------
                        % Evaluate the spline
                        %--------------------------------------------------
                        Y = obj.evalSpline( X, Name, C( R,: ), K( R,: ) );
                        Ok( R ) = all( ( Y >= Lo ) & ( Y <= Hi ) );
                    end
                    Idx( :, Q ) = Ok;
                end
            end
            %--------------------------------------------------------------
            % Now select the valid design points
            %--------------------------------------------------------------
            Idx = all( Idx, 2 );
            Des = Des( Idx,: );
            NumFeasible = sum( Idx );
            if ( NumFeasible > Sz )
                Des = Des( 1:Sz,: );
            end
            obj.Design = Des;
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
            Ok = contains( Name, obj.Factors.Name );
            assert( Ok, 'Parameter "%s" not defined', Name );
            %--------------------------------------------------------------
            % Check that parameter is of type "distributed"
            %--------------------------------------------------------------
            Idx = contains( obj.Factors.Name, Name );
            Ok = ~obj.Factors{ Idx, "Fixed" };
            assert( Ok, 'Parameter "%s" cannot be of type "Fixed"', Name );
            %--------------------------------------------------------------
            % Clip axial tube dimension
            %--------------------------------------------------------------
            X( X > obj.TubeLength ) = obj.TubeLength;
            X( X < 0 ) = 0;
            %--------------------------------------------------------------
            % Capture the relevant B-spline object
            %--------------------------------------------------------------
            B = obj.Bspline{ Name, "Object"};
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

        function obj = addFactor( obj, S, Opts )
            %--------------------------------------------------------------
            % Add a factor to the factor table property. If factor is
            % already defined, details are overwritten with the new
            % information supplied.
            %
            % obj.addFactor( S, Name, Value );
            %
            % Input Arguments:
            %
            % S      --> (struct) Multidimensional Structure defining factor 
            %                     properties with fields:
            %
            %             Name  - (string) Name of factor
            %             Units - (string) Factor units
            %             Fixed - (logical) True if fixed factor. False
            %                     if distributed factor.
            %             Lo    - (double) Low natural limit for factor
            %             Hi    - (double) High natural limit for factor
            %             Sz    - (int64) Size of corresponding lookup
            %                     table (only for distributed parameters).
            %                     Set to 1 if the factor is fixed.
            %             Type  - (string) Set to "Parameter" to denote am 
            %                     identifable coefficient of "Boundary" to
            %                     denote a boundary condition.
            %
            % Note each dimension of S must define a different factor.
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
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   
                S       (1,:)   struct   { mustBeNonempty( S ) }
                Opts.M  (:,1)   int8 = 4
                Opts.K  (:,1)   int8 = 2
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
            T = obj.parseFactor( S );            
            T = T( :, obj.Expected );
            %--------------------------------------------------------------
            % Add factors to current list
            %--------------------------------------------------------------
            if ~isempty( obj.Factors )
                Idx = contains(  [ obj.Factors.Name ], T.Name,...
                    'IgnoreCase', true );
                obj.Factors( Idx,: ) = T( Idx,: );
                obj.Factors = vertcat(  obj.Factors, T( ~Idx,: ) );
            else
                obj.Factors = T;
            end
            %--------------------------------------------------------------
            % Duplicate knot or order data if a scalar
            %--------------------------------------------------------------
            if isscalar( Opts.M )
                Opts.M = repmat( Opts.M, obj.NumDist, 1 );
            end
            if isscalar( Opts.K )
                Opts.K = repmat( Opts.K, obj.NumDist, 1 );
            end            
            %--------------------------------------------------------------
            % Create B-spline array for distributed factors
            %--------------------------------------------------------------
            obj = obj.createBsplineTable( Opts.M, Opts.K );
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
                Name = obj.Factors.Name( Q );
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
                    [ D( Q,: ), Finish ] = obj.parseDistributed( obj.Factors.Name( Q ) ,...
                                                         Finish );
                end
            end
            %--------------------------------------------------------------
            % Update the design info table property with the pointers
            %--------------------------------------------------------------
            D = cell2table( D );
            D.Properties.VariableNames = [ "Coefficients", "Knots" ];
            D.Properties.RowNames = obj.Factors.Name;
            obj.DesignInfo = D;
        end % genDesignInfo

        function obj = createBsplineTable( obj, M, K )
            %--------------------------------------------------------------
            % Create a B-spline representation for each distributed
            % parameter
            %
            % obj = obj.createBsplineArray( M, K, Names );
            %
            % Input Arguments:
            %
            % M     --> (int8) Spline order (1 <= M <= 4)
            % K     --> (int8) Number of knots (1 <= K <= 7)
            %
            % Note both M and K may be scalers or vectors. If vecctors,
            % they must be vectors of the same length. This permits the
            % user to modulate the complexity of the axial distribution
            % of each distributed parameter as required. If M and K are
            % both scalars then every distributed parameter is assigned  
            % order M and number of knots K.
            %--------------------------------------------------------------
            arguments
                obj (1,1)                   { mustBeNonempty( obj ) }
                M   (:,1) int8  { mustBeGreaterThan(M,0), ...
                                  mustBeLessThan(M,5)} = 4;
                K   (:,1) int8  { mustBeGreaterThan(K,0), ...
                                  mustBeLessThan(K,8)} = 2;
            end
            %--------------------------------------------------------------
            % Fetch names of distributed factors
            %--------------------------------------------------------------
            D = ~obj.Factors{ :, "Fixed" };                                 % Point to distributed factors
            RowNames = obj.Factors.Name( D );                               % Names of distributed factors
            N = sum( D );
            %--------------------------------------------------------------
            % Parse the spline data and create the necessary Bspline
            % objects.
            %--------------------------------------------------------------
            if ( N > 0 )
                B( N, 1 ) = bSplineTools();
                for Q = 1:N
                    DK = linspace( 0, obj.TubeLength, K( Q ) + 2 ).';
                    DK = DK( 2:end-1 );
                    B( Q ) = bSplineTools( ( M( Q ) - 1 ), DK,...
                        0, obj.TubeLength );
                end
                NumBasis = double( [B.nb].' );
                NumKnots = double( [B.k].' );
                NumPar = NumBasis + NumKnots;
                obj.Bspline = table( B, NumBasis, NumKnots, NumPar );
                obj.Bspline.Properties.RowNames = RowNames;
                obj.Bspline.Properties.VariableNames = [ "Object",...
                                 "NumBasis", "NumKnots", "NumPar"];
            end
        end % createBsplineTable

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
            T = struct2table( S );
        end % parseFactor

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
            Ok = false( size( obj.Expected ) );
            F = fieldnames( S );
            N = numel( obj.Expected );
            for Q = 1:N
                Ok( Q ) = contains( obj.Expected( Q ), F, "IgnoreCase",...
                    true );
            end
        end % fieldCheck
    end % protected methods

    methods ( Access = private )
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
            Idx = contains( obj.Factors.Name, Name );
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