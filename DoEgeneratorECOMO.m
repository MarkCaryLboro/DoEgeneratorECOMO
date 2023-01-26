classdef DoEgeneratorECOMO < handle
    % A class to handle mixed fixed and distributed parameter DoE
    % generation for the ECOMO model

    events
        DESIGN_AVAILABLE
    end % events

    properties ( Constant = true, Access = protected)
        Expected    string = [ "Name", "Units", "Fixed", "Lo", "Hi" ];
    end % constant properties

    properties ( Access = protected )
        DesignInfo  (:,:)    table                                          % Table of pointers to make it easy to populate the design table
        Design      (:,:)    double                                         % Design array
        NumPoints_  (1,1)   double                                          % Number of points in the design
        NumColDes_  (1,1)   double                                          % Number of colums in the design matrix
    end % protected properties

    properties ( SetAccess = protected )
        Bspline     (:,4)    table                                          % Table of bSplineTools objects (one row for each distributed parameter)
        Factors     (:,:)    table                                          % Factor details and type
        Scramble    (1,1)    logical = false                                % Set to true to apply scramble to design
        TubeLength  (1,1)    double  = 185.00                               % Length of the tube [mm]
        TubeIntDia  (1,1)    double  = 4.5                                  % Clean inner diameter of tube [mm]
    end % SetAccess protected

    properties ( Access = private )
    end % private properties

    properties ( SetAccess = protected, Dependent = true )
        NumPoints           int64                                           % Number of points in the design
        NumFactors          int64                                           % Number of factors
        NumFixed            int64                                           % Number of fixed factors
        NumDist             int64                                           % Number of B-spline factors
    end % Accessible dependent properties

    methods ( Abstract = true )
        obj = generate( obj, varargin )
    end % abstract method signatures

    methods
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
                Xc      (:,1) double { mustBeGreaterThanOrEqual( Xc, 0 ),...
                                       mustBeLessThanOrEqual( Xc, 1 ) }
                Name    (1,1) string { mustBeNonempty( Name ) }
            end
            Idx = contains( obj.Factors.Name, Name );
            Ok = obj.Factors{ Idx, "Fixed" };
            assert( Ok, 'Factor "%s" is not a fixed parameter', Name);
            A = obj.Factors.Lo( Idx );
            B = obj.Factors.Hi( Idx );
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
            A = obj.Factors.Lo( Idx );
            B = obj.Factors.Hi( Idx );
            Xc = ( X - A ) ./ ( B - A );
        end % code

        function obj = applyConstraints( obj, Des )
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
            %--------------------------------------------------------------
            arguments
                obj (1,1)         { mustBeNonempty( obj ) }
                Des (:,:)  double                           = obj.Design
            end
            N = 1:obj.NumFactors;
            Idx = true( size( Des, 1 ), 1 );
            for Q = 1:N
                %----------------------------------------------------------
                % Evaluate the splines
                %----------------------------------------------------------

            end
        end % applyConstraints

        function Y = evalSpline( obj, X, Name, Coeff, Knots )
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
            % Knots --> (double) Strictly increasing knot sequence
            %--------------------------------------------------------------
            arguments
                obj   (1,1)          { mustBeNonempty( obj ) }
                X     (:,1)  double  { mustBeNonempty( X ) }
                Name  (1,1)  string  { mustBeNonempty( Name ) }
                Coeff (:,1)  double  { mustBeNonempty( Coeff ) }
                Knots (:,1)  double  { mustBeNonempty( Knots ) }
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
            Ok = ( numel( Knots ) == B.k );
            assert( Ok, '"%s" spline must have %2.0f knots', Name, B.k);
            B.n = Knots;
            %--------------------------------------------------------------
            % Calculate the response variables 
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
            end            %--------------------------------------------------------------
            % Create B-spline array for distributed factors
            %--------------------------------------------------------------
            obj = obj.createBsplineTable( Opts.M, Opts.K );
            %--------------------------------------------------------------
            % Generate the design information linking columns of the design
            % matrix to factor values
            %--------------------------------------------------------------
            obj = obj.genDesignInfo();
        end % addFactor
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
                    Finish = Finish + 1;
                    D( Q,: ) = { Finish, nan };
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
    end % private methods

    methods ( Static = true, Access = protected )
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
    end % set/get methods

end % DoEgeneratorECOMO