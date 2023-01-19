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
        Bspline     (:,1)    bSplineTools                                   % Array of bSplineTools objects (one for each distributed parameter)
    end % protected properties

    properties ( SetAccess = protected )
        Factors     (:,:)    table                                          % Factor details and type
        Scramble    (1,1)    logical = false                                % Set to true to apply scramble to design
        Design      (:,:)    double                                         % Design array
    end % SetAccess protected

    properties ( Access = private )
        NumPoints_  (1,1)   double                                          % Number of points in the design
    end % private properties

    properties ( SetAccess = protected, Dependent = true )
        NumPoints           int64                                           % Number of points in the design
        NumFactors          int64                                           % Number of factors
        NumFixed            int64                                           % Number of fixed factors
        NumDist             int64                                           % Number of B-spline factors
    end % Accessible dependent properties

    methods ( Abstract = true )
        obj = generate( obj, varargin )
        obj = scrambleDesign( obj, varargin )
    end % abstract method signatures

    methods
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
            % obj.addFactor( S );
            %
            % Input Arguments:
            %
            % S     --> (struct) Multidimensional Structure defining factor 
            %                    properties with fields:
            %
            %            Name  - (string) Name of factor
            %            Units - (string) Factor units
            %            Fixed - (logical) True if fixed factor. False
            %                    if distributed factor.
            %            Lo    - (double) Low natural limit for factor
            %            Hi    - (double) High natural limit for factor
            %
            % Note each dimension of S must define a different factor
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   
                S       (1,:)   struct   { mustBeNonempty( S ) }
                Opts.M  (1,1)   int8 = 4
                Opts.K  (1,1)   int8 = 2
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
            % Create B-spline array for distributed factors
            %--------------------------------------------------------------
            obj = obj.createBsplineArray( Opts.M, Opts.K );
        end % addFactor
    end % ordinary methods

    methods ( Access = protected )
        function obj = createBsplineArray( obj, M, K )
            %--------------------------------------------------------------
            % Create a B-spline representation for each distributed
            % parameter
            %
            % obj = obj.createBsplineArray( M, K );
            %
            % Input Arguments:
            %
            % M     --> (int8) Spline order (1 <= M <= 4)
            % K     --> (int8) Number of knots (1 <= K <= 7)
            %--------------------------------------------------------------
            arguments
                obj (1,1)                   { mustBeNonempty( obj ) }
                M   (1,1) int8  { mustBeGreaterThan(M,0), ...
                                  mustBeLessThan(M,5)} = 4;
                K   (1,1) int8  { mustBeGreaterThan(K,0), ...
                                  mustBeLessThan(K,8)} = 2;
            end
            D = ~obj.Factors{ :, "Fixed" };
            N = sum( D );
            if ( N > 0 )
                obj.Bspline( N ) = bSplineTools.empty;
                Lo = obj.Factors{ :, "Lo" }( D );
                Hi = obj.Factors{ :, "Hi" }( D );
                for Q = 1:N
                    DK = linspace( Lo( Q ), Hi( Q ), K + 2 ).';
                    DK = DK( 2:end-1 );
                    obj.Bspline( Q ) = bSplineTools( ( M - 1 ), DK,...
                        Lo( Q ), Hi( Q ) );
                end
            end
        end % createBsplineArray

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