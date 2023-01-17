classdef DoEgeneratorECOMO < handle
    % A class to handle mixed fixed and distributed parameter DoE
    % generation for the ECOMO model

    events
        DESIGN_AVAILABLE
    end % events

    properties ( Access = protected )
        Bspline     ( :, 1 )    bSplineTools 
    end % protected properties

    properties ( SetAccess = protected )
    end % SetAccess protected

    properties ( Access = private )
    end % private properties

    properties ( SetAccess = protected, Dependent = true )
    end % Accessible dependent properties

    methods ( Abstract = true )
        obj = generate( obj, varargin )
        obj = scrambleDesign( obj, varargin )
    end % abstract method signatures

    methods
    end % Constructor method

    methods
    end % ordinary methods



end % DoEgeneratorECOMO