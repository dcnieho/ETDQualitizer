classdef ScreenConfiguration
    %SCREENCONFIGURATION Screen and Viewing Geometry Configuration
    %
    %   Provides methods for converting between pixel, millimeter, and degree units.
    %
    %   Example:
    %       sc = ScreenConfiguration(500, 300, 1920, 1080, 600);
    %       [azi, ele] = sc.pix_to_deg(960, 540)
    properties (SetAccess=private)
        screen_size_x_mm       % Screen width in millimeters
        screen_size_y_mm       % Screen height in millimeters
        screen_res_x_pix       % Horizontal screen resolution in pixels
        screen_res_y_pix       % Vertical screen resolution in pixels
        viewing_distance_mm    % Viewing distance in millimeters
    end

    methods
        function obj = ScreenConfiguration(screen_size_x_mm, screen_size_y_mm, screen_res_x_pix, screen_res_y_pix, viewing_distance_mm)
            %SCREENCONFIGURATION Construct a ScreenConfiguration object
            %
            %   Inputs:
            %       screen_size_x_mm     - Screen width in mm
            %       screen_size_y_mm     - Screen height in mm
            %       screen_res_x_pix     - Horizontal resolution in pixels
            %       screen_res_y_pix     - Vertical resolution in pixels
            %       viewing_distance_mm  - Viewing distance in mm
            %
            %   Output:
            %       obj - ScreenConfiguration object
            arguments
                screen_size_x_mm     (1,1) {mustBeNumeric}
                screen_size_y_mm     (1,1) {mustBeNumeric}
                screen_res_x_pix     (1,1) {mustBeNumeric}
                screen_res_y_pix     (1,1) {mustBeNumeric}
                viewing_distance_mm  (1,1) {mustBeNumeric}
            end
            obj.screen_size_x_mm = screen_size_x_mm;
            obj.screen_size_y_mm = screen_size_y_mm;
            obj.screen_res_x_pix = screen_res_x_pix;
            obj.screen_res_y_pix = screen_res_y_pix;
            obj.viewing_distance_mm = viewing_distance_mm;
        end

        function [x_mm, y_mm] = pix_to_mm(obj, x, y)
            %PIX_TO_MM Convert pixel coordinates to millimeters
            %
            %   Inputs:
            %       x - Horizontal pixel coordinate
            %       y - Vertical pixel coordinate
            %
            %   Outputs:
            %       x_mm - Horizontal position in mm
            %       y_mm - Vertical position in mm
            x_mm = x/obj.screen_res_x_pix*obj.screen_size_x_mm;
            y_mm = y/obj.screen_res_y_pix*obj.screen_size_y_mm;
        end

        function [azi, ele] = pix_to_deg(obj, x, y)
            %PIX_TO_DEG Convert pixel coordinates to angular gaze direction in degrees
            %
            %   Inputs:
            %       x - Horizontal pixel coordinate
            %       y - Vertical pixel coordinate
            %
            %   Outputs:
            %       azi - Azimuth in degrees (Fick angle)
            %       ele - Elevation in degrees (Fick angle)
            [x_mm , y_mm ] = obj.pix_to_mm(x, y);
            [azi, ele] = obj.mm_to_deg(x_mm, y_mm);
        end

        function [azi, ele] = mm_to_deg(obj, x, y)
            %MM_TO_DEG Convert millimeter coordinates to angular gaze direction in degrees
            %
            %   Inputs:
            %       x - Horizontal position in mm
            %       y - Vertical position in mm
            %
            %   Outputs:
            %       azi - Azimuth in degrees (Fick angle)
            %       ele - Elevation in degrees (Fick angle)
            azi = atan2(x,obj.viewing_distance_mm)*180/pi;          % azimuth
            ele = atan2(y,hypot(obj.viewing_distance_mm,x))*180/pi; % elevation
        end

        function [x_pix, y_pix] = mm_to_pix(obj, x, y)
            %MM_TO_PIX Convert millimeter coordinates to pixel coordinates
            %
            %   Inputs:
            %       x - Horizontal position in mm
            %       y - Vertical position in mm
            %
            %   Outputs:
            %       x_pix - Horizontal pixel coordinate
            %       y_pix - Vertical pixel coordinate
            x_pix = x/obj.screen_size_x_mm*obj.screen_res_x_pix;
            y_pix = y/obj.screen_size_y_mm*obj.screen_res_y_pix;
        end

        function [x_pix, y_pix] = deg_to_pix(obj, azi, ele)
            %DEG_TO_MM Convert angular gaze direction in degrees to millimeter coordinates
            %
            %   Inputs:
            %       azi - Azimuth in degrees (Fick angle)
            %       ele - Elevation in degrees (Fick angle)
            %
            %   Outputs:
            %       x_mm - Horizontal position in mm
            %       y_mm - Vertical position in mm
            [x_mm , y_mm ] = obj.deg_to_mm(azi, ele);
            [x_pix, y_pix] = obj.mm_to_pix(x_mm, y_mm);
        end

        function [x_mm, y_mm] = deg_to_mm(obj, azi, ele)
            %DEG_TO_MM Convert angular gaze direction in degrees to millimeter coordinates
            %
            %   Inputs:
            %       azi - Azimuth in degrees (Fick angle)
            %       ele - Elevation in degrees (Fick angle)
            %
            %   Outputs:
            %       x_mm - Horizontal position in mm
            %       y_mm - Vertical position in mm
            x_mm = obj.viewing_distance_mm.*tand(azi);
            y_mm = obj.viewing_distance_mm.*tand(ele)./cosd(azi);
        end

        function [x_deg, y_deg] = screen_extents(obj)
            %SCREEN_EXTENTS Compute screen extents in degrees
            %
            %   Outputs:
            %       x_deg - Horizontal extent in degrees
            %       y_deg - Vertical extent in degrees
            x_deg = obj.mm_to_deg(obj.screen_size_x_mm/2,0)*2;
            y_deg = obj.mm_to_deg(obj.screen_size_y_mm/2,0)*2;
        end
    end
end
