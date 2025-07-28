classdef ScreenConfiguration
    properties (SetAccess=private)
        screen_size_x_mm
        screen_size_y_mm
        screen_res_x_pix
        screen_res_y_pix
        viewing_distance_mm
    end

    methods
        function obj = ScreenConfiguration(screen_size_x_mm, screen_size_y_mm, screen_res_x_pix, screen_res_y_pix, viewing_distance_mm)
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
            x_mm = x/obj.screen_res_x_pix*obj.screen_size_x_mm;
            y_mm = y/obj.screen_res_y_pix*obj.screen_size_y_mm;
        end

        function [x_deg, y_deg] = pix_to_deg(obj, x, y)
            % N.B.: output is in Fick angles
            [x_mm , y_mm ] = obj.pix_to_mm(x, y);
            [x_deg, y_deg] = obj.mm_to_deg(x_mm, y_mm);
        end

        function [x_deg, y_deg] = mm_to_deg(obj, x, y)
            % N.B.: output is in Fick angles
            x_deg = atan2(x,obj.viewing_distance_mm)*180/pi;            % azimuth
            y_deg = atan2(y,hypot(obj.viewing_distance_mm,x))*180/pi;   % elevation
        end

        function [x_pix, y_pix] = mm_to_pix(obj, x, y)
            x_pix = x/obj.screen_size_x_mm*obj.screen_res_x_pix;
            y_pix = y/obj.screen_size_y_mm*obj.screen_res_y_pix;
        end

        function [x_pix, y_pix] = deg_to_pix(obj, x, y)
            % N.B.: input is in Fick angles
            [x_mm , y_mm ] = obj.deg_to_mm(x, y);
            [x_pix, y_pix] = obj.mm_to_pix(x_mm, y_mm);
        end

        function [x_mm, y_mm] = deg_to_mm(obj, x, y)
            % N.B.: input is in Fick angles
            [x,y,z] = Fick_to_cartesian(x, y);
            x_mm = x./z.*obj.viewing_distance_mm;
            y_mm = y./z.*obj.viewing_distance_mm;
        end

        function [x_deg, y_deg] = screen_extents(obj)
            x_deg = obj.mm_to_deg(obj.screen_size_x_mm/2,0)*2;
            y_deg = obj.mm_to_deg(obj.screen_size_y_mm/2,0)*2;
        end
    end
end
