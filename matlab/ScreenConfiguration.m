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
            x_deg = np.atan2(x,obj.viewing_distance_mm)*180/pi;             % azimuth
            y_deg = np.atan2(y,np.hypot(obj.viewing_distance_mm,x))*180/pi; % elevation
        end
    end
end
