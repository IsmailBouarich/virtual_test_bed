hex_size = '${fparse 1.775 * 2.54 *0.01}'

A_hexa = '${fparse 6 * hex_size * hex_size /sqrt(3)}'
A_fuel = '${fparse pi * 0.007777 * 0.007777}'
heated_length = 1.982
n_heater = 42
radius_fuel = ${units 12.7 mm -> m}
section_fuel_channel = '${fparse pi * radius_fuel * radius_fuel}'
volume_fuel_channel = '${fparse section_fuel_channel * heated_length}'
R_l = 0.009
# power =3.802e+04 # W

[Functions]
    [rho_G348]
        type = ParsedFunction
        value = '-0.033488 * t + 1899.3'
    []
    [cp_G348]
        type = PiecewiseLinear
        x = '280 374.15 472.45 574.75 674.75 774.75 874.75 974.85 1074.45 1173.95 1274.05'
        y = '749.40 961.86 1188.22 1378.01 1525.18 1641.99 1734.99 1809.46 1869.08 1917.20 1956.54'
        # scaling = 0.01
    []
    [k_G348]
        type = PiecewiseLinear
        x = '280 374.15 472.45 574.75 674.75 774.75 874.75 974.85 1074.45 1173.95 1274.05'
        y = '131.510 126.985 116.145 105.030 95.635 87.830 81.295 75.770 70.850 67.280 63.685'
    []
    [cp_greencast_fn]
        type = PiecewiseLinear
        x = '335.45 422.05 506.35 589.25 671.25 752.35 832.75 912.45 991.35 1069.75 1147.45 1224.55 1301.25 1377.65 1453.85 1529.95 1606.05 1682.25 1758.15 1834.25'
        y = '940 1100 1250 1260 1180 1180 1190 1200 1210 1250 1200 1200 1180 1190 1220 1290 1290 1290 3050 1270'
        # scale_factor = 0.01
    []
    [k_greencast_half_fn]
        type = PiecewiseLinear
        x = '478.15 698.15 923.15 1143.15 1368.15'
        y = '2.625	1.79	1.415	1.245	1.235'
    []
    [rho_greencast]
        type = PiecewiseLinear
        x = '0 1800'
        y = '2912 2912'
    []
    [k_greencast_fn]
        type = PiecewiseLinear
        x = '478.15 698.15 923.15 1143.15 1368.15'
        y = '5.25 3.58 2.83 2.49 2.47'
    []
    # [power_fn]
    #     type = PiecewiseLinear
    #     data_file = 'power.txt'
    #     x_index_in_file = 0
    #     y_index_in_file = 1
    #     format = columns
    # []
    [power_fn]
        type = ConstantFunction
        value = 2.2e6
    []
    [power_density_sub]
        type = ParsedFunction
        expression = 'power/(210 * A_fuel * heated_length)'
        symbol_names = 'power A_fuel heated_length'
        symbol_values = 'power_fn ${A_fuel} 1.982'
    []
[]    