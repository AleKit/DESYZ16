import clastromoon

mod = input('graph, text or both? in quotation marks: ')
mymoon = clastromoon.MoonInterpolate()
mymoon.config_mod(mod)
mymoon.config_moon_plan('berlin', '2016-10-21 00:00:00', 0, 3, 60, 6000)
mymoon.get_interpolation()
mymoon.config_moon_graph_1()

