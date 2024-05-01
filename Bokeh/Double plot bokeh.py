from bokeh.layouts import row
from bokeh.models import BoxAnnotation, CustomJS
from bokeh.plotting import figure, show
import numpy as np
from bokeh.models import WheelPanTool, WheelZoomTool
from bokeh.models.tickers import FixedTicker

Peak_positions = [129.995153, 128.98897799999997, 127.98280299999998, 126.97662799999998, 125.97045299999998,
                  124.96427799999998, 131.993108, 130.986933, 130.98693300000002, 129.980758, 128.974583, 127.968408,
                  126.962233, 128.991798, 127.98562299999998, 126.97944799999998, 125.97327299999998,
                  124.96709799999998, 123.96092299999998, 130.989753, 129.983578, 128.977403, 127.971228, 126.965053,
                  125.958878, 127.98844299999998, 126.98226799999998, 125.97609299999998, 124.96991799999998,
                  123.96374299999998, 122.95756799999998, 129.98639799999998, 128.980223, 127.97404799999998,
                  126.96787299999998, 125.96169799999998, 124.95552299999999, 128.998119, 127.99194399999999,
                  126.98576899999999, 125.97959399999999, 124.97341899999999, 123.967244, 130.996074,
                  129.98989899999998, 128.983724, 127.97754899999998, 126.97137399999998, 125.96519899999998,
                  127.99476399999998, 126.98858899999998, 125.98241399999998, 124.97623899999998, 123.97006399999998,
                  122.96388899999998, 129.992719, 128.98654399999998, 127.98036899999998, 126.97419399999998,
                  125.96801899999998, 124.96184399999999, 126.99140899999998, 125.98523399999998, 124.97905899999998,
                  123.97288399999998, 122.96670899999998, 121.96053399999998, 128.989364, 127.98318899999998,
                  126.97701399999998, 125.97083899999998, 124.96466399999998, 123.95848899999999, 123.95248899999999]
Peak_intensities = [1.7233190184375e-26, 5.743535068615782e-22, 7.656898051140651e-18, 5.103833010955323e-14,
                    1.701022478667893e-10, 2.2676897670614572e-07, 1.6764028565625e-26, 1.1174342640893439e-22,
                    4.469737056357375e-22, 7.448444326331537e-18, 4.964884706455059e-14, 1.6547133245830305e-10,
                    2.2059534901124576e-07, 3.09884092588125e-24, 1.0327920332474559e-19, 1.3768494859232918e-15,
                    9.177619723336025e-12, 3.058747693459174e-08, 4.077718508406875e-05, 3.0144771366187495e-24,
                    1.0046749883494187e-19, 1.3393657161348891e-15, 8.927765408516457e-12, 2.9754754145683945e-08,
                    3.966705457674946e-05, 1.3930698525893435e-22, 4.642869640371518e-18, 6.189564279900617e-14,
                    4.125757230172421e-10, 1.3750461221959655e-06, 0.0018331198203701813, 1.3551444945981563e-22,
                    4.5164707430798886e-18, 6.02105769662457e-14, 4.0134363586467184e-10, 1.3376114477309736e-06,
                    0.0017832144080184187, 4.640385778565624e-24, 1.5465632402329466e-19, 2.0617750076625486e-15,
                    1.3743104942742665e-11, 4.580347825667084e-08, 6.106214364657646e-05, 4.514054502684374e-24,
                    1.5044591315029905e-19, 2.0056446168443536e-15, 1.336895846767885e-11, 4.455651041302899e-08,
                    5.9399769248622725e-05, 8.344257336366185e-22, 2.781001899255245e-17, 3.707446331960474e-13,
                    2.4712601433404533e-09, 8.236298181063176e-06, 0.010980083648448024, 8.117090733008809e-22,
                    2.7052910564662867e-17, 3.6065136837437575e-13, 2.4039818044607964e-09, 8.01207069063376e-06,
                    0.010681158506706888, 3.7511229571209806e-20, 1.2501867628924708e-15, 1.6666656465040492e-11,
                    1.1109437644380494e-07, 0.0003702590409577946, 0.493604669468868, 3.6490012431571427e-20,
                    1.2161512976568899e-15, 1.6212918332829887e-11, 1.0806990930053309e-07, 0.00036017899604712673,
                    0.48016662559695955, 0.0017276163431410381]

# precision parameter
x = 0.004

mass_range = np.linspace(min(Peak_positions) - 1, max(Peak_positions) + 1, 1000)

intensity = np.zeros_like(mass_range)

for peak_position, peak_intensity in zip(Peak_positions, Peak_intensities):
    peak_shape = peak_intensity * np.exp(-((mass_range - peak_position) ** 2) / (2 * x ** 2))  # Gaussian example

    intensity += peak_shape

ticked_peaks = []
for i in range(len(Peak_positions)):
    if Peak_intensities[i] > 0.0001:
        ticked_peaks.append(Peak_positions[i])

print(ticked_peaks)

# Create a new plot with a title and axis labels
p1 = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [Th]', y_axis_label='Intensity')
p1 = figure(width=700, title=f'Mass spectrum of molecule')
p1.height = 500
p1.xaxis.ticker = FixedTicker(ticks=ticked_peaks)
p1.add_tools(WheelPanTool(dimension="height"))
p1.add_tools(WheelZoomTool(dimensions="height"))
p1.line(mass_range, intensity, legend_label="Intensity", line_width=1)
p1.xaxis.major_label_orientation = "vertical"
def overview_plot(Mass_range, intensity, p1):
    p2 = figure(title="Simulated Mass Spectrum", x_axis_label='Mass [Th]', y_axis_label='Intensity')
    p2 = figure(width=300, title=f'Mass spectrum of molecule')
    p2 = figure(toolbar_location=None)
    p2.height = 300
    p2.line(mass_range, intensity, legend_label="Intensity", line_width=1)

    box = BoxAnnotation(left=0, right=0, bottom=0, top=0,
    fill_alpha=0.1, line_color='black', fill_color='black')

    jscode = """
        box[%r] = cb_obj.start
        box[%r] = cb_obj.end
    """

    xcb = CustomJS(args=dict(box=box), code=jscode % ('left', 'right'))
    ycb = CustomJS(args=dict(box=box), code=jscode % ('bottom', 'top'))

    p1.x_range.js_on_change('start', xcb)
    p1.x_range.js_on_change('end', xcb)
    p1.y_range.js_on_change('start', ycb)
    p1.y_range.js_on_change('end', ycb)

    p2.add_layout(box)

    return p2

layout = row(p1, overview_plot(mass_range, intensity, p1))
show(layout)

