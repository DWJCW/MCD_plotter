import matplotlib.pyplot as plt
import seaborn as sns


def set_plot_config():
    sns.set_context("paper", font_scale=0.75)
    sns.set_palette(sns.color_palette("bright", 8))
    sns.set_style('white')

    SMALL_SIZE = 6
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 10

    rc = {
        'xtick.major.size': 2,
        'xtick.major.width': 0.5,
        'ytick.major.size': 2,
        'ytick.major.width': 0.5,
        'xtick.bottom': True,
        'ytick.left': True,
        'font.size': MEDIUM_SIZE,
        'axes.titlesize': MEDIUM_SIZE,
        'axes.labelsize': MEDIUM_SIZE,
        'xtick.labelsize': SMALL_SIZE,
        'ytick.labelsize': SMALL_SIZE,
        'legend.fontsize': SMALL_SIZE,
        'figure.titlesize': BIGGER_SIZE,
        'savefig.dpi': 300,
        'figure.dpi': 300,
        "font.family": "serif",
        "font.serif": ["Liberation Serif", "DejaVu Serif", "Nimbus Roman No9 L", "Times"]
    }

    plt.rcParams.update(rc)
