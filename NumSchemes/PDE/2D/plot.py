import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def round_up(number: float, decimals:int=1) -> float:
    """ Rounding of a number

    :param number: A number to round
    :type number: float
    :param decimals: Number of decimals defaults to 1
    :type decimals: int, optional
    :return: The rounded number
    :rtype: float
    """

    if number != 0.0:
        power = int(np.log10(number))
        digit = number / 10**power * 10**decimals
        return np.ceil(digit) * 10**(power - decimals)
    else:
        return 1.0

def plot_ax_scalar(fig, ax, X, Y, field, title, cmap='RdBu', 
        field_ticks=None, max_value=None, cbar=True):
    """ Plot a 2D field on mesh X and Y with contourf and contour. Automatic
    handling of maximum values for the colorbar with an up rounding to a certain
    number of decimals. """
    # If not specified find the maximum value and round this up to 1 decimal
    if max_value is None:
        max_value = round_up(np.max(np.abs(field)), decimals=1)
    else:
        max_value = round_up(max_value, decimals=1)

    # Depending on the scale (log is typically for streamers) the treatment
    # is not the same
    if cmap == 'Blues':
        field_ticks = np.linspace(0, max_value, 5)
        levels = np.linspace(0, max_value, 101)
    else:
        field_ticks = np.linspace(-max_value, max_value, 5)
        levels = np.linspace(-max_value, max_value, 101)
    cs1 = ax.contourf(X, Y, field, levels, cmap=cmap)

    # Contours
    clevels = np.linspace(np.min(field), np.max(field), 6)[1:-1]
    cs2 = ax.contour(X, Y, field, levels=clevels, colors='k', linewidths=0.9)
    ax.clabel(cs2, fmt='%.2f')

    # Put colorbar if specified
    xmax, ymax = np.max(X), np.max(Y)
    if cbar:
        # Adjust the size of the colorbar
        xmax, ymax = np.max(X), np.max(Y)
        fraction_cbar = 0.1
        aspect = 0.85 * ymax / fraction_cbar / xmax
        # Set the colorbar in scientific notation
        sfmt = mpl.ticker.ScalarFormatter(useMathText=True) 
        sfmt.set_powerlimits((0, 0))
        fig.colorbar(cs1, ax=ax, pad=0.05, fraction=fraction_cbar, aspect=aspect, 
            ticks=field_ticks, format=sfmt)

    # Apply same formatting to x and y axis with scientific notation
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    ax.set_aspect("equal")
    ax.set_title(title)
