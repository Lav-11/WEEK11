#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import sys
import argparse

# Parameters
defLW = 1.2  # Default line width (balanced for readability)
defMS = 5  # Marker size (reduced for more subtle symbols)

# Support for up to 12 models
dashes = ['-', '--', '-.', ':', (0, (5, 1)), (0, (3, 5, 1, 5)), (0, (1, 1)), '-', '--', (0, (4, 2, 1, 2)), (0, (2, 2)), (0, (1, 2))]
markers = ['o', 'x', 's', '^', '*', '+', 'd', 'p', 'h', 'v', '<', '>']  # Extended marker list
colors = ['r', 'g', 'b', 'orange', 'purple', 'cyan', 'brown', 'pink', 'gray', 'lime', 'teal', 'gold']  # Extended color palette


class CmdLineParser(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(usage='usage: python3 perfprof.py [options] csvfile.csv outputfile.pdf')
        # Default options
        self.parser.add_argument("-D", "--delimiter", dest="delimiter", default=None, help="delimiter for input files")
        self.parser.add_argument("-M", "--maxratio", dest="maxratio", default=4, type=float, help="max ratio for the performance profile")
        self.parser.add_argument("-S", "--shift", dest="shift", default=0, type=float, help="shift for data")
        self.parser.add_argument("-L", "--logplot", dest="logplot", action="store_true", default=False, help="logarithmic scale for x-axis")
        self.parser.add_argument("-T", "--timelimit", dest="timelimit", default=1e99, type=float, help="time limit for runs")
        self.parser.add_argument("-P", "--plot-title", dest="plottitle", default=None, help="title for the plot")
        self.parser.add_argument("-X", "--x-label", dest="xlabel", default='Cost Ratio', help="label for x-axis")
        self.parser.add_argument("-B", "--bw", dest="bw", action="store_true", default=False, help="create black-and-white plots")
        
        # File input and output should be positional arguments
        self.parser.add_argument("input", help="CSV file with performance data")
        self.parser.add_argument("output", help="Output PDF file")

    def parseArgs(self):
        return self.parser.parse_args()


def readTable(fp, delimiter):
    """
    Read a CSV file with the performance profile data
    The format is as follows:
    ncols algo1 algo2 ...
    instance_name time(algo1) time(algo2) ...
    ...
    """
    firstline = fp.readline().strip().split(delimiter)
    ncols = int(firstline[0])
    assert ncols <= len(markers), f"Too many columns for the defined markers ({len(markers)} maximum allowed)"
    cnames = firstline[1:]  # Algorithm names
    rnames = []  # Instance names
    rows = []  # Data for all rows
    for row in fp:
        row = row.strip().split(delimiter)
        rnames.append(row[0])  # Add the name of the instance
        rdata = np.empty(ncols)  # Initialize the row data
        for j in range(ncols):
            rdata[j] = float(row[j + 1])  # Convert data to float
        rows.append(rdata)
    data = np.array(rows)  # Convert the list of rows into a numpy array
    return (rnames, cnames, data)


def main():
    parser = CmdLineParser()
    opt = parser.parseArgs()
    print(opt)
    
    # Read data from the CSV file
    with open(opt.input, 'r') as fp:
        rnames, cnames, data = readTable(fp, opt.delimiter)

    nrows, ncols = data.shape
    # Add shift to the data
    data = data + opt.shift
    # Compute the ratio of each value compared to the minimum
    minima = data.min(axis=1)
    ratio = data
    for j in range(ncols):
        ratio[:, j] = data[:, j] / minima
    # Compute the maximum ratio for scaling
    if opt.maxratio == -1:
        opt.maxratio = ratio.max()
    # Any time >= timelimit will count as maxratio + large value (so it does not show up in plots)
    for i in range(nrows):
        for j in range(ncols):
            if data[i, j] >= opt.timelimit:
                ratio[i, j] = opt.maxratio + 1e6
    # Sort the ratios to prepare for plotting
    ratio.sort(axis=0)
    
    # Plot the performance profile
    y = np.arange(nrows, dtype=np.float64) / nrows  # Proportion of instances solved
    for j in range(ncols):
        options = dict(label=cnames[j],  # Algorithm name for legend
                       linewidth=defLW, linestyle=dashes[j % len(dashes)],  # Line style
                       marker=markers[j % len(markers)], markeredgewidth=defLW, markersize=defMS)  # Marker style
        if opt.bw:
            options['markerfacecolor'] = 'w'
            options['markeredgecolor'] = 'k'
            options['color'] = 'k'  # Use black for black-and-white plots
        else:
            options['color'] = colors[j % len(colors)]  # Use predefined colors for each algorithm
        
        # Plot using a linear scale for the x-axis
        plt.plot(ratio[:, j], y, **options)
    
    # Adjust the x-axis to zoom into the range of interest
    plt.axis([1, 1.05, 0, 1])
    
    # Add a legend to the plot
    plt.legend(loc='lower right', fontsize=10)  # Improve the size of the legend
    
    # Add title and axis labels
    if opt.plottitle is not None:
        plt.title(opt.plottitle, fontsize=14)
    plt.xlabel(opt.xlabel, fontsize=12)
    plt.ylabel('Proportion', fontsize=12)
    plt.savefig(opt.output)  # Save the plot as a PDF


if __name__ == '__main__':
    main()