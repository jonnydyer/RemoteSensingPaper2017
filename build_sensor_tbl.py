#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 18:28:15 2017

@author: jonny
"""
import pandas as pd
sensors = pd.read_csv('figures/sensors.csv')

tbl = '''
\\begin{table}[t]
\\centering
\\resizebox{0.96 \\linewidth}{!}{%
\\begin{tabular}{cccccccccc}
\\toprule
\\textbf{Model} & \\textbf{Type} & \\textbf{Shutter} & \\textbf{Width} 
& \\textbf{Height} & \\textbf{FPS} & \\textbf{Pixel Size (um)} & 
\\textbf{QE (525 nm)} & \\textbf{$N_{rd}$}  & \\textbf{$N_{e^-}^{FWC}$} \\\ \n
\\midrule \n                        
'''

for s in sensors.itertuples():
    tbl += '%s & %s & %s & %d & %d & %d & %.1f & %.0f & %.1f & %d' % (s.Model, 
                                                               s.Type,
                                                               s.Shutter,
                                                               s.Width,
                                                               s.Height,
                                                               s.FPS,
                                                               s._8,
                                                               s._9,
                                                               s._10,
                                                               s.FWC)
    tbl += '\\\ \n'

tbl += '''
\\end{tabular}%
}
\\caption{All sensors considered in analyses}
\\label{table:all_sensors}
\\end{table}
'''
with open('sensor_table.tex', 'w') as f:
    f.write(tbl)