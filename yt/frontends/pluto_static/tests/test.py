#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 19:49:22 2023

@author: adutt
"""

import yt

ds = yt.load("/u/adutt/jan23/pluto-clcrush/mhdRT/output/data.%04d.dbl.h5"%5)
s = yt.SlicePlot(ds, "z", ("gas", "density"))