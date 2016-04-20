"""
IPython Imports
Nick Eubank
"""
import pandas as pd
import numpy as np

import re
import os
import sys
import time
from pandas_profiling import ProfileReport

from see import see

pd.options.display.max_rows = 200
pd.options.display.max_columns = 90
pd.set_option('display.width', 150)
pd.set_option('mode.chained_assignment','raise')
pd.set_option('io.hdf.default_format', 'table')

from pandas.util.testing import assert_frame_equal, assert_series_equal
import pandas.util.testing as tm
from pandas import DataFrame, Series, Panel, Index, date_range


