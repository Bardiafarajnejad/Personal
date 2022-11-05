#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import datetime
import matplotlib.pyplot as plt
import statsmodels.formula.api as sm
from matplotlib.ticker import FormatStrFormatter

def plot_func_p1(p1_i,host,fig):
    par1 = host.twinx()
    par2 = host.twinx()

    host.set_xlabel("Date (10/27/2021) and Time (in US/Eastern 9:30-16:00)", size=17)
    host.set_ylabel("Avg Spread (Ask-Bid)", size=17)
    par1.set_ylabel("Avg Total Bin Volume (% of Full Day's Volume)", size=15)
    par2.set_ylabel("Volume-weighted average price (VWAP)", size=17)

    p1, = host.plot(list(p1_i.index), list(p1_i[' avg_spread'].values),    color=plt.cm.viridis(0), label="Avg Spread (Ask-Bid)")
    p2, = par1.plot(list(p1_i.index), list(p1_i[' avg_total_bin_volume'].values),    color=plt.cm.viridis(0.3), label="Avg Total Bin Volume (% of Full Day's Volume)")
    p3, = par2.plot(list(p1_i.index), list(p1_i[' vwap'].values), color=plt.cm.viridis(0.9), label="Volume-weighted average price (VWAP)")

    lns = [p1, p2, p3]
    host.legend(handles=lns, loc='upper center', fontsize=16)

    # right, left, top, bottom
    par2.spines['right'].set_position(('outward', 60))

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())
    host.tick_params(axis='y', colors=p1.get_color())
    par1.tick_params(axis='y', colors=p2.get_color())
    par2.tick_params(axis='y', colors=p3.get_color())

    host.tick_params(axis='both', labelsize=15)
    par1.tick_params(axis='both', labelsize=15)
    par2.tick_params(axis='both', labelsize=15)
    host.tick_params(axis='x', labelsize=14)

    par1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    fig.tight_layout()
    plt.title(p1_i.index.name, size=20)
    plt.grid()


def ols_coef(x,formula):
    return sm.ols(formula,data=x).fit().params

# read stock and convert 'timestamp' from mili-seconds to 15 seconds
def get_vol(ticker, lookback_freq):
    stock_i = pd.read_csv(ticker+'.csv')
    stock_i['timestamp'] = round(stock_i['timestamp']/1000)

    # sort timestamp, create a date index (localize in UTC, convert to ET), drop 'timestamp', & keep necessary cols
    stock_i = stock_i.sort_values(by='timestamp')
    stock_i.index = [datetime.datetime.fromtimestamp(i, datetime.timezone.utc) for i in stock_i['timestamp']]
    stock_i.index = stock_i.index.tz_convert('US/Eastern')
    stock_i = stock_i.drop('timestamp', axis=1)
    stock_i = stock_i[[' bid_price',' ask_price',' nbb_agg_size',' nbo_agg_size',\
                       ' trade_price',' trade_size',' volume',' vwap']].ffill()

    # throw away observations before 9:30 and after 16:00 eastern standard time
    stock_i = stock_i[stock_i.index>='2021-10-27 09:30:00-04:00']
    stock_i = stock_i[stock_i.index<='2021-10-27 16:00:00-04:00']

    # midpoint price:
    stock_i[' mid price'] = (stock_i[' bid_price'] + stock_i[' ask_price'])/2
    #########################################
    #########################################
    # returns: create 15 sec interval and get last
    p2 = stock_i[[' mid price']].resample('15s').last().ffill()

    p2 = p2.pct_change(1).dropna()

    #########################################
    #########################################
    vol_15sec = np.sqrt(p2.resample(lookback_freq).var())
    mult_ann = np.sqrt(4*60*6.5*252) # 4 15secs in 1min, 60min per hour, 6.5hrs in 1 trading day, 252 days in a trading yr
    vol_ann = mult_ann*vol_15sec
    vol_ann = vol_ann.rename(columns={' mid price':'Ann Vol'})
    #vol_ann.plot(title=ticker+' Ann Vol')
    #stock_i[[' mid price']].plot(title=ticker + ' Mid Price')
    
    return(vol_ann, stock_i[[' mid price']])
    
    
def plot_func_p3(p1_i,p1_c,host,fig):
    par1 = host.twinx()

    host.set_xlabel("Date (10/27/2021) and Time (in US/Eastern 9:30-16:00)", size=18)
    host.set_ylabel("Annualized Volatility", size=18)
    par1.set_ylabel("Mid Price", size=18)

    p1, = host.plot(list(p1_i.index), list(p1_i['Ann Vol'].values),    color=plt.cm.viridis(0), label="Ann Volatility")
    p2, = par1.plot(list(p1_c.index), list(p1_c[' mid price'].values),    color=plt.cm.viridis(0.3), label="Mid Price")

    lns = [p1, p2]
    host.legend(handles=lns, loc='upper right', fontsize=16)
    
    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    host.tick_params(axis='y', colors=p1.get_color())
    par1.tick_params(axis='y', colors=p2.get_color())

    host.tick_params(axis='both', labelsize=15)
    par1.tick_params(axis='both', labelsize=15)

    fig.tight_layout()
    plt.title(p1_i.index.name, size=20)

