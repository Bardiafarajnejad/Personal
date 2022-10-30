#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
import statsmodels.api as sm
path="/Users/bardiafarajnejad/Desktop/1st Half MFE/AFP/Final Delivery/MFE Group 16 AFP Code"
os.chdir(path)
print(os.getcwd())  #check if correct wd


# In[2]:


def get_data(flag=False): #run this function once to get the data!
    if(flag):
        data = pd.read_csv('afp_data_sample_1995_onward_v2.csv')
        data

        data.columns.values #view column names


        ## First, generate daily returns for universe

        ### Note:
        ### pricepctchgd is the raw ret on adj close
        ### so, "ret" OR mkt cap growth = pricepctchgd

        #### Should use unadjustedpriceclose for pricepctchgd, which is what is going on above ^ (back when it was checked)

        ## Next, generate 1 week, 2 week, 3 week, and 4 week returns for universe
        ## (RECALL: MAX HOLDING PERIOD O'NEIL TEAM WANTS IS 10 DAYS, ie 1-2 weeks)

        data['date'] = pd.to_datetime(data['tradedate'], format='%Y%m%d', errors='ignore')
        data['week_count'] = (data['date']-datetime.datetime(1900,1,1)).dt.days // 7 +1 #acutally, we wont use this since we want rolling weekly returns, but it is useful to have just in case
        data['2week_count'] = (data['date']-datetime.datetime(1900,1,1)).dt.days // 14 +1 #acutally, we wont use this since we want rolling weekly returns, but it is useful to have just in case
        data['3week_count'] = (data['date']-datetime.datetime(1900,1,1)).dt.days // 21 +1 #acutally, we wont use this since we want rolling weekly returns, but it is useful to have just in case
        data['4week_count'] = (data['date']-datetime.datetime(1900,1,1)).dt.days // 28 +1 #acutally, we wont use this since we want rolling weekly returns, but it is useful to have just in case

        data = data.sort_values(['osid','tradedate']).reset_index(drop=True).copy()

        data


        #Next, compute rolling/moving weekly returns:
        data['pricepctchgd'] = data['pricepctchgd']/100 #convert percentage returns to decimal returns


        #rolling weekly returns (1 week = 5 trading days backward)
        #NOTE: we must avoid 5 days turning into many months because of missing data in between day 1 and day 5
        #for example, IF (day 1 = November 1st, 1999) and (day 2 = December 1st, 1999), 
        #THEN we should not bridge the gap and compute rolling return...
        #so we will restrict the rolling return only to the times when week.diff(1) = 1 OR week.diff(1) = 0
        #and, of course, we will group by OSID so to not compute rolling ret over 2 different stocks...

        data['cumret']=0
        data['cumalpha']=0
        for i in range(0,20): #0-19: 20 trading days ~ 4 weeks
            print(i,end=",") #to keep track of slow loop
            #FIRST, CHECK IF LAG WEEK == WEEK 
            data['r_shifted'] = data[['osid','pricepctchgd']].groupby(['osid'])['pricepctchgd'].shift(i) #group by OSID
            data['alpha_shifted'] = data[['osid','alpha']].groupby(['osid'])['alpha'].shift(i) #group by OSID

            data['cumret'] = (1+data['cumret'])*(1+data['r_shifted']) -1
            data['cumalpha'] = (1+data['cumalpha'])*(1+data['alpha_shifted']) -1
            #Note: taking avg is like weekly regression except for Beta differences throughout the week

            if i in [1-1,5-1,10-1,15-1,20-1]:
                data['ret'+str(i+1)+'d'] = data['cumret']
                data['alpha'+str(i+1)+'d'] = data['cumalpha']

        data.drop(['cumret','r_shifted','cumalpha','alpha_shifted'], axis=1, inplace=True)
        data


        #### NOW, restrict the rolling return only to the times when week.diff(1) = 0 OR week.diff(1) = 1
        #sort the data first
        data = data.sort_values(['osid','tradedate']).reset_index(drop=True).copy()
        #now that the data is sorted, we can just check at which times week.diff(1) is NOT <=1
        xx = list((data['date'].diff(1)).dt.days >20)
        yy = list((data['osid'].diff(1)) ==0)
        zz = [(a and b) for a, b in zip(xx, yy)]
        zz

        print('\nThere are ', sum(zz), ' problem rows') #these 110 dates are problematic

        problem_rows = list(np.where(zz)[0]) #gets the row index of the rows for which zz is True
        #problem_rows

        problem_rows2 = [i-1 for i in problem_rows]

        problem_rows3 = np.sort(problem_rows + problem_rows2)
        problem_rows3 #this is the indices we want to check to confirm that there is a large gap in the date (for example, instead of being the next day, it is instead the next year)


        ###########################################################################################
        ###########################################################################################
        ###########################################################################################
        #this one will take a long time to run


        data['index'] = list(data.index) #create a col called index
        data[data.index.isin(problem_rows3)] 
        #Now, investigate why df['tradedate'].diff(1) IS NOT < 20 days and make sure all 120 exceptions are valid
        data.drop(['index'], axis=1, inplace=True)


        #after investigating, looks like all exceptions are valid and that ret5d,ret10d, ret15d, and ret20d should be NaN'd out
        for i in problem_rows:
            
            for j in range(1): #0
                if( (data['osid'].diff(1))[i+j] ==0 ): #make sure we are still in one unique OSID
                    data.loc[i+j,'ret1d'] = np.nan
                    data.loc[i+j,'alpha1d'] = np.nan
        #j iterates 0 days forward, 
        #i+j days forward could potentially have a corrupted return that used the return on spot i-1


        
            for j in range(4): #0,1,2,3
                if( (data['osid'].diff(1))[i+j] ==0 ): #make sure we are still in one unique OSID
                    data.loc[i+j,'ret5d'] = np.nan
                    data.loc[i+j,'alpha5d'] = np.nan
        #j iterates 4 days forward, 
        #i+j days forward could potentially have a corrupted return that used the return on spot i-1


            for j in range(9): #0,1,2,3,...,8
                if( (data['osid'].diff(1))[i+j] ==0 ): #make sure we are still in one unique OSID
                    data.loc[i+j,'ret10d'] = np.nan
                    data.loc[i+j,'alpha10d'] = np.nan
        #j iterates 9 days forward, 
        #i+j days forward could potentially have a corrupted return that used the return on spot i-1

            for j in range(14): #0,1,2,3,...,13
                if( (data['osid'].diff(1))[i+j] ==0 ): #make sure we are still in one unique OSID
                    data.loc[i+j,'ret15d'] = np.nan
                    data.loc[i+j,'alpha15d'] = np.nan
        #j iterates 14 days forward, 
        #i+j days forward could potentially have a corrupted return that used the return on spot i-1


            for j in range(19): #0,1,2,3,...,19
                if( (data['osid'].diff(1))[i+j] ==0 ): #make sure we are still in one unique OSID
                    data.loc[i+j,'ret20d'] = np.nan
                    data.loc[i+j,'alpha20d'] = np.nan
        #j iterates 19 days forward, 
        #i+j days forward could potentially have a corrupted return that used the return on spot i-1



        #data[data.index.isin(problem_rows3)] #this is to visually check 

        ###########################################################################################
        ###########################################################################################
        ###########################################################################################

        #for example, we have the following for ret1d, ret5d, ret10d, ret15d, and ret20d
        example_row = problem_rows3[-1].copy()
        data[example_row-1:example_row+22]  #run this line to see what the result looks like if interested
        
        
        ##############################################################################################################
        ##################### (BELOW) GET lagged returns and alphas for 5,10,15,&20 day horizons #####################
        ##############################################################################################################
        #Idea:
        #Lets lag ret5d by 5 days so it represents forward looking returns instead of backward looking returns
        #that is, right now 'ret5d' (at time t) is the trailing 5 day cumulitive return (from t-4,t-3,t-2,t-1,t)
        #we want to create another vairable called 'ret5d_lag'
        #that is, 'ret5d_lag' (at time t) will represent the 5 day cumulitive return (on t+1,t+2,t+3,t+4,t+5)

        for i in [1,5,10,15,20]:
            data['ret'+str(i)+'d_lag'] = data[['osid','ret'+str(i)+'d']].groupby(['osid'])['ret'+str(i)+'d'].shift(-i) #group by OSID, shift 5d by 5, 10d by 10, etc...
            data['alpha'+str(i)+'d_lag'] = data[['osid','alpha'+str(i)+'d']].groupby(['osid'])['alpha'+str(i)+'d'].shift(-i) #group by OSID

        data

        #print(data['ret5d_lag'][0:10])
        #print(data['ret5d'][0:15])  #use this line and the one above to see the results

        ##############################################################################################################
        ##################### (ABOVE) GOT lagged returns and alphas for 5,10,15,&20 day horizons #####################
        ##############################################################################################################

        #save data to use later without running all the code above again every time
        data.to_pickle('data_with_rets.pkl')
        ### Now we have ret5d, ret10d, ret15d, and ret20d for valid indices only
        ### Now we have alpha5d, alpha10d, alpha15d, and alpha20d for valid indices only


# In[3]:


def trade_on_valid_lags_only(data):
    
    #VALID LAGS ONLY!
    #same company, but longer than 20 day difference between observations --> don't trade even if 'trade_tomorrow'==1

    xx = list(abs((data[['date','osid']].groupby(['osid'])['date'].diff(-1)).dt.days) > 20) # abs(today - tomorrow)>20
    yy = list(data['buy_tomorrow']>0) #signal is saying buy, can be =1, or can be =2 or =3 if it is moved on 2 or 3 indices at once
    zz = list(data['sell_tomorrow']>0) #signal is saying sell

    zz_buy = [(a and b) for a, b in zip(xx, yy)] #this is where you should not trade!
    zz_sell = [(a and b) for a, b in zip(xx, zz)] #this is where you should not trade!

    print('\nThere are ', sum(zz_buy), ' problems for buying')
    print('\nThere are ', sum(zz_sell), ' problems for selling')


    problem_rowsa = list(np.where(zz_buy)[0]) #gets the row index of the rows for which zz is True
    problem_rowsb = list(np.where(zz_sell)[0]) #gets the row index of the rows for which zz is True
    problem_rows = problem_rowsa + problem_rowsb

    problem_rows2 = [i+1 for i in problem_rows]

    problem_rows3 = np.sort(problem_rows + problem_rows2)
    problem_rows3 #this is the indices we want to check to confirm that there is a large gap in the date (for example, instead of being the next day, it is instead the next year)


    data['index'] = list(data.index) #create a col called index
    data[data.index.isin(problem_rows3)][['tradedate','symbol', 'coname','volume','avgvol50d','buy_tomorrow', 'sell_tomorrow']]  #THESE ROWS ARE PROBLEMATIC!

    #Now, investigate why df['tradedate'].diff(1) IS NOT < 20 days and make sure all 120 exceptions are valid
    #data.drop(['index'], axis=1, inplace=True)

# WE SHOULDN'T TRADE ON ANY OF THE EVENTS ABOVE!!! 
#Because the dates jump from one year to another

#after investigating, looks like all exceptions are valid 
#and that we shouldn't trade on the events from above

    data.loc[data.index.isin(problem_rowsa),'buy_tomorrow'] = 0 #erase signal and dont trade on it
    data.loc[data.index.isin(problem_rowsb),'sell_tomorrow'] = 0 #erase signal and dont trade on it

    data[data.index.isin(problem_rows3)][['tradedate','symbol', 'coname','sp100f','sp500f','nasdaq100f','buy_tomorrow', 'sell_tomorrow']]  #THESE ROWS ARE PROBLEMATIC!
    #data.drop(['index'], axis=1, inplace=True)

    print('Total number of Buys: ',np.count_nonzero(data['buy_tomorrow']))
    print('Total number of Sells: ',np.count_nonzero(data['sell_tomorrow']))
    
    return(data)


# In[4]:


def pull_signal_forward(data):
    data['ones'] = 1
    data_smaller = data.copy()
    data_smaller = data_smaller.loc[data_smaller['ret1d_lag'].notna(),:].reset_index(drop=True).copy()  #DROP MISSING RETURNS
    data_smaller = data_smaller.loc[data_smaller['alpha1d_lag'].notna(),:].reset_index(drop=True).copy()  #DROP MISSING ALPHAS

    data_smaller['index'] = list(data_smaller.index) #create a col called index

    ################################################################################################################
    ################################################################################################################
    #Goal: make this part much faster!
    #VALID LAGS ONLY! THIS WILL TAKE A LONG TIME TO RUN
    #same company, AND smaller than 20 day difference in-between observations --> set 'trade_tomorrow'=1 in the next 4,9,14,19 days forward

    xx = list(abs((data_smaller[['date','osid']].groupby(['osid'])['date'].diff(-1)).dt.days) < 20) # abs(today - tomorrow)<20
    yy = list(data_smaller['buy_tomorrow']>0) #signal is saying buy, can be =1, or can be =2 or =3 if it is moved on 2 or 3 indices at once
    zz = list(data_smaller['sell_tomorrow']>0) #signal is saying sell

    zz_buy = [(a and b) for a, b in zip(xx, yy)] #this is where you should bring buy_signal forward!
    zz_sell = [(a and b) for a, b in zip(xx, zz)] #this is where you should bring sell_signal forward!

    rowsa = list(np.where(zz_buy)[0]) #gets the row index of the rows for which we should bring buy_signal forward
    rowsb = list(np.where(zz_sell)[0]) #gets the row index of the rows for which we should bring sell_signal forward

    ################################################################################################################

    #now, ALSO create data_smaller5, data_smaller10, data_smaller15, data_smaller20: holding periods of 5,10,15,20 days
    data_smaller5 = data_smaller.copy()
    for i in range(1,5): #lag signal i=1,2,3,4 days forward
        data_smaller5['buy_tomorrow'] += np.where(data_smaller5[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsa), 1, 0)
        data_smaller5['sell_tomorrow'] += np.where(data_smaller5[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsb), 1, 0)
        #if index (i spots above) is in rowsa/OR/rowsb, then set 'buy_tomorrow'/'sell_tomorrow'=1

    data_smaller10 = data_smaller5.copy()
    for i in range(5,10): #lag signal i=5,...,9 days forward
        data_smaller10['buy_tomorrow'] += np.where(data_smaller10[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsa), 1, 0)
        data_smaller10['sell_tomorrow'] += np.where(data_smaller10[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsb), 1, 0)
        #if index (i spots above) is in rowsa/OR/rowsb, then set 'buy_tomorrow'/'sell_tomorrow'=1

    data_smaller15 = data_smaller10.copy()
    for i in range(10,15): #lag signal i=10,...,14 days forward
        data_smaller15['buy_tomorrow'] += np.where(data_smaller15[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsa), 1, 0)
        data_smaller15['sell_tomorrow'] += np.where(data_smaller15[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsb), 1, 0)
        #if index (i spots above) is in rowsa/OR/rowsb, then set 'buy_tomorrow'/'sell_tomorrow'=1

    data_smaller20 = data_smaller15.copy()
    for i in range(15,20): #lag signal i=15,...,19 days forward
        data_smaller20['buy_tomorrow'] += np.where(data_smaller20[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsa), 1, 0)
        data_smaller20['sell_tomorrow'] += np.where(data_smaller20[['osid','index']].groupby(['osid'])['index'].shift(i).isin(rowsb), 1, 0)
        #if index (i spots above) is in rowsa/OR/rowsb, then set 'buy_tomorrow'/'sell_tomorrow'=1

    return(data_smaller,data_smaller5,data_smaller10,data_smaller15,data_smaller20)


# In[5]:


def mini_summ(ew_alphas,ew_returns, horizon_days):
    #horizon_days=1,5,10,15,20
    answers = np.zeros(8)
    answers[0] = (252/horizon_days)*100*ew_alphas.mean()
    answers[1] = np.sqrt((252/horizon_days)) * (100) * np.sqrt( np.sum(ew_alphas**2)/(len(ew_alphas)-1) )
    answers[2] = answers[0]/answers[1]
    answers[3] = ew_alphas.skew()

    answers[4] = (252/horizon_days)*100*ew_returns.mean()
    answers[5] = np.sqrt((252/horizon_days)) * (100) * np.sqrt( np.sum((ew_returns-np.mean(ew_returns))**2)/(len(ew_returns)-1) )
    answers[6] = answers[4]/answers[5]
    answers[7] = ew_returns.skew()
                 
    return(answers)

############################################################################################################
############################################################################################################
############################################################################################################

def get_portfolio_returns(data_smaller):
    data_buy = data_smaller[data_smaller['buy_tomorrow']>0]

    #without weights
    ew_returns_buy2 = data_buy.groupby(['date'])['ret1d_lag'].mean().to_frame().reset_index().rename(columns={'ret1d_lag':'ew_return'})
    ew_returns_buy2['ew_return'] = round(ew_returns_buy2['ew_return'],4) #round to 4 decimal digits
    ew_returns_buy2

    ew_alphas_buy2 = data_buy.groupby(['date'])['alpha1d_lag'].mean().to_frame().reset_index().rename(columns={'alpha1d_lag':'ew_alpha'})
    ew_alphas_buy2['ew_alpha'] = round(ew_alphas_buy2['ew_alpha'],4) #round to 4 decimal digits
    ew_alphas_buy2

    ##with weights
    w = data_buy.groupby(['date'])['ones'].sum().to_frame().reset_index().rename(columns={'ones':'N'})
    w['weight'] = 1/w['N']
    data_buy = data_buy.merge(w, how='left', on='date').copy()
    data_buy['ew_return_contr'] = data_buy['ret1d_lag'] * data_buy['weight'] #return contribtuion !
    ew_returns_buy = data_buy.groupby(['date'])['ew_return_contr'].sum().to_frame().reset_index().rename(columns={'ew_return_contr':'ew_return'})
    ew_returns_buy['ew_return'] = round(ew_returns_buy['ew_return'],4) #round to 4 decimal digits
    #print(ew_returns_buy)

    data_buy['ew_alpha_contr'] = data_buy['alpha1d_lag'] * data_buy['weight'] #return contribtuion !
    ew_alphas_buy = data_buy.groupby(['date'])['ew_alpha_contr'].sum().to_frame().reset_index().rename(columns={'ew_alpha_contr':'ew_alpha'})
    ew_alphas_buy['ew_alpha'] = round(ew_alphas_buy['ew_alpha'],4) #round to 4 decimal digits
    #print(ew_returns_buy)

    print('\nThere are ',sum(ew_returns_buy['ew_return'] != ew_returns_buy2['ew_return']),' returns that dont match using .mean() vs. using equal weights')
    print('\nThere are ',sum(ew_alphas_buy['ew_alpha'] != ew_alphas_buy2['ew_alpha']),' alphas that dont match using .mean() vs. using equal weights')

    #for i in range(len(ew_returns_buy['ew_return'])):
    #    if(ew_returns_buy['ew_return'][i] != ew_returns_buy2['ew_return'][i]):
    #        print(ew_returns_buy['ew_return'][i], ' != ', ew_returns_buy2['ew_return'][i])
    #these numbers actually match and are no problem !

    ew_returns_buy = ew_returns_buy.sort_values(by='date')

    ew_alphas_buy = ew_alphas_buy.sort_values(by='date')

    ew_returns_buy['ew_alpha'] = ew_alphas_buy['ew_alpha']
    ew_returns_buy

    ####################################################################################################################################
    ####################################################################################################################################
    ####################################################################################################################################

    data_sell = data_smaller[data_smaller['sell_tomorrow']>0]

    #without weights
    ew_returns_sell2 = data_sell.groupby(['date'])['ret1d_lag'].mean().to_frame().reset_index().rename(columns={'ret1d_lag':'ew_return'})
    ew_returns_sell2['ew_return'] = round(ew_returns_sell2['ew_return'],4) #round to 4 decimal digits
    ew_returns_sell2

    ew_alphas_sell2 = data_sell.groupby(['date'])['alpha1d_lag'].mean().to_frame().reset_index().rename(columns={'alpha1d_lag':'ew_alpha'})
    ew_alphas_sell2['ew_alpha'] = round(ew_alphas_sell2['ew_alpha'],4) #round to 4 decimal digits
    ew_alphas_sell2


    ##with weights
    w = data_sell.groupby(['date'])['ones'].sum().to_frame().reset_index().rename(columns={'ones':'N'})
    w['weight'] = 1/w['N']
    data_sell = data_sell.merge(w, how='left', on='date').copy()
    data_sell['ew_return_contr'] = data_sell['ret1d_lag'] * data_sell['weight'] #return contribtuion !
    ew_returns_sell = data_sell.groupby(['date'])['ew_return_contr'].sum().to_frame().reset_index().rename(columns={'ew_return_contr':'ew_return'})
    ew_returns_sell['ew_return'] = round(ew_returns_sell['ew_return'],4) #round to 4 decimal digits
    #print(ew_returns_sell)

    data_sell['ew_alpha_contr'] = data_sell['alpha1d_lag'] * data_sell['weight'] #return contribtuion !
    ew_alphas_sell = data_sell.groupby(['date'])['ew_alpha_contr'].sum().to_frame().reset_index().rename(columns={'ew_alpha_contr':'ew_alpha'})
    ew_alphas_sell['ew_alpha'] = round(ew_alphas_sell['ew_alpha'],4) #round to 4 decimal digits
    #print(ew_returns_buy)


    print('\nThere are ',sum(ew_returns_sell['ew_return'] != ew_returns_sell2['ew_return']),' returns that dont match using .mean() vs. using equal weights')
    print('\nThere are ',sum(ew_alphas_sell['ew_alpha'] != ew_alphas_sell2['ew_alpha']),' alphas that dont match using .mean() vs. using equal weights')
    #for i in range(len(ew_returns_sell['ew_return'])):
    #    if(ew_returns_sell['ew_return'][i] != ew_returns_sell2['ew_return'][i]):
    #        print(ew_returns_sell['ew_return'][i], ' != ', ew_returns_sell2['ew_return'][i])
    #these numbers actually match and are no problem !

    ew_returns_sell = ew_returns_sell.sort_values(by='date')
    ew_returns_sell

    ew_alphas_sell = ew_alphas_sell.sort_values(by='date')
    ew_alphas_sell

    ew_returns_sell['ew_alpha'] = ew_alphas_sell['ew_alpha']
    ew_returns_sell

    ####################################################################################################################################
    ####################################################################################################################################
    ####################################################################################################################################

    #make sure weights sum to 1 on each date that we want to trade!
    x = sum(data_buy.groupby(['date'])['weight'].sum() < 1.0000001) #do this because weights may be 0.99999999999999999 and not =1.0
    y = sum(data_buy.groupby(['date'])['weight'].sum() > 0.9999999) 
    total = len(data_buy.groupby(['date'])['weight'].sum())
    if(x==total & y==total):
        print('\nAll weights sum to one for Buys')
    else:
        print('\n!!!Weights dont sum to one for Buys!!!')


    x = sum(data_sell.groupby(['date'])['weight'].sum() < 1.0000001) #do this because weights may be 0.99999999999999999 and not =1.0
    y = sum(data_sell.groupby(['date'])['weight'].sum() > 0.9999999) 
    total = len(data_sell.groupby(['date'])['weight'].sum())
    if(x==total & y==total):
        print('\nAll weights sum to one for Sells')
    else:
        print('\n!!!Weights dont sum to one for Sells!!!')

    ####################################################################################################################################
    ####################################################################################################################################
    ####################################################################################################################################

    ew_returns_buy['cum_ret'] = (ew_returns_buy['ew_return']+1).cumprod()
    ew_returns_buy['cum_alpha'] = (ew_returns_buy['ew_alpha']+1).cumprod()

    ew_returns_sell['cum_ret'] = (ew_returns_sell['ew_return']+1).cumprod()
    ew_returns_sell['cum_alpha'] = (ew_returns_sell['ew_alpha']+1).cumprod()

    ####################################################################################################################################
    ####################################################################################################################################
    ####################################################################################################################################

    summ_stats_buy = mini_summ(ew_returns_buy['ew_alpha'], ew_returns_buy['ew_return'], horizon_days=1)
    summ_stats_sell = mini_summ(ew_returns_sell['ew_alpha'], ew_returns_sell['ew_return'], horizon_days=1)

    return(ew_returns_buy, ew_returns_sell, data_buy, data_sell, summ_stats_buy, summ_stats_sell)


# In[6]:


def turnover(data_buy): #can also input 'data_sell' to get 'turnover_sell'!!
    turnover_buy = data_buy[['date','osid','symbol','ret1d_lag','alpha1d_lag','buy_tomorrow','N','weight',                         'ew_return_contr','ew_alpha_contr']].sort_values(by=['osid','date']).reset_index(drop=True)
    turnover_buy['weight_end'] = turnover_buy['weight'] * (1+turnover_buy['ret1d_lag'])

    turnover_buy = turnover_buy.merge(turnover_buy[['date','weight_end']]                             .groupby(['date']).sum().reset_index().rename(columns={"weight_end":"weight_end_sum"}),                             on=['date'],how='left')
    #group by date and sum up all the 'weight_end's to get total weight of Portfolio (now no longer =1)

    turnover_buy.loc[turnover_buy['weight_end_sum'] <=0,'weight_end_sum'] = np.nan
    turnover_buy['weight_in_portfolio_end'] = turnover_buy['weight_end'] / turnover_buy['weight_end_sum']

    turnover_buy['lag_weight_in_portfolio_end'] = turnover_buy[['osid','weight_in_portfolio_end']].groupby('osid').shift(1)['weight_in_portfolio_end']
    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################
    ### Note: IF a stock is traded today, and then not traded again for more than 15 days, 
    #THEN we should take it out of the portfolio before putting it back in 
    #(ie weight should go from + number to 0, then from 0 to a + number again)

    ### That is, if a stock is traded discretely once today and never again until 1 year later, 
    #then turnover for that stock should be 100% when it is traded today
    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################

    #BEFORE WE COMPUTE CHANGE IN WEIGHTS: find the problem rows when change in date > 5:
    xx = list((turnover_buy[['date','osid']].groupby(['osid'])['date'].diff(1)).dt.days >5) #(today - yesterday) > 5
    #IF the difference between trade dates for one stock is > 5 days
    #THEN, we should change 'lag_weight_end' to 0 !

    print('\nThere are ', sum(xx), ' problems for turnover using lag_weight_in_portfolio_end')

    problem_rows = list(np.where(xx)[0]) #gets the row index of the rows for which xx is True

    problem_rows2 = [i-1 for i in problem_rows]

    problem_rows3 = np.sort(problem_rows + problem_rows2)
    problem_rows3 #this is the indices we want to check to confirm that there is a large gap in the date (for example, instead of being the next day, it is instead the next year)


    turnover_buy['index'] = list(turnover_buy.index) #create a col called index
    turnover_buy[turnover_buy.index.isin(problem_rows3)]
    #THESE ROWS ARE PROBLEMATIC AND WE SHOULD NOT use their 'lag_weight_in_portfolio_end' to get T-cost!

    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################

    #after investigating, looks like all exceptions are valid and that we shouldn't use 'lag_weight_in_portfolio_end'!
    #DONT DO THE FOLLOWING 3 LINES BELOW
    #for i in problem_rows:
    #    if( (turnover_buy['osid'].diff(1))[i] ==0 ): #make sure we are still in one unique OSID
    #        turnover_buy.loc[i,'lag_weight_in_portfolio_end'] = 0

    turnover_buy.loc[turnover_buy.index.isin(problem_rows),'lag_weight_in_portfolio_end'] = 0
    turnover_buy[turnover_buy.index.isin(problem_rows3)]
    #THESE ROWS ARE PROBLEMATIC AND WE SHOULD NOT use their 'lag_weight_in_portfolio_end' to get T-cost!
    #turnover_buy.drop(['index'], axis=1, inplace=True)

    ##############################################################################################################
    ##############################################################################################################
    ##############################################################################################################

    ## Now we can compute the change in weights using 'lag_weight_in_portfolio_end'
    #computing change in portfolio weights
    turnover_buy['delta_share'] = turnover_buy['weight'] - turnover_buy['lag_weight_in_portfolio_end'].fillna(0)
    turnover_buy['delta_share'] = turnover_buy['delta_share'].abs()

    turnover_buy = turnover_buy.sort_values(by=['date']).reset_index(drop=True)
    turnover_buy2 = turnover_buy[['date','delta_share']].groupby(['date']).sum().reset_index().rename(columns={'delta_share':'turnover'})
    turnover_buy2['turnover'] = round(turnover_buy2['turnover'],4)
    turnover_buy2

    print('Total number of days when turnover is not equal to 100%: ',sum(list(turnover_buy2['turnover']!=1)))
    #listt = []
    #for i in range(len(turnover_buy2)):
    #    if(turnover_buy2['turnover'][i]!=1):
    #        #print('Row ',i,' has ',turnover_buy2['turnover'][i],'% turnover')
    #        listt.append(i)

    ##### IF output is: "Row  122  has  0.0 % turnover" --> 
    ##### check row 122 via the following 2 lines below to make sure we are trading the same stock 
    ##### on 2 consecutive days and that turnover should truly be =0

    #turnover_buy[turnover_buy['date']==turnover_buy2['date'][122-1]]
    #turnover_buy[turnover_buy['date']==turnover_buy2['date'][122]]

    ##### ^^^ This confirms that we dont have any turnover on that day ^^^
    ##### since we ONLY (N=1 on both days) trade the SAME (symbol is the same on both days) stock on both days
    
    return(turnover_buy2)


# ##### NOTE: Now we have the following 2 ew_returns series, grouped by date:
# #### 1) ew_returns_buy #(includes both ret and alpha, cum_ret and cum_alpha, AND TURNOVER NOW!)
# #### 2) ew_returns_sell #(includes both ret and alpha, cum_ret and cum_alpha, AND TURNOVER NOW!)
# 
# ##### NOTE: The weights (and the return contributions) are in the following dataframes:
# #### 1) data_buy
# #### 2) data_sell
# 
# ##### NOTE: The summary statistics are in the following lists:
# #### 1) summ_stats_buy
# #### 2) summ_stats_sell

# In[7]:


def T_cost_and_capital(ew_returns_buy): #can also input 'ew_returns_sell' to get 'T_cost_sell'
    capital = [1]
    T_cost_buy = 0
    alpha_capital = [1]

    for i in range(len(ew_returns_buy['ew_return'])):
        starting_capital_afterTcost = capital[i] - capital[i]*0.0005*ew_returns_buy['turnover'][i] 
        #^^5 bp when we buy * turnover when we buy^^^^^^^
        returns = 1+ew_returns_buy['ew_return'][i]
        ending_capital = (starting_capital_afterTcost*returns)

        if(i != len(ew_returns_buy['ew_return']) - 1):
            ending_capital_afterTcost = ending_capital - ending_capital*0.0005*ew_returns_buy['turnover'][i+1] 
            #5 bp when we sell * turnover when we sell ^^^^^^^^^^^^^^^^^^^^^^^
        else: #assumption: 100% turnover on the last day of sample when we sell and close the position!
            ending_capital_afterTcost = ending_capital - ending_capital*0.0005*1
            #5 bp when we sell*(turnover=100%) when we sell (on the last day of sample)^^^^^^^^^^^^^^^^^^^^^^^

        capital.append(ending_capital_afterTcost)
        T_cost_buy += (capital[i]-starting_capital_afterTcost) + (ending_capital-ending_capital_afterTcost)
        #take off 5bp when we buy^^^^^^^^^^^^^^^^^^^^^^^^, then take off 5bp when we sell^^^^^^^^^^^^^^^

    for i in range(len(ew_returns_buy['ew_alpha'])):
        starting_capital_afterTcost = alpha_capital[i] - alpha_capital[i]*0.0005*ew_returns_buy['turnover'][i] 
        #^^5 bp when we buy * turnover when we buy^^^^^^^
        returns = 1+ew_returns_buy['ew_alpha'][i]
        ending_capital = (starting_capital_afterTcost*returns)

        if(i != len(ew_returns_buy['ew_alpha']) - 1):
            ending_capital_afterTcost = ending_capital - ending_capital*0.0005*ew_returns_buy['turnover'][i+1] 
            #5 bp when we sell * turnover when we sell ^^^^^^^^^^^^^^^^^^^^^^^
        else: #assumption: 100% turnover on the last day of sample when we sell and close the position!
            ending_capital_afterTcost = ending_capital - ending_capital*0.0005*1
            #5 bp when we sell*(turnover=100%) when we sell (on the last day of sample)^^^^^^^^^^^^^^^^^^^^^^^

        alpha_capital.append(ending_capital_afterTcost)


    del capital[0] #remove the initial capital 1 so that the size of capital fits in df 'ew_returns_buy'
    del alpha_capital[0] #remove the initial capital 1 so that the size of capital fits in df 'ew_returns_buy'
    ew_returns_buy['capital'] = capital
    ew_returns_buy['alpha_capital'] = alpha_capital

    return(ew_returns_buy,T_cost_buy)


# In[8]:


def print_summary(ew_returns_buy, ew_returns_sell, trade_days_summ_stats, non_trade_days_summ_stats, trade_days_T_cost, non_trade_days_T_cost):
    print('Summary Statistics:')
    print(' ')
    SS = pd.DataFrame(columns=[' ','Buy', 'Sell'])
    SS.loc[0, ' '] = 'Ann. Mean Alpha (without T-costs)'
    SS.loc[0, 'Buy'] = trade_days_summ_stats[0]
    SS.loc[0, 'Sell'] = non_trade_days_summ_stats[0]
    SS.loc[1, ' '] = 'Ann. Tracking Error'
    SS.loc[1, 'Buy'] = trade_days_summ_stats[1]
    SS.loc[1, 'Sell'] = non_trade_days_summ_stats[1]
    SS.loc[2, ' '] = 'Ann. Information Ratio (without T-costs)'
    SS.loc[2, 'Buy'] = trade_days_summ_stats[2]
    SS.loc[2, 'Sell'] = non_trade_days_summ_stats[2]
    SS.loc[3, ' '] = 'Skewness of Alpha'
    SS.loc[3, 'Buy'] = trade_days_summ_stats[3]
    SS.loc[3, 'Sell'] = non_trade_days_summ_stats[3]
    SS.loc[4, ' '] = '--------------------------'
    SS.loc[4, 'Buy'] = '--'
    SS.loc[4, 'Sell'] = '--'
    SS.loc[5, ' '] = 'Ann. Mean Return (without T-costs)'
    SS.loc[5, 'Buy'] = trade_days_summ_stats[4]
    SS.loc[5, 'Sell'] = non_trade_days_summ_stats[4]
    SS.loc[6, ' '] = 'Ann. Volatility'
    SS.loc[6, 'Buy'] = trade_days_summ_stats[5]
    SS.loc[6, 'Sell'] = non_trade_days_summ_stats[5]
    SS.loc[7, ' '] = 'Ann. Sharpe Ratio (without T-costs)'
    SS.loc[7, 'Buy'] = trade_days_summ_stats[6]
    SS.loc[7, 'Sell'] = non_trade_days_summ_stats[6]
    SS.loc[8, ' '] = 'Skewness of Return'
    SS.loc[8, 'Buy'] = trade_days_summ_stats[7]
    SS.loc[8, 'Sell'] = non_trade_days_summ_stats[7]
    SS.loc[9, ' '] = 'Average Turnover'
    SS.loc[9, 'Buy'] = np.mean(ew_returns_buy['turnover'])
    SS.loc[9, 'Sell'] = np.mean(ew_returns_sell['turnover'])
    SS.loc[10, ' '] = 'T-costs (on 1 unit of capital)'
    SS.loc[10, 'Buy'] = trade_days_T_cost
    SS.loc[10, 'Sell'] = non_trade_days_T_cost
    SS = SS.set_index(' ').copy()
    pd.set_option('display.expand_frame_repr', False, 'display.float_format', '{:,.2f}'.format)
    print(SS)
    print(" ")


# In[9]:


def get_plots_and_summary(ew_returns_buy,ew_returns_sell, summ_stats_buy, summ_stats_sell, T_cost_buy, T_cost_sell):    # New start: get performance plots
    fig, axs = plt.subplots(2,2,figsize=(9,9))
    
    axs[0, 0].plot(ew_returns_buy['date'],ew_returns_buy['capital']) #cumulative_ret
    axs[0, 0].set_xlabel('Date')
    axs[0, 0].set_ylabel("Cumulative Return")
    axs[0, 0].set_title("Cum Return (net of T-cost) on buy days")

    axs[1, 0].plot(ew_returns_buy['date'],ew_returns_buy['alpha_capital']) #cumulative_ret
    axs[1, 0].set_xlabel('Date')
    axs[1, 0].set_ylabel("Cumulative Alpha")
    axs[1, 0].set_title("Cum Alpha (net of T-cost) on buy days")
    print('Number of buy periods', len(ew_returns_buy['capital']))
    print('Check number of buy periods', len(ew_returns_buy['alpha_capital']))


    axs[0, 1].plot(ew_returns_sell['date'],ew_returns_sell['capital'], color="r") #cumulative_ret
    axs[0, 1].set_xlabel('Date')
    axs[0, 1].set_ylabel("Cumulative Return")
    axs[0, 1].set_title("Cum Return (net of T-cost) on sell days")

    axs[1, 1].plot(ew_returns_sell['date'],ew_returns_sell['alpha_capital'], color="r") #cumulative_ret
    axs[1, 1].set_xlabel('Date')
    axs[1, 1].set_ylabel("Cumulative Alpha")
    axs[1, 1].set_title("Cum Alpha (net of T-cost) on sell days")
    print('Number of sell periods', len(ew_returns_sell['capital']))
    print('Check number of sell periods', len(ew_returns_sell['alpha_capital']))

    fig.tight_layout()
    print_summary(ew_returns_buy, ew_returns_sell, summ_stats_buy, summ_stats_sell, T_cost_buy, T_cost_sell)


# In[10]:


def FF_regression(regression_input):
    y = regression_input['ew_return-RF']
    X = regression_input[['Mkt-RF','SMB','HML','RMW','CMA','MOM']]
    X = sm.add_constant(X) #intercept
    model = sm.OLS(y,X)
    results = model.fit()
    print(results.summary())


# In[11]:


#one function to do it all, instead of repeating all these lines 5 times for each strategy we want to test
def get_returns_turnover_Tcost_summary_plusPlots(data_smaller): #will input data_smaller, data_smaller5, data_smaller10, data_smaller15, and data_smaller20 in here
    
#FIRST GET FF6 factors    
    FF_MOM = pd.read_csv('F-F_Momentum_Factor_daily.csv')[['date','MOM']]
    for i in range(len(FF_MOM)):
        FF_MOM.loc[i,'date'] = FF_MOM.loc[i,'date'].replace(',', '') #remove the comma from the date
    FF_MOM['date'] = pd.to_datetime(FF_MOM['date'], format='%Y%m%d', errors='ignore')

    FF5 = pd.read_csv('F-F_Research_Data_5_Factors_2x3_daily.csv').rename(columns={'Unnamed: 0':'date'})
    FF5['date'] = pd.to_datetime(FF5['date'], format='%Y%m%d', errors='ignore')

    FF6 = FF5.merge(FF_MOM, how='left', on='date').copy()
    FF6[['Mkt-RF','SMB','HML','RMW','CMA','RF','MOM']] = FF6[['Mkt-RF','SMB','HML','RMW','CMA','RF','MOM']]/100
    FF6
    
    
#THEN GET PORTFOLIO RETURNS, TURNOVER, T_COST, AND SUMMARY STATS + PLOTS
    ew_returns_buy, ew_returns_sell, data_buy, data_sell, summ_stats_buy, summ_stats_sell = get_portfolio_returns(data_smaller)

    buy_regresion = ew_returns_buy[['date','ew_return']].merge(FF6, how='left', on='date').copy()
    sell_regression = ew_returns_sell[['date','ew_return']].merge(FF6, how='left', on='date').copy()
    buy_regresion['ew_return-RF'] = buy_regresion['ew_return'] - buy_regresion['RF']
    sell_regression['ew_return-RF'] = sell_regression['ew_return'] - sell_regression['RF']
    
    print("\033[1m" + 'Regression for Buy Portfolio' + "\033[0m")
    FF_regression(buy_regresion)
    print("\033[1m" + 'Regression for Sell Portfolio' + "\033[0m")
    FF_regression(sell_regression)
    
    turnover_buy2 = turnover(data_buy) #get turnover
    turnover_sell2 = turnover(data_sell) #get turnover
    ew_returns_buy['turnover'] = turnover_buy2['turnover']
    ew_returns_sell['turnover'] = turnover_sell2['turnover']
    ew_returns_sell

    ew_returns_buy,T_cost_buy = T_cost_and_capital(ew_returns_buy)
    ew_returns_buy
    ew_returns_sell,T_cost_sell = T_cost_and_capital(ew_returns_sell)
    ew_returns_sell

    get_plots_and_summary(ew_returns_buy,ew_returns_sell, summ_stats_buy, summ_stats_sell, T_cost_buy, T_cost_sell)


# In[ ]:





# In[ ]:




