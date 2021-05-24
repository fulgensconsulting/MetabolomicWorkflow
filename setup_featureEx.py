""" Setup for feature creation and extraction 

    Fulgens Consulting 
""" 

def extendWindowsDuo(df, lowCol, highCol, sort=True):
    """ Create bins of overlapping values without joins (via iteration), in the
        event self-joins would be memory-prohibitive. Join the bins onto the
        index of the original df """
    
    if sort is True:
        df.sort_values(lowCol, inplace=True)
        df.reset_index(inplace=True, drop=True)

    ixCollect = []
    ixLow, ixHigh = dict(), dict()
    for i, low, high in zip(df.index, df[lowCol], df[highCol]):
        if i == 0:
            ixCollect.append(i)
            actLow, actHigh = low, high 
        
        if low > actHigh:
            for ix in ixCollect:
                ixLow[ix] = actLow 
                ixHigh[ix] = actHigh
            ixCollect = [i]
            actLow, actHigh = low, high 
        
        elif low <= actHigh and i != (len(df)-1):
            ixCollect.append(i)
            actHigh = high

        elif i == (len(df)-1):
            for ix in ixCollect + [i]:
                ixLow[ix] = actLow 
                ixHigh[ix] = high

    if len(df)-1 not in ixLow and len(df)-1 not in ixHigh:
        ixLow[len(df)-1] = actLow 
        ixHigh[len(df)-1] = high

    df[lowCol] = df.index.map(ixLow)
    df[highCol] = df.index.map(ixHigh)

    if any(df[lowCol] > df[highCol]):
        raise Exception('Window extends yielded "low" values g.t. "high"')

    return df 


def prepMZRTSort(df, mzCol='mz', rtCol='rt'):
    """ Prep dataframe with sort values on mass and rt, employ rounding to get 
        approximation """
        
    df['mzround'] = df[mzCol].round(4)
    df['rtround'] = df[rtCol].divide(100).round(1)
    df.sort_values(['mzround', 'rtround'], inplace=True)
    df.reset_index(inplace=True, drop=True)
    df.drop(['mzround', 'rtround'], axis='columns', inplace=True)

    return df 
