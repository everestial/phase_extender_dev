import matplotlib.pyplot as plt
import pandas as pd

''' function that computes several stats from haplotype block before and after phase extension.
These data can be used to compare the improvements in phase extension, and plotting histogram.
This can include the size based on genome co-ordinate of the "POS" in each block.
This can also include the average number of variants per haplotype before/after phase extension.
and make a plot out of it - using pyPlot, MatlibPlot or R. '''
def compute_haplotype_stats(hap_data, soi, prefix, outputdir) :

    print()
    ''' this function computes the stats of the haplotype file both before and after phase extension.
        - we pass in "prefix" as "initial" vs. "final" for the haplotype block before and after
          extension to compute statistics for respective phased file. '''

    ''' Step 09 - A : write this haplotype information in a more concise format.
        - in this step we are restructuring the data so proper statistics can be handled properly. '''
    hap_stats = hap_data.groupby(['CHROM', soi +':PI'], sort=False)\
        .size().rename('num_Vars_by_PI')\
        .reset_index(level=1).astype(str)\
        .groupby(level=0).agg(','.join).reset_index()

    ''' add other required columns '''
    hap_stats['range_of_PI'] = hap_data.groupby(['CHROM', soi +':PI'], sort=False)['POS']\
        .apply(lambda g: g.max() - g.min()).rename('range').reset_index(level=1)\
        .astype(str).groupby(level=0).agg(','.join).reset_index()['range']

    hap_stats['total_haplotypes'] = hap_stats\
        .apply(lambda row: len(row[soi +':PI'].split(',')), axis=1)

    hap_stats['total_Vars'] = hap_stats \
        .apply(lambda row: sum([int (i) for i in row.num_Vars_by_PI.split(',')]), axis=1)

    ''' write this concise statistics to a file. '''
    pd.DataFrame.to_csv(hap_stats, outputdir + '/' + prefix+'_haplotype_stats_' + soi + '.txt',
                        sep='\t', header=True, index=False)


    ## ** for future - plot haplotype STATs - may need improvement.
    ''' Step 09 - B : Now, use the above dataframe to plot several statistical plots from here on.
        - plot the haplotype statistics using matplotlib, pylab or ggplot (python).
        - start with "with open() .... " to open a file, write the plot and then close. '''
    # some good links for hist-plots: http://pandas.pydata.org/pandas-docs/version/0.18/visualization.html

    # plots total number of variants per-chromosome
    with open(outputdir + '/' + 'total_vars_'+ soi +'_'+ prefix+'.png', 'wb') as fig_initial:
        hap_stats.plot(x='CHROM', y='total_Vars', kind='bar')
        plt.xlabel('chromosomes')
        plt.ylabel('number of variants')
        plt.suptitle('number of variants for each chromosome')
        plt.savefig(fig_initial)
        # plt.show()  # optional: to show the plot to the console

    # plots total number of haplotypes per-chromosome
    with open(outputdir + '/' + 'total_haps_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:
        hap_stats.plot(x='CHROM', y='total_haplotypes', kind='bar')
        plt.xlabel('chromosomes')
        plt.ylabel('number of haplotypes')
        plt.suptitle('number of haplotypes for each chromosome')
        plt.savefig(fig_initial)


    # plot histogram of haplotype size (by number of variants) per-chromosome
    # make different plots for each chromosome, but X-axis is shared.
    with open(outputdir + '/' + 'hap_size_byVar_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:

        ## making histogram plot.
        # raising if-else when number of chromosome is 1 vs. more than 1. ** This can be optimize in the future.
        fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)

        # when we have only one chromosome
        if len(hap_stats) == 1:
            # fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True, squeeze=False)
            data_i = hap_stats['num_Vars_by_PI'].str.split(',')
            data_i = pd.Series([int(x) for x in data_i[0]])
            data_i.plot(kind='hist', label=str(hap_stats['CHROM']), alpha=0.5)
            plt.ylabel('frequency of the haplotypes')

        elif len(hap_stats) > 1:
            for i, data in hap_stats.iterrows():
                # first convert data to list of integers
                data_i = [int(x) for x in data['num_Vars_by_PI'].split(',')]
                # print(data)
                # print(data_i)
                ax[i].hist(data_i, label=str(data['CHROM']), alpha=0.5)
                ax[i].legend()

            fig.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')

        # add a common x-label
        plt.xlabel('size of the haplotype (by number of variants)')
        # plt.ylabel('frequency of the haplotypes')
        plt.suptitle(prefix + ' histogram of size of the haplotype (by number of variants) \n'
                     'for each chromosome')
        plt.savefig(fig_initial)



    # plot histogram of haplotype size (by genomic ranges of haplotype) per-chromosome
    # make different plots for each chromosome, but X-axis is shared.
    with open(outputdir + '/' + 'hap_size_byGenomicRange_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:

        ## making histogram plot.
        # raising if-else when number of chromosome is 1 vs. more than 1. ** This can be optimize in the future.
        fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)

        # when we have only one chromosome
        if len(hap_stats) == 1:
            # fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True, squeeze=False)
            data_i = hap_stats['range_of_PI'].str.split(',')
            data_i = pd.Series([int(x) for x in data_i[0]])
            data_i.plot(kind='hist', label=str(hap_stats['CHROM']), alpha=0.5)
            plt.ylabel('frequency of the haplotypes')

        # when the dataframe has more than one chromosome
        elif len(hap_stats) > 1:
            for i, data in hap_stats.iterrows():
                # first convert data to list of integers
                data_i = [int(x) for x in data['range_of_PI'].split(',')]
                # print(data)
                # print(data_i)
                ax[i].hist(data_i, label=str(data['CHROM']), alpha=0.5)
                ax[i].legend()

            # adding y-label separately if there are multiple chromosomes
            fig.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')

        # writing a common x-label
        plt.xlabel('size of the haplotype (by Genomic Distance)')
        #plt.ylabel('frequency of the haplotypes')
        plt.suptitle(prefix + ' histogram of size of the haplotype (by Genomic Distance)\n'
                     'for each chromosome')
        plt.savefig(fig_initial)


    #print('Global maximum memory usage: %.2f (mb)' % current_mem_usage())

    # clear memory
    del hap_data, hap_stats


