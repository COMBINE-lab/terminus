import pandas as pd
import numpy as np
from collections import defaultdict

def relDiff(c1, c2, DF, cutoff=0.01, verbose=False):
    import pandas as pd
    """
    Computes the relative difference between the values
    in columns c1 and c2 of DF.
    c1 and c2 are column names and DF is a Pandas data frame.
    Values less than cutoff will be treated as 0.

    The relative difference is defined as

    d(x_i, y_i) =
        0.0 if x_i and y_i are 0
        (x_i - y_i) / (|x_i + y_i|) otherwise

    This function returns two values.

    rd is a DataFrame where the "relDiff" column contains all
    of the computed relative differences.

    nz is a set of indexes into rd for data points where
    x_i and y_i were not *both* zero.
    """
    import numpy as np
    rd = pd.DataFrame(data = {"Name" : DF.index, "relDiff" : np.zeros(len(DF.index))*np.nan})
    rd.set_index("Name", inplace=True)
    bothZero = DF.loc[(DF[c1] < cutoff) & (DF[c2] < cutoff)].index
    nonZero = DF.index.difference(bothZero)
    if (verbose):
        print("Zero occurs in both columns {} times".format(len(rd.loc[bothZero])))
        print("Nonzero occurs in at least 1 column {} times".format(len(rd.loc[nonZero])))
    allRD = 2.0 * ((DF[c1] - DF[c2]) / (DF[c1] + DF[c2]).abs())
    assert(len(rd.loc[nonZero]["relDiff"]) == len(allRD[nonZero]))
    rd["relDiff"][nonZero] = allRD[nonZero]
    if len(bothZero) > 0:
        rd["relDiff"][bothZero] = 0.0
    return rd, nonZero

def group_map(gf, tnames, dup_dict = None):
    group_map = defaultdict(list)
    groups = {}
    if dup_dict is None:
        with open(gf) as fp:
            for line in fp:
                names = line.strip().split(',')
                for n in names[1:]:
                    groups[n] = names[0]
                    group_map[names[0]].append(n)
    else:
        with open(gf) as fp:
            for line in fp:
                names = line.strip().split(',')
                for n in names[1:]:
                    if n in dup_dict:
                        groups[dup_dict[n]] = names[0]
                    groups[n] = names[0]
    for n in tnames:
        if not n in groups:
            groups[n] = n
    return groups,group_map

def get_swim(
    mercury_dir
):

    # Truth + Salmon
    samples = [1,2,3,4]
    conditions = [1,2]
    quant_df = {}
    quant_mercury_df = {}
    for s in samples :
        df = {}
        mdf = {}
        vdf = {}
        for c in conditions:
            df[c] = pd.read_csv(
                '/mnt/scratch1/hirak/ASE_Analysis/simulation/quant_witohut_decoy/out_{}_sample_0{}/quant.sf'.format(s,c),
                 sep = '\t'
             )
            mdf[c] =  pd.read_csv(
                '{}/out_{}_sample_0{}/quant.sf'.format(mercury_dir,s,c),
                 sep = '\t'
             )
        quant_df[s] = df.copy()
        quant_mercury_df[s] = mdf.copy()

    true_df = {}
    for s in samples:
        column_name = ['sample_01', 'sample_02']
        file_name = '/mnt/scratch1/hirak/ASE_Analysis/simulation/fastq/out_{}/sim_counts_matrix.csv'.format(s)
        true_df_tmp = pd.read_csv(file_name)
        for c in conditions:
            column_name += ['salmon_'+str(c)]
            true_df_tmp = true_df_tmp.join(quant_df[s][c].set_index('Name')[['NumReads']], how='outer', rsuffix='_'+str(c))
        true_df[s] = true_df_tmp.fillna(0.0).copy()
        true_df[s].columns = column_name

    # True group
    true_group_df = {}
    gr,gmap = group_map(
        '{}/out_1_sample_01/clusters.txt'.format(mercury_dir)
        , quant_df[1][1]['Name'].values
    )

    for s in samples:
        true_group_tmp = {}
        for c in conditions:
            # tg_temp = true_df[s]['sample_0'+str(c)].join
            tg = pd.DataFrame(true_df[s][['sample_0'+str(c)]].groupby(gr)['sample_0'+str(c)].sum())
            true_group_tmp[c] = tg.join(quant_mercury_df[s][c].set_index('Name')[['NumReads']])
            true_group_tmp[c].columns = ['sample_0'+str(c), 'mercury_'+str(c)]

        true_group_df[s] = true_group_tmp

    return {
        'salmon_with_truth': true_df,
        'mercury_with_truth': true_group_df,
        'groups':gr,
        'gmap':gmap
    }


def get_swim_mmcollapse(
    true_df
):
    samples = [1,2,3,4]
    conditions = [1,2]

    mmcollapse_df = pd.read_csv(
                '/mnt/scratch1/hirak/ASE_Analysis/simulation/mmseq_quant/out_1_sample_01.collapsed.mmseq',
                skiprows = 1,
                sep = '\t'
            )
    mmseq_df = {}
    for s in samples:
        mmcollapse_df_tmp = {}
        mmseq_df_tmp = {}
        for c in conditions:
            with open(
                '/mnt/scratch1/hirak/ASE_Analysis/simulation/mmseq_quant/out_{}_sample_0{}.mmseq'.format(
                    s,
                    c
                )
            ) as fp:
                numreads = int(fp.readline().strip().split(' ')[-1])
            d = pd.read_csv(
                '/mnt/scratch1/hirak/ASE_Analysis/simulation/mmseq_quant/out_{}_sample_0{}.mmseq'.format(
                    s,
                    c
                ),
                skiprows = 1,
                sep = '\t'
            )

            d['mmseq_numreads'] = np.exp(d.log_mu_em)
            d.mmseq_numreads = d.mmseq_numreads * d.effective_length
            d.mmseq_numreads = d.mmseq_numreads / 1e9
            d.mmseq_numreads = d.mmseq_numreads * numreads
            mmseq_df_tmp[c] = d.copy()
        mmseq_df[s] = mmseq_df_tmp

    import re
    mmcollapse_groups = {}
    group_names = []
    grouped_transcripts_mmcollapse = []
    for t in mmcollapse_df.feature_id.values:
        if ('+' in t) or ('*' in t):
            cands = re.split('\+|\*', t)
            for c in cands:
                mmcollapse_groups[c] = t
                grouped_transcripts_mmcollapse += [c]
            group_names += [t]
        else:
            mmcollapse_groups[t] = t

    true_mmcollapse_df = {}
    for s in samples:
        true_mmcollapse_tmp = {}
        for c in conditions:
            tg = pd.DataFrame(true_df[s][['sample_0'+str(c)]].groupby(mmcollapse_groups)['sample_0'+str(c)].sum())
            mmg = pd.DataFrame(mmseq_df[s][c].set_index('feature_id')[['mmseq_numreads']].groupby(
                mmcollapse_groups)['mmseq_numreads'].sum())
            true_mmcollapse_tmp[c] = tg.join(mmg).fillna(0)
            true_mmcollapse_tmp[c].columns = ['sample_0'+str(c), 'mmcollapse_'+str(c)]
        true_mmcollapse_df[s] = true_mmcollapse_tmp

    return {
        'mmcollapse_df' : true_mmcollapse_df,
        'groups' : mmcollapse_groups
    }


def calculate_spearman_mard(
    true_df,
    true_group_df,
    mmcollapse_df
):
    samples = [1,2,3,4]
    conditions = [1,2]
    results = []
    for s in samples:
        results += [["salmon",s,1,
                  true_df[s][['sample_01','salmon_1']].corr(method='spearman').values[0][1],
                  relDiff('sample_01','salmon_1',true_df[s])[0]['relDiff'].abs().mean()
                       ]]
        results += [["salmon",s,2 ,
              true_df[s][['sample_02','salmon_1']].corr(method='spearman').values[0][1],
              relDiff('sample_02','salmon_2',true_df[s])[0]['relDiff'].abs().mean()
                   ]]

        for c in conditions:
            results += [["mercury",s, c,
                  true_group_df[s][c][['sample_0'+str(c),'mercury_'+str(c)]].corr(method='spearman').values[0][1],
                  relDiff('sample_0'+str(c),'mercury_'+str(c),true_group_df[s][c])[0]['relDiff'].abs().mean()
                       ]]
            results += [["mmcollapse",s, c,
                  mmcollapse_df[s][c][['sample_0'+str(c),'mmcollapse_'+str(c)]].corr(method='spearman').values[0][1],
                  relDiff('sample_0'+str(c),'mmcollapse_'+str(c),mmcollapse_df[s][c])[0]['relDiff'].abs().mean()
                       ]]
    result_df = pd.DataFrame(results)
    result_df.columns = ['method','sample','condition','correlation','mard']
    return {
        'acc_df':result_df
    }


def get_joined_df_mouse_allele():
    true_df = pd.read_csv('/mnt/scratch1/hirak/ASE_Analysis/simulation/fastq/mouse_diploid/sim_counts.tsv',
           header = None, sep = '\t', names = ['transcript','count'])

    with open(
        '/mnt/scratch1/hirak/ASE_Analysis/simulation/mmseq_quant_mouse_feynman/sample_01.mmseq'
    ) as fp:
        numreads = int(fp.readline().strip().split(' ')[-1])
    d = pd.read_csv(
        '/mnt/scratch1/hirak/ASE_Analysis/simulation/mmseq_quant_mouse_feynman/sample_01.mmseq',
        skiprows = 1,
        sep = '\t'
    )

    d['mmseq_numreads'] = np.exp(d.log_mu_em)
    d.mmseq_numreads = d.mmseq_numreads * d.effective_length
    d.mmseq_numreads = d.mmseq_numreads / 1e9
    d.mmseq_numreads = d.mmseq_numreads * numreads

    mmseq_df = d.copy()

    salmon_df = pd.read_csv(
                '/mnt/scratch1/hirak/ASE_Analysis/simulation/quant_witohut_decoy/mouse_emase/quant.sf',
                 sep = '\t'
             )

    joint_df = true_df.set_index('transcript').join(salmon_df.set_index('Name')[['NumReads']], how = 'outer').join(
        mmseq_df.set_index('feature_id')['mmseq_numreads']
    ).fillna(0)


    # terminus
    gr,gmap = group_map(
        '/mnt/scratch1/hirak/ASE_Analysis/simulation/terminus_result_without_decoy_mouse/mouse_emase/clusters.txt'
        ,joint_df.index.values
    )

    terminus_tg = pd.DataFrame(joint_df[['count']].groupby(gr)['count'].sum())
    terminus_df = terminus_tg.join(pd.read_csv(
        '/mnt/scratch1/hirak/ASE_Analysis/simulation/terminus_result_without_decoy_mouse/mouse_emase/quant.sf',
        sep = '\t'
    ).set_index('Name')['NumReads'])

    # mmcollapse
    mmcollapse_df = pd.read_csv(
            '/mnt/scratch1/hirak/ASE_Analysis/simulation/mmseq_quant_mouse_feynman/sample_01.collapsed.mmseq',
            skiprows = 1,
            sep = '\t'
    )
    import re
    mmcollapse_groups = {}
    group_names = []
    grouped_transcripts_mmcollapse = []
    for t in joint_df.index.values:
        mmcollapse_groups[t] = t
    for t in mmcollapse_df.feature_id.values:
        if ('+' in t) or ('*' in t):
            cands = re.split('\+|\*', t)
            for c in cands:
                mmcollapse_groups[c] = t
                grouped_transcripts_mmcollapse += [c]
            group_names += [t]
        else:
            mmcollapse_groups[t] = t

    mmcollapse_df = joint_df[['count','mmseq_numreads']].groupby(mmcollapse_groups).sum()

    return(joint_df, terminus_df, mmcollapse_df)


def get_spearman_mard(t,m):
    print('Correlation')
    print('Terminus', t.corr('spearman')['count']['NumReads'])
    print('mmcollapse', m.corr('spearman')['count']['mmseq_numreads'])

    print('MARD')
    print('Terminus', relDiff('count','NumReads',t)[0]['relDiff'].abs().mean())
    print('mmcollapse', relDiff('count','mmseq_numreads',m)[0]['relDiff'].abs().mean())


def write_dataframes(
    true_df,
    true_group_df,
    mmcollapse_df,
    j,
    t,
    m,
):
    from datetime import datetime
    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
    #print("date and time =", dt_string)

    samples = [1,2,3,4]
    conditions = [1,2]

    for s in samples:
        true_df[s].to_csv(
            '/mnt/scratch1/hirak/ASE_Analysis/simulation/dataframes/salmon_truth_{}_{}'.format(s, dt_string),
        )
        for c in conditions:
            true_group_df[s][c].to_csv(
                '/mnt/scratch1/hirak/ASE_Analysis/simulation/dataframes/merc_truth_{}_{}_{}'.format(s,c,dt_string)
            )
            mmcollapse_df[s][c].to_csv(
                '/mnt/scratch1/hirak/ASE_Analysis/simulation/dataframes/mmcollapse_truth_{}_{}_{}'.format(s,c,dt_string)
            )
    j.to_csv('/mnt/scratch1/hirak/ASE_Analysis/simulation/dataframes/allelic_joint_{}'.format(dt_string))
    t.to_csv('/mnt/scratch1/hirak/ASE_Analysis/simulation/dataframes/allelic_terminus_{}'.format(dt_string))
    m.to_csv('/mnt/scratch1/hirak/ASE_Analysis/simulation/dataframes/allelic_mmseq_{}'.format(dt_string))


def main():
    print('Parsing swimming downstream data')
    inf_res = get_swim('/mnt/scratch1/hirak/ASE_Analysis/simulation/terminus_result_without_decoy/')
    res_mmcollapse_swim = get_swim_mmcollapse(inf_res['salmon_with_truth'])
    inf_acc = calculate_spearman_mard(inf_res['salmon_with_truth'],inf_res['mercury_with_truth'],
                                 res_mmcollapse_swim['mmcollapse_df'])

    print(inf_acc['acc_df'].groupby(['method'])['correlation','mard'].mean())

    print('Parsing mouse allelic data')
    j,t,m = get_joined_df_mouse_allele()
    print(j.corr('spearman'))
    get_spearman_mard(t,m)

    # write_dataframes(inf_res['salmon_with_truth'], inf_res['mercury_with_truth'],res_mmcollapse_swim['mmcollapse_df'],j,t,m)


if __name__ == '__main__':
    main()
