from feat_select import *
from util import *


def feature_select(dname, switch=True):

    # dat_fs = pd.read_csv("sel_df_tf.csv", delimiter=" ", header=None)
    # print dat_fs

    dat = pd.read_csv(dname, delimiter=",")
    print "data imported..."

    dat_libsvm = tolibsvm(dat, dname, gen_file=2)
    # print "data converted to libsvm format..."

    dat_fs = fselect_optimizer(dat, dname, switch)
    print "fselect done, new feature selected data created..."
    return dat_fs


def svm(dat_fs, dname):
    # if type(dat_fs) == str:
    #     dat_fs = pd.read_table(dat_fs, delimiter=" ", header=None)
    print "initiating IFS..."
    ifs_trainer(dat_fs, markers=(500, 1000))

    dfname, c, g = gen_final_df(dname+".libsvm.cvr", dat_fs, dname)
    df_final = svm_scale(dfname, -1, 1)

    # if sys.platform != "win32":
    svm_train(c, g, 5, df_final)

    # print "plotting CV rates..."
    # plot_cvrate(dname+".libsvm.cvr")


def post_fs_svm(svmr_name, finfscore_name, num_fts):
    finfs_df = pd.read_csv(finfscore_name)
    red_ind = finfs_df.iloc[:num_fts, 0]

    # print red_ind

    print "transforming dataset..."
    red_df, red_ft_name = kmer_selector(svmr_name, red_ind)
    red_name = svmr_name + '_' + str(num_fts)
    # print red_df
    red_df.to_csv(red_name, index=None)

    print 'converting to libsvm format...'
    libsvm_name = tolibsvm(red_name, red_name)

    print "df scaling..."
    df_scaled = svm_scale(libsvm_name, -1, 1)

    print "grid searching..."
    c, g, cvrate = para_search(df_scaled)

    print 'svming...'
    svm_train(c, g, 5, df_scaled)


if __name__ == "__main__":
    # dname = "bp_sp_df9merfs.csv" # kmer file
    # fs_dat = feature_select(dname, switch=False)

    #
    fs = "bp_sp_df7mer.csv.finfscore"
    dat = "bp_sp_df7mer.csv.libsvm"
    threshold = range(3600, 3890, 50)

