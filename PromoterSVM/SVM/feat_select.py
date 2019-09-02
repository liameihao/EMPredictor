from util import *
import pandas as pd
import sys


def fselect(input_fname):
    """
    input filename, output .select file, .fscore file
    """
    if sys.platform != 'win32':
        fselect_exe = './PromoterSVM/libsvm/tools/fselect.py'
    else:
        fselect_exe = r'D:\LiYuan\"Bioinformatic Research"\IMU\zuo\Tasks\"20160206-Promoter SVM"\biopromoter_script\PromoterSVM\libsvm\tools\fselect.py'

    #dat_libsvm = "dat_libsvm.csv"
    cmd = "%s %s" % (fselect_exe, input_fname)
    print cmd
    Popen(cmd, shell=True, stdout=PIPE).communicate()


def feature_postprocess(fselect_res_fname, tar_fts):
    """
    input filename, output max_acc, num of optimal feats, feat list and respective fscores as a csv file
    :param fselect_res_fname:
    :return:
    """
    select_file = open(fselect_res_fname+".select", "r")
    feat_acc_dict = {}
    select_list = [line.split() for line in select_file]
    for i in range(len(select_list)):
        if "%#Feat" in select_list[i]:
            for j in range(i+1, len(select_list)):
                if "max" in select_list[j]:
                    break
                feat_acc_dict[int(select_list[j][0].rstrip(":"))] = float(select_list[j][1])
    select_file.close()
    max_acc = max(feat_acc_dict.values())
    opti_feat = min([key for key in feat_acc_dict.keys() if feat_acc_dict[key] == max(feat_acc_dict.values())])
    max_feat = max(feat_acc_dict.keys())

    fscore_df = pd.read_csv(fselect_res_fname+".fscore", delimiter=r": \t", header=None, engine="python")
    fs_ind_full = fscore_df.ix[:, 0].tolist()
    fs_ind_opti = fs_ind_full[:opti_feat]
    fs_ind_tar = fs_ind_full[:tar_fts]
    fs_score_full = fscore_df.ix[:, 1].tolist()
    fs_score_opti = fs_score_full[:opti_feat]
    fs_score_tar = fs_score_full[:tar_fts]

    return feat_acc_dict, max_feat, opti_feat, max_acc, fs_ind_opti, fs_score_opti, fs_ind_tar, fs_score_tar, fs_ind_full, fs_score_full


def gen_fscore_doc(sel_ind, sel_kmer, sel_fscores, dname, l=None):
    # if l:
    #     file_out = open("final_fscore"+str(l), "w")
    # else:
    #     file_out = open("final_fscore", "w")
    file_out = open(dname+".finfscore", "w")
    for ind, kmer, fs in zip(sel_ind, sel_kmer, sel_fscores):
        print >>file_out, "{}, {}, {}".format(ind, kmer, fs)
    file_out.close()


def kmer_selector(df_svm_header, sel_ind):
    """
    transforms the input df into a new df with the optimal features and rearranges the feature columns in
    decreasing fscore
    :param df_svm_header: pandas dataframe/data file with kmers as the header
    :param sel_ind: list of selected kmer indices
    :return: transformed dataframe with the new set of kmer columns
    """
    if type(df_svm_header) == str:
        df_svm_header = pd.read_csv(df_svm_header)
    kmer_list = df_svm_header.columns.values.tolist()[1:]
    sel_kmer = []
    for ind in sel_ind:
        sel_kmer.append(kmer_list[ind-1])
    sel_df = pd.DataFrame()
    sel_df["Label"] = df_svm_header["Label"]
    # sel_df["Label"] = sel_df["Label"].replace(to_replace=0, value=?)
    for kmer in sel_kmer:
        sel_df[kmer] = df_svm_header[kmer]

    return sel_df, sel_kmer


def fselect_main(svmr_fname):
    # if type(svmr_fname) is not str:
    inp = tolibsvm(svmr_fname, svmr_fname)
    print "data converted to libsvm format..."

    df_scaled = svm_scale(inp, -1, 1)
    print "df scaled..."

    fselect(df_scaled)
    print "fselect done..."


def gen_finfscore(svmr_name):
    raw_df = pd.read_csv(svmr_name)
    fs_name_full_unordered = raw_df.columns.values.tolist()[1:]
    fscore_df = pd.read_csv(svmr_name+".libsvm.scale.fscore", delimiter=r": \t", header=None, engine="python")
    fs_ind_full = fscore_df.ix[:, 0].tolist()
    fs_score_full = fscore_df.ix[:, 1].tolist()
    fs_name_full = []

    for i in fs_ind_full:
        fs_name_full.append(fs_name_full_unordered[i-1])

    finfscores_df = pd.DataFrame({"Features": fs_name_full, "Feature_Index": fs_ind_full, "Fscores": fs_score_full})
    finfscores_df.to_csv(svmr_name + '.finfscore', index=None)



def fselect_optimizer(inp, pre_dfname, tar_fts):
    """
    :param inp: input df filename (libsvm format, unscaled, no header), pre_df filename (header)
    :param num_of_loop:
    :return: fselect process data in libsvm format, unscaled
    """
    if type(inp) is not str:
        inp = tolibsvm(inp, pre_dfname)
        print "data converted to libsvm format..."


    df_scaled = svm_scale(inp, -1, 1)
    print "df scaled..."
    fselect(df_scaled)
    print "fselect done..."


    feat_dict, max_feat, opti_feat, acc, sel_inds, sel_scores, tar_inds, tar_scores, ful_inds, ful_scores = \
        feature_postprocess(df_scaled, tar_fts)
    print "optimal parameters selected...", max_feat, opti_feat, acc
    # print feat_dict, max_feat, opti_feat, acc, sel_inds, sel_scores
    df_tar, tar_kmers = kmer_selector(pre_dfname, tar_inds)
    print "dataframe transformed..."
    df_sel_svm, df_sel_svm_file = tolibsvm(df_tar, pre_dfname, gen_file=2)

    gen_fscore_doc(tar_inds, tar_kmers, tar_scores, pre_dfname)

    return df_sel_svm_file


