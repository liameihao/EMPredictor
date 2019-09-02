from subprocess import *
import pandas as pd
import sys


def step_modifier(num_feat, markers):
    fir = markers[0]
    sec = markers[1]
    step_list = []
    i = 100
    while i <= num_feat:
        # if i <= 500:
        step_list.append(i)
        i += 100
        # elif i > 500 and i <= 1000:
        #     step_list.append(i)
        #     i += 200
        # elif i > 1000:
        #     step_list.append(i)
        #     i += 250
    return step_list


def svm_train(c, g, v, dname):
    """
    Train the model using SVM via libsvm library
    Inputs: SVM parameters c, g, v, and the dataset file
    Output: generate SVM results
    """
    if sys.platform != 'win32':
        # grid_py = '../tools/grid.py'
        # gnuplot_exe = '/usr/bin/gnuplot'
        svmtrain_exe = './PromoterSVM/libsvm/windows/my-svm-train'
    else:
        # grid_py = r"D:\LiYuan\Data_Science\libsvm-3.21\windows\grid.py"
        # gnuplot_exe = r"D:\Program Files\gnuplot\bin\gnuplot.exe"
        svmtrain_exe = r'D:\LiYuan\"Bioinformatic Research"\IMU\zuo\Tasks\"20160206-Promoter SVM"\biopromoter_script\PromoterSVM\libsvm\windows\my-svm-train.exe'

    df_model = dname + ".model"
    cmd = '{0} -c {1} -g {2} -v {3} "{4}" "{5}"'.format(svmtrain_exe, c, g, v, dname, df_model)
    Popen(cmd, shell=True, stdout=PIPE).communicate()


def gen_final_df(res_fname, dat_sel, dname):
    """
    Generate reduced (feature selected) dataset
    Inputs: SVM result files, selected features, and output name
    Output: Desired final output
    """
    res_file = open(res_fname, "r")
    rdict = {}
    pdict = {}
    for p in res_file:
        i, c, g, rate = p.split(",")
        rdict[i] = rate
        pdict[i] = [c, g]
    res_file.close()

    opti_feat = min([key for key in rdict.keys() if rdict[key] == max(rdict.values())])
    opti_c, opti_g = pdict[opti_feat][0], pdict[opti_feat][1]
    opti_df = df_transformer(dat_sel, opti_feat)
    opti_df = opti_df.replace(to_replace=-1, value=2)
    opti_df.to_csv(dname+".res", sep=" ", header=None, index=None)
    return dname+".res", opti_c, opti_g


def para_search(dname):
    """
    Find optimal parameter(c,g) for svm via grid.py and calculate the CV rate on data (dname)
    :param dname: dataframe modified (scaled, libsvm format) and ready for process (grid.py, svm-trian)
    :return: c, g, cv rate
    """
    if sys.platform != 'win32':
        grid_py = './PromoterSVM/libsvm/tools/grid.py'
        # gnuplot_exe = '/usr/bin/gnuplot'
        # svmtrain_exe = '../my-svm-train'
    else:
        grid_py = r'D:\LiYuan\"Bioinformatic Research"\IMU\zuo\Tasks\"20160206-Promoter SVM"\biopromoter_script\PromoterSVM\libsvm\tools\grid.py'
        # gnuplot_exe = r"D:\Program Files\gnuplot\bin\gnuplot.exe"
        # svmtrain_exe = r'D:\LiYuan\"Bioinformatic Research"\IMU\zuo\Tasks\"20160206-Promoter SVM"\biopromoter_script\PromoterSVM\libsvm\windows\my-svm-train.exe'

    ###################################
    # grid.py: find the best parameter(c,g), and generates a model file
    ###################################
    # cg_results = []
    cmd = "{0} {1}".format(grid_py, dname)
    print "Cross Validating..."
    grid_out = Popen(cmd, shell=True, stdout=PIPE).stdout

    print cmd
    print grid_out.readline()

    line = ""
    while True:
        last_line = line
        line = grid_out.readline()
        if not line:
            break
    c, g, cvrate = map(float, last_line.split())
    # cg_results.append('Best c={0}, g={1} CV rate={2}'.format(c,g,rate))

    print('Best c={0}, g={1} CV rate={2}'.format(c, g, cvrate))

    return c, g, cvrate


def ifs_trainer(dname, markers=(500, 1000)):
    """
    Incremental feature selection (IFS), adds n features one at time, finds the best svm parameters (c,g) using grid.py,
    train the input dat via svm algorithm(svm-train) to obtain accuracy
    :param dname: output datafile from fselect (non-scaled, libsvm format) and ready for process (grid.py, svm-train)
    :param markers: step sizes of feature increment
    :return:
    """

    # import data
    dat = pd.read_csv(dname, delimiter=" ", header=None)
    num_feat = len(dat.columns) - 1
    # labels = dname
    print "number of features: ", num_feat
    # print dat
    # print num_feat
    # if single is True:
    #     step_list = [num_feat]
    #     labels = None
    # else:
    step_list = step_modifier(num_feat, markers)
    para_list = []
    # i = 1
    for i in step_list:
        print i
        # transform data to specified number of features
        dat_tf = df_transformer(dat, i, dname=dname+str(i)+".tf")

        # scale data
        scale_dat = svm_scale(dat_tf, -1, 1)

        # find svm parameter using grid.py
        c, g, cvrate = para_search(scale_dat)

        para_list.append([i, c, g, cvrate])

        # i += step

    res_file = open(dname+".cvr", "w")
    # print >> res_file, "n, c, g, CV rate"
    for i, c, g, rate in para_list:
        print >> res_file, "{},{},{},{}".format(i, c, g, rate)
    res_file.close()


def svm_scale(df_libsvm_fname, lower, upper):
    """
    :param df_libsvm_fname: input scale required filename
    :param lower: lower bound
    :param upper: upper bound
    :return: output scaled filename
    """
    if sys.platform != 'win32':
        scale_exe = './PromoterSVM/libsvm/windows/svm-scale'
    else:
        scale_exe = r'D:\LiYuan\"Bioinformatic Research"\IMU\zuo\Tasks\"20160206-Promoter SVM"\biopromoter_script\PromoterSVM\libsvm\windows\svm-scale'
    # range_file = df_libsvm_fname + ".range"
    scale_file = df_libsvm_fname + ".scale"
    print scale_file
    # cmd = '{0} -l {1} -u {2} -s "{3}" "{4}" > "{5}"'\
    #     .format(scale_exe, lower, upper, range_file, df_libsvm_fname, scale_file)
    cmd = '{0} -l {1} -u {2} "{3}" > "{4}"'\
        .format(scale_exe, lower, upper, df_libsvm_fname, scale_file)
    print cmd
    Popen(cmd, shell=True, stdout=PIPE).communicate()
    return scale_file


def df_transformer(dat_sel, num_of_feat_req, dname=None):
    if type(dat_sel) == str:
        dat_sel = pd.read_csv(dat_sel, header=None, delimiter=" ")
        # dat_sel = open(dat_sel, "r")
    dat_tf = dat_sel.ix[:, range(int(num_of_feat_req)+1)]
    if dname:
        # if ".tf" in dname:
        #     dname = dname
        #     dat_tf.to_csv(dname, sep=" ", header=None, index=None)
        # else:
        #     dname = dname+".tf"
        #     dat_tf.to_csv(dname, sep=" ", header=None, index=None)
        dat_tf.to_csv(dname, sep=" ", header=None, index=None)
        return dname
    else:
        return dat_tf


def tolibsvm(df_in, fname_out, gen_file=1):
    """
    Convert a pandas dataframe to libsvm format, given the first column of the dataframe is the label
    :param df_in: pandas dataframe or a csv file with header
    :return: df file in libsvm format
    """
    if type(df_in) is str:
        df_in = pd.read_csv(df_in)
    df_libsvm = []
    for ind, row in df_in.iterrows():
        row_res = []
        col_ind = 0
        for i in row:
            if col_ind == 0:
                ele = str(i)
                row_res.append(ele)
                col_ind += 1
                continue
            else:
                ele = "{}:{}".format(col_ind, i)
                row_res.append(ele)
                col_ind += 1
        df_libsvm.append(" ".join(row_res))
    if gen_file == 2:
        if ".libsvm" in fname_out:
            fname = fname_out
        else:
            fname = str(fname_out)+".libsvm"
        file_out = open(fname, "w")
        for i in df_libsvm:
            print >>file_out, i
        file_out.close()
        return df_libsvm, fname
    elif gen_file == 0:
        return df_libsvm
    elif gen_file == 1:
        if ".libsvm" in fname_out:
            fname = fname_out
        else:
            fname = str(fname_out)+".libsvm"
        file_out = open(fname, "w")
        for i in df_libsvm:
            print >>file_out, i
        file_out.close()
        return fname

