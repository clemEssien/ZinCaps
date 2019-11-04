import sys
import pandas as pd
from keras.models import load_model, clone_model
import numpy as np
import argparse
import csv
from DProcess import convertRawToXY
from EXtractfragment_sort import extractFragforPredict
from capsulenet import Capsnet_main
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import mean_squared_error

import os

os.environ["CUDA_VISIBLE_DEVICES"] = "0";
import tensorflow as tf

config = tf.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.8


# load model weights from file
def load_model_weights(inputweights, model_arch):
    all_weights = list()
    n_models = 0
    filename = inputweights + str(n_models)
    while os.path.exists(filename):
        # define filename for this ensemble
        model_arch[0].load_weights(filename)
        # add to list of members
        all_weights.append(model_arch[0].get_weights())
        print('>loaded %s' % filename)
        n_models += 1
        filename = inputweights + str(n_models)

    return all_weights


# evaluate a specific number of members in an ensemble
def polynomial_decay(decay_steps, init_value, end_value, power=1):
    weight_list = []
    weight = init_value
    for i in range(decay_steps):
        decayed_weight = (weight - end_value) * (1 - i / decay_steps) ** (power) + end_value
        weight = decayed_weight
        weight_list.append(decayed_weight)

    return weight_list


def predict_by_avg_members(members, model_arch, testX):
    n_models = len(members)
    print("length", n_models)
    increase_steps = int(n_models / 2)
    decrease_steps = n_models - increase_steps
    temp_weight = polynomial_decay(increase_steps, 1, 0.5, power=1)
    temp_weight.reverse()
    # weights = temp_weight+polynomial_decay(decrease_steps,1,0.5,power=1)
    # weights = polynomial_decay(n_models,1,0.5,power=1)
    # weights = [1.0 / n_models for i in range(1, n_models + 1)]
    weights = list(np.ones(increase_steps)) + list(np.ones(decrease_steps) * 0.5)
    # determine how many layers need to be averaged
    n_layers = len(members[0])

    # create an set of average model weights
    avg_model_weights = list()
    for layer in range(n_layers):
        # collect this layer from each model
        layer_weights = np.array([x[layer] for x in members])
        # weighted average of weights for this layer
        avg_layer_weights = np.average(layer_weights, axis=0, weights=weights)
        # store average layer weights
        avg_model_weights.append(avg_layer_weights)

    # set the weights in the new
    model_arch[0].set_weights(avg_model_weights)
    predict_proba = model_arch[1].predict(testX, verbose=0)[0]
    return predict_proba


def predict_by_snapshot(members, model_arch, testX):
    n_models = len(members)
    increase_steps = int(n_models / 2)
    decrease_steps = n_models - increase_steps
    temp_weight = polynomial_decay(increase_steps, 1, 0.5, power=1)
    temp_weight.reverse()
    # weights = temp_weight+polynomial_decay(decrease_steps,1,0.5,power=1)
    # weights = [1.0 / n_models for i in range(1, n_models + 1)]
    weights = list(np.ones(increase_steps)) + list(np.ones(decrease_steps) * 0.5)
    # weights = polynomial_decay(n_models,1,0.5,power=1)
    # weights = [1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    predict_list = []
    for member_weight in members:
        model_arch[0].set_weights(member_weight)
        predict_list.append(model_arch[1].predict(testX, verbose=0)[0])

    predict_list = np.asarray(predict_list)
    avg_predict_results = np.average(predict_list, axis=0, weights=weights)
    return avg_predict_results


def write_output(outputfile, predict_proba, ids, poses, focuses):
    poses = poses + 1  # start from 1
    results = np.column_stack((ids, poses, focuses, predict_proba[:, 1]))
    result = pd.DataFrame(results)
    result.to_csv(outputfile + ".txt", index=False, header=None, sep='\t',
                  quoting=csv.QUOTE_NONNUMERIC)


def evaluate(predict_proba, testY):
    from sklearn.metrics import roc_auc_score, average_precision_score
    true_label = [np.argmax(x) for x in testY]
    false_label = [np.argmin(x) for x in testY]

    np.savetxt("predict_proba.txt", predict_proba, delimiter=",")
    np.savetxt("true_label.txt", true_label, delimiter=",")
    np.savetxt("false_label.txt", false_label, delimiter=",")
    np.savetxt("testY.txt", testY, delimiter=",")

    roc_score = roc_auc_score(true_label, predict_proba[:, 1])
    pr_score = average_precision_score(true_label, predict_proba[:, 1])

    predict_proba[:, 1] = [(x >= 0.5).astype('int') for x in predict_proba[:, 1]]
    predict_proba[:, 0] = [(x >= 0.5).astype('int') for x in predict_proba[:, 0]]

    accuracy = accuracy_score(true_label, predict_proba[:, 1])
    sensitivity = recall_score(true_label, predict_proba[:, 1])
    specificity = recall_score(false_label, predict_proba[:, 0])
    f1 = f1_score(true_label, predict_proba[:, 1])
    mcc = matthews_corrcoef(true_label, predict_proba[:, 1])

    return roc_score, pr_score, accuracy, sensitivity, specificity, f1, mcc


def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument('-input',  dest='inputfile', type=str, help='Protein sequences to be predicted in fasta format.', required=True)
    parser.add_argument('-output', dest='outputfile', type=str, help='prefix of the prediction results.', required=True)
    parser.add_argument('-window', dest='window', type=int, help='specify the window size', required=True)
    parser.add_argument('-model-prefix', dest='modelprefix', type=str,
                        help='prefix of custom model used for prediciton. If you do not have one, please run train_models.py to train a model.',
                        required=False, default=None)
    parser.add_argument('-residue-types', dest='residues', type=str,
                        help='Residue types that to be predicted. For multiple residues, seperate each with \',\'',
                        required=False, default="C,H,E,D")
    parser.add_argument('-codingMode',  dest='codingMode', type=int, help='Set the input sequence encoding mode.', required=False, default=0)

    args = parser.parse_args()

    # inputfile=args.inputfile;
    outputfile = args.outputfile;
    residues = args.residues.split(",")
    modelprefix = args.modelprefix;
    window = args.window;
    codemode = args.codingMode
    print(outputfile, residues, modelprefix, window)
    # outputfile = r'/home/ucexbw/ZinCaps/ActiveSitePrediction/data/output/'
    # fp = open(outputfile+"eval_by_AUC_precision_scores_polynomial_decay_increase_decrease_1_0.5_1",'w')
    # fp = open(outputfile+"eval_by_AUC_precision_scores_polynomial_decay_1_0.5_1",'w')
    # fp = open(outputfile+"eval_by_AUC_precision_scores_10fold",'w')
    fp = open(outputfile + "eval_by_AUC_precision_scores_10fold_constantweight1_0.5_25.txt", 'w')

    model_arch = Capsnet_main(np.zeros([3, 2 * window + 1, 6]), [], nb_epoch=1, compiletimes=0, lr=0.001,
                              batch_size=500, lam_recon=0, routings=3, modeltype='nogradientstop', nb_classes=2,
                              predict=True)
    # model_arch=Capsnet_main(np.zeros([3,2*16+1,21]),[],nb_epoch=1,compiletimes=0,lr=0.001,batch_size=500,lam_recon=0,routings=3,modeltype='nogradientstop',nb_classes=2,predict=True)

    roc_average_weight = np.zeros(5)
    roc_average_predict = np.zeros(5)
    roc_average_last_predict = np.zeros(5)

    accuracy_average_last_predict = np.zeros(5)
    sensitivity_average_last_predict = np.zeros(5)
    specificity_average_last_predict = np.zeros(5)
    f1_score_average_last_predict = np.zeros(5)
    mcc_average_last_predict = np.zeros(5)

    pr_average_weight = np.zeros(5)
    pr_average_predict = np.zeros(5)
    pr_average_last_predict = np.zeros(5)

    for time in range(5):
        fp.write("############################" + str(time) + "\n")
        inputfile = '/scratch/ucexbw/ZinCaps25/ActiveSitePrediction/lib/K-Fold/annotated_sequence.fasta_training_annotated_' + str(
            time) + '.fasta'
        # if os.path.exists(outputfile+"eval_by_AUC_precision_scores"):
        #   os.rm(outputfile+"eval_by_AUC_precision_scores")

        checkpointweights = '/scratch/ucexbw/ZinCaps25/ActiveSitePrediction/data/weights/Zinc_' + str(time) + '_weights'
        modelprefix = '/scratch/ucexbw/ZinCaps25/ActiveSitePrediction/data/models/Zinc_' + str(time) + '_model'
        eval_type = 'average_last_predict'  # all evaluate by all method
        # average_weight
        # average_predict
        # average_last_predict

        if modelprefix is None:
            # print ("Please specify the prefix for an existing custom model by "
            #        "-model-prefix!\n\
            # It indicates two files [-model-prefix]_HDF5model and [-model-prefix]_parameters.\n \
            # If you don't have such files, please run train_models.py to get the "
            #        "custom model first!\n")
            exit()
        else:  # custom prediction
            model = modelprefix + str("_HDF5model")
            parameter = modelprefix + str("_parameters")
            try:
                f = open(parameter, 'r')
            except IOError:
                print('cannot open ' + parameter + " ! check if the model exists. "
                                                   "please run train_general.py or train_kinase.py to get the custom model first!\n")
            else:
                f = open(parameter, 'r')
                parameters = f.read()
                f.close()

            nclass = int(parameters.split("\t")[0])
            window = int(parameters.split("\t")[1])
            residues = parameters.split("\t")[2]
            residues = residues.split(",")
            codemode = int(parameters.split("\t")[4])
            modeltype = str(parameters.split("\t")[5])
            nb_classes = int(parameters.split("\t")[6])

        testfrag, ids, poses, focuses = extractFragforPredict(inputfile, window, '-', focus=residues)

        testX, testY = convertRawToXY(testfrag.as_matrix(), codingMode=codemode)
        if len(testX.shape) > 3:
            testX.shape = (testX.shape[0], testX.shape[2], testX.shape[3])

        predict_average_weight = np.zeros((testX.shape[0], 2))
        predict_average_predict = np.zeros((testX.shape[0], 2))
        predict_average_last_predict = np.zeros((testX.shape[0], 2))

        for bt in range(nclass):  # 0 648 bt=2 len(tf.trainable_variables())=1530
            # load all involving mode weights
            # sess = tf.Session()
            inputweights = checkpointweights + "_nclass" + str(bt) + "_iteration"
            model_members = load_model_weights(inputweights, model_arch)
            if eval_type == "all" or eval_type == "average_weight":
                predict_temp = predict_by_avg_members(model_members, model_arch, testX)
                predict_average_weight += predict_temp
                auc_score, pr_score, accuracy, sensitivity, specificity, f1_score, mcc = evaluate(predict_temp, testY)
                roc_average_last_predict[time] = auc_score
                pr_average_last_predict[time] = pr_score
                accuracy_average_last_predict[time]  = accuracy
                sensitivity_average_last_predict[time]  = sensitivity
                specificity_average_last_predict[time]  = specificity
                f1_score_average_last_predict[time]  = f1_score
                mcc_average_last_predict[time]  = mcc
                # fp.write(
                #     "average_weight_results_bt" + str(bt) + "\t" + str(auc_score) + "\t" + str(pr_score) + "\t" + str(
                #         accuracy) + "\t" + str(sensitivity) + "\t" + str(specificity) + "\t" + str(
                #         f1_score) + "\t" + str(mcc) + "\n")

            if eval_type == "all" or eval_type == "average_predict":
                predict_temp = predict_by_snapshot(model_members, model_arch, testX)
                predict_average_predict += predict_temp
                auc_score, pr_score, accuracy, sensitivity, specificity, f1_score, mcc = evaluate(predict_temp, testY)
                roc_average_last_predict[time] = auc_score
                pr_average_last_predict[time] = pr_score
                accuracy_average_last_predict[time]  = accuracy
                sensitivity_average_last_predict[time]  = sensitivity
                specificity_average_last_predict[time]  = specificity
                f1_score_average_last_predict[time]  = f1_score
                mcc_average_last_predict[time]  = mcc
                print("average_predict results:")
                # fp.write(
                #     "average_predict_results_bt" + str(bt) + "\t" + str(auc_score) + "\t" + str(pr_score) + "\t" + str(
                #         accuracy) + "\t" + str(sensitivity) + "\t" + str(specificity) + "\t" + str(
                #         f1_score) + "\t" + str(mcc) + "\n")

            del model_members
            # sess.close()

        if eval_type == "all" or eval_type == "average_weight1":
            predict_average_weight = predict_average_weight / float(nclass)
            auc_score, pr_score, accuracy, sensitivity, specificity, f1_score, mcc = evaluate(predict_average_weight,
                                                                                             testY)
            print("average_weight1")
            roc_average_last_predict[time] = auc_score
            pr_average_last_predict[time] = pr_score
            accuracy_average_last_predict[time]  = accuracy
            sensitivity_average_last_predict[time]  = sensitivity
            specificity_average_last_predict[time]  = specificity
            f1_score_average_last_predict[time]  = f1_score
            mcc_average_last_predict[time]  = mcc
            # fp.write(
            #     "average_weight_results\t" + str(auc_score) + "\t" + str(pr_score) + "\t" + str(accuracy) + "\t" + str(
            #         sensitivity) + "\t" + str(specificity) + "\t" + str(f1_score) + "\t" + str(mcc) + "\n")
            # roc_average_weight[time] = auc_score
            # pr_average_weight[time] = pr_score
            # write_output(outputfile + "average_weight_results_fold"+str(time)+".txt",predict_average_weight,ids,poses,focuses)

        if eval_type == "all" or eval_type == "average_predict":
            predict_average_predict = predict_average_predict / float(nclass)
            auc_score, pr_score, accuracy, sensitivity, specificity, f1_score, mcc = evaluate(predict_average_predict,
                                                                                              testY)
            roc_average_last_predict[time] = auc_score
            pr_average_last_predict[time] = pr_score
            accuracy_average_last_predict[time]  = accuracy
            sensitivity_average_last_predict[time]  = sensitivity
            specificity_average_last_predict[time]  = specificity
            f1_score_average_last_predict[time]  = f1_score
            mcc_average_last_predict[time]  = mcc
            # fp.write("average_predict_results:\t" + str(auc_score) + "\t" + str(pr_score) + "\t" + str(
            #     accuracy) + "\t" + str(sensitivity) + "\t" + str(specificity) + "\t" + str(f1_score) + "\t" + str(
            #     mcc) + "\n")
            # roc_average_predict[time] = auc_score
            # pr_average_predict[time] = pr_score
            # write_output(outputfile + "average_predict_results_fold"+str(time)+".txt",predict_average_predict,ids,poses,focuses)

        if eval_type == "all" or eval_type == "average_last_predict":
            nclass_ini = 1
            for bt in range(nclass):
                model_arch[0].load_weights(model + "_class" + str(bt))
                predict_temp = model_arch[1].predict(testX)[0]
                predict_average_last_predict += predict_temp
                auc_score, pr_score, accuracy, sensitivity, specificity, f1_score, mcc = evaluate(predict_temp, testY)
                # fp.write("average_last_predict_results_bt" + str(bt) + "\t" + str(auc_score) + "\t" + str(
                #     pr_score) + "\t" + str(accuracy) + "\t" + str(sensitivity) + "\t" + str(specificity) + "\t" + str(
                #     f1_score) + "\t" + str(mcc) + "\n")

            predict_average_last_predict = predict_average_last_predict / (nclass * nclass_ini)
            auc_score, pr_score, accuracy, sensitivity, specificity, f1_score, mcc = evaluate(
                predict_average_last_predict, testY)
            # fp.write("average_last_predict_results\t" + str(auc_score) + "\t" + str(pr_score) + "\t" + str(
            #     accuracy) + "\t" + str(sensitivity) + "\t" + str(specificity) + "\t" + str(f1_score) + "\t" + str(
            #     mcc) + "\n")
            roc_average_last_predict[time] = auc_score
            pr_average_last_predict[time] = pr_score
            accuracy_average_last_predict[time]  = accuracy
            sensitivity_average_last_predict[time]  = sensitivity
            specificity_average_last_predict[time]  = specificity
            f1_score_average_last_predict[time]  = f1_score
            mcc_average_last_predict[time]  = mcc
            # write_output(outputfile + "average_last_predict_results_fold"+str(time)+".txt",predict_average_last_predict,ids,poses,focuses)
            print("Successfully predicted from custom models !\n")

    fp.write("!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    # fp.write("average_weight_results\t" + ",".join([str(x) for x in roc_average_weight]) + "\t" + ",".join(
    #     [str(x) for x in pr_average_weight]) + "\t" + str(np.mean(roc_average_weight)) + "," + str(
    #     np.std(roc_average_weight)) + "\t" + str(np.mean(pr_average_weight)) + "," + str(
    #     np.std(pr_average_weight)) + "\n")
    # fp.write("average_predict_results\t" + ",".join([str(x) for x in roc_average_predict]) + "\t" + ",".join(
    #     [str(x) for x in pr_average_predict]) + "\t" + str(np.mean(roc_average_predict)) + "," + str(
    #     np.std(roc_average_predict)) + "\t" + str(np.mean(pr_average_predict)) + "," + str(
    #     np.std(pr_average_predict)) + "\n")
    # fp.write("average_last_predict_results\t" + ",".join([str(x) for x in roc_average_last_predict]) + "\t" + ",".join(
    #     [str(x) for x in pr_average_last_predict]) + "\t" + str(np.mean(roc_average_last_predict)) + "," + str(
    #     np.std(roc_average_last_predict)) + "\t" + str(np.mean(pr_average_last_predict)) + "," + str(
    #     np.std(pr_average_last_predict)) + "\n")
    #

    print("roc: \n")
    print(roc_average_last_predict)
    fp.write("average_last_predict_results: \t" + "\t" + str(np.mean(roc_average_last_predict)) + ","  + "\t" + str(np.mean(pr_average_last_predict)) + "," +str(np.mean(accuracy_average_last_predict))+","  + "\t" +
             str(np.mean(sensitivity_average_last_predict)) +","  + "\t" +str(np.mean(specificity_average_last_predict)) +","  + "\t"  +str(np.mean(f1_score_average_last_predict)) +","  + "\t"  +str(np.mean(mcc_average_last_predict)) 
         + "\n")
    fp.close()


if __name__ == "__main__":
    main()
